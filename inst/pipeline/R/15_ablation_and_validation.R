# =============================================================================
# 15_ablation_and_validation.R
# MosqEdit-R Manuscript 1
# STEP 15 â€” Evidence ablation, rank concordance, and validation analyses
# PURE R
# =============================================================================

source("R/helpers.R")

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

log_step("15", "Starting evidence ablation and validation analysis")

INTEGRATED_FILE <- "data_processed/09_targeted_integrated_prioritization.csv"
UNCERTAINTY_FILE <- "data_processed/12_bootstrap_rank_uncertainty.csv"
METRICS_FILE <- "data_processed/08_pu_metrics_summary.csv"
POP_FINAL_FILE <- "data_processed/14_final_population_validated_prioritization.csv"

GENE_ABLATION_FILE <- "data_processed/15_gene_level_ablation.csv"
SUMMARY_FILE <- "data_processed/15_ablation_summary.csv"
BENCHMARK_FILE <- "data_processed/15_benchmark_ablation.csv"
QC_FILE <- "data_processed/15_ablation_qc.csv"
PROVENANCE_FILE <- "data_processed/15_ablation_provenance.csv"
SESSION_INFO_FILE <- "logs/15_ablation_sessionInfo.txt"
CHECKSUM_FILE <- "logs/15_checksums.tsv"

FIG_DIR <- "figures"
FIG15A_PNG <- file.path(FIG_DIR, "15A_ablation_rank_concordance.png")
FIG15A_PDF <- file.path(FIG_DIR, "15A_ablation_rank_concordance.pdf")
FIG15B_PNG <- file.path(FIG_DIR, "15B_top50_overlap_by_ablation.png")
FIG15B_PDF <- file.path(FIG_DIR, "15B_top50_overlap_by_ablation.pdf")
FIG15C_PNG <- file.path(FIG_DIR, "15C_model_set_performance.png")
FIG15C_PDF <- file.path(FIG_DIR, "15C_model_set_performance.pdf")

dir.create("data_processed", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

for (path in c(INTEGRATED_FILE, UNCERTAINTY_FILE, METRICS_FILE)) {
    if (!file.exists(path)) {
        stop(paste0("Required file is missing:\n", path), call. = FALSE)
    }
}

x <- fread(INTEGRATED_FILE)
u <- fread(UNCERTAINTY_FILE)
metrics <- fread(METRICS_FILE)

if (uniqueN(x$gene_id) != nrow(x)) {
    stop("Step 09 integrated table contains duplicate gene IDs.", call. = FALSE)
}

score_available <- function(biology, pu, crispr, keep) {
    comp <- cbind(
        biology = biology,
        pu = pu,
        crispr = crispr
    )

    comp <- comp[, keep, drop = FALSE]

    w <- c(
        biology = 0.50,
        pu = 0.30,
        crispr = 0.20
    )[keep]

    observed <- is.finite(comp)
    wmat <- matrix(w, nrow = nrow(comp), ncol = length(w), byrow = TRUE)

    denom <- rowSums(observed * wmat)
    num <- rowSums(ifelse(observed, comp, 0) * wmat)

    raw <- num / denom
    raw[denom <= 0] <- NA_real_

    coverage <- rowMeans(observed)

    raw * (0.80 + 0.20 * coverage)
}

rank_score <- function(score, gene_id) {
    z <- data.table(
        idx = seq_along(gene_id),
        gene_id = gene_id,
        score = score
    )
    z[, sort_score := fifelse(is.finite(score), score, -Inf)]
    setorder(z, -sort_score, gene_id)
    z[, r := seq_len(.N)]

    out <- rep(NA_integer_, nrow(z))
    out[z$idx] <- z$r
    out[!is.finite(score)] <- NA_integer_
    out
}

biology <- x$biological_discovery_percentile
pu <- x$pu_consensus_percentile
crispr <- x$crispr_tractability_percentile

ablations <- list(
    full = c("biology", "pu", "crispr"),
    no_biology = c("pu", "crispr"),
    no_pu = c("biology", "crispr"),
    no_crispr = c("biology", "pu"),
    biology_only = c("biology"),
    pu_only = c("pu"),
    crispr_only = c("crispr")
)

gene_results <- list()
summary_results <- list()

full_rank <- x$final_prepopulation_rank

topk_overlap <- function(rank_a, rank_b, k) {
    a <- which(rank_a <= k)
    b <- which(rank_b <= k)

    if (length(a) == 0L || length(b) == 0L) {
        return(NA_real_)
    }

    length(intersect(a, b)) / k
}

counter <- 0L

for (name in names(ablations)) {
    keep <- ablations[[name]]

    score <- score_available(
        biology = biology,
        pu = pu,
        crispr = crispr,
        keep = keep
    )

    rank <- rank_score(score, x$gene_id)

    counter <- counter + 1L

    gene_results[[counter]] <- data.table(
        gene_id = x$gene_id,
        is_external_benchmark = x$is_external_benchmark,
        benchmark_name =
            if ("benchmark_name" %in% names(x)) {
                x$benchmark_name
            } else {
                NA_character_
            },
        ablation = name,
        ablation_score = score,
        ablation_rank = rank,
        frozen_final_rank = full_rank,
        rank_shift_vs_full = rank - full_rank
    )

    summary_results[[counter]] <- data.table(
        ablation = name,
        retained_domains = paste(keep, collapse = "+"),
        spearman_vs_full = suppressWarnings(
            cor(
                full_rank,
                rank,
                method = "spearman",
                use = "complete.obs"
            )
        ),
        top10_overlap = topk_overlap(full_rank, rank, 10L),
        top20_overlap = topk_overlap(full_rank, rank, 20L),
        top50_overlap = topk_overlap(full_rank, rank, 50L),
        top100_overlap = topk_overlap(full_rank, rank, 100L),
        benchmark_median_rank = median(
            rank[x$is_external_benchmark == TRUE],
            na.rm = TRUE
        )
    )
}

gene_ablation <- rbindlist(gene_results)
ablation_summary <- rbindlist(summary_results)

# Optional population second-stage comparison when Step 14 exists.
population_status <- "NOT_AVAILABLE"

if (file.exists(POP_FINAL_FILE)) {
    pop_final <- fread(POP_FINAL_FILE)
    population_status <- "AVAILABLE"

    pop_compare <- pop_final[
        ,
        .(
            gene_id,
            population_robustness_class,
            final_population_validated_rank
        )
    ]

    pop_compare <- merge(
        x[
            gene_id %in% pop_compare$gene_id,
            .(
                gene_id,
                final_prepopulation_rank
            )
        ],
        pop_compare,
        by = "gene_id",
        all.x = TRUE
    )

    population_summary <- data.table(
        ablation = "population_second_stage",
        retained_domains = "biology+PU+CRISPR then population robustness filter",
        spearman_vs_full = suppressWarnings(
            cor(
                pop_compare$final_prepopulation_rank,
                pop_compare$final_population_validated_rank,
                method = "spearman",
                use = "complete.obs"
            )
        ),
        top10_overlap = topk_overlap(
            pop_compare$final_prepopulation_rank,
            pop_compare$final_population_validated_rank,
            10L
        ),
        top20_overlap = topk_overlap(
            pop_compare$final_prepopulation_rank,
            pop_compare$final_population_validated_rank,
            20L
        ),
        top50_overlap = topk_overlap(
            pop_compare$final_prepopulation_rank,
            pop_compare$final_population_validated_rank,
            50L
        ),
        top100_overlap = topk_overlap(
            pop_compare$final_prepopulation_rank,
            pop_compare$final_population_validated_rank,
            100L
        ),
        benchmark_median_rank = NA_real_
    )

    ablation_summary <- rbind(
        ablation_summary,
        population_summary,
        fill = TRUE
    )
}

benchmark_ablation <- gene_ablation[
    is_external_benchmark == TRUE
]

fwrite(gene_ablation, GENE_ABLATION_FILE, na = "NA")
fwrite(ablation_summary, SUMMARY_FILE, na = "NA")
fwrite(benchmark_ablation, BENCHMARK_FILE, na = "NA")

# Figure 15A.
plot_a <- ablation_summary[
    ablation != "full"
]

p15a <- ggplot(
    plot_a,
    aes(
        x = reorder(ablation, spearman_vs_full),
        y = spearman_vs_full
    )
) +
    geom_col() +
    coord_flip() +
    labs(
        x = NULL,
        y = "Spearman correlation with frozen full rank",
        title = "Rank concordance after evidence-domain ablation"
    ) +
    theme_classic(base_size = 11)

ggsave(FIG15A_PNG, p15a, width = 8, height = 5.5, dpi = 400, bg = "white")
ggsave(FIG15A_PDF, p15a, width = 8, height = 5.5, device = cairo_pdf)

# Figure 15B.
overlap_long <- melt(
    plot_a,
    id.vars = c("ablation", "retained_domains"),
    measure.vars = c(
        "top10_overlap",
        "top20_overlap",
        "top50_overlap",
        "top100_overlap"
    ),
    variable.name = "top_k",
    value.name = "overlap"
)

p15b <- ggplot(
    overlap_long,
    aes(
        x = top_k,
        y = overlap,
        group = ablation,
        linetype = ablation
    )
) +
    geom_line() +
    geom_point() +
    labs(
        x = NULL,
        y = "Fraction retained from full ranking",
        title = "Top-k overlap under evidence-domain ablation"
    ) +
    theme_classic(base_size = 11)

ggsave(FIG15B_PNG, p15b, width = 8, height = 5.8, dpi = 400, bg = "white")
ggsave(FIG15B_PDF, p15b, width = 8, height = 5.8, device = cairo_pdf)

# Figure 15C: M1/M2/M3 performance.
p15c <- ggplot(
    metrics,
    aes(
        x = apparent_auroc_mean,
        y = apparent_auprc_mean,
        shape = algorithm
    )
) +
    geom_point(size = 3) +
    geom_text(
        aes(label = model_set),
        nudge_y = 0.002,
        check_overlap = TRUE,
        size = 3
    ) +
    labs(
        x = "Mean apparent AUROC",
        y = "Mean apparent AUPRC",
        title = "PU model performance across M1â€“M3 feature sets"
    ) +
    theme_classic(base_size = 11)

ggsave(FIG15C_PNG, p15c, width = 8, height = 6, dpi = 400, bg = "white")
ggsave(FIG15C_PDF, p15c, width = 8, height = 6, device = cairo_pdf)

qc <- data.table(
    metric = c(
        "Ablation scenarios",
        "Gene-level ablation rows",
        "Benchmark ablation rows",
        "Population second-stage available",
        "Figure 15A written",
        "Figure 15B written",
        "Figure 15C written"
    ),
    value = c(
        length(ablations),
        nrow(gene_ablation),
        nrow(benchmark_ablation),
        as.integer(population_status == "AVAILABLE"),
        as.integer(file.exists(FIG15A_PNG) && file.exists(FIG15A_PDF)),
        as.integer(file.exists(FIG15B_PNG) && file.exists(FIG15B_PDF)),
        as.integer(file.exists(FIG15C_PNG) && file.exists(FIG15C_PDF))
    )
)

provenance <- data.table(
    component = c(
        "Full model",
        "Ablation principle",
        "No-biology",
        "No-PU",
        "No-CRISPR",
        "Single-domain baselines",
        "Population analysis"
    ),
    specification = c(
        "Frozen Step 09 biology + PU + CRISPR ranking.",
        "Remove one evidence domain and renormalize over the remaining observed domains using the same completeness rule.",
        "PU + CRISPR only.",
        "Biology + CRISPR only.",
        "Biology + PU only.",
        "Biology-only, PU-only, and CRISPR-only rankings are included.",
        paste0(
            "Population second-stage comparison status: ",
            population_status,
            ". Population metrics are never back-fed into PU training."
        )
    )
)

fwrite(qc, QC_FILE, na = "NA")
fwrite(provenance, PROVENANCE_FILE, na = "NA")
capture.output(sessionInfo(), file = SESSION_INFO_FILE)

checksum_files <- c(
    INTEGRATED_FILE,
    UNCERTAINTY_FILE,
    METRICS_FILE,
    GENE_ABLATION_FILE,
    SUMMARY_FILE,
    BENCHMARK_FILE,
    QC_FILE,
    PROVENANCE_FILE,
    SESSION_INFO_FILE,
    FIG15A_PNG,
    FIG15A_PDF,
    FIG15B_PNG,
    FIG15B_PDF,
    FIG15C_PNG,
    FIG15C_PDF
)
if (file.exists(POP_FINAL_FILE)) {
    checksum_files <- c(checksum_files, POP_FINAL_FILE)
}
checksum_files <- unique(checksum_files[file.exists(checksum_files)])
write_checksum(checksum_files, CHECKSUM_FILE)

cat(
    "\n============================================================\n",
    "MOSQEDIT-R STEP 15 COMPLETED SUCCESSFULLY\n",
    "Evidence ablation and validation\n",
    "============================================================\n",
    "Ablation scenarios:                       ", length(ablations), "\n",
    "Population second-stage available:        ", population_status, "\n",
    "Publication figures written:              3 PNG + 3 PDF\n",
    "============================================================\n",
    sep = ""
)

print(ablation_summary)

log_step("15", "Evidence ablation and validation completed successfully")

