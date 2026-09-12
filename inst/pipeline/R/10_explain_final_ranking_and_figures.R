# =============================================================================
# 10_explain_final_ranking_and_figures.R
#
# MosqEdit-R Manuscript 1
#
# STEP 10
# Explain final pre-population ranking + publication-ready evidence figures
#
# PURE R
#
# PURPOSE
# -------
# Explain the Step 09 ranking without fitting any new predictive model.
#
# For every targeted gene this script:
#   - decomposes the primary integrated score into biology / PU / CRISPR
#     contributions that sum back to the final score;
#   - quantifies rank movement from Step 04 to Step 09;
#   - summarizes weighting-scenario sensitivity;
#   - generates publication-ready figures for model performance, evidence-layer
#     profiles, rank movement, and benchmark recovery.
#
# IMPORTANT
# ---------
# This remains PRE-POPULATION-GENOMICS.
# Ag1000G is still pending and is not inferred or imputed.
#
# INPUTS
# ------
# data_processed/09_targeted_integrated_prioritization.csv
# data_processed/09_novel_candidate_ranking.csv
# data_processed/09_benchmark_recovery.csv
# data_processed/09_selected_pu_models.csv
# data_processed/08_pu_metrics_summary.csv
#
# OUTPUTS
# -------
# data_processed/10_integrated_score_contributions.csv
# data_processed/10_top20_candidate_evidence.csv
# data_processed/10_benchmark_evidence.csv
# data_processed/10_rank_shift_summary.csv
# data_processed/10_explainability_qc.csv
# data_processed/10_explainability_provenance.csv
#
# figures/10A_pu_model_performance.png/.pdf
# figures/10B_top20_evidence_heatmap.png/.pdf
# figures/10C_rank_shift_top50.png/.pdf
# figures/10D_benchmark_evidence_heatmap.png/.pdf
#
# logs/10_explainability_sessionInfo.txt
# logs/10_checksums.tsv
#
# =============================================================================

source("R/helpers.R")

required_packages <- c("data.table", "ggplot2")
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
    stop(
        paste0(
            "Missing required package(s): ",
            paste(missing_packages, collapse = ", "),
            "\nInstall with:\ninstall.packages(c(",
            paste0('"', missing_packages, '"', collapse = ", "),
            "))"
        ),
        call. = FALSE
    )
}

suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
})

log_step(
    "10",
    "Starting integrated-ranking explainability and publication figures"
)

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------

INTEGRATED_FILE <- "data_processed/09_targeted_integrated_prioritization.csv"
NOVEL_FILE <- "data_processed/09_novel_candidate_ranking.csv"
BENCHMARK_FILE <- "data_processed/09_benchmark_recovery.csv"
SELECTED_MODELS_FILE <- "data_processed/09_selected_pu_models.csv"
METRICS_FILE <- "data_processed/08_pu_metrics_summary.csv"

CONTRIBUTIONS_FILE <- "data_processed/10_integrated_score_contributions.csv"
TOP20_FILE <- "data_processed/10_top20_candidate_evidence.csv"
BENCHMARK_EVIDENCE_FILE <- "data_processed/10_benchmark_evidence.csv"
RANK_SHIFT_FILE <- "data_processed/10_rank_shift_summary.csv"
QC_FILE <- "data_processed/10_explainability_qc.csv"
PROVENANCE_FILE <- "data_processed/10_explainability_provenance.csv"
SESSION_INFO_FILE <- "logs/10_explainability_sessionInfo.txt"
CHECKSUM_FILE <- "logs/10_checksums.tsv"

FIG_DIR <- "figures"

FIG10A_PNG <- file.path(FIG_DIR, "10A_pu_model_performance.png")
FIG10A_PDF <- file.path(FIG_DIR, "10A_pu_model_performance.pdf")

FIG10B_PNG <- file.path(FIG_DIR, "10B_top20_evidence_heatmap.png")
FIG10B_PDF <- file.path(FIG_DIR, "10B_top20_evidence_heatmap.pdf")

FIG10C_PNG <- file.path(FIG_DIR, "10C_rank_shift_top50.png")
FIG10C_PDF <- file.path(FIG_DIR, "10C_rank_shift_top50.pdf")

FIG10D_PNG <- file.path(FIG_DIR, "10D_benchmark_evidence_heatmap.png")
FIG10D_PDF <- file.path(FIG_DIR, "10D_benchmark_evidence_heatmap.pdf")

PRIMARY_WEIGHTS <- c(
    biology = 0.50,
    pu = 0.30,
    crispr = 0.20
)

COMPLETENESS_FLOOR <- 0.80

TOP20_N <- 20L
TOP50_N <- 50L

dir.create("data_processed", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Helpers
# -----------------------------------------------------------------------------

assert_file <- function(path) {
    if (!file.exists(path)) {
        stop(
            paste0("Required input file is missing:\n", path),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

assert_unique_gene_ids <- function(dt, label) {
    if (!"gene_id" %in% names(dt)) {
        stop(paste0(label, " lacks gene_id."), call. = FALSE)
    }

    if (uniqueN(dt$gene_id) != nrow(dt)) {
        stop(
            paste0(label, " contains duplicate gene_id values."),
            call. = FALSE
        )
    }

    invisible(TRUE)
}

safe_label <- function(gene_id, gene_name = NULL) {
    if (is.null(gene_name)) {
        return(gene_id)
    }

    out <- gene_id
    ok <- !is.na(gene_name) & nzchar(gene_name)

    out[ok] <- paste0(
        gene_id[ok],
        " (",
        gene_name[ok],
        ")"
    )

    out
}

save_plot_pair <- function(
    plot_object,
    png_file,
    pdf_file,
    width,
    height
) {
    ggsave(
        filename = png_file,
        plot = plot_object,
        width = width,
        height = height,
        units = "in",
        dpi = 400,
        bg = "white"
    )

    ggsave(
        filename = pdf_file,
        plot = plot_object,
        width = width,
        height = height,
        units = "in",
        device = cairo_pdf
    )
}

# -----------------------------------------------------------------------------
# 3. Load inputs
# -----------------------------------------------------------------------------

for (path in c(
    INTEGRATED_FILE,
    NOVEL_FILE,
    BENCHMARK_FILE,
    SELECTED_MODELS_FILE,
    METRICS_FILE
)) {
    assert_file(path)
}

integrated <- fread(INTEGRATED_FILE)
novel <- fread(NOVEL_FILE)
benchmarks <- fread(BENCHMARK_FILE)
selected_models <- fread(SELECTED_MODELS_FILE)
metrics <- fread(METRICS_FILE)

assert_unique_gene_ids(integrated, "Step 09 integrated prioritization")
assert_unique_gene_ids(novel, "Step 09 novel ranking")
assert_unique_gene_ids(benchmarks, "Step 09 benchmark recovery")

# -----------------------------------------------------------------------------
# 4. Explain the primary integrated score
# -----------------------------------------------------------------------------

required_domains <- c(
    "biological_discovery_percentile",
    "pu_consensus_percentile",
    "crispr_tractability_percentile",
    "final_integrated_score",
    "final_evidence_domain_coverage",
    "final_prepopulation_rank"
)

missing_domains <- setdiff(
    required_domains,
    names(integrated)
)

if (length(missing_domains) > 0L) {
    stop(
        paste0(
            "Integrated table is missing required column(s): ",
            paste(missing_domains, collapse = ", ")
        ),
        call. = FALSE
    )
}

explain <- copy(integrated)

components <- as.matrix(
    explain[
        ,
        .(
            biological_discovery_percentile,
            pu_consensus_percentile,
            crispr_tractability_percentile
        )
    ]
)

colnames(components) <- c(
    "biology",
    "pu",
    "crispr"
)

observed <- is.finite(components)

weight_matrix <- matrix(
    PRIMARY_WEIGHTS[
        colnames(components)
    ],
    nrow = nrow(components),
    ncol = ncol(components),
    byrow = TRUE
)

observed_weight <- rowSums(
    observed * weight_matrix
)

domain_coverage <- rowMeans(
    observed
)

completeness_factor <-
    COMPLETENESS_FLOOR +
    (1 - COMPLETENESS_FLOOR) *
    domain_coverage

contribution_matrix <- matrix(
    0,
    nrow = nrow(components),
    ncol = ncol(components)
)

colnames(contribution_matrix) <- colnames(components)

for (j in seq_len(ncol(components))) {
    contribution_matrix[, j] <- ifelse(
        observed[, j] &
            observed_weight > 0,
        components[, j] *
            weight_matrix[, j] /
            observed_weight *
            completeness_factor,
        0
    )
}

explain[
    ,
    contribution_biology :=
        contribution_matrix[, "biology"]
]

explain[
    ,
    contribution_pu :=
        contribution_matrix[, "pu"]
]

explain[
    ,
    contribution_crispr :=
        contribution_matrix[, "crispr"]
]

explain[
    ,
    reconstructed_final_score :=
        contribution_biology +
        contribution_pu +
        contribution_crispr
]

explain[
    ,
    score_reconstruction_error :=
        abs(
            reconstructed_final_score -
            final_integrated_score
        )
]

explain[
    ,
    dominant_final_driver :=
        {
            mat <- cbind(
                contribution_biology,
                contribution_pu,
                contribution_crispr
            )

            labels <- c(
                "biology",
                "PU",
                "CRISPR"
            )

            apply(
                mat,
                1,
                function(z) {
                    if (!any(is.finite(z))) {
                        return(NA_character_)
                    }
                    labels[[which.max(z)]]
                }
            )
        }
]

# -----------------------------------------------------------------------------
# 5. Rank movement and robustness
# -----------------------------------------------------------------------------

if ("preliminary_rank" %in% names(explain)) {
    explain[
        ,
        rank_change_from_step04 :=
            preliminary_rank -
            final_prepopulation_rank
    ]
} else {
    explain[
        ,
        rank_change_from_step04 :=
            NA_integer_
    ]
}

if ("final_rank_sensitivity_range" %in% names(explain)) {
    explain[
        ,
        weighting_robustness :=
            fcase(
                is.na(final_rank_sensitivity_range),
                "not_estimable",

                final_rank_sensitivity_range == 0,
                "invariant",

                final_rank_sensitivity_range <= 5,
                "highly_stable",

                final_rank_sensitivity_range <= 15,
                "moderately_stable",

                default = "sensitive"
            )
    ]
} else {
    explain[
        ,
        weighting_robustness :=
            "not_estimable"
    ]
}

if ("pu_repeat_sd_percentile" %in% names(explain)) {
    explain[
        ,
        pu_stability_class :=
            fcase(
                is.na(pu_repeat_sd_percentile),
                "not_estimable",

                pu_repeat_sd_percentile <= 0.01,
                "very_stable",

                pu_repeat_sd_percentile <= 0.03,
                "stable",

                pu_repeat_sd_percentile <= 0.06,
                "moderately_stable",

                default = "variable"
            )
    ]
} else {
    explain[
        ,
        pu_stability_class :=
            "not_estimable"
    ]
}

# -----------------------------------------------------------------------------
# 6. Output tables
# -----------------------------------------------------------------------------

setorder(
    explain,
    final_prepopulation_rank,
    gene_id,
    na.last = TRUE
)

top20 <- explain[
    is_external_benchmark == FALSE
][
    1:min(TOP20_N, .N)
]

benchmark_evidence <- explain[
    is_external_benchmark == TRUE
]

rank_shift <- explain[
    is_external_benchmark == FALSE &
        !is.na(preliminary_rank) &
        !is.na(final_prepopulation_rank)
][
    ,
    .(
        gene_id,
        gene_name =
            if ("gene_name" %in% names(explain)) {
                gene_name
            } else {
                NA_character_
            },
        preliminary_rank,
        final_prepopulation_rank,
        rank_change_from_step04,
        final_integrated_score,
        final_rank_sensitivity_range,
        weighting_robustness,
        pu_repeat_sd_percentile,
        pu_stability_class
    )
]

setorder(
    rank_shift,
    final_prepopulation_rank
)

fwrite(
    explain,
    CONTRIBUTIONS_FILE,
    na = "NA"
)

fwrite(
    top20,
    TOP20_FILE,
    na = "NA"
)

fwrite(
    benchmark_evidence,
    BENCHMARK_EVIDENCE_FILE,
    na = "NA"
)

fwrite(
    rank_shift,
    RANK_SHIFT_FILE,
    na = "NA"
)

# -----------------------------------------------------------------------------
# 7. Figure 10A â€” Step 08 model performance
# -----------------------------------------------------------------------------

metrics_plot <- copy(metrics)

metrics_plot[
    ,
    model_label :=
        paste0(
            model_set,
            "\n",
            algorithm
        )
]

selected_keys <- paste(
    selected_models$model_set,
    selected_models$algorithm,
    sep = "::"
)

metrics_plot[
    ,
    selected_for_pu_consensus :=
        paste(
            model_set,
            algorithm,
            sep = "::"
        ) %in%
        selected_keys
]

p10a <- ggplot(
    metrics_plot,
    aes(
        x = apparent_auroc_mean,
        y = apparent_auprc_mean
    )
) +
    geom_point(
        aes(
            shape = selected_for_pu_consensus,
            size = ndcg_at_100_mean
        )
    ) +
    geom_text(
        aes(
            label = paste0(
                algorithm,
                "\n",
                sub(
                    "_expression.*$",
                    "",
                    model_set
                )
            )
        ),
        nudge_y = 0.002,
        check_overlap = TRUE,
        size = 3
    ) +
    labs(
        x = "Mean apparent AUROC",
        y = "Mean apparent AUPRC",
        shape = "PU consensus model",
        size = "Mean NDCG@100",
        title = "Cross-validated positiveâ€“unlabeled model performance",
        subtitle = "Observed-positive versus unlabeled discrimination"
    ) +
    theme_classic(base_size = 11) +
    theme(
        legend.position = "right"
    )

save_plot_pair(
    p10a,
    FIG10A_PNG,
    FIG10A_PDF,
    width = 8.2,
    height = 6.2
)

# -----------------------------------------------------------------------------
# 8. Figure 10B â€” top-20 evidence heatmap
# -----------------------------------------------------------------------------

top20_plot <- copy(top20)

top20_plot[
    ,
    display_label :=
        safe_label(
            gene_id,
            if ("gene_name" %in% names(top20_plot)) {
                gene_name
            } else {
                NULL
            }
        )
]

top20_plot[
    ,
    display_label :=
        factor(
            display_label,
            levels = rev(display_label)
        )
]

top20_long <- melt(
    top20_plot,
    id.vars = c(
        "gene_id",
        "display_label",
        "final_prepopulation_rank"
    ),
    measure.vars = c(
        "biological_discovery_percentile",
        "pu_consensus_percentile",
        "crispr_tractability_percentile"
    ),
    variable.name = "evidence_domain",
    value.name = "percentile"
)

top20_long[
    ,
    evidence_domain :=
        factor(
            evidence_domain,
            levels = c(
                "biological_discovery_percentile",
                "pu_consensus_percentile",
                "crispr_tractability_percentile"
            ),
            labels = c(
                "Biological discovery",
                "PU consensus",
                "CRISPR tractability"
            )
        )
]

p10b <- ggplot(
    top20_long,
    aes(
        x = evidence_domain,
        y = display_label,
        fill = percentile
    )
) +
    geom_tile() +
    geom_text(
        aes(
            label = sprintf(
                "%.2f",
                percentile
            )
        ),
        size = 3
    ) +
    scale_fill_viridis_c(
        limits = c(0, 1),
        option = "C"
    ) +
    labs(
        x = NULL,
        y = NULL,
        fill = "Percentile",
        title = "Evidence profile of the top 20 novel MosqEdit-R candidates",
        subtitle = "Final pre-population-genomics ranking"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = 20,
            hjust = 1
        )
    )

save_plot_pair(
    p10b,
    FIG10B_PNG,
    FIG10B_PDF,
    width = 8.5,
    height = 8.5
)

# -----------------------------------------------------------------------------
# 9. Figure 10C â€” rank movement among top-50 novel candidates
# -----------------------------------------------------------------------------

rank_top50 <- explain[
    is_external_benchmark == FALSE &
        !is.na(final_prepopulation_rank) &
        final_prepopulation_rank <= TOP50_N &
        !is.na(preliminary_rank)
]

rank_top50[
    ,
    display_label :=
        safe_label(
            gene_id,
            if ("gene_name" %in% names(rank_top50)) {
                gene_name
            } else {
                NULL
            }
        )
]

p10c <- ggplot(
    rank_top50,
    aes(
        x = final_prepopulation_rank,
        y = rank_change_from_step04
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = 2
    ) +
    geom_point(
        aes(
            size = final_rank_sensitivity_range
        )
    ) +
    geom_text(
        data = rank_top50[
            final_prepopulation_rank <= 15
        ],
        aes(
            label = gene_id
        ),
        nudge_y = 4,
        check_overlap = TRUE,
        size = 3
    ) +
    labs(
        x = "Final pre-population rank",
        y = "Rank change from Step 04\n(positive = promoted after PU/CRISPR integration)",
        size = "Weighting\nsensitivity range",
        title = "How PU support and CRISPR tractability changed candidate priority"
    ) +
    theme_classic(base_size = 11)

save_plot_pair(
    p10c,
    FIG10C_PNG,
    FIG10C_PDF,
    width = 8.3,
    height = 6.3
)

# -----------------------------------------------------------------------------
# 10. Figure 10D â€” external benchmark evidence heatmap
# -----------------------------------------------------------------------------

bench_plot <- copy(benchmark_evidence)

bench_plot[
    ,
    display_label :=
        if ("benchmark_name" %in% names(bench_plot)) {
            paste0(
                gene_id,
                "\n",
                benchmark_name
            )
        } else {
            gene_id
        }
]

bench_plot[
    ,
    display_label :=
        factor(
            display_label,
            levels = rev(display_label)
        )
]

bench_long <- melt(
    bench_plot,
    id.vars = c(
        "gene_id",
        "display_label",
        "final_prepopulation_rank"
    ),
    measure.vars = c(
        "biological_discovery_percentile",
        "pu_consensus_percentile",
        "crispr_tractability_percentile"
    ),
    variable.name = "evidence_domain",
    value.name = "percentile"
)

bench_long[
    ,
    evidence_domain :=
        factor(
            evidence_domain,
            levels = c(
                "biological_discovery_percentile",
                "pu_consensus_percentile",
                "crispr_tractability_percentile"
            ),
            labels = c(
                "Biological discovery",
                "PU consensus",
                "CRISPR tractability"
            )
        )
]

p10d <- ggplot(
    bench_long,
    aes(
        x = evidence_domain,
        y = display_label,
        fill = percentile
    )
) +
    geom_tile() +
    geom_text(
        aes(
            label = ifelse(
                is.na(percentile),
                "NA",
                sprintf("%.2f", percentile)
            )
        ),
        size = 3.2
    ) +
    scale_fill_viridis_c(
        limits = c(0, 1),
        option = "C",
        na.value = "grey90"
    ) +
    labs(
        x = NULL,
        y = NULL,
        fill = "Percentile",
        title = "External benchmark recovery across evidence layers",
        subtitle = "Benchmarks were not used for model selection or weight tuning"
    ) +
    theme_minimal(base_size = 11) +
    theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(
            angle = 20,
            hjust = 1
        )
    )

save_plot_pair(
    p10d,
    FIG10D_PNG,
    FIG10D_PDF,
    width = 8.7,
    height = 5.8
)

# -----------------------------------------------------------------------------
# 11. QC
# -----------------------------------------------------------------------------

max_reconstruction_error <- max(
    explain$score_reconstruction_error,
    na.rm = TRUE
)

n_invariant <- explain[
    weighting_robustness == "invariant",
    .N
]

n_highly_stable <- explain[
    weighting_robustness %in% c(
        "invariant",
        "highly_stable"
    ),
    .N
]

n_top20_three_domains <- top20[
    final_evidence_domain_coverage == 1,
    .N
]

qc <- data.table(
    metric = c(
        "Integrated targeted genes",
        "External benchmark genes",
        "Novel genes",
        "Top-20 novel genes",
        "Top-20 with all three evidence domains",
        "Maximum score reconstruction error",
        "Genes invariant across weighting scenarios",
        "Genes invariant or highly stable across weighting scenarios",
        "Figure 10A written",
        "Figure 10B written",
        "Figure 10C written",
        "Figure 10D written",
        "Population-genomics information added"
    ),

    value = c(
        nrow(explain),
        sum(explain$is_external_benchmark, na.rm = TRUE),
        sum(!explain$is_external_benchmark, na.rm = TRUE),
        nrow(top20),
        n_top20_three_domains,
        max_reconstruction_error,
        n_invariant,
        n_highly_stable,
        as.integer(
            file.exists(FIG10A_PNG) &
                file.exists(FIG10A_PDF)
        ),
        as.integer(
            file.exists(FIG10B_PNG) &
                file.exists(FIG10B_PDF)
        ),
        as.integer(
            file.exists(FIG10C_PNG) &
                file.exists(FIG10C_PDF)
        ),
        as.integer(
            file.exists(FIG10D_PNG) &
                file.exists(FIG10D_PDF)
        ),
        0
    )
)

if (
    is.finite(max_reconstruction_error) &&
    max_reconstruction_error > 1e-8
) {
    stop(
        paste0(
            "Integrated-score reconstruction failed. Maximum absolute error = ",
            signif(max_reconstruction_error, 6)
        ),
        call. = FALSE
    )
}

# -----------------------------------------------------------------------------
# 12. Provenance
# -----------------------------------------------------------------------------

provenance <- data.table(
    component = c(
        "Score explanation",
        "Biology contribution",
        "PU contribution",
        "CRISPR contribution",
        "Completeness adjustment",
        "Rank movement",
        "Weighting robustness",
        "PU stability",
        "Figure 10A",
        "Figure 10B",
        "Figure 10C",
        "Figure 10D",
        "Population genomics"
    ),

    specification = c(
        paste0(
            "Step 09 primary integrated score is algebraically decomposed into ",
            "domain-specific contributions that sum to the stored final score."
        ),

        "50% nominal primary weight before available-domain renormalization.",

        "30% nominal primary weight before available-domain renormalization.",

        "20% nominal primary weight before available-domain renormalization.",

        paste0(
            "Available-domain score multiplied by [",
            COMPLETENESS_FLOOR,
            " + ",
            1 - COMPLETENESS_FLOOR,
            " * domain coverage]."
        ),

        paste0(
            "rank_change_from_step04 = preliminary_rank - ",
            "final_prepopulation_rank; positive values indicate promotion."
        ),

        paste0(
            "Sensitivity classes are based on final_rank_sensitivity_range: ",
            "0=invariant; <=5=highly stable; <=15=moderately stable; >15=sensitive."
        ),

        paste0(
            "PU stability classes use repeat-level PU percentile SD: <=0.01 very ",
            "stable; <=0.03 stable; <=0.06 moderately stable; >0.06 variable."
        ),

        "Scatter plot of Step 08 apparent AUPRC versus apparent AUROC.",

        "Heatmap of biology, PU, and CRISPR percentiles for top-20 novel candidates.",

        "Rank movement from Step 04 to Step 09 among top-50 novel candidates.",

        "Heatmap of external benchmark evidence profiles.",

        "Ag1000G population-genomic information remains pending and is not added."
    )
)

fwrite(
    qc,
    QC_FILE,
    na = "NA"
)

fwrite(
    provenance,
    PROVENANCE_FILE,
    na = "NA"
)

capture.output(
    sessionInfo(),
    file = SESSION_INFO_FILE
)

checksum_files <- c(
    INTEGRATED_FILE,
    NOVEL_FILE,
    BENCHMARK_FILE,
    SELECTED_MODELS_FILE,
    METRICS_FILE,
    CONTRIBUTIONS_FILE,
    TOP20_FILE,
    BENCHMARK_EVIDENCE_FILE,
    RANK_SHIFT_FILE,
    QC_FILE,
    PROVENANCE_FILE,
    SESSION_INFO_FILE,
    FIG10A_PNG,
    FIG10A_PDF,
    FIG10B_PNG,
    FIG10B_PDF,
    FIG10C_PNG,
    FIG10C_PDF,
    FIG10D_PNG,
    FIG10D_PDF
)

checksum_files <- checksum_files[
    file.exists(checksum_files)
]

write_checksum(
    checksum_files,
    CHECKSUM_FILE
)

# -----------------------------------------------------------------------------
# 13. Console summary
# -----------------------------------------------------------------------------

cat(
    "\n",
    "============================================================\n",
    "MOSQEDIT-R STEP 10 COMPLETED SUCCESSFULLY\n",
    "Integrated ranking explainability + publication figures\n",
    "============================================================\n",
    "Integrated targeted genes:                 ",
    nrow(explain),
    "\n",
    "Novel genes:                               ",
    sum(!explain$is_external_benchmark, na.rm = TRUE),
    "\n",
    "External benchmark genes:                  ",
    sum(explain$is_external_benchmark, na.rm = TRUE),
    "\n",
    "Top-20 novel genes explained:              ",
    nrow(top20),
    "\n",
    "Top-20 with all three evidence domains:    ",
    n_top20_three_domains,
    "\n",
    "Max score reconstruction error:            ",
    format(
        max_reconstruction_error,
        scientific = TRUE
    ),
    "\n",
    "Invariant weighting-scenario ranks:         ",
    n_invariant,
    "\n",
    "Invariant/highly-stable ranks:              ",
    n_highly_stable,
    "\n",
    "Publication figures written:               4 PNG + 4 PDF\n",
    "Population genomics integrated:            NO\n",
    "Python used:                               NO\n",
    "============================================================\n",
    sep = ""
)

cat(
    "\nTop 20 novel candidates with evidence contributions:\n"
)

display_cols <- c(
    "final_prepopulation_rank",
    "gene_id",
    "gene_name",
    "preliminary_rank",
    "rank_change_from_step04",
    "biological_discovery_percentile",
    "pu_consensus_percentile",
    "crispr_tractability_percentile",
    "contribution_biology",
    "contribution_pu",
    "contribution_crispr",
    "final_integrated_score",
    "final_rank_sensitivity_range",
    "weighting_robustness",
    "pu_repeat_sd_percentile",
    "pu_stability_class"
)

display_cols <- display_cols[
    display_cols %in% names(top20)
]

print(
    top20[
        ,
        ..display_cols
    ]
)

cat(
    "\nExternal benchmark explainability:\n"
)

benchmark_display <- c(
    "gene_id",
    "benchmark_name",
    "preliminary_rank",
    "final_prepopulation_rank",
    "rank_change_from_step04",
    "biological_discovery_percentile",
    "pu_consensus_percentile",
    "crispr_tractability_percentile",
    "contribution_biology",
    "contribution_pu",
    "contribution_crispr",
    "final_integrated_score",
    "final_rank_sensitivity_range",
    "weighting_robustness"
)

benchmark_display <- benchmark_display[
    benchmark_display %in% names(benchmark_evidence)
]

print(
    benchmark_evidence[
        ,
        ..benchmark_display
    ]
)

cat(
    "\nQC summary:\n"
)

print(qc)

log_step(
    "10",
    "Integrated-ranking explainability and publication figures completed successfully"
)

