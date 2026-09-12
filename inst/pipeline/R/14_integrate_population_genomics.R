# =============================================================================
# 14_integrate_population_genomics.R
# MosqEdit-R Manuscript 1
# STEP 14 â€” Integrate REAL targeted population-genomic validation
# PURE R
# =============================================================================
#
# This script intentionally refuses to fabricate population results.
# It runs only after the Step 13 population-site template has been populated
# from real population-genomic analysis.
#
# Population genomics is treated as a TARGETED SECOND-STAGE ROBUSTNESS FILTER,
# not as a genome-wide PU predictor.
# =============================================================================

source("R/helpers.R")

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
}
suppressPackageStartupMessages(library(data.table))

log_step("14", "Integrating targeted population-genomic robustness")

UNCERTAINTY_FILE <- "data_processed/12_bootstrap_rank_uncertainty.csv"
PANEL_FILE <- "data_processed/13_population_validation_panel.csv"

# User should overwrite/copy the Step 13 template to this filename only after
# real population metrics have been populated.
POPULATION_METRICS_FILE <- "data_processed/13_population_site_metrics.csv"

GENE_POP_FILE <- "data_processed/14_gene_population_robustness.csv"
FINAL_FILE <- "data_processed/14_final_population_validated_prioritization.csv"
BENCHMARK_FILE <- "data_processed/14_population_validated_benchmarks.csv"
TOP50_FILE <- "data_processed/14_final_top50_novel.csv"
QC_FILE <- "data_processed/14_population_integration_qc.csv"
PROVENANCE_FILE <- "data_processed/14_population_integration_provenance.csv"
SESSION_INFO_FILE <- "logs/14_population_integration_sessionInfo.txt"
CHECKSUM_FILE <- "logs/14_checksums.tsv"

dir.create("data_processed", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

for (path in c(UNCERTAINTY_FILE, PANEL_FILE, POPULATION_METRICS_FILE)) {
    if (!file.exists(path)) {
        stop(
            paste0(
                "Required file is missing:\n",
                path,
                "\n\nStep 14 must not run until real population metrics exist."
            ),
            call. = FALSE
        )
    }
}

u <- fread(UNCERTAINTY_FILE)
panel <- fread(PANEL_FILE)
pop <- fread(POPULATION_METRICS_FILE)

required_pop <- c(
    "population_site_id",
    "gene_id",
    "population_callable_fraction",
    "pam_intact_fraction",
    "protospacer_exact_match_fraction",
    "target_23bp_exact_match_fraction",
    "max_target_variant_alt_af",
    "population_target_site_robustness_score"
)

missing_pop <- setdiff(required_pop, names(pop))
if (length(missing_pop) > 0L) {
    stop(
        paste0(
            "Population metrics file is missing: ",
            paste(missing_pop, collapse = ", ")
        ),
        call. = FALSE
    )
}

# Require genuine populated metrics.
n_populated <- pop[
    is.finite(population_target_site_robustness_score),
    .N
]

if (n_populated == 0L) {
    stop(
        paste0(
            "No real population robustness scores are populated. ",
            "Step 14 refuses to substitute or fabricate them."
        ),
        call. = FALSE
    )
}

# Validate range.
fraction_cols <- intersect(
    c(
        "population_callable_fraction",
        "pam_intact_fraction",
        "protospacer_exact_match_fraction",
        "target_23bp_exact_match_fraction",
        "max_target_variant_alt_af",
        "population_target_site_robustness_score"
    ),
    names(pop)
)

for (col in fraction_cols) {
    bad <- pop[
        !is.na(get(col)) &
        (!is.finite(get(col)) | get(col) < 0 | get(col) > 1),
        .N
    ]

    if (bad > 0L) {
        stop(
            paste0(
                col,
                " contains ",
                bad,
                " values outside [0,1]."
            ),
            call. = FALSE
        )
    }
}

# -------------------------------------------------------------------------
# Aggregate site robustness to gene level.
# Editing feasibility is driven by the availability of at least one robust
# target site, while multiple robust sites provide redundancy.
# -------------------------------------------------------------------------

pop_valid <- pop[
    is.finite(population_target_site_robustness_score)
]

setorder(
    pop_valid,
    gene_id,
    -population_target_site_robustness_score,
    population_site_id
)

pop_valid[
    ,
    population_site_rank :=
        seq_len(.N),
    by = gene_id
]

gene_pop <- pop_valid[
    ,
    {
        scores <- population_target_site_robustness_score
        scores <- scores[is.finite(scores)]

        best <- if (length(scores) >= 1L) scores[[1]] else NA_real_
        second <- if (length(scores) >= 2L) scores[[2]] else NA_real_
        top3 <- head(scores, 3L)

        redundancy_component <- if (length(top3) >= 2L) {
            median(top3)
        } else {
            best
        }

        gene_score <- if (is.finite(best)) {
            0.80 * best + 0.20 * redundancy_component
        } else {
            NA_real_
        }

        list(
            n_population_sites_assessed = length(scores),
            best_population_site_robustness = best,
            second_best_population_site_robustness = second,
            median_top3_population_site_robustness =
                if (length(top3) > 0L) median(top3) else NA_real_,
            gene_population_robustness_score = gene_score,
            best_population_site_id =
                population_site_id[[1]],
            best_target_23bp_exact_match_fraction =
                target_23bp_exact_match_fraction[[1]],
            best_pam_intact_fraction =
                pam_intact_fraction[[1]],
            best_max_target_variant_alt_af =
                max_target_variant_alt_af[[1]],
            best_population_callable_fraction =
                population_callable_fraction[[1]]
        )
    },
    by = gene_id
]

gene_pop[
    ,
    population_robustness_class :=
        fcase(
            gene_population_robustness_score >= 0.95 &
                best_target_23bp_exact_match_fraction >= 0.95 &
                best_pam_intact_fraction >= 0.99 &
                best_max_target_variant_alt_af <= 0.05,
            "ROBUST",

            gene_population_robustness_score >= 0.80 &
                best_target_23bp_exact_match_fraction >= 0.80,
            "INTERMEDIATE",

            default = "CONCERN"
        )
]

# -------------------------------------------------------------------------
# Final targeted population-validated ordering.
# Population genomics acts first as a robustness class, then as a within-class
# tie-breaker. The frozen Step 09 rank remains the biological priority anchor.
# -------------------------------------------------------------------------

final <- merge(
    panel,
    gene_pop,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)

class_order <- c(
    "ROBUST" = 1L,
    "INTERMEDIATE" = 2L,
    "CONCERN" = 3L
)

final[
    ,
    population_class_order :=
        unname(
            class_order[
                population_robustness_class
            ]
        )
]

final[
    is.na(population_class_order),
    population_class_order := 4L
]

setorder(
    final,
    population_class_order,
    final_prepopulation_rank,
    -gene_population_robustness_score,
    gene_id
)

final[
    ,
    final_population_validated_rank :=
        seq_len(.N)
]

final[
    ,
    population_rank_shift :=
        final_prepopulation_rank -
        final_population_validated_rank
]

benchmark_final <- final[
    is_external_benchmark == TRUE
]

top50 <- final[
    is_external_benchmark == FALSE
][
    1:min(50L, .N)
]

fwrite(gene_pop, GENE_POP_FILE, na = "NA")
fwrite(final, FINAL_FILE, na = "NA")
fwrite(benchmark_final, BENCHMARK_FILE, na = "NA")
fwrite(top50, TOP50_FILE, na = "NA")

qc <- data.table(
    metric = c(
        "Population panel genes",
        "Genes with population metrics",
        "Genes classified ROBUST",
        "Genes classified INTERMEDIATE",
        "Genes classified CONCERN",
        "Genes lacking population assessment",
        "External benchmarks",
        "Final top-50 novel genes"
    ),
    value = c(
        nrow(final),
        sum(is.finite(final$gene_population_robustness_score)),
        final[population_robustness_class == "ROBUST", .N],
        final[population_robustness_class == "INTERMEDIATE", .N],
        final[population_robustness_class == "CONCERN", .N],
        final[is.na(gene_population_robustness_score), .N],
        nrow(benchmark_final),
        nrow(top50)
    )
)

provenance <- data.table(
    component = c(
        "Population scope",
        "Site-level input",
        "Gene-level aggregation",
        "Population robustness class",
        "Final ordering",
        "PU model leakage",
        "Interpretation"
    ),
    specification = c(
        "Targeted validation panel only; not a genome-wide predictor.",
        "Only real populated Step 13 site-level population metrics are accepted.",
        "80% best-site robustness + 20% median of up to the top three sites.",
        paste0(
            "ROBUST requires gene score >=0.95, best 23-bp exact-match fraction >=0.95, ",
            "PAM intact fraction >=0.99, and max target-site alternate AF <=0.05. ",
            "INTERMEDIATE requires gene score >=0.80 and exact-match fraction >=0.80."
        ),
        paste0(
            "Population robustness class is applied as a second-stage filter; ",
            "frozen pre-population rank remains the within-class priority anchor."
        ),
        "Population metrics are never fed back into Step 08 PU model training.",
        "Final_population_validated_rank is valid only for the targeted Step 13 panel."
    )
)

fwrite(qc, QC_FILE, na = "NA")
fwrite(provenance, PROVENANCE_FILE, na = "NA")
capture.output(sessionInfo(), file = SESSION_INFO_FILE)

checksum_files <- c(
    UNCERTAINTY_FILE,
    PANEL_FILE,
    POPULATION_METRICS_FILE,
    GENE_POP_FILE,
    FINAL_FILE,
    BENCHMARK_FILE,
    TOP50_FILE,
    QC_FILE,
    PROVENANCE_FILE,
    SESSION_INFO_FILE
)
checksum_files <- checksum_files[file.exists(checksum_files)]
write_checksum(checksum_files, CHECKSUM_FILE)

cat(
    "\n============================================================\n",
    "MOSQEDIT-R STEP 14 COMPLETED SUCCESSFULLY\n",
    "Targeted population-genomic validation integrated\n",
    "============================================================\n",
    "Population panel genes:                   ", nrow(final), "\n",
    "Genes with population metrics:            ",
    sum(is.finite(final$gene_population_robustness_score)), "\n",
    "ROBUST:                                   ",
    final[population_robustness_class == "ROBUST", .N], "\n",
    "INTERMEDIATE:                             ",
    final[population_robustness_class == "INTERMEDIATE", .N], "\n",
    "CONCERN:                                  ",
    final[population_robustness_class == "CONCERN", .N], "\n",
    "Population metrics used in PU training:   NO\n",
    "============================================================\n",
    sep = ""
)

print(qc)

log_step("14", "Population-genomic integration completed successfully")

