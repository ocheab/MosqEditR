# =============================================================================
# 16_build_publication_tables_and_figures.R
# MosqEdit-R Manuscript 1
# STEP 16 â€” Build manuscript-ready tables, figure manifest, and result extracts
# PURE R
# =============================================================================

source("R/helpers.R")

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
}
suppressPackageStartupMessages(library(data.table))

log_step("16", "Building manuscript-ready publication outputs")

OUTDIR <- "publication"
TABLE_DIR <- file.path(OUTDIR, "tables")
FIGURE_DIR <- file.path(OUTDIR, "figures")
SUPP_DIR <- file.path(OUTDIR, "supplement")

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

FILES <- list(
    master = "data_processed/master_gene_feature_matrix.csv",
    step04 = "data_processed/04_preliminary_genomewide_ranking.csv",
    step08_metrics = "data_processed/08_pu_metrics_summary.csv",
    step09 = "data_processed/09_targeted_integrated_prioritization.csv",
    step10_top20 = "data_processed/10_top20_candidate_evidence.csv",
    step11_importance = "data_processed/11_consensus_feature_importance.csv",
    step12 = "data_processed/12_bootstrap_rank_uncertainty.csv",
    step15_ablation = "data_processed/15_ablation_summary.csv"
)

required <- c(
    "master",
    "step08_metrics",
    "step09",
    "step10_top20",
    "step11_importance",
    "step12",
    "step15_ablation"
)

for (nm in required) {
    path <- FILES[[nm]]
    if (!file.exists(path)) {
        stop(
            paste0(
                "Required file is missing for publication assembly:\n",
                path
            ),
            call. = FALSE
        )
    }
}

master <- fread(FILES$master)
metrics <- fread(FILES$step08_metrics)
step09 <- fread(FILES$step09)
top20 <- fread(FILES$step10_top20)
importance <- fread(FILES$step11_importance)
step12 <- fread(FILES$step12)
ablation <- fread(FILES$step15_ablation)

pop_file <- "data_processed/14_final_population_validated_prioritization.csv"
population_available <- file.exists(pop_file)

if (population_available) {
    final_pop <- fread(pop_file)
}

# -------------------------------------------------------------------------
# Main tables
# -------------------------------------------------------------------------

table1 <- data.table(
    evidence_layer = c(
        "Genome annotation",
        "Mosquito tissue expression",
        "Drosophila orthology/phenotype",
        "PU modelling",
        "CRISPR reference editability",
        "Rank uncertainty",
        "Population genomics"
    ),
    pipeline_stage = c(
        "Step 01",
        "Step 02",
        "Step 03",
        "Steps 07â€“08",
        "Step 05",
        "Step 12",
        if (population_available) "Steps 13â€“14" else "Pending"
    ),
    analysis_scope = c(
        "13,845-gene AgamP4 universe",
        "Genome-wide partial coverage",
        "Genome-wide/partial genome-wide",
        "Leakage-safe genome-wide PU training pools",
        "205-gene targeted panel",
        "205-gene targeted panel",
        if (population_available) "Targeted validation panel" else "Not integrated"
    ),
    role = c(
        "Universe",
        "Biological discovery evidence",
        "Functional/orthology evidence",
        "Independent predictive support",
        "Reference editability validation",
        "Ranking robustness",
        "Second-stage population robustness"
    )
)

table2 <- metrics[
    order(
        -apparent_auprc_mean,
        -apparent_auroc_mean
    )
]

top20_cols <- intersect(
    c(
        "final_prepopulation_rank",
        "gene_id",
        "gene_name",
        "preliminary_rank",
        "rank_change_from_step04",
        "biological_discovery_percentile",
        "pu_consensus_percentile",
        "crispr_tractability_percentile",
        "final_integrated_score",
        "final_rank_sensitivity_range",
        "weighting_robustness",
        "pu_repeat_sd_percentile",
        "pu_stability_class"
    ),
    names(top20)
)

table3 <- top20[, ..top20_cols]

benchmark_cols <- intersect(
    c(
        "gene_id",
        "benchmark_name",
        "benchmark_class",
        "preliminary_rank",
        "final_prepopulation_rank",
        "biological_discovery_percentile",
        "pu_consensus_percentile",
        "crispr_tractability_percentile",
        "final_integrated_score",
        "final_rank_sensitivity_range"
    ),
    names(step09)
)

table4 <- step09[
    is_external_benchmark == TRUE,
    ..benchmark_cols
]

table5 <- ablation

table6 <- importance[
    1:min(15L, .N)
]

uncertainty_cols <- intersect(
    c(
        "final_prepopulation_rank",
        "gene_id",
        "bootstrap_rank_median",
        "bootstrap_rank_q025",
        "bootstrap_rank_q975",
        "probability_top10",
        "probability_top20",
        "probability_top50",
        "probability_top100",
        "bootstrap_stability_class"
    ),
    names(step12)
)

table7 <- step12[
    is_external_benchmark == FALSE &
        final_prepopulation_rank <= 20,
    ..uncertainty_cols
]

if (population_available) {
    pop_cols <- intersect(
        c(
            "final_population_validated_rank",
            "gene_id",
            "final_prepopulation_rank",
            "population_robustness_class",
            "gene_population_robustness_score",
            "best_population_site_robustness",
            "best_target_23bp_exact_match_fraction",
            "best_max_target_variant_alt_af"
        ),
        names(final_pop)
    )

    table8 <- final_pop[
        is_external_benchmark == FALSE
    ][
        1:min(50L, .N),
        ..pop_cols
    ]
}

fwrite(table1, file.path(TABLE_DIR, "Table1_data_sources_and_roles.csv"), na = "NA")
fwrite(table2, file.path(TABLE_DIR, "Table2_PU_model_performance.csv"), na = "NA")
fwrite(table3, file.path(TABLE_DIR, "Table3_top20_novel_candidates.csv"), na = "NA")
fwrite(table4, file.path(TABLE_DIR, "Table4_external_benchmark_recovery.csv"), na = "NA")
fwrite(table5, file.path(TABLE_DIR, "Table5_ablation_analysis.csv"), na = "NA")
fwrite(table6, file.path(TABLE_DIR, "Table6_top_global_predictors.csv"), na = "NA")
fwrite(table7, file.path(TABLE_DIR, "Table7_top20_rank_uncertainty.csv"), na = "NA")

if (population_available) {
    fwrite(
        table8,
        file.path(TABLE_DIR, "Table8_population_validated_top50.csv"),
        na = "NA"
    )
}

# -------------------------------------------------------------------------
# Supplementary tables
# -------------------------------------------------------------------------

supp_files <- c(
    "data_processed/04_preliminary_genomewide_ranking.csv",
    "data_processed/05_crispr_editability.csv",
    "data_processed/08_pu_gene_oof_scores.csv",
    "data_processed/09_targeted_integrated_prioritization.csv",
    "data_processed/10_integrated_score_contributions.csv",
    "data_processed/11_xgboost_global_shap_importance.csv",
    "data_processed/11_ranger_permutation_importance.csv",
    "data_processed/12_bootstrap_rank_uncertainty.csv",
    "data_processed/15_gene_level_ablation.csv"
)

for (src in supp_files[file.exists(supp_files)]) {
    file.copy(
        src,
        file.path(SUPP_DIR, basename(src)),
        overwrite = TRUE
    )
}

if (population_available) {
    file.copy(
        pop_file,
        file.path(SUPP_DIR, basename(pop_file)),
        overwrite = TRUE
    )
}

# -------------------------------------------------------------------------
# Figure manifest and copies
# -------------------------------------------------------------------------

figure_candidates <- list.files(
    "figures",
    pattern = "^(10|11|12|15).+\\.(png|pdf)$",
    full.names = TRUE
)

for (src in figure_candidates) {
    file.copy(
        src,
        file.path(FIGURE_DIR, basename(src)),
        overwrite = TRUE
    )
}

figure_manifest <- data.table(
    filename = basename(figure_candidates),
    source_path = figure_candidates,
    publication_path = file.path(
        FIGURE_DIR,
        basename(figure_candidates)
    )
)

fwrite(
    figure_manifest,
    file.path(OUTDIR, "figure_manifest.csv"),
    na = "NA"
)

# -------------------------------------------------------------------------
# Results snapshot
# -------------------------------------------------------------------------

best_model <- metrics[
    order(-apparent_auprc_mean, -apparent_auroc_mean)
][1]

top_candidate <- top20[
    order(final_prepopulation_rank)
][1]

rho <- suppressWarnings(
    cor(
        step12[
            is_external_benchmark == FALSE,
            final_prepopulation_rank
        ],
        step12[
            is_external_benchmark == FALSE,
            bootstrap_rank_median
        ],
        method = "spearman",
        use = "complete.obs"
    )
)

results_snapshot <- data.table(
    result = c(
        "Master gene universe",
        "Best PU model set",
        "Best PU algorithm",
        "Best mean apparent AUPRC",
        "Best mean apparent AUROC",
        "Top novel candidate",
        "Top candidate final rank",
        "Top candidate final integrated score",
        "Frozen-vs-bootstrap Spearman rho",
        "Population genomics integrated"
    ),
    value = c(
        as.character(nrow(master)),
        as.character(best_model$model_set),
        as.character(best_model$algorithm),
        format(best_model$apparent_auprc_mean, digits = 6),
        format(best_model$apparent_auroc_mean, digits = 6),
        as.character(top_candidate$gene_id),
        as.character(top_candidate$final_prepopulation_rank),
        format(top_candidate$final_integrated_score, digits = 6),
        format(rho, digits = 6),
        if (population_available) "YES" else "NO"
    )
)

fwrite(
    results_snapshot,
    file.path(OUTDIR, "results_snapshot.csv"),
    na = "NA"
)

# -------------------------------------------------------------------------
# Publication manifest
# -------------------------------------------------------------------------

publication_files <- list.files(
    OUTDIR,
    recursive = TRUE,
    full.names = TRUE
)

publication_manifest <- data.table(
    file = publication_files,
    size_bytes = file.info(publication_files)$size
)

fwrite(
    publication_manifest,
    file.path(OUTDIR, "publication_manifest.csv"),
    na = "NA"
)

capture.output(
    sessionInfo(),
    file = "logs/16_publication_sessionInfo.txt"
)

cat(
    "\n============================================================\n",
    "MOSQEDIT-R STEP 16 COMPLETED SUCCESSFULLY\n",
    "Manuscript-ready publication package assembled\n",
    "============================================================\n",
    "Main tables written:                      ",
    7L + as.integer(population_available),
    "\n",
    "Figures copied:                           ",
    length(figure_candidates),
    "\n",
    "Population genomics integrated:           ",
    if (population_available) "YES" else "NO",
    "\n",
    "Publication directory:                    ",
    OUTDIR,
    "\n",
    "============================================================\n",
    sep = ""
)

print(results_snapshot)

log_step("16", "Publication outputs assembled successfully")

