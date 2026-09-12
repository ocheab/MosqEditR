# =============================================================================
# 09_integrate_pu_crispr_final_prioritization.R
#
# MosqEdit-R Manuscript 1
#
# STEP 09
# Robust PU consensus + CRISPR-aware final pre-population prioritization
#
# PURE R
#
# PURPOSE
# -------
# Integrate three conceptually distinct evidence layers on the Step 05 targeted
# panel:
#
#   1. Step 04 biological discovery evidence
#   2. Step 08 leakage-safe PU support
#   3. Step 05 exon-aware CRISPR tractability
#
# IMPORTANT
# ---------
# This is a PRE-POPULATION-GENOMICS prioritization.
# Ag1000G conservation remains pending and is NOT imputed.
#
# External benchmark genes remain evaluation-only objects. Their recovery is
# reported, but they are NEVER used to select models or tune integration weights.
#
# PU CONSENSUS
# ------------
# The top three CV-performing model/set combinations from Step 08 are selected
# by apparent AUPRC. Their scores are transformed to within-model percentiles
# before aggregation, avoiding direct averaging of incomparable score scales.
#
# PRIMARY INTEGRATION
# -------------------
# Biological discovery evidence: 50%
# PU consensus support:           30%
# CRISPR tractability:            20%
#
# Completeness adjustment:
#   adjusted_score = raw_score * [0.80 + 0.20 * domain_coverage]
#
# SENSITIVITY SCENARIOS
# ---------------------
# biology_led:    0.60 / 0.25 / 0.15
# primary:        0.50 / 0.30 / 0.20
# validation_led: 0.40 / 0.35 / 0.25
#
# Weights are prespecified and NOT optimized against benchmark recovery.
#
# INPUTS
# ------
# data_processed/04_preliminary_genomewide_ranking.csv
# data_processed/07_targeted_crispr_validation_panel.csv
# data_processed/08_pu_gene_oof_scores.csv
# data_processed/08_external_benchmark_scores.csv
# data_processed/08_pu_metrics_summary.csv
# data_processed/08_pu_cv_predictions.csv
#
# OUTPUTS
# -------
# data_processed/09_selected_pu_models.csv
# data_processed/09_targeted_integrated_prioritization.csv
# data_processed/09_novel_candidate_ranking.csv
# data_processed/09_benchmark_recovery.csv
# data_processed/09_top50_prepopulation_shortlist.csv
# data_processed/09_top100_prepopulation_shortlist.csv
# data_processed/09_pu_repeat_stability.csv
# data_processed/09_prioritization_qc.csv
# data_processed/09_prioritization_provenance.csv
# logs/09_prioritization_sessionInfo.txt
# logs/09_checksums.tsv
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Helpers and packages
# -----------------------------------------------------------------------------

source("R/helpers.R")


if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
        "Package 'data.table' is required. Install with install.packages('data.table').",
        call. = FALSE
    )
}


suppressPackageStartupMessages({
    library(data.table)
})


log_step(
    "09",
    "Starting robust PU-consensus and CRISPR-aware final prioritization"
)


# -----------------------------------------------------------------------------
# 2. Configuration
# -----------------------------------------------------------------------------

STEP04_FILE <-
    "data_processed/04_preliminary_genomewide_ranking.csv"


TARGETED_FILE <-
    "data_processed/07_targeted_crispr_validation_panel.csv"


OOF_FILE <-
    "data_processed/08_pu_gene_oof_scores.csv"


BENCHMARK_SCORE_FILE <-
    "data_processed/08_external_benchmark_scores.csv"


METRICS_FILE <-
    "data_processed/08_pu_metrics_summary.csv"


CV_PREDICTIONS_FILE <-
    "data_processed/08_pu_cv_predictions.csv"


SELECTED_MODELS_FILE <-
    "data_processed/09_selected_pu_models.csv"


INTEGRATED_FILE <-
    "data_processed/09_targeted_integrated_prioritization.csv"


NOVEL_FILE <-
    "data_processed/09_novel_candidate_ranking.csv"


BENCHMARK_RECOVERY_FILE <-
    "data_processed/09_benchmark_recovery.csv"


TOP50_FILE <-
    "data_processed/09_top50_prepopulation_shortlist.csv"


TOP100_FILE <-
    "data_processed/09_top100_prepopulation_shortlist.csv"


PU_STABILITY_FILE <-
    "data_processed/09_pu_repeat_stability.csv"


QC_FILE <-
    "data_processed/09_prioritization_qc.csv"


PROVENANCE_FILE <-
    "data_processed/09_prioritization_provenance.csv"


SESSION_INFO_FILE <-
    "logs/09_prioritization_sessionInfo.txt"


CHECKSUM_FILE <-
    "logs/09_checksums.tsv"


N_PU_MODELS_FOR_CONSENSUS <- 3L

TOP50_N <- 50L
TOP100_N <- 100L

COMPLETENESS_FLOOR <- 0.80


# Primary weights.
PRIMARY_WEIGHTS <- c(
    biology = 0.50,
    pu = 0.30,
    crispr = 0.20
)


# Sensitivity scenarios.
SCENARIO_WEIGHTS <- list(
    biology_led = c(
        biology = 0.60,
        pu = 0.25,
        crispr = 0.15
    ),

    primary = PRIMARY_WEIGHTS,

    validation_led = c(
        biology = 0.40,
        pu = 0.35,
        crispr = 0.25
    )
)


dir.create(
    "data_processed",
    recursive = TRUE,
    showWarnings = FALSE
)


dir.create(
    "logs",
    recursive = TRUE,
    showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# 3. Utility functions
# -----------------------------------------------------------------------------

assert_file <- function(path) {

    if (!file.exists(path)) {

        stop(
            paste0(
                "Required input file is missing:\n",
                path
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


assert_unique_gene_ids <- function(
        dt,
        label
) {

    if (!"gene_id" %in% names(dt)) {

        stop(
            paste0(
                label,
                " lacks gene_id."
            ),
            call. = FALSE
        )
    }


    if (
        uniqueN(
            dt$gene_id
        ) !=
            nrow(dt)
    ) {

        stop(
            paste0(
                label,
                " contains duplicate gene_id values."
            ),
            call. = FALSE
        )
    }


    invisible(TRUE)
}


rank01 <- function(
        x,
        higher_is_better = TRUE
) {

    x <- suppressWarnings(
        as.numeric(
            x
        )
    )


    out <- rep(
        NA_real_,
        length(
            x
        )
    )


    ok <- is.finite(
        x
    )


    n_ok <- sum(
        ok
    )


    if (n_ok == 0L) {

        return(
            out
        )
    }


    if (n_ok == 1L) {

        out[ok] <- 0.5

        return(
            out
        )
    }


    z <- x[ok]


    if (!higher_is_better) {

        z <- -z
    }


    r <- rank(
        z,
        ties.method = "average"
    )


    out[ok] <- (
        r -
        1
    ) /
        (
            n_ok -
            1
        )


    out
}


weighted_available_score <- function(
        biology,
        pu,
        crispr,
        weights,
        completeness_floor
) {

    components <- cbind(
        biology,
        pu,
        crispr
    )


    colnames(
        components
    ) <- c(
        "biology",
        "pu",
        "crispr"
    )


    w <- weights[
        colnames(
            components
        )
    ]


    observed <- is.finite(
        components
    )


    numerator <- rowSums(
        sweep(
            ifelse(
                observed,
                components,
                0
            ),
            2,
            w,
            `*`
        )
    )


    observed_weight <- rowSums(
        sweep(
            observed * 1,
            2,
            w,
            `*`
        )
    )


    raw <- numerator /
        observed_weight


    raw[
        observed_weight <=
            0
    ] <- NA_real_


    domain_coverage <- rowMeans(
        observed
    )


    adjusted <- raw *
        (
            completeness_floor +
            (
                1 -
                completeness_floor
            ) *
            domain_coverage
        )


    list(
        raw = raw,
        adjusted = adjusted,
        domain_coverage = domain_coverage,
        observed_weight = observed_weight
    )
}


deterministic_rank <- function(
        score,
        gene_id
) {

    z <- data.table(
        gene_id =
            gene_id,

        score =
            suppressWarnings(
                as.numeric(
                    score
                )
            ),

        row_index =
            seq_along(
                gene_id
            )
    )


    z[
        ,
        score_for_sort :=
            fifelse(
                is.finite(
                    score
                ),
                score,
                -Inf
            )
    ]


    data.table::setorder(
        z,
        -score_for_sort,
        gene_id
    )


    z[
        ,
        rank :=
            seq_len(
                .N
            )
    ]


    out <- integer(
        nrow(
            z
        )
    )


    out[
        z$row_index
    ] <- z$rank


    out[
        !is.finite(
            score
        )
    ] <- NA_integer_


    out
}


# -----------------------------------------------------------------------------
# 4. Load inputs
# -----------------------------------------------------------------------------

for (
    path in c(
        STEP04_FILE,
        TARGETED_FILE,
        OOF_FILE,
        BENCHMARK_SCORE_FILE,
        METRICS_FILE,
        CV_PREDICTIONS_FILE
    )
) {

    assert_file(
        path
    )
}


step04 <- fread(
    STEP04_FILE
)


targeted <- fread(
    TARGETED_FILE
)


oof <- fread(
    OOF_FILE
)


benchmark_scores <- fread(
    BENCHMARK_SCORE_FILE
)


metrics <- fread(
    METRICS_FILE
)


cv_predictions <- fread(
    CV_PREDICTIONS_FILE
)


assert_unique_gene_ids(
    step04,
    "Step 04 ranking"
)


assert_unique_gene_ids(
    targeted,
    "Step 07 targeted CRISPR panel"
)


# -----------------------------------------------------------------------------
# 5. Select top PU models by apparent AUPRC
# -----------------------------------------------------------------------------
#
# This selection uses CV metrics ONLY.
# Benchmark recovery is not inspected or used.

required_metric_columns <- c(
    "model_set",
    "algorithm",
    "apparent_auprc_mean",
    "apparent_auroc_mean"
)


missing_metric_columns <- setdiff(
    required_metric_columns,
    names(
        metrics
    )
)


if (
    length(
        missing_metric_columns
    ) >
        0L
) {

    stop(
        paste0(
            "Metrics table is missing required column(s): ",
            paste(
                missing_metric_columns,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}


selected_models <- copy(
    metrics
)


data.table::setorder(
    selected_models,
    -apparent_auprc_mean,
    -apparent_auroc_mean,
    model_set,
    algorithm
)


selected_models <- selected_models[
    seq_len(
        min(
            N_PU_MODELS_FOR_CONSENSUS,
            .N
        )
    )
]


selected_models[
    ,
    pu_consensus_model_rank :=
        seq_len(
            .N
        )
]


selected_models[
    ,
    model_key :=
        paste(
            model_set,
            algorithm,
            sep = "::"
        )
]


if (
    nrow(
        selected_models
    ) <
        2L
) {

    stop(
        "Fewer than two PU models were available for consensus.",
        call. = FALSE
    )
}


# -----------------------------------------------------------------------------
# 6. Convert selected OOF model scores to percentiles
# -----------------------------------------------------------------------------

selected_keys <- selected_models$model_key


oof[
    ,
    model_key :=
        paste(
            model_set,
            algorithm,
            sep = "::"
        )
]


selected_oof <- oof[
    model_key %in%
        selected_keys
]


if (
    nrow(
        selected_oof
    ) ==
        0L
) {

    stop(
        "No OOF scores matched the selected PU models.",
        call. = FALSE
    )
}


selected_oof[
    ,
    pu_score_percentile :=
        rank01(
            mean_oof_pu_score,
            higher_is_better = TRUE
        ),
    by = model_key
]


# Wide percentile table for non-holdout genes.
pu_wide <- dcast(
    selected_oof,
    gene_id ~ model_key,
    value.var =
        "pu_score_percentile"
)


selected_percentile_columns <- setdiff(
    names(
        pu_wide
    ),
    "gene_id"
)


pu_wide[
    ,
    pu_consensus_percentile :=
        apply(
            .SD,
            1,
            function(z) {

                z <- as.numeric(
                    z
                )


                z <- z[
                    is.finite(
                        z
                    )
                ]


                if (
                    length(
                        z
                    ) ==
                        0L
                ) {

                    return(
                        NA_real_
                    )
                }


                median(
                    z
                )
            }
        ),
    .SDcols =
        selected_percentile_columns
]


pu_wide[
    ,
    pu_models_available :=
        rowSums(
            !is.na(
                .SD
            )
        ),
    .SDcols =
        selected_percentile_columns
]


# -----------------------------------------------------------------------------
# 7. Add external benchmark PU percentiles
# -----------------------------------------------------------------------------

benchmark_scores[
    ,
    model_key :=
        paste(
            model_set,
            algorithm,
            sep = "::"
        )
]


selected_benchmark <- benchmark_scores[
    model_key %in%
        selected_keys
]


benchmark_pu_wide <- dcast(
    selected_benchmark,
    gene_id ~ model_key,
    value.var =
        "percentile_vs_crossfitted_nonholdout_genes"
)


benchmark_percentile_columns <- setdiff(
    names(
        benchmark_pu_wide
    ),
    "gene_id"
)


benchmark_pu_wide[
    ,
    pu_consensus_percentile :=
        apply(
            .SD,
            1,
            function(z) {

                z <- as.numeric(
                    z
                )


                z <- z[
                    is.finite(
                        z
                    )
                ]


                if (
                    length(
                        z
                    ) ==
                        0L
                ) {

                    return(
                        NA_real_
                    )
                }


                median(
                    z
                )
            }
        ),
    .SDcols =
        benchmark_percentile_columns
]


benchmark_pu_wide[
    ,
    pu_models_available :=
        rowSums(
            !is.na(
                .SD
            )
        ),
    .SDcols =
        benchmark_percentile_columns
]


# Combine non-holdout OOF percentiles and benchmark external percentiles.
all_pu <- rbindlist(
    list(
        pu_wide,
        benchmark_pu_wide
    ),
    fill = TRUE,
    use.names = TRUE
)


# Benchmark rows should override any accidental duplicate.
all_pu[
    ,
    is_benchmark_pu_row :=
        gene_id %in%
        benchmark_pu_wide$gene_id
]


data.table::setorder(
    all_pu,
    gene_id,
    -is_benchmark_pu_row
)


all_pu <- unique(
    all_pu,
    by = "gene_id"
)


all_pu[
    ,
    is_benchmark_pu_row :=
        NULL
]


# -----------------------------------------------------------------------------
# 8. PU repeat stability for non-holdout genes
# -----------------------------------------------------------------------------

cv_predictions[
    ,
    model_key :=
        paste(
            model_set,
            algorithm,
            sep = "::"
        )
]


selected_cv <- cv_predictions[
    model_key %in%
        selected_keys
]


# Each gene has one OOF prediction per repeat for a given model.
selected_cv[
    ,
    repeat_model_percentile :=
        rank01(
            pu_score,
            higher_is_better = TRUE
        ),
    by = .(
        model_key,
        repeat_id
    )
]


repeat_consensus <- selected_cv[
    ,
    .(
        pu_repeat_consensus_percentile =
            median(
                repeat_model_percentile,
                na.rm = TRUE
            ),

        pu_models_available_in_repeat =
            sum(
                is.finite(
                    repeat_model_percentile
                )
            )
    ),
    by = .(
        gene_id,
        repeat_id
    )
]


pu_stability <- repeat_consensus[
    ,
    .(
        pu_repeat_mean_percentile =
            mean(
                pu_repeat_consensus_percentile,
                na.rm = TRUE
            ),

        pu_repeat_sd_percentile =
            stats::sd(
                pu_repeat_consensus_percentile,
                na.rm = TRUE
            ),

        pu_repeat_min_percentile =
            min(
                pu_repeat_consensus_percentile,
                na.rm = TRUE
            ),

        pu_repeat_max_percentile =
            max(
                pu_repeat_consensus_percentile,
                na.rm = TRUE
            ),

        n_pu_repeats =
            .N
    ),
    by = gene_id
]


pu_stability[
    !is.finite(
        pu_repeat_sd_percentile
    ),
    pu_repeat_sd_percentile :=
        NA_real_
]


# -----------------------------------------------------------------------------
# 9. Assemble targeted panel
# -----------------------------------------------------------------------------

integrated <- copy(
    targeted
)


integrated <- merge(
    integrated,
    all_pu[
        ,
        .(
            gene_id,
            pu_consensus_percentile,
            pu_models_available
        )
    ],
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


integrated <- merge(
    integrated,
    pu_stability,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


# -----------------------------------------------------------------------------
# 10. Biological discovery percentile
# -----------------------------------------------------------------------------
#
# Use Step 04 primary score percentile across ALL genes with Step 04 primary
# score, not only the targeted panel.

if (
    !"score_primary" %in%
    names(
        step04
    )
) {

    stop(
        "Step 04 table lacks score_primary.",
        call. = FALSE
    )
}


step04[
    ,
    biological_discovery_percentile :=
        rank01(
            score_primary,
            higher_is_better = TRUE
        )
]


biology_map <- step04[
    ,
    .(
        gene_id,
        biological_discovery_percentile,
        preliminary_rank
    )
]


if (
    "preliminary_rank" %in%
    names(
        integrated
    )
) {

    biology_map[
        ,
        preliminary_rank :=
            NULL
    ]
}


integrated <- merge(
    integrated,
    biology_map,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


# -----------------------------------------------------------------------------
# 11. CRISPR percentile within targeted validation panel
# -----------------------------------------------------------------------------

if (
    !"crispr_exon_aware_reference_editability_score" %in%
    names(
        integrated
    )
) {

    stop(
        "Targeted panel lacks crispr_exon_aware_reference_editability_score.",
        call. = FALSE
    )
}


integrated[
    ,
    crispr_tractability_percentile :=
        rank01(
            crispr_exon_aware_reference_editability_score,
            higher_is_better = TRUE
        )
]


# -----------------------------------------------------------------------------
# 12. Integrated primary score
# -----------------------------------------------------------------------------

primary <- weighted_available_score(
    biology =
        integrated$
            biological_discovery_percentile,

    pu =
        integrated$
            pu_consensus_percentile,

    crispr =
        integrated$
            crispr_tractability_percentile,

    weights =
        PRIMARY_WEIGHTS,

    completeness_floor =
        COMPLETENESS_FLOOR
)


integrated[
    ,
    final_integrated_raw :=
        primary$raw
]


integrated[
    ,
    final_evidence_domain_coverage :=
        primary$domain_coverage
]


integrated[
    ,
    final_integrated_score :=
        primary$adjusted
]


integrated[
    ,
    final_integrated_rank :=
        deterministic_rank(
            final_integrated_score,
            gene_id
        )
]


# -----------------------------------------------------------------------------
# 13. Sensitivity scenarios
# -----------------------------------------------------------------------------

scenario_rank_columns <- character()


for (
    scenario_name in names(
        SCENARIO_WEIGHTS
    )
) {

    w <- SCENARIO_WEIGHTS[[scenario_name]]


    result <- weighted_available_score(
        biology =
            integrated$
                biological_discovery_percentile,

        pu =
            integrated$
                pu_consensus_percentile,

        crispr =
            integrated$
                crispr_tractability_percentile,

        weights =
            w,

        completeness_floor =
            COMPLETENESS_FLOOR
    )


    score_col <- paste0(
        "score_",
        scenario_name
    )


    rank_col <- paste0(
        "rank_",
        scenario_name
    )


    integrated[
        ,
        (score_col) :=
            result$adjusted
    ]


    integrated[
        ,
        (rank_col) :=
            deterministic_rank(
                get(
                    score_col
                ),
                gene_id
            )
    ]


    scenario_rank_columns <- c(
        scenario_rank_columns,
        rank_col
    )
}


integrated[
    ,
    final_consensus_rank_median :=
        apply(
            .SD,
            1,
            function(z) {

                z <- as.numeric(
                    z
                )


                z <- z[
                    is.finite(
                        z
                    )
                ]


                if (
                    length(
                        z
                    ) ==
                        0L
                ) {

                    return(
                        NA_real_
                    )
                }


                median(
                    z
                )
            }
        ),
    .SDcols =
        scenario_rank_columns
]


integrated[
    ,
    final_rank_sensitivity_range :=
        apply(
            .SD,
            1,
            function(z) {

                z <- as.numeric(
                    z
                )


                z <- z[
                    is.finite(
                        z
                    )
                ]


                if (
                    length(
                        z
                    ) <=
                        1L
                ) {

                    return(
                        NA_real_
                    )
                }


                max(
                    z
                ) -
                min(
                    z
                )
            }
        ),
    .SDcols =
        scenario_rank_columns
]


# Primary final ordering:
# consensus median rank first, then primary adjusted score.
integrated[
    ,
    final_consensus_missing :=
        is.na(
            final_consensus_rank_median
        )
]


data.table::setorder(
    integrated,
    final_consensus_missing,
    final_consensus_rank_median,
    -final_integrated_score,
    gene_id
)


integrated[
    ,
    final_prepopulation_rank :=
        data.table::fifelse(
            !final_consensus_missing,
            seq_len(
                .N
            ),
            NA_integer_
        )
]


integrated[
    ,
    final_consensus_missing :=
        NULL
]


# Reassign rank sequentially ONLY among scoreable genes, by key.
rankable <- integrated[
    !is.na(
        final_consensus_rank_median
    ),
    .(
        gene_id,
        final_consensus_rank_median,
        final_integrated_score
    )
]


data.table::setorder(
    rankable,
    final_consensus_rank_median,
    -final_integrated_score,
    gene_id
)


rankable[
    ,
    final_prepopulation_rank :=
        seq_len(
            .N
        )
]


integrated[
    ,
    final_prepopulation_rank :=
        NA_integer_
]


integrated[
    rankable,
    on = "gene_id",
    final_prepopulation_rank :=
        i.final_prepopulation_rank
]


# -----------------------------------------------------------------------------
# 14. Benchmark status
# -----------------------------------------------------------------------------

if (
    "is_external_landmark_validation" %in%
    names(
        integrated
    )
) {

    integrated[
        ,
        is_external_benchmark :=
            !is.na(
                is_external_landmark_validation
            ) &
            as.logical(
                is_external_landmark_validation
            )
    ]

} else if (
    "benchmark_name" %in%
    names(
        integrated
    )
) {

    integrated[
        ,
        is_external_benchmark :=
            !is.na(
                benchmark_name
            )
    ]

} else {

    integrated[
        ,
        is_external_benchmark :=
            FALSE
    ]
}


# -----------------------------------------------------------------------------
# 15. Novel candidate table and top-N flags
# -----------------------------------------------------------------------------

novel <- integrated[
    is_external_benchmark ==
        FALSE
]


data.table::setorder(
    novel,
    final_prepopulation_rank,
    gene_id,
    na.last = TRUE
)


novel[
    ,
    novel_candidate_rank :=
        seq_len(
            .N
        )
]


integrated <- merge(
    integrated,
    novel[
        ,
        .(
            gene_id,
            novel_candidate_rank
        )
    ],
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


integrated[
    ,
    top50_novel_prepopulation :=
        !is.na(
            novel_candidate_rank
        ) &
        novel_candidate_rank <=
            TOP50_N
]


integrated[
    ,
    top100_novel_prepopulation :=
        !is.na(
            novel_candidate_rank
        ) &
        novel_candidate_rank <=
            TOP100_N
]


# Restore final overall order.
integrated[
    ,
    final_rank_missing :=
        is.na(
            final_prepopulation_rank
        )
]


data.table::setorder(
    integrated,
    final_rank_missing,
    final_prepopulation_rank,
    gene_id
)


integrated[
    ,
    final_rank_missing :=
        NULL
]


novel <- integrated[
    is_external_benchmark ==
        FALSE
]


benchmark_recovery <- integrated[
    is_external_benchmark ==
        TRUE
]


top50 <- integrated[
    top50_novel_prepopulation ==
        TRUE
]


top100 <- integrated[
    top100_novel_prepopulation ==
        TRUE
]


# -----------------------------------------------------------------------------
# 16. Integrity checks
# -----------------------------------------------------------------------------

if (
    uniqueN(
        integrated$gene_id
    ) !=
        nrow(
            integrated
        )
) {

    stop(
        "Integrated targeted table contains duplicate gene_id values.",
        call. = FALSE
    )
}


if (
    nrow(
        integrated
    ) !=
        nrow(
            targeted
        )
) {

    stop(
        paste0(
            "Integrated panel row count changed from ",
            nrow(
                targeted
            ),
            " to ",
            nrow(
                integrated
            ),
            "."
        ),
        call. = FALSE
    )
}


ranked <- integrated[
    !is.na(
        final_prepopulation_rank
    )
]


if (
    uniqueN(
        ranked$
            final_prepopulation_rank
    ) !=
        nrow(
            ranked
        )
) {

    stop(
        "Final pre-population ranks are not unique.",
        call. = FALSE
    )
}


if (
    nrow(
        top50
    ) !=
        min(
            TOP50_N,
            nrow(
                novel
            )
        )
) {

    stop(
        "Top-50 novel shortlist size is incorrect.",
        call. = FALSE
    )
}


if (
    nrow(
        top100
    ) !=
        min(
            TOP100_N,
            nrow(
                novel
            )
        )
) {

    stop(
        "Top-100 novel shortlist size is incorrect.",
        call. = FALSE
    )
}


# -----------------------------------------------------------------------------
# 17. QC
# -----------------------------------------------------------------------------

qc <- data.table(
    metric = c(
        "Targeted panel genes",
        "External benchmark genes",
        "Novel candidate genes",
        "PU models selected for consensus",
        "Genes with biological discovery percentile",
        "Genes with PU consensus percentile",
        "Genes with CRISPR tractability percentile",
        "Genes with all three integration domains",
        "Genes with final pre-population rank",
        "Top-50 novel candidates",
        "Top-100 novel candidates",
        "Population-genomics domains integrated",
        "Benchmark genes used to choose PU models",
        "Benchmark genes used to tune integration weights"
    ),

    value = c(
        nrow(
            integrated
        ),

        sum(
            integrated$
                is_external_benchmark,
            na.rm = TRUE
        ),

        sum(
            !integrated$
                is_external_benchmark,
            na.rm = TRUE
        ),

        nrow(
            selected_models
        ),

        sum(
            is.finite(
                integrated$
                    biological_discovery_percentile
            )
        ),

        sum(
            is.finite(
                integrated$
                    pu_consensus_percentile
            )
        ),

        sum(
            is.finite(
                integrated$
                    crispr_tractability_percentile
            )
        ),

        sum(
            integrated$
                final_evidence_domain_coverage ==
                1,
            na.rm = TRUE
        ),

        nrow(
            ranked
        ),

        nrow(
            top50
        ),

        nrow(
            top100
        ),

        0,
        0,
        0
    )
)


# -----------------------------------------------------------------------------
# 18. Provenance
# -----------------------------------------------------------------------------

selected_model_text <- paste(
    selected_models$model_key,
    collapse = "; "
)


provenance <- data.table(
    component = c(
        "PU model selection",
        "Selected PU models",
        "PU score harmonization",
        "Benchmark PU scoring",
        "Biological discovery domain",
        "CRISPR domain",
        "Primary integration weights",
        "Completeness adjustment",
        "Sensitivity scenarios",
        "Benchmark independence",
        "Population genomics",
        "Interpretation"
    ),

    specification = c(
        paste0(
            "Top ",
            nrow(
                selected_models
            ),
            " model/set combinations selected exclusively from Step 08 ",
            "cross-validated apparent AUPRC; benchmark recovery was not used."
        ),

        selected_model_text,

        paste0(
            "Each selected model's OOF score is transformed to a within-model ",
            "percentile; gene-level PU support is the median percentile across ",
            "selected models."
        ),

        paste0(
            "External benchmark percentiles come only from Step 08 models in ",
            "which benchmark genes were completely held out from fitting."
        ),

        paste0(
            "Step 04 score_primary transformed to a percentile across all genes ",
            "with a primary Step 04 score."
        ),

        paste0(
            "Step 05 exon-aware reference editability transformed to a percentile ",
            "within the targeted panel."
        ),

        paste0(
            "biology=",
            PRIMARY_WEIGHTS[["biology"]],
            ", PU=",
            PRIMARY_WEIGHTS[["pu"]],
            ", CRISPR=",
            PRIMARY_WEIGHTS[["crispr"]],
            "."
        ),

        paste0(
            "Available-domain weighted score multiplied by [",
            COMPLETENESS_FLOOR,
            " + ",
            1 -
                COMPLETENESS_FLOOR,
            " * domain coverage]."
        ),

        paste0(
            "biology-led 0.60/0.25/0.15; primary 0.50/0.30/0.20; ",
            "validation-led 0.40/0.35/0.25. Final order uses median scenario rank ",
            "then primary adjusted score."
        ),

        paste0(
            "Benchmark recovery is reported after ranking and is not used to ",
            "choose models, weights, thresholds, or sensitivity scenarios."
        ),

        paste0(
            "Ag1000G population conservation remains pending and contributes ",
            "no value to Step 09."
        ),

        paste0(
            "Step 09 ranks reference-biologically-supported, PU-supported, ",
            "CRISPR-tractable targets prior to population-genomic validation."
        )
    )
)


# -----------------------------------------------------------------------------
# 19. Write outputs
# -----------------------------------------------------------------------------

fwrite(
    selected_models,
    SELECTED_MODELS_FILE,
    na = "NA"
)


fwrite(
    integrated,
    INTEGRATED_FILE,
    na = "NA"
)


fwrite(
    novel,
    NOVEL_FILE,
    na = "NA"
)


fwrite(
    benchmark_recovery,
    BENCHMARK_RECOVERY_FILE,
    na = "NA"
)


fwrite(
    top50,
    TOP50_FILE,
    na = "NA"
)


fwrite(
    top100,
    TOP100_FILE,
    na = "NA"
)


fwrite(
    pu_stability,
    PU_STABILITY_FILE,
    na = "NA"
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
    STEP04_FILE,
    TARGETED_FILE,
    OOF_FILE,
    BENCHMARK_SCORE_FILE,
    METRICS_FILE,
    CV_PREDICTIONS_FILE,
    SELECTED_MODELS_FILE,
    INTEGRATED_FILE,
    NOVEL_FILE,
    BENCHMARK_RECOVERY_FILE,
    TOP50_FILE,
    TOP100_FILE,
    PU_STABILITY_FILE,
    QC_FILE,
    PROVENANCE_FILE,
    SESSION_INFO_FILE
)


checksum_files <- checksum_files[
    file.exists(
        checksum_files
    )
]


write_checksum(
    checksum_files,
    CHECKSUM_FILE
)


# -----------------------------------------------------------------------------
# 20. Console summary
# -----------------------------------------------------------------------------

cat(
    "\n",
    "============================================================\n",
    "MOSQEDIT-R STEP 09 COMPLETED SUCCESSFULLY\n",
    "Robust PU consensus + CRISPR-aware pre-population ranking\n",
    "============================================================\n",
    "Targeted panel genes:                     ",
    nrow(
        integrated
    ),
    "\n",
    "External benchmarks:                      ",
    nrow(
        benchmark_recovery
    ),
    "\n",
    "Novel candidates:                         ",
    nrow(
        novel
    ),
    "\n",
    "PU models in consensus:                   ",
    nrow(
        selected_models
    ),
    "\n",
    "Genes with all three evidence domains:    ",
    sum(
        integrated$
            final_evidence_domain_coverage ==
            1,
        na.rm = TRUE
    ),
    "\n",
    "Top-50 novel shortlist:                   ",
    nrow(
        top50
    ),
    "\n",
    "Top-100 novel shortlist:                  ",
    nrow(
        top100
    ),
    "\n",
    "Population genomics integrated:           NO\n",
    "Benchmark-informed weight tuning:         NO\n",
    "Python used:                              NO\n",
    "============================================================\n",
    sep = ""
)


cat(
    "\nPU models selected by cross-validated apparent AUPRC:\n"
)


print(
    selected_models[
        ,
        .(
            pu_consensus_model_rank,
            model_set,
            algorithm,
            apparent_auprc_mean,
            apparent_auprc_sd,
            apparent_auroc_mean
        )
    ]
)


cat(
    "\nTop 20 NOVEL pre-population candidates:\n"
)


top_display_columns <- c(
    "novel_candidate_rank",
    "final_prepopulation_rank",
    "gene_id",
    "gene_name",
    "preliminary_rank",
    "biological_discovery_percentile",
    "pu_consensus_percentile",
    "pu_repeat_sd_percentile",
    "crispr_tractability_percentile",
    "final_integrated_score",
    "final_consensus_rank_median",
    "final_rank_sensitivity_range"
)


top_display_columns <- top_display_columns[
    top_display_columns %in%
        names(
            novel
        )
]


print(
    novel[
        1:min(
            20L,
            .N
        ),
        ..top_display_columns
    ]
)


cat(
    "\nExternal benchmark recovery after integration:\n"
)


benchmark_display_columns <- c(
    "gene_id",
    "benchmark_name",
    "benchmark_class",
    "preliminary_rank",
    "biological_discovery_percentile",
    "pu_consensus_percentile",
    "crispr_tractability_percentile",
    "final_integrated_score",
    "final_prepopulation_rank",
    "final_rank_sensitivity_range"
)


benchmark_display_columns <- benchmark_display_columns[
    benchmark_display_columns %in%
        names(
            benchmark_recovery
        )
]


print(
    benchmark_recovery[
        ,
        ..benchmark_display_columns
    ]
)


cat(
    "\nQC summary:\n"
)


print(
    qc
)


log_step(
    "09",
    "Robust PU consensus and CRISPR-aware prioritization completed successfully"
)

