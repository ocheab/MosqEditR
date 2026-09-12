# =============================================================================
# 12_bootstrap_rank_uncertainty.R
#
# MosqEdit-R Manuscript 1
#
# STEP 12
# Resampling uncertainty and rank-stability analysis
#
# PURE R
#
# PURPOSE
# -------
# Quantify uncertainty in the Step 09 final pre-population ranking arising from
# the PU-learning component, while preserving the already frozen biological and
# CRISPR evidence layers.
#
# For each bootstrap/resampling replicate:
#   1. sample the selected PU model keys with replacement;
#   2. for each sampled model, sample one Step 08 CV repeat;
#   3. use within-model/repeat PU score percentiles for non-benchmark genes;
#   4. use fully external Step 08 benchmark percentiles for benchmark genes;
#   5. aggregate sampled PU percentiles by the median;
#   6. recompute the Step 09 primary integrated score;
#   7. rerank the complete 205-gene targeted panel.
#
# OUTPUTS include:
#   - median bootstrap rank
#   - 2.5% and 97.5% empirical rank interval
#   - rank SD and interval width
#   - P(top 10), P(top 20), P(top 50), P(top 100)
#   - score interval
#   - robust rank envelope combining bootstrap rank uncertainty with the three
#     deterministic Step 09 weighting scenarios
#
# IMPORTANT INTERPRETATION
# ------------------------
# The empirical intervals quantify uncertainty from PU model/repeat resampling.
# They are NOT full biological confidence intervals because Step 04 and Step 05
# evidence values are held fixed.
#
# LANDMARK VALIDATION
# -------------------
# All six external landmarks were already excluded from every Step 08/11 model
# fit. A leave-one-landmark-out analysis would therefore be weaker/redundant and
# is not performed. Benchmark stability is assessed from their fully external
# PU predictions within the same resampling framework.
#
# POPULATION GENOMICS
# -------------------
# Ag1000G remains pending and is not used here.
#
# INPUTS
# ------
# data_processed/09_targeted_integrated_prioritization.csv
# data_processed/09_selected_pu_models.csv
# data_processed/08_pu_cv_predictions.csv
# data_processed/08_external_benchmark_scores.csv
#
# OUTPUTS
# -------
# data_processed/12_bootstrap_rank_uncertainty.csv
# data_processed/12_top20_rank_uncertainty.csv
# data_processed/12_benchmark_rank_uncertainty.csv
# data_processed/12_bootstrap_qc.csv
# data_processed/12_bootstrap_provenance.csv
# figures/12A_top50_rank_intervals.png/.pdf
# figures/12B_topk_inclusion_probabilities.png/.pdf
# figures/12C_final_vs_bootstrap_rank.png/.pdf
# figures/12D_benchmark_rank_intervals.png/.pdf
# logs/12_bootstrap_sessionInfo.txt
# logs/12_checksums.tsv
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
    "12",
    "Starting PU-resampling rank-uncertainty analysis"
)

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------

INTEGRATED_FILE <- "data_processed/09_targeted_integrated_prioritization.csv"
SELECTED_MODELS_FILE <- "data_processed/09_selected_pu_models.csv"
CV_PREDICTIONS_FILE <- "data_processed/08_pu_cv_predictions.csv"
BENCHMARK_SCORES_FILE <- "data_processed/08_external_benchmark_scores.csv"

UNCERTAINTY_FILE <- "data_processed/12_bootstrap_rank_uncertainty.csv"
TOP20_FILE <- "data_processed/12_top20_rank_uncertainty.csv"
BENCHMARK_FILE <- "data_processed/12_benchmark_rank_uncertainty.csv"
QC_FILE <- "data_processed/12_bootstrap_qc.csv"
PROVENANCE_FILE <- "data_processed/12_bootstrap_provenance.csv"
SESSION_INFO_FILE <- "logs/12_bootstrap_sessionInfo.txt"
CHECKSUM_FILE <- "logs/12_checksums.tsv"

FIG_DIR <- "figures"

FIG12A_PNG <- file.path(FIG_DIR, "12A_top50_rank_intervals.png")
FIG12A_PDF <- file.path(FIG_DIR, "12A_top50_rank_intervals.pdf")

FIG12B_PNG <- file.path(FIG_DIR, "12B_topk_inclusion_probabilities.png")
FIG12B_PDF <- file.path(FIG_DIR, "12B_topk_inclusion_probabilities.pdf")

FIG12C_PNG <- file.path(FIG_DIR, "12C_final_vs_bootstrap_rank.png")
FIG12C_PDF <- file.path(FIG_DIR, "12C_final_vs_bootstrap_rank.pdf")

FIG12D_PNG <- file.path(FIG_DIR, "12D_benchmark_rank_intervals.png")
FIG12D_PDF <- file.path(FIG_DIR, "12D_benchmark_rank_intervals.pdf")

BOOTSTRAP_SEED <- 20260912L
N_BOOTSTRAP <- 2000L

PRIMARY_WEIGHTS <- c(
    biology = 0.50,
    pu = 0.30,
    crispr = 0.20
)

COMPLETENESS_FLOOR <- 0.80

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

rank01 <- function(x) {
    x <- suppressWarnings(as.numeric(x))

    out <- rep(NA_real_, length(x))
    ok <- is.finite(x)
    n_ok <- sum(ok)

    if (n_ok == 0L) {
        return(out)
    }

    if (n_ok == 1L) {
        out[ok] <- 0.5
        return(out)
    }

    r <- rank(
        x[ok],
        ties.method = "average"
    )

    out[ok] <- (r - 1) / (n_ok - 1)
    out
}

weighted_available_score <- function(
    biology,
    pu,
    crispr,
    weights = PRIMARY_WEIGHTS,
    completeness_floor = COMPLETENESS_FLOOR
) {
    components <- cbind(
        biology = biology,
        pu = pu,
        crispr = crispr
    )

    observed <- is.finite(components)

    weight_matrix <- matrix(
        weights[colnames(components)],
        nrow = nrow(components),
        ncol = ncol(components),
        byrow = TRUE
    )

    observed_weight <- rowSums(
        observed * weight_matrix
    )

    numerator <- rowSums(
        ifelse(observed, components, 0) *
            weight_matrix
    )

    raw <- numerator / observed_weight
    raw[observed_weight <= 0] <- NA_real_

    domain_coverage <- rowMeans(observed)

    adjusted <- raw *
        (
            completeness_floor +
            (1 - completeness_floor) *
            domain_coverage
        )

    adjusted
}

deterministic_rank <- function(score, gene_id) {
    z <- data.table(
        row_index = seq_along(gene_id),
        gene_id = gene_id,
        score = suppressWarnings(as.numeric(score))
    )

    z[
        ,
        score_for_sort :=
            fifelse(
                is.finite(score),
                score,
                -Inf
            )
    ]

    setorder(
        z,
        -score_for_sort,
        gene_id
    )

    z[
        ,
        rank_value := seq_len(.N)
    ]

    out <- rep(NA_integer_, nrow(z))
    out[z$row_index] <- z$rank_value
    out[!is.finite(score)] <- NA_integer_

    out
}

quantile_safe <- function(x, prob) {
    x <- x[is.finite(x)]

    if (length(x) == 0L) {
        return(NA_real_)
    }

    as.numeric(
        stats::quantile(
            x,
            probs = prob,
            names = FALSE,
            type = 7
        )
    )
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
# 3. Load data
# -----------------------------------------------------------------------------

for (path in c(
    INTEGRATED_FILE,
    SELECTED_MODELS_FILE,
    CV_PREDICTIONS_FILE,
    BENCHMARK_SCORES_FILE
)) {
    assert_file(path)
}

integrated <- fread(INTEGRATED_FILE)
selected_models <- fread(SELECTED_MODELS_FILE)
cv_predictions <- fread(CV_PREDICTIONS_FILE)
benchmark_scores <- fread(BENCHMARK_SCORES_FILE)

assert_unique_gene_ids(
    integrated,
    "Step 09 integrated targeted panel"
)

required_integrated <- c(
    "gene_id",
    "biological_discovery_percentile",
    "crispr_tractability_percentile",
    "final_prepopulation_rank",
    "is_external_benchmark"
)

missing_integrated <- setdiff(
    required_integrated,
    names(integrated)
)

if (length(missing_integrated) > 0L) {
    stop(
        paste0(
            "Integrated table is missing required column(s): ",
            paste(missing_integrated, collapse = ", ")
        ),
        call. = FALSE
    )
}

if (!"model_key" %in% names(selected_models)) {
    selected_models[
        ,
        model_key :=
            paste(
                model_set,
                algorithm,
                sep = "::"
            )
    ]
}

selected_keys <- selected_models$model_key

if (length(selected_keys) < 2L) {
    stop(
        "At least two selected PU models are required.",
        call. = FALSE
    )
}

# -----------------------------------------------------------------------------
# 4. Build repeat-level PU percentile table for non-benchmark genes
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
    model_key %in% selected_keys
]

if (nrow(selected_cv) == 0L) {
    stop(
        "No Step 08 CV predictions matched the selected Step 09 models.",
        call. = FALSE
    )
}

selected_cv[
    ,
    repeat_model_percentile :=
        rank01(pu_score),
    by = .(
        model_key,
        repeat_id
    )
]

# Every gene should have one OOF prediction per model/repeat.
duplicate_cv <- selected_cv[
    ,
    .N,
    by = .(
        model_key,
        repeat_id,
        gene_id
    )
][
    N != 1L
]

if (nrow(duplicate_cv) > 0L) {
    stop(
        "Selected CV predictions are not unique by model/repeat/gene.",
        call. = FALSE
    )
}

repeat_ids <- sort(
    unique(
        selected_cv$repeat_id
    )
)

if (length(repeat_ids) < 2L) {
    stop(
        "At least two Step 08 repeats are required for resampling.",
        call. = FALSE
    )
}

# -----------------------------------------------------------------------------
# 5. Build external benchmark PU percentile table
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
    model_key %in% selected_keys
]

if (
    !"percentile_vs_crossfitted_nonholdout_genes" %in%
        names(selected_benchmark)
) {
    stop(
        paste0(
            "Benchmark score table lacks ",
            "percentile_vs_crossfitted_nonholdout_genes."
        ),
        call. = FALSE
    )
}

# -----------------------------------------------------------------------------
# 6. Prebuild model/repeat PU vectors on the exact 205-gene panel
# -----------------------------------------------------------------------------

panel_gene_ids <- integrated$gene_id
n_panel <- length(panel_gene_ids)

benchmark_mask <- as.logical(
    integrated$is_external_benchmark
)

benchmark_gene_ids <- integrated[
    benchmark_mask == TRUE,
    gene_id
]

novel_gene_ids <- integrated[
    benchmark_mask == FALSE,
    gene_id
]

pu_vectors <- list()

for (key in selected_keys) {
    pu_vectors[[key]] <- list()

    benchmark_key <- selected_benchmark[
        model_key == key,
        .(
            gene_id,
            external_percentile =
                percentile_vs_crossfitted_nonholdout_genes
        )
    ]

    for (rep_id in repeat_ids) {
        z <- selected_cv[
            model_key == key &
                repeat_id == rep_id,
            .(
                gene_id,
                pu_percentile =
                    repeat_model_percentile
            )
        ]

        vec <- z$pu_percentile[
            match(
                panel_gene_ids,
                z$gene_id
            )
        ]

        # External benchmarks are never in CV training pools, so replace their
        # missing OOF positions with fully external benchmark percentiles.
        if (nrow(benchmark_key) > 0L) {
            benchmark_positions <- match(
                benchmark_key$gene_id,
                panel_gene_ids
            )

            valid_positions <- !is.na(
                benchmark_positions
            )

            vec[
                benchmark_positions[
                    valid_positions
                ]
            ] <- benchmark_key$
                external_percentile[
                    valid_positions
                ]
        }

        pu_vectors[[key]][[as.character(rep_id)]] <- vec
    }
}

# Coverage QC before resampling.
coverage_matrix <- do.call(
    cbind,
    lapply(
        selected_keys,
        function(key) {
            pu_vectors[[key]][[as.character(repeat_ids[[1]])]]
        }
    )
)

n_no_pu_support <- sum(
    rowSums(
        is.finite(
            coverage_matrix
        )
    ) ==
        0L
)

if (n_no_pu_support > 0L) {
    stop(
        paste0(
            n_no_pu_support,
            " targeted genes have no PU support from any selected model."
        ),
        call. = FALSE
    )
}

# -----------------------------------------------------------------------------
# 7. Bootstrap / resampling analysis
# -----------------------------------------------------------------------------

biology <- integrated$biological_discovery_percentile
crispr <- integrated$crispr_tractability_percentile

rank_matrix <- matrix(
    NA_integer_,
    nrow = n_panel,
    ncol = N_BOOTSTRAP
)

score_matrix <- matrix(
    NA_real_,
    nrow = n_panel,
    ncol = N_BOOTSTRAP
)

set.seed(BOOTSTRAP_SEED)

for (b in seq_len(N_BOOTSTRAP)) {
    sampled_models <- sample(
        selected_keys,
        size = length(selected_keys),
        replace = TRUE
    )

    sampled_pu <- matrix(
        NA_real_,
        nrow = n_panel,
        ncol = length(sampled_models)
    )

    for (j in seq_along(sampled_models)) {
        key <- sampled_models[[j]]

        sampled_repeat <- sample(
            repeat_ids,
            size = 1L,
            replace = TRUE
        )

        sampled_pu[, j] <-
            pu_vectors[[key]][[as.character(sampled_repeat)]]
    }

    pu_boot <- apply(
        sampled_pu,
        1,
        function(z) {
            z <- z[is.finite(z)]

            if (length(z) == 0L) {
                return(NA_real_)
            }

            median(z)
        }
    )

    score_boot <- weighted_available_score(
        biology = biology,
        pu = pu_boot,
        crispr = crispr
    )

    rank_boot <- deterministic_rank(
        score = score_boot,
        gene_id = panel_gene_ids
    )

    score_matrix[, b] <- score_boot
    rank_matrix[, b] <- rank_boot

    if (
        b %% 250L == 0L ||
        b == N_BOOTSTRAP
    ) {
        log_step(
            "12",
            paste0(
                "Completed ",
                format(b, big.mark = ","),
                "/",
                format(N_BOOTSTRAP, big.mark = ","),
                " resampling replicates"
            )
        )
    }
}

# -----------------------------------------------------------------------------
# 8. Summarize empirical uncertainty
# -----------------------------------------------------------------------------

uncertainty <- copy(integrated)

uncertainty[
    ,
    bootstrap_rank_mean :=
        rowMeans(
            rank_matrix,
            na.rm = TRUE
        )
]

uncertainty[
    ,
    bootstrap_rank_median :=
        apply(
            rank_matrix,
            1,
            stats::median,
            na.rm = TRUE
        )
]

uncertainty[
    ,
    bootstrap_rank_sd :=
        apply(
            rank_matrix,
            1,
            stats::sd,
            na.rm = TRUE
        )
]

uncertainty[
    ,
    bootstrap_rank_q025 :=
        apply(
            rank_matrix,
            1,
            quantile_safe,
            prob = 0.025
        )
]

uncertainty[
    ,
    bootstrap_rank_q975 :=
        apply(
            rank_matrix,
            1,
            quantile_safe,
            prob = 0.975
        )
]

uncertainty[
    ,
    bootstrap_rank_interval_width :=
        bootstrap_rank_q975 -
        bootstrap_rank_q025
]

uncertainty[
    ,
    bootstrap_score_mean :=
        rowMeans(
            score_matrix,
            na.rm = TRUE
        )
]

uncertainty[
    ,
    bootstrap_score_median :=
        apply(
            score_matrix,
            1,
            stats::median,
            na.rm = TRUE
        )
]

uncertainty[
    ,
    bootstrap_score_q025 :=
        apply(
            score_matrix,
            1,
            quantile_safe,
            prob = 0.025
        )
]

uncertainty[
    ,
    bootstrap_score_q975 :=
        apply(
            score_matrix,
            1,
            quantile_safe,
            prob = 0.975
        )
]

uncertainty[
    ,
    probability_top10 :=
        rowMeans(
            rank_matrix <= 10,
            na.rm = TRUE
        )
]

uncertainty[
    ,
    probability_top20 :=
        rowMeans(
            rank_matrix <= 20,
            na.rm = TRUE
        )
]

uncertainty[
    ,
    probability_top50 :=
        rowMeans(
            rank_matrix <= 50,
            na.rm = TRUE
        )
]

uncertainty[
    ,
    probability_top100 :=
        rowMeans(
            rank_matrix <= 100,
            na.rm = TRUE
        )
]

uncertainty[
    ,
    bootstrap_median_rank_shift_from_step09 :=
        bootstrap_rank_median -
        final_prepopulation_rank
]

# -----------------------------------------------------------------------------
# 9. Robust rank envelope with deterministic Step 09 weighting scenarios
# -----------------------------------------------------------------------------

scenario_rank_candidates <- intersect(
    c(
        "rank_biology_led",
        "rank_primary",
        "rank_validation_led"
    ),
    names(uncertainty)
)

if (length(scenario_rank_candidates) > 0L) {
    uncertainty[
        ,
        scenario_rank_min :=
            apply(
                .SD,
                1,
                function(z) {
                    z <- as.numeric(z)
                    z <- z[is.finite(z)]

                    if (length(z) == 0L) {
                        return(NA_real_)
                    }

                    min(z)
                }
            ),
        .SDcols = scenario_rank_candidates
    ]

    uncertainty[
        ,
        scenario_rank_max :=
            apply(
                .SD,
                1,
                function(z) {
                    z <- as.numeric(z)
                    z <- z[is.finite(z)]

                    if (length(z) == 0L) {
                        return(NA_real_)
                    }

                    max(z)
                }
            ),
        .SDcols = scenario_rank_candidates
    ]
} else {
    uncertainty[
        ,
        `:=`(
            scenario_rank_min = NA_real_,
            scenario_rank_max = NA_real_
        )
    ]
}

uncertainty[
    ,
    robust_rank_envelope_min :=
        pmin(
            bootstrap_rank_q025,
            scenario_rank_min,
            na.rm = TRUE
        )
]

uncertainty[
    ,
    robust_rank_envelope_max :=
        pmax(
            bootstrap_rank_q975,
            scenario_rank_max,
            na.rm = TRUE
        )
]

uncertainty[
    ,
    robust_rank_envelope_width :=
        robust_rank_envelope_max -
        robust_rank_envelope_min
]

# pmin/pmax with all NA can return Inf/-Inf; normalize.
uncertainty[
    !is.finite(robust_rank_envelope_min),
    robust_rank_envelope_min := NA_real_
]

uncertainty[
    !is.finite(robust_rank_envelope_max),
    robust_rank_envelope_max := NA_real_
]

uncertainty[
    !is.finite(robust_rank_envelope_width),
    robust_rank_envelope_width := NA_real_
]

# Descriptive stability class.
uncertainty[
    ,
    bootstrap_stability_class :=
        fcase(
            probability_top20 >= 0.95 &
                bootstrap_rank_interval_width <= 10,
            "very_high",

            probability_top50 >= 0.95 &
                bootstrap_rank_interval_width <= 20,
            "high",

            probability_top50 >= 0.75,
            "moderate",

            default = "variable"
        )
]

# Stable ordering by frozen Step 09 rank.
setorder(
    uncertainty,
    final_prepopulation_rank,
    gene_id,
    na.last = TRUE
)

top20_uncertainty <- uncertainty[
    is_external_benchmark == FALSE
][
    1:min(20L, .N)
]

benchmark_uncertainty <- uncertainty[
    is_external_benchmark == TRUE
]

# -----------------------------------------------------------------------------
# 10. Integrity checks
# -----------------------------------------------------------------------------

if (nrow(uncertainty) != n_panel) {
    stop(
        "Uncertainty output row count differs from the Step 09 targeted panel.",
        call. = FALSE
    )
}

if (uniqueN(uncertainty$gene_id) != nrow(uncertainty)) {
    stop(
        "Uncertainty output contains duplicate gene IDs.",
        call. = FALSE
    )
}

if (
    any(
        uncertainty$probability_top10 < 0 |
            uncertainty$probability_top10 > 1,
        na.rm = TRUE
    ) ||
    any(
        uncertainty$probability_top20 < 0 |
            uncertainty$probability_top20 > 1,
        na.rm = TRUE
    ) ||
    any(
        uncertainty$probability_top50 < 0 |
            uncertainty$probability_top50 > 1,
        na.rm = TRUE
    ) ||
    any(
        uncertainty$probability_top100 < 0 |
            uncertainty$probability_top100 > 1,
        na.rm = TRUE
    )
) {
    stop(
        "One or more top-k probabilities fall outside [0,1].",
        call. = FALSE
    )
}

# -----------------------------------------------------------------------------
# 11. Publication figures
# -----------------------------------------------------------------------------

# Figure 12A â€” top-50 novel rank intervals
plot_top50 <- uncertainty[
    is_external_benchmark == FALSE &
        final_prepopulation_rank <= 50
]

plot_top50[
    ,
    display_label :=
        factor(
            gene_id,
            levels = rev(
                gene_id[
                    order(
                        final_prepopulation_rank
                    )
                ]
            )
        )
]

p12a <- ggplot(
    plot_top50,
    aes(
        y = display_label
    )
) +
    geom_segment(
        aes(
            x = bootstrap_rank_q025,
            xend = bootstrap_rank_q975,
            yend = display_label
        )
    ) +
    geom_point(
        aes(
            x = bootstrap_rank_median
        )
    ) +
    geom_point(
        aes(
            x = final_prepopulation_rank
        ),
        shape = 4
    ) +
    labs(
        x = "Rank",
        y = NULL,
        title = "PU-resampling uncertainty for the top 50 novel candidates",
        subtitle = "Point = bootstrap median; Ã— = frozen Step 09 rank; line = empirical 95% rank interval"
    ) +
    theme_classic(base_size = 10)

save_plot_pair(
    p12a,
    FIG12A_PNG,
    FIG12A_PDF,
    width = 9,
    height = 12
)

# Figure 12B â€” top-k probabilities for top-30 novel genes
plot_topk <- uncertainty[
    is_external_benchmark == FALSE &
        final_prepopulation_rank <= 30,
    .(
        gene_id,
        final_prepopulation_rank,
        probability_top10,
        probability_top20,
        probability_top50,
        probability_top100
    )
]

plot_topk_long <- melt(
    plot_topk,
    id.vars = c(
        "gene_id",
        "final_prepopulation_rank"
    ),
    measure.vars = c(
        "probability_top10",
        "probability_top20",
        "probability_top50",
        "probability_top100"
    ),
    variable.name = "threshold",
    value.name = "probability"
)

plot_topk_long[
    ,
    threshold :=
        factor(
            threshold,
            levels = c(
                "probability_top10",
                "probability_top20",
                "probability_top50",
                "probability_top100"
            ),
            labels = c(
                "Top 10",
                "Top 20",
                "Top 50",
                "Top 100"
            )
        )
]

gene_order_top30 <- plot_topk[
    order(final_prepopulation_rank),
    gene_id
]

plot_topk_long[
    ,
    gene_id :=
        factor(
            gene_id,
            levels = rev(gene_order_top30)
        )
]

p12b <- ggplot(
    plot_topk_long,
    aes(
        x = threshold,
        y = gene_id,
        fill = probability
    )
) +
    geom_tile() +
    geom_text(
        aes(
            label = sprintf("%.2f", probability)
        ),
        size = 2.8
    ) +
    scale_fill_viridis_c(
        limits = c(0, 1),
        option = "C"
    ) +
    labs(
        x = NULL,
        y = NULL,
        fill = "Probability",
        title = "Probability of retaining high priority under PU resampling"
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
        panel.grid = element_blank()
    )

save_plot_pair(
    p12b,
    FIG12B_PNG,
    FIG12B_PDF,
    width = 8.2,
    height = 9.5
)

# Figure 12C â€” frozen final rank vs bootstrap median rank
plot_compare <- uncertainty[
    is_external_benchmark == FALSE
]

p12c <- ggplot(
    plot_compare,
    aes(
        x = final_prepopulation_rank,
        y = bootstrap_rank_median
    )
) +
    geom_abline(
        slope = 1,
        intercept = 0,
        linetype = 2
    ) +
    geom_point(
        aes(
            size = bootstrap_rank_interval_width
        )
    ) +
    labs(
        x = "Frozen Step 09 final rank",
        y = "Bootstrap median rank",
        size = "95% rank\ninterval width",
        title = "Concordance between frozen and resampled candidate ranks"
    ) +
    theme_classic(base_size = 11)

save_plot_pair(
    p12c,
    FIG12C_PNG,
    FIG12C_PDF,
    width = 7.5,
    height = 6.5
)

# Figure 12D â€” benchmark uncertainty
benchmark_plot <- copy(
    benchmark_uncertainty
)

benchmark_plot[
    ,
    display_label :=
        if ("benchmark_name" %in% names(benchmark_plot)) {
            paste0(
                gene_id,
                "\n",
                benchmark_name
            )
        } else {
            gene_id
        }
]

benchmark_plot[
    ,
    display_label :=
        factor(
            display_label,
            levels = rev(display_label)
        )
]

p12d <- ggplot(
    benchmark_plot,
    aes(
        y = display_label
    )
) +
    geom_segment(
        aes(
            x = bootstrap_rank_q025,
            xend = bootstrap_rank_q975,
            yend = display_label
        )
    ) +
    geom_point(
        aes(
            x = bootstrap_rank_median
        )
    ) +
    geom_point(
        aes(
            x = final_prepopulation_rank
        ),
        shape = 4
    ) +
    labs(
        x = "Rank among 205 targeted genes",
        y = NULL,
        title = "External benchmark rank stability under PU model resampling",
        subtitle = "Benchmarks remained excluded from every fitted model"
    ) +
    theme_classic(base_size = 10.5)

save_plot_pair(
    p12d,
    FIG12D_PNG,
    FIG12D_PDF,
    width = 8.5,
    height = 5.8
)

# -----------------------------------------------------------------------------
# 12. QC
# -----------------------------------------------------------------------------

top20_novel <- uncertainty[
    is_external_benchmark == FALSE &
        final_prepopulation_rank <= 20
]

n_top20_prob_ge_075 <- top20_novel[
    probability_top20 >= 0.75,
    .N
]

n_top20_prob_ge_095 <- top20_novel[
    probability_top20 >= 0.95,
    .N
]

spearman_rank <- suppressWarnings(
    stats::cor(
        uncertainty[
            is_external_benchmark == FALSE,
            final_prepopulation_rank
        ],
        uncertainty[
            is_external_benchmark == FALSE,
            bootstrap_rank_median
        ],
        method = "spearman",
        use = "complete.obs"
    )
)

qc <- data.table(
    metric = c(
        "Targeted panel genes",
        "Novel genes",
        "External benchmark genes",
        "Selected PU models",
        "Step 08 repeat IDs available",
        "Bootstrap/resampling replicates",
        "Genes lacking all selected-model PU support",
        "Top-20 novel genes with P(top20) >= 0.75",
        "Top-20 novel genes with P(top20) >= 0.95",
        "Spearman frozen-rank vs bootstrap-median-rank",
        "Leave-one-landmark-out models fitted",
        "Benchmark genes used for fitting",
        "Population-genomics predictors used",
        "Figure 12A written",
        "Figure 12B written",
        "Figure 12C written",
        "Figure 12D written"
    ),

    value = c(
        nrow(uncertainty),
        sum(!uncertainty$is_external_benchmark, na.rm = TRUE),
        sum(uncertainty$is_external_benchmark, na.rm = TRUE),
        length(selected_keys),
        length(repeat_ids),
        N_BOOTSTRAP,
        n_no_pu_support,
        n_top20_prob_ge_075,
        n_top20_prob_ge_095,
        spearman_rank,
        0,
        0,
        0,
        as.integer(
            file.exists(FIG12A_PNG) &&
                file.exists(FIG12A_PDF)
        ),
        as.integer(
            file.exists(FIG12B_PNG) &&
                file.exists(FIG12B_PDF)
        ),
        as.integer(
            file.exists(FIG12C_PNG) &&
                file.exists(FIG12C_PDF)
        ),
        as.integer(
            file.exists(FIG12D_PNG) &&
                file.exists(FIG12D_PDF)
        )
    )
)

# -----------------------------------------------------------------------------
# 13. Provenance
# -----------------------------------------------------------------------------

provenance <- data.table(
    component = c(
        "Uncertainty target",
        "Resampling unit",
        "PU model uncertainty",
        "PU repeat uncertainty",
        "Benchmark treatment",
        "Biological evidence",
        "CRISPR evidence",
        "Integration formula",
        "Empirical rank interval",
        "Top-k probabilities",
        "Weighting robustness",
        "Leave-one-landmark-out",
        "Population genomics"
    ),

    specification = c(
        "Step 09 pre-population integrated rank on the 205-gene targeted panel.",

        paste0(
            N_BOOTSTRAP,
            " resampling replicates with selected PU model keys sampled with replacement."
        ),

        paste0(
            "Each replicate resamples the ",
            length(selected_keys),
            " selected Step 09 PU models with replacement."
        ),

        paste0(
            "For each sampled PU model, one of ",
            length(repeat_ids),
            " Step 08 repeated-CV OOF prediction sets is sampled."
        ),

        paste0(
            "Benchmarks use the fully external Step 08 model-specific percentiles; ",
            "they are never inserted into OOF training data."
        ),

        "Step 04 biological discovery percentiles are held fixed.",

        "Step 05 exon-aware CRISPR tractability percentiles are held fixed.",

        paste0(
            "Primary Step 09 integration retained: biology=0.50, PU=0.30, ",
            "CRISPR=0.20 with completeness floor ",
            COMPLETENESS_FLOOR,
            "."
        ),

        "2.5th and 97.5th empirical quantiles of resampled ranks.",

        "Fraction of resampling replicates in which each gene appears in top 10/20/50/100.",

        paste0(
            "A descriptive robust rank envelope combines the PU-resampling rank interval ",
            "with Step 09 biology-led, primary, and validation-led scenario ranks. ",
            "This envelope is not interpreted as a formal confidence interval."
        ),

        paste0(
            "Not performed because all six landmarks were already excluded from every ",
            "Step 08 and Step 11 model fit; the existing complete landmark holdout is ",
            "stronger than leave-one-landmark-out."
        ),

        "Ag1000G remains pending and is not used in Step 12."
    )
)

# -----------------------------------------------------------------------------
# 14. Write outputs
# -----------------------------------------------------------------------------

fwrite(
    uncertainty,
    UNCERTAINTY_FILE,
    na = "NA"
)

fwrite(
    top20_uncertainty,
    TOP20_FILE,
    na = "NA"
)

fwrite(
    benchmark_uncertainty,
    BENCHMARK_FILE,
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
    INTEGRATED_FILE,
    SELECTED_MODELS_FILE,
    CV_PREDICTIONS_FILE,
    BENCHMARK_SCORES_FILE,
    UNCERTAINTY_FILE,
    TOP20_FILE,
    BENCHMARK_FILE,
    QC_FILE,
    PROVENANCE_FILE,
    SESSION_INFO_FILE,
    FIG12A_PNG,
    FIG12A_PDF,
    FIG12B_PNG,
    FIG12B_PDF,
    FIG12C_PNG,
    FIG12C_PDF,
    FIG12D_PNG,
    FIG12D_PDF
)

checksum_files <- checksum_files[
    file.exists(checksum_files)
]

write_checksum(
    checksum_files,
    CHECKSUM_FILE
)

# -----------------------------------------------------------------------------
# 15. Console summary
# -----------------------------------------------------------------------------

cat(
    "\n",
    "============================================================\n",
    "MOSQEDIT-R STEP 12 COMPLETED SUCCESSFULLY\n",
    "PU-resampling rank uncertainty\n",
    "============================================================\n",
    "Targeted panel genes:                     ",
    nrow(uncertainty),
    "\n",
    "Novel genes:                              ",
    sum(!uncertainty$is_external_benchmark, na.rm = TRUE),
    "\n",
    "External benchmark genes:                 ",
    sum(uncertainty$is_external_benchmark, na.rm = TRUE),
    "\n",
    "Selected PU models:                       ",
    length(selected_keys),
    "\n",
    "Step 08 CV repeats available:             ",
    length(repeat_ids),
    "\n",
    "Resampling replicates:                    ",
    N_BOOTSTRAP,
    "\n",
    "Top-20 with P(top20) >= 0.75:             ",
    n_top20_prob_ge_075,
    "\n",
    "Top-20 with P(top20) >= 0.95:             ",
    n_top20_prob_ge_095,
    "\n",
    "Frozen-vs-bootstrap Spearman rho:          ",
    format(
        spearman_rank,
        digits = 4
    ),
    "\n",
    "Leave-one-landmark-out fitted:            NO (redundant)\n",
    "Benchmark genes used for fitting:         NO\n",
    "Population genomics integrated:           NO\n",
    "Publication figures written:              4 PNG + 4 PDF\n",
    "Python used:                              NO\n",
    "============================================================\n",
    sep = ""
)

cat(
    "\nTop 20 novel candidates with rank uncertainty:\n"
)

display_cols <- c(
    "final_prepopulation_rank",
    "gene_id",
    "preliminary_rank",
    "bootstrap_rank_median",
    "bootstrap_rank_q025",
    "bootstrap_rank_q975",
    "bootstrap_rank_interval_width",
    "probability_top10",
    "probability_top20",
    "probability_top50",
    "robust_rank_envelope_min",
    "robust_rank_envelope_max",
    "bootstrap_stability_class"
)

display_cols <- display_cols[
    display_cols %in% names(top20_uncertainty)
]

print(
    top20_uncertainty[
        ,
        ..display_cols
    ]
)

cat(
    "\nExternal benchmark rank uncertainty:\n"
)

benchmark_display <- c(
    "gene_id",
    "benchmark_name",
    "final_prepopulation_rank",
    "bootstrap_rank_median",
    "bootstrap_rank_q025",
    "bootstrap_rank_q975",
    "probability_top20",
    "probability_top50",
    "robust_rank_envelope_min",
    "robust_rank_envelope_max"
)

benchmark_display <- benchmark_display[
    benchmark_display %in% names(benchmark_uncertainty)
]

print(
    benchmark_uncertainty[
        ,
        ..benchmark_display
    ]
)

cat(
    "\nQC summary:\n"
)

print(qc)

log_step(
    "12",
    "PU-resampling rank-uncertainty analysis completed successfully"
)

