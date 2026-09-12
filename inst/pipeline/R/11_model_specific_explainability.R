# =============================================================================
# 11_model_specific_explainability.R
#
# MosqEdit-R Manuscript 1
#
# STEP 11
# Model-specific explainability for selected PU models
#
# PURE R
#
# PURPOSE
# -------
# Reconstruct the selected Step 08 PU models using the same:
#   - modelling datasets
#   - predictors
#   - fold assignment
#   - repeat structure
#   - bagging seeds
#   - hyperparameters
#   - fold-safe preprocessing
#
# Then compute:
#
# 1. Cross-fitted XGBoost SHAP importance
#    - global mean absolute SHAP across out-of-fold genes
#    - local SHAP for Step 10 top-20 novel candidates
#    - fully external SHAP for six benchmark genes
#
# 2. Cross-fitted bag-averaged Ranger permutation importance
#
# 3. Cross-model consensus feature importance
#
# IMPORTANT
# ---------
# For non-benchmark genes, SHAP values are calculated only from models in which
# the gene was in the held-out fold. This avoids in-sample local explanations.
#
# Benchmark genes remain absent from all model fitting and are therefore fully
# external explanations.
#
# FlyBase phenotype variables, Step 04 ranking variables, Step 05 CRISPR
# variables, and Ag1000G variables remain excluded from model fitting.
#
# INPUTS
# ------
# data_processed/07_pu_M2_expression_orthology.csv
# data_processed/07_pu_M3_expression_orthology_annotation.csv
# data_processed/07_predictor_manifest.csv
# data_processed/master_gene_feature_matrix.csv
# data_processed/08_model_configuration.csv
# data_processed/09_selected_pu_models.csv
# data_processed/10_top20_candidate_evidence.csv
# data_processed/10_benchmark_evidence.csv
#
# OUTPUTS
# -------
# data_processed/11_xgboost_global_shap_importance.csv
# data_processed/11_xgboost_local_shap_top20.csv
# data_processed/11_xgboost_local_shap_benchmarks.csv
# data_processed/11_ranger_permutation_importance.csv
# data_processed/11_consensus_feature_importance.csv
# data_processed/11_refit_summary.csv
# data_processed/11_explainability_qc.csv
# data_processed/11_explainability_provenance.csv
#
# figures/11A_consensus_global_importance.png/.pdf
# figures/11B_M3_xgboost_global_shap.png/.pdf
# figures/11C_top20_local_shap_heatmap.png/.pdf
# figures/11D_benchmark_local_shap_heatmap.png/.pdf
#
# logs/11_explainability_sessionInfo.txt
# logs/11_checksums.tsv
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Helpers and packages
# -----------------------------------------------------------------------------

source("R/helpers.R")


required_packages <- c(
    "data.table",
    "ggplot2",
    "xgboost",
    "ranger"
)


missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]


if (length(missing_packages) > 0L) {

    stop(
        paste0(
            "Missing required package(s): ",
            paste(
                missing_packages,
                collapse = ", "
            ),
            "\nInstall with:\ninstall.packages(c(",
            paste0(
                '"',
                missing_packages,
                '"',
                collapse = ", "
            ),
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
    "11",
    "Starting cross-fitted model-specific explainability"
)


# -----------------------------------------------------------------------------
# 2. Files
# -----------------------------------------------------------------------------

M2_FILE <-
    "data_processed/07_pu_M2_expression_orthology.csv"


M3_FILE <-
    "data_processed/07_pu_M3_expression_orthology_annotation.csv"


PREDICTOR_MANIFEST_FILE <-
    "data_processed/07_predictor_manifest.csv"


MASTER_FILE <-
    "data_processed/master_gene_feature_matrix.csv"


STEP08_CONFIG_FILE <-
    "data_processed/08_model_configuration.csv"


SELECTED_MODELS_FILE <-
    "data_processed/09_selected_pu_models.csv"


TOP20_FILE <-
    "data_processed/10_top20_candidate_evidence.csv"


BENCHMARK_FILE <-
    "data_processed/10_benchmark_evidence.csv"


XGB_GLOBAL_FILE <-
    "data_processed/11_xgboost_global_shap_importance.csv"


XGB_TOP20_LOCAL_FILE <-
    "data_processed/11_xgboost_local_shap_top20.csv"


XGB_BENCH_LOCAL_FILE <-
    "data_processed/11_xgboost_local_shap_benchmarks.csv"


RANGER_IMPORTANCE_FILE <-
    "data_processed/11_ranger_permutation_importance.csv"


CONSENSUS_IMPORTANCE_FILE <-
    "data_processed/11_consensus_feature_importance.csv"


REFIT_SUMMARY_FILE <-
    "data_processed/11_refit_summary.csv"


QC_FILE <-
    "data_processed/11_explainability_qc.csv"


PROVENANCE_FILE <-
    "data_processed/11_explainability_provenance.csv"


SESSION_INFO_FILE <-
    "logs/11_explainability_sessionInfo.txt"


CHECKSUM_FILE <-
    "logs/11_checksums.tsv"


FIG_DIR <- "figures"


FIG11A_PNG <-
    file.path(
        FIG_DIR,
        "11A_consensus_global_importance.png"
    )


FIG11A_PDF <-
    file.path(
        FIG_DIR,
        "11A_consensus_global_importance.pdf"
    )


FIG11B_PNG <-
    file.path(
        FIG_DIR,
        "11B_M3_xgboost_global_shap.png"
    )


FIG11B_PDF <-
    file.path(
        FIG_DIR,
        "11B_M3_xgboost_global_shap.pdf"
    )


FIG11C_PNG <-
    file.path(
        FIG_DIR,
        "11C_top20_local_shap_heatmap.png"
    )


FIG11C_PDF <-
    file.path(
        FIG_DIR,
        "11C_top20_local_shap_heatmap.pdf"
    )


FIG11D_PNG <-
    file.path(
        FIG_DIR,
        "11D_benchmark_local_shap_heatmap.png"
    )


FIG11D_PDF <-
    file.path(
        FIG_DIR,
        "11D_benchmark_local_shap_heatmap.pdf"
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


dir.create(
    FIG_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# 3. Basic helpers
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


safe_numeric <- function(x) {

    suppressWarnings(
        as.numeric(
            x
        )
    )
}


rank01 <- function(x) {

    x <- safe_numeric(
        x
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


    r <- rank(
        x[ok],
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


save_plot_pair <- function(
        plot_object,
        png_file,
        pdf_file,
        width,
        height
) {

    ggsave(
        filename =
            png_file,

        plot =
            plot_object,

        width =
            width,

        height =
            height,

        units =
            "in",

        dpi =
            400,

        bg =
            "white"
    )


    ggsave(
        filename =
            pdf_file,

        plot =
            plot_object,

        width =
            width,

        height =
            height,

        units =
            "in",

        device =
            cairo_pdf
    )
}


# -----------------------------------------------------------------------------
# 4. Read Step 08 configuration
# -----------------------------------------------------------------------------

for (
    path in c(
        M2_FILE,
        M3_FILE,
        PREDICTOR_MANIFEST_FILE,
        MASTER_FILE,
        STEP08_CONFIG_FILE,
        SELECTED_MODELS_FILE,
        TOP20_FILE,
        BENCHMARK_FILE
    )
) {

    assert_file(
        path
    )
}


config <- fread(
    STEP08_CONFIG_FILE
)


if (
    !all(
        c(
            "parameter",
            "value"
        ) %in%
            names(
                config
            )
    )
) {

    stop(
        "Step 08 configuration file lacks parameter/value columns.",
        call. = FALSE
    )
}


get_cfg <- function(
        name,
        default = NA_real_
) {

    z <- config[
        parameter ==
            name,
        value
    ]


    if (
        length(
            z
        ) ==
            0L
    ) {

        return(
            default
        )
    }


    as.numeric(
        z[[1]]
    )
}


MASTER_SEED <-
    as.integer(
        get_cfg(
            "MASTER_SEED",
            20260911
        )
    )


N_FOLDS <-
    as.integer(
        get_cfg(
            "N_FOLDS",
            5
        )
    )


N_REPEATS <-
    as.integer(
        get_cfg(
            "N_REPEATS",
            3
        )
    )


N_BAGS <-
    as.integer(
        get_cfg(
            "N_BAGS",
            10
        )
    )


UNLABELED_TO_POSITIVE_RATIO <-
    get_cfg(
        "UNLABELED_TO_POSITIVE_RATIO",
        1
    )


RANGER_NUM_TREES <-
    as.integer(
        get_cfg(
            "RANGER_NUM_TREES",
            400
        )
    )


RANGER_MIN_NODE_SIZE <-
    as.integer(
        get_cfg(
            "RANGER_MIN_NODE_SIZE",
            5
        )
    )


XGB_NROUNDS <-
    as.integer(
        get_cfg(
            "XGB_NROUNDS",
            150
        )
    )


XGB_MAX_DEPTH <-
    as.integer(
        get_cfg(
            "XGB_MAX_DEPTH",
            3
        )
    )


XGB_ETA <-
    get_cfg(
        "XGB_ETA",
        0.05
    )


XGB_SUBSAMPLE <-
    get_cfg(
        "XGB_SUBSAMPLE",
        0.80
    )


XGB_COLSAMPLE <-
    get_cfg(
        "XGB_COLSAMPLE",
        0.80
    )


XGB_MIN_CHILD_WEIGHT <-
    get_cfg(
        "XGB_MIN_CHILD_WEIGHT",
        2
    )


MODEL_SETS_STEP08 <- c(
    "M1_expression",
    "M2_expression_orthology",
    "M3_expression_orthology_annotation"
)


ALGORITHMS_STEP08 <- c(
    "elastic_net",
    "ranger",
    "xgboost",
    "svm_radial"
)


# -----------------------------------------------------------------------------
# 5. Load modelling and explanation targets
# -----------------------------------------------------------------------------

datasets <- list(

    M2_expression_orthology =
        fread(
            M2_FILE
        ),

    M3_expression_orthology_annotation =
        fread(
            M3_FILE
        )
)


predictor_manifest <- fread(
    PREDICTOR_MANIFEST_FILE
)


master <- fread(
    MASTER_FILE
)


selected_models <- fread(
    SELECTED_MODELS_FILE
)


top20 <- fread(
    TOP20_FILE
)


benchmarks <- fread(
    BENCHMARK_FILE
)


for (
    nm in names(
        datasets
    )
) {

    assert_unique_gene_ids(
        datasets[[nm]],
        nm
    )
}


assert_unique_gene_ids(
    master,
    "Master feature matrix"
)


assert_unique_gene_ids(
    top20,
    "Step 10 top-20 table"
)


assert_unique_gene_ids(
    benchmarks,
    "Step 10 benchmark table"
)


# Recreate gene_length_bp exactly as in Step 07 if needed.
if (
    !"gene_length_bp" %in%
        names(
            master
        ) &&
    all(
        c(
            "gene_start",
            "gene_end"
        ) %in%
            names(
                master
            )
    )
) {

    master[
        ,
        gene_length_bp :=
            as.numeric(
                gene_end
            ) -
            as.numeric(
                gene_start
            ) +
            1
    ]
}


# External benchmark master rows.
if (
    "is_external_landmark_validation" %in%
        names(
            master
        )
) {

    master[
        ,
        external_validation_holdout :=
            !is.na(
                is_external_landmark_validation
            ) &
            as.logical(
                is_external_landmark_validation
            )
    ]

} else {

    stop(
        "Master matrix lacks is_external_landmark_validation.",
        call. = FALSE
    )
}


benchmark_master <- master[
    external_validation_holdout ==
        TRUE
]


top20_gene_ids <- top20$gene_id


# -----------------------------------------------------------------------------
# 6. Restrict to selected explainable algorithms
# -----------------------------------------------------------------------------

if (
    !"model_key" %in%
        names(
            selected_models
        )
) {

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


supported_selected <- selected_models[
    algorithm %in%
        c(
            "xgboost",
            "ranger"
        )
]


if (
    nrow(
        supported_selected
    ) ==
        0L
) {

    stop(
        paste0(
            "None of the selected PU models uses XGBoost or Ranger; ",
            "Step 11 has nothing to explain."
        ),
        call. = FALSE
    )
}


selected_xgb <- supported_selected[
    algorithm ==
        "xgboost"
]


selected_ranger <- supported_selected[
    algorithm ==
        "ranger"
]


# All selected model sets must be available.
missing_model_sets <- setdiff(
    unique(
        supported_selected$
            model_set
    ),
    names(
        datasets
    )
)


if (
    length(
        missing_model_sets
    ) >
        0L
) {

    stop(
        paste0(
            "Selected model set(s) not available to Step 11: ",
            paste(
                missing_model_sets,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}


# -----------------------------------------------------------------------------
# 7. Predictor sets
# -----------------------------------------------------------------------------

get_predictors <- function(
        model_set_name
) {

    unique(
        predictor_manifest[
            modelling_set ==
                model_set_name &
                allowed_for_pu_training ==
                TRUE,
            predictor
        ]
    )
}


predictor_sets <- lapply(
    unique(
        supported_selected$
            model_set
    ),
    get_predictors
)


names(
    predictor_sets
) <- unique(
    supported_selected$
        model_set
)


for (
    nm in names(
        predictor_sets
    )
) {

    missing_predictors <- setdiff(
        predictor_sets[[nm]],
        names(
            datasets[[nm]]
        )
    )


    if (
        length(
            missing_predictors
        ) >
            0L
    ) {

        stop(
            paste0(
                nm,
                " is missing predictor(s): ",
                paste(
                    missing_predictors,
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
    }
}


# -----------------------------------------------------------------------------
# 8. Stratified folds â€” identical logic to Step 08
# -----------------------------------------------------------------------------

make_stratified_folds <- function(
        y,
        k,
        seed
) {

    set.seed(
        seed
    )


    fold <- integer(
        length(
            y
        )
    )


    for (
        cls in c(
            0L,
            1L
        )
    ) {

        idx <- which(
            y ==
                cls
        )


        idx <- sample(
            idx,
            length(
                idx
            ),
            replace = FALSE
        )


        fold[
            idx
        ] <- rep(
            seq_len(
                k
            ),
            length.out =
                length(
                    idx
                )
        )
    }


    fold
}


# -----------------------------------------------------------------------------
# 9. Fold-safe preprocessing â€” identical design to Step 08
# -----------------------------------------------------------------------------

fit_preprocessor <- function(
        train_dt,
        predictors
) {

    specs <- list()

    pieces <- list()


    for (
        col in predictors
    ) {

        x <- train_dt[[col]]


        if (
            is.numeric(
                x
            ) ||
            is.integer(
                x
            ) ||
            is.logical(
                x
            )
        ) {

            raw <- if (
                is.logical(
                    x
                )
            ) {

                as.numeric(
                    x
                )

            } else {

                safe_numeric(
                    x
                )
            }


            missing <- is.na(
                raw
            ) |
                !is.finite(
                    raw
                )


            finite_values <- raw[
                !missing
            ]


            med <- if (
                length(
                    finite_values
                ) >
                    0L
            ) {

                median(
                    finite_values,
                    na.rm = TRUE
                )

            } else {

                0
            }


            imp <- raw

            imp[
                missing
            ] <- med


            mu <- mean(
                imp
            )


            sdv <- stats::sd(
                imp
            )


            if (
                !is.finite(
                    sdv
                ) ||
                sdv ==
                    0
            ) {

                sdv <- 1
            }


            mat <- cbind(
                (
                    imp -
                    mu
                ) /
                    sdv,

                as.numeric(
                    missing
                )
            )


            colnames(
                mat
            ) <- c(
                paste0(
                    "num__",
                    col
                ),
                paste0(
                    "missing__",
                    col
                )
            )


            pieces[[length( pieces ) + 1L]] <- mat


            specs[[col]] <- list(
                type =
                    "numeric",

                median =
                    med,

                mean =
                    mu,

                sd =
                    sdv
            )

        } else {

            raw <- as.character(
                x
            )


            raw[
                is.na(
                    raw
                ) |
                    !nzchar(
                        raw
                    )
            ] <- "__MISSING__"


            levels_train <- sort(
                unique(
                    raw
                )
            )


            levels_all <- unique(
                c(
                    levels_train,
                    "__MISSING__",
                    "__OTHER__"
                )
            )


            mat <- matrix(
                0,
                nrow =
                    length(
                        raw
                    ),
                ncol =
                    length(
                        levels_all
                    )
            )


            colnames(
                mat
            ) <- paste0(
                "cat__",
                col,
                "__",
                make.names(
                    levels_all,
                    unique = TRUE
                )
            )


            for (
                j in seq_along(
                    levels_all
                )
            ) {

                mat[
                    ,
                    j
                ] <- as.numeric(
                    raw ==
                        levels_all[[j]]
                )
            }


            pieces[[length( pieces ) + 1L]] <- mat


            specs[[col]] <- list(
                type =
                    "categorical",

                levels =
                    levels_all
            )
        }
    }


    X <- do.call(
        cbind,
        pieces
    )


    X <- as.matrix(
        X
    )


    storage.mode(
        X
    ) <- "double"


    colnames(
        X
    ) <- make.names(
        colnames(
            X
        ),
        unique = TRUE
    )


    keep <- vapply(
        seq_len(
            ncol(
                X
            )
        ),
        function(j) {

            z <- X[
                ,
                j
            ]


            all(
                is.finite(
                    z
                )
            ) &&
                is.finite(
                    stats::sd(
                        z
                    )
                ) &&
                stats::sd(
                    z
                ) >
                0
        },
        logical(
            1
        )
    )


    if (!any(keep)) {

        stop(
            "Preprocessing produced no nonconstant predictors.",
            call. = FALSE
        )
    }


    X <- X[
        ,
        keep,
        drop = FALSE
    ]


    list(
        specs =
            specs,

        matrix_columns =
            colnames(
                X
            ),

        X =
            X
    )
}


apply_preprocessor <- function(
        dt,
        prep
) {

    pieces <- list()


    for (
        col in names(
            prep$specs
        )
    ) {

        spec <- prep$specs[[col]]


        x <- dt[[col]]


        if (
            spec$type ==
                "numeric"
        ) {

            raw <- if (
                is.logical(
                    x
                )
            ) {

                as.numeric(
                    x
                )

            } else {

                safe_numeric(
                    x
                )
            }


            missing <- is.na(
                raw
            ) |
                !is.finite(
                    raw
                )


            imp <- raw

            imp[
                missing
            ] <- spec$median


            mat <- cbind(
                (
                    imp -
                    spec$mean
                ) /
                    spec$sd,

                as.numeric(
                    missing
                )
            )


            colnames(
                mat
            ) <- c(
                paste0(
                    "num__",
                    col
                ),
                paste0(
                    "missing__",
                    col
                )
            )


            pieces[[length( pieces ) + 1L]] <- mat

        } else {

            raw <- as.character(
                x
            )


            raw[
                is.na(
                    raw
                ) |
                    !nzchar(
                        raw
                    )
            ] <- "__MISSING__"


            raw[
                !raw %in%
                    spec$levels
            ] <- "__OTHER__"


            mat <- matrix(
                0,
                nrow =
                    length(
                        raw
                    ),
                ncol =
                    length(
                        spec$levels
                    )
            )


            colnames(
                mat
            ) <- paste0(
                "cat__",
                col,
                "__",
                make.names(
                    spec$levels,
                    unique = TRUE
                )
            )


            for (
                j in seq_along(
                    spec$levels
                )
            ) {

                mat[
                    ,
                    j
                ] <- as.numeric(
                    raw ==
                        spec$levels[[j]]
                )
            }


            pieces[[length( pieces ) + 1L]] <- mat
        }
    }


    X <- do.call(
        cbind,
        pieces
    )


    X <- as.matrix(
        X
    )


    storage.mode(
        X
    ) <- "double"


    colnames(
        X
    ) <- make.names(
        colnames(
            X
        ),
        unique = TRUE
    )


    absent <- setdiff(
        prep$matrix_columns,
        colnames(
            X
        )
    )


    if (
        length(
            absent
        ) >
            0L
    ) {

        extra <- matrix(
            0,
            nrow =
                nrow(
                    X
                ),
            ncol =
                length(
                    absent
                )
        )


        colnames(
            extra
        ) <- absent


        X <- cbind(
            X,
            extra
        )
    }


    X[
        ,
        prep$matrix_columns,
        drop = FALSE
    ]
}


# -----------------------------------------------------------------------------
# 10. Map engineered matrix columns back to original predictors
# -----------------------------------------------------------------------------

source_feature_from_engineered <- function(
        engineered_name
) {

    parts <- strsplit(
        engineered_name,
        "__",
        fixed = TRUE
    )[[1]]


    if (
        length(
            parts
        ) >=
            2L &&
        parts[[1]] %in%
            c(
                "num",
                "missing",
                "cat"
            )
    ) {

        return(
            parts[[2]]
        )
    }


    engineered_name
}


source_feature_map <- function(
        engineered_columns
) {

    data.table(
        engineered_feature =
            engineered_columns,

        source_feature =
            vapply(
                engineered_columns,
                source_feature_from_engineered,
                character(
                    1
                )
            )
    )
}


collapse_matrix_to_source_features <- function(
        mat,
        source_map
) {

    source_features <- unique(
        source_map$
            source_feature
    )


    out <- matrix(
        0,
        nrow =
            nrow(
                mat
            ),
        ncol =
            length(
                source_features
            )
    )


    colnames(
        out
    ) <- source_features


    for (
        feature in source_features
    ) {

        cols <- source_map[
            source_feature ==
                feature,
            engineered_feature
        ]


        out[
            ,
            feature
        ] <- rowSums(
            mat[
                ,
                cols,
                drop = FALSE
            ]
        )
    }


    out
}


# -----------------------------------------------------------------------------
# 11. Base XGBoost and Ranger fitters
# -----------------------------------------------------------------------------

fit_xgboost_model <- function(
        X,
        y,
        seed
) {

    dtrain <- xgboost::xgb.DMatrix(
        data =
            X,

        label =
            y,

        nthread =
            1
    )


    params <- list(
        objective =
            "binary:logistic",

        eval_metric =
            "logloss",

        max_depth =
            XGB_MAX_DEPTH,

        eta =
            XGB_ETA,

        subsample =
            XGB_SUBSAMPLE,

        colsample_bytree =
            XGB_COLSAMPLE,

        min_child_weight =
            XGB_MIN_CHILD_WEIGHT,

        nthread =
            1,

        seed =
            seed
    )


    xgboost::xgb.train(
        params =
            params,

        data =
            dtrain,

        nrounds =
            XGB_NROUNDS,

        verbose =
            0
    )
}


fit_ranger_model <- function(
        X,
        y,
        seed
) {

    df <- as.data.frame(
        X
    )


    df$.y <- factor(
        y,
        levels = c(
            0,
            1
        )
    )


    ranger::ranger(
        dependent.variable.name =
            ".y",

        data =
            df,

        probability =
            TRUE,

        num.trees =
            RANGER_NUM_TREES,

        mtry =
            max(
                1L,
                floor(
                    sqrt(
                        ncol(
                            X
                        )
                    )
                )
            ),

        min.node.size =
            RANGER_MIN_NODE_SIZE,

        importance =
            "permutation",

        seed =
            seed,

        num.threads =
            1
    )
}


# -----------------------------------------------------------------------------
# 12. SHAP helper
# -----------------------------------------------------------------------------

xgb_shap_matrix <- function(
        model,
        X
) {

    d <- xgboost::xgb.DMatrix(
        data =
            X,

        nthread =
            1
    )


    shap <- predict(
        model,
        d,
        predcontrib =
            TRUE
    )


    shap <- as.matrix(
        shap
    )


    # Last column is the bias/intercept contribution.
    bias_index <- ncol(
        shap
    )


    if (
        bias_index <=
            1L
    ) {

        stop(
            "XGBoost SHAP output contains no feature columns.",
            call. = FALSE
        )
    }


    shap_feature <- shap[
        ,
        seq_len(
            bias_index -
                1L
        ),
        drop = FALSE
    ]


    colnames(
        shap_feature
    ) <- colnames(
        X
    )


    shap_feature
}


# -----------------------------------------------------------------------------
# 13. Accumulators
# -----------------------------------------------------------------------------

xgb_global_records <- list()

xgb_global_counter <- 0L


top20_local_records <- list()

top20_counter <- 0L


benchmark_local_records <- list()

benchmark_counter <- 0L


ranger_records <- list()

ranger_counter <- 0L


refit_records <- list()

refit_counter <- 0L


# -----------------------------------------------------------------------------
# 14. Explain selected XGBoost models
# -----------------------------------------------------------------------------

for (
    selected_row in seq_len(
        nrow(
            selected_xgb
        )
    )
) {

    model_set_name <-
        selected_xgb$
            model_set[[selected_row]]


    dt <- copy(
        datasets[[model_set_name]]
    )


    predictors <- predictor_sets[[model_set_name]]


    y <- as.integer(
        dt$pu_label
    )


    # Benchmark predictor coverage for this model set.
    benchmark_predictor_coverage <- rowSums(
        !is.na(
            benchmark_master[
                ,
                ..predictors
            ]
        )
    ) /
        length(
            predictors
        )


    for (
        repeat_id in seq_len(
            N_REPEATS
        )
    ) {

        repeat_seed <-
            MASTER_SEED +
            match(
                model_set_name,
                MODEL_SETS_STEP08
            ) *
            100000L +
            repeat_id *
            10000L


        fold_id <- make_stratified_folds(
            y =
                y,

            k =
                N_FOLDS,

            seed =
                repeat_seed
        )


        for (
            fold in seq_len(
                N_FOLDS
            )
        ) {

            train_idx <- which(
                fold_id !=
                    fold
            )


            test_idx <- which(
                fold_id ==
                    fold
            )


            train_dt <- dt[
                train_idx
            ]


            test_dt <- dt[
                test_idx
            ]


            prep <- fit_preprocessor(
                train_dt =
                    train_dt,

                predictors =
                    predictors
            )


            X_train <- prep$X


            X_test <- apply_preprocessor(
                dt =
                    test_dt,

                prep =
                    prep
            )


            X_benchmark <- apply_preprocessor(
                dt =
                    benchmark_master,

                prep =
                    prep
            )


            y_train <- as.integer(
                train_dt$pu_label
            )


            positive_idx <- which(
                y_train ==
                    1L
            )


            unlabeled_idx <- which(
                y_train ==
                    0L
            )


            n_unlabeled_bag <- min(
                length(
                    unlabeled_idx
                ),
                max(
                    length(
                        positive_idx
                    ),
                    as.integer(
                        round(
                            UNLABELED_TO_POSITIVE_RATIO *
                            length(
                                positive_idx
                            )
                        )
                    )
                )
            )


            source_map <- source_feature_map(
                colnames(
                    X_train
                )
            )


            global_abs_sum <- setNames(
                rep(
                    0,
                    length(
                        unique(
                            source_map$
                                source_feature
                        )
                    )
                ),
                unique(
                    source_map$
                        source_feature
                )
            )


            global_signed_sum <- global_abs_sum


            global_n <-
                0L


            top20_in_test <- which(
                test_dt$gene_id %in%
                    top20_gene_ids
            )


            for (
                bag in seq_len(
                    N_BAGS
                )
            ) {

                fit_seed <-
                    repeat_seed +
                    fold *
                    1000L +
                    match(
                        "xgboost",
                        ALGORITHMS_STEP08
                    ) *
                    100L


                bag_seed <-
                    fit_seed +
                    bag *
                    1009L


                set.seed(
                    bag_seed
                )


                sampled_unlabeled <- sample(
                    unlabeled_idx,
                    size =
                        n_unlabeled_bag,
                    replace =
                        FALSE
                )


                bag_idx <- c(
                    positive_idx,
                    sampled_unlabeled
                )


                bag_y <- c(
                    rep(
                        1L,
                        length(
                            positive_idx
                        )
                    ),
                    rep(
                        0L,
                        length(
                            sampled_unlabeled
                        )
                    )
                )


                model <- fit_xgboost_model(
                    X =
                        X_train[
                            bag_idx,
                            ,
                            drop = FALSE
                        ],

                    y =
                        bag_y,

                    seed =
                        bag_seed
                )


                # -------------------------------------------------------------
                # OOF global SHAP
                # -------------------------------------------------------------

                shap_test_engineered <-
                    xgb_shap_matrix(
                        model,
                        X_test
                    )


                shap_test <-
                    collapse_matrix_to_source_features(
                        shap_test_engineered,
                        source_map
                    )


                global_abs_sum[
                    colnames(
                        shap_test
                    )
                ] <-
                    global_abs_sum[
                        colnames(
                            shap_test
                        )
                    ] +
                    colSums(
                        abs(
                            shap_test
                        )
                    )


                global_signed_sum[
                    colnames(
                        shap_test
                    )
                ] <-
                    global_signed_sum[
                        colnames(
                            shap_test
                        )
                    ] +
                    colSums(
                        shap_test
                    )


                global_n <-
                    global_n +
                    nrow(
                        shap_test
                    )


                # -------------------------------------------------------------
                # OOF local SHAP for top-20 genes
                # -------------------------------------------------------------

                if (
                    length(
                        top20_in_test
                    ) >
                        0L
                ) {

                    local_mat <-
                        shap_test[
                            top20_in_test,
                            ,
                            drop = FALSE
                        ]


                    local_dt <- as.data.table(
                        local_mat
                    )


                    local_dt[
                        ,
                        gene_id :=
                            test_dt$
                                gene_id[
                                    top20_in_test
                                ]
                    ]


                    local_long <- melt(
                        local_dt,
                        id.vars =
                            "gene_id",
                        variable.name =
                            "feature",
                        value.name =
                            "shap_value"
                    )


                    local_long[
                        ,
                        `:=`(
                            model_set =
                                model_set_name,

                            algorithm =
                                "xgboost",

                            repeat_id =
                                repeat_id,

                            fold =
                                fold,

                            bag =
                                bag
                        )
                    ]


                    top20_counter <-
                        top20_counter +
                        1L


                    top20_local_records[[top20_counter]] <- local_long
                }


                # -------------------------------------------------------------
                # External benchmark SHAP
                # -------------------------------------------------------------

                shap_bench_engineered <-
                    xgb_shap_matrix(
                        model,
                        X_benchmark
                    )


                shap_bench <-
                    collapse_matrix_to_source_features(
                        shap_bench_engineered,
                        source_map
                    )


                # A benchmark with zero measured predictor coverage should not
                # receive an imputation-only explanation.
                zero_coverage <-
                    benchmark_predictor_coverage ==
                    0


                if (
                    any(
                        zero_coverage
                    )
                ) {

                    shap_bench[
                        zero_coverage,
                    ] <- NA_real_
                }


                bench_dt <- as.data.table(
                    shap_bench
                )


                bench_dt[
                    ,
                    gene_id :=
                        benchmark_master$
                            gene_id
                ]


                bench_long <- melt(
                    bench_dt,
                    id.vars =
                        "gene_id",
                    variable.name =
                        "feature",
                    value.name =
                        "shap_value"
                )


                bench_long[
                    ,
                    `:=`(
                        model_set =
                            model_set_name,

                        algorithm =
                            "xgboost",

                        repeat_id =
                            repeat_id,

                        fold =
                            fold,

                        bag =
                            bag
                    )
                ]


                benchmark_counter <-
                    benchmark_counter +
                    1L


                benchmark_local_records[[benchmark_counter]] <- bench_long
            }


            # Store fold/repeat aggregate global SHAP.
            xgb_global_counter <-
                xgb_global_counter +
                1L


            xgb_global_records[[xgb_global_counter]] <- data.table(
                model_set =
                    model_set_name,

                algorithm =
                    "xgboost",

                repeat_id =
                    repeat_id,

                fold =
                    fold,

                feature =
                    names(
                        global_abs_sum
                    ),

                sum_abs_shap =
                    as.numeric(
                        global_abs_sum
                    ),

                sum_signed_shap =
                    as.numeric(
                        global_signed_sum
                    ),

                n_gene_model_explanations =
                    global_n
            )


            refit_counter <-
                refit_counter +
                1L


            refit_records[[refit_counter]] <- data.table(
                model_set =
                    model_set_name,

                algorithm =
                    "xgboost",

                repeat_id =
                    repeat_id,

                fold =
                    fold,

                n_training_genes =
                    length(
                        train_idx
                    ),

                n_test_genes =
                    length(
                        test_idx
                    ),

                n_observed_positive_training =
                    length(
                        positive_idx
                    ),

                n_unlabeled_per_bag =
                    n_unlabeled_bag,

                n_bags =
                    N_BAGS,

                engineered_predictors =
                    ncol(
                        X_train
                    )
            )
        }
    }
}


# -----------------------------------------------------------------------------
# 15. Explain selected Ranger models with permutation importance
# -----------------------------------------------------------------------------

for (
    selected_row in seq_len(
        nrow(
            selected_ranger
        )
    )
) {

    model_set_name <-
        selected_ranger$
            model_set[[selected_row]]


    dt <- copy(
        datasets[[model_set_name]]
    )


    predictors <- predictor_sets[[model_set_name]]


    y <- as.integer(
        dt$pu_label
    )


    for (
        repeat_id in seq_len(
            N_REPEATS
        )
    ) {

        repeat_seed <-
            MASTER_SEED +
            match(
                model_set_name,
                MODEL_SETS_STEP08
            ) *
            100000L +
            repeat_id *
            10000L


        fold_id <- make_stratified_folds(
            y =
                y,

            k =
                N_FOLDS,

            seed =
                repeat_seed
        )


        for (
            fold in seq_len(
                N_FOLDS
            )
        ) {

            train_idx <- which(
                fold_id !=
                    fold
            )


            test_idx <- which(
                fold_id ==
                    fold
            )


            train_dt <- dt[
                train_idx
            ]


            prep <- fit_preprocessor(
                train_dt =
                    train_dt,

                predictors =
                    predictors
            )


            X_train <- prep$X


            y_train <- as.integer(
                train_dt$pu_label
            )


            positive_idx <- which(
                y_train ==
                    1L
            )


            unlabeled_idx <- which(
                y_train ==
                    0L
            )


            n_unlabeled_bag <- min(
                length(
                    unlabeled_idx
                ),
                max(
                    length(
                        positive_idx
                    ),
                    as.integer(
                        round(
                            UNLABELED_TO_POSITIVE_RATIO *
                            length(
                                positive_idx
                            )
                        )
                    )
                )
            )


            source_map <- source_feature_map(
                colnames(
                    X_train
                )
            )


            for (
                bag in seq_len(
                    N_BAGS
                )
            ) {

                fit_seed <-
                    repeat_seed +
                    fold *
                    1000L +
                    match(
                        "ranger",
                        ALGORITHMS_STEP08
                    ) *
                    100L


                bag_seed <-
                    fit_seed +
                    bag *
                    1009L


                set.seed(
                    bag_seed
                )


                sampled_unlabeled <- sample(
                    unlabeled_idx,
                    size =
                        n_unlabeled_bag,
                    replace =
                        FALSE
                )


                bag_idx <- c(
                    positive_idx,
                    sampled_unlabeled
                )


                bag_y <- c(
                    rep(
                        1L,
                        length(
                            positive_idx
                        )
                    ),
                    rep(
                        0L,
                        length(
                            sampled_unlabeled
                        )
                    )
                )


                fit <- fit_ranger_model(
                    X =
                        X_train[
                            bag_idx,
                            ,
                            drop = FALSE
                        ],

                    y =
                        bag_y,

                    seed =
                        bag_seed
                )


                engineered_importance <-
                    fit$
                        variable.importance


                if (
                    is.null(
                        engineered_importance
                    )
                ) {

                    stop(
                        "Ranger returned no permutation importance.",
                        call. = FALSE
                    )
                }


                importance_dt <- data.table(
                    engineered_feature =
                        names(
                            engineered_importance
                        ),

                    permutation_importance =
                        as.numeric(
                            engineered_importance
                        )
                )


                importance_dt <- merge(
                    importance_dt,
                    source_map,
                    by =
                        "engineered_feature",
                    all.x =
                        TRUE,
                    sort =
                        FALSE
                )


                collapsed <- importance_dt[
                    ,
                    .(
                        permutation_importance =
                            sum(
                                permutation_importance,
                                na.rm = TRUE
                            )
                    ),
                    by =
                        source_feature
                ]


                collapsed[
                    ,
                    `:=`(
                        model_set =
                            model_set_name,

                        algorithm =
                            "ranger",

                        repeat_id =
                            repeat_id,

                        fold =
                            fold,

                        bag =
                            bag
                    )
                ]


                ranger_counter <-
                    ranger_counter +
                    1L


                ranger_records[[ranger_counter]] <- collapsed
            }


            refit_counter <-
                refit_counter +
                1L


            refit_records[[refit_counter]] <- data.table(
                model_set =
                    model_set_name,

                algorithm =
                    "ranger",

                repeat_id =
                    repeat_id,

                fold =
                    fold,

                n_training_genes =
                    length(
                        train_idx
                    ),

                n_test_genes =
                    length(
                        test_idx
                    ),

                n_observed_positive_training =
                    length(
                        positive_idx
                    ),

                n_unlabeled_per_bag =
                    n_unlabeled_bag,

                n_bags =
                    N_BAGS,

                engineered_predictors =
                    ncol(
                        X_train
                    )
            )
        }
    }
}


# -----------------------------------------------------------------------------
# 16. Aggregate XGBoost global SHAP
# -----------------------------------------------------------------------------

xgb_global_raw <- if (
    length(
        xgb_global_records
    ) >
        0L
) {

    rbindlist(
        xgb_global_records,
        fill = TRUE,
        use.names = TRUE
    )

} else {

    data.table()
}


if (
    nrow(
        xgb_global_raw
    ) >
        0L
) {

    xgb_global <- xgb_global_raw[
        ,
        .(
            mean_abs_shap =
                sum(
                    sum_abs_shap,
                    na.rm = TRUE
                ) /
                sum(
                    n_gene_model_explanations,
                    na.rm = TRUE
                ),

            mean_signed_shap =
                sum(
                    sum_signed_shap,
                    na.rm = TRUE
                ) /
                sum(
                    n_gene_model_explanations,
                    na.rm = TRUE
                ),

            total_gene_model_explanations =
                sum(
                    n_gene_model_explanations,
                    na.rm = TRUE
                )
        ),
        by = .(
            model_set,
            algorithm,
            feature
        )
    ]


    xgb_global[
        ,
        importance_percentile :=
            rank01(
                mean_abs_shap
            ),
        by = .(
            model_set,
            algorithm
        )
    ]


    xgb_global[
        ,
        feature_rank :=
            frank(
                -mean_abs_shap,
                ties.method =
                    "min"
            ),
        by = .(
            model_set,
            algorithm
        )
    ]


    data.table::setorder(
        xgb_global,
        model_set,
        feature_rank,
        feature
    )

} else {

    xgb_global <- data.table()
}


# -----------------------------------------------------------------------------
# 17. Aggregate local XGBoost SHAP
# -----------------------------------------------------------------------------

top20_local_raw <- if (
    length(
        top20_local_records
    ) >
        0L
) {

    rbindlist(
        top20_local_records,
        fill = TRUE,
        use.names = TRUE
    )

} else {

    data.table()
}


if (
    nrow(
        top20_local_raw
    ) >
        0L
) {

    top20_local <- top20_local_raw[
        ,
        .(
            mean_shap =
                mean(
                    shap_value,
                    na.rm = TRUE
                ),

            mean_abs_shap =
                mean(
                    abs(
                        shap_value
                    ),
                    na.rm = TRUE
                ),

            sd_shap =
                stats::sd(
                    shap_value,
                    na.rm = TRUE
                ),

            n_oof_model_explanations =
                sum(
                    is.finite(
                        shap_value
                    )
                )
        ),
        by = .(
            model_set,
            algorithm,
            gene_id,
            feature
        )
    ]


    top20_local[
        ,
        local_abs_importance_rank :=
            frank(
                -mean_abs_shap,
                ties.method =
                    "min"
            ),
        by = .(
            model_set,
            algorithm,
            gene_id
        )
    ]


    top20_local <- merge(
        top20_local,
        top20[
            ,
            .(
                gene_id,
                final_prepopulation_rank,
                novel_candidate_rank =
                    if (
                        "novel_candidate_rank" %in%
                            names(
                                top20
                            )
                    ) {
                        novel_candidate_rank
                    } else {
                        NA_integer_
                    }
            )
        ],
        by =
            "gene_id",
        all.x =
            TRUE,
        sort =
            FALSE
    )

} else {

    top20_local <- data.table()
}


benchmark_local_raw <- if (
    length(
        benchmark_local_records
    ) >
        0L
) {

    rbindlist(
        benchmark_local_records,
        fill = TRUE,
        use.names = TRUE
    )

} else {

    data.table()
}


if (
    nrow(
        benchmark_local_raw
    ) >
        0L
) {

    benchmark_local <- benchmark_local_raw[
        ,
        .(
            mean_shap =
                mean(
                    shap_value,
                    na.rm = TRUE
                ),

            mean_abs_shap =
                mean(
                    abs(
                        shap_value
                    ),
                    na.rm = TRUE
                ),

            sd_shap =
                stats::sd(
                    shap_value,
                    na.rm = TRUE
                ),

            n_external_model_explanations =
                sum(
                    is.finite(
                        shap_value
                    )
                )
        ),
        by = .(
            model_set,
            algorithm,
            gene_id,
            feature
        )
    ]


    benchmark_local[
        ,
        local_abs_importance_rank :=
            frank(
                -mean_abs_shap,
                ties.method =
                    "min"
            ),
        by = .(
            model_set,
            algorithm,
            gene_id
        )
    ]


    benchmark_local <- merge(
        benchmark_local,
        benchmarks[
            ,
            .(
                gene_id,
                benchmark_name =
                    if (
                        "benchmark_name" %in%
                            names(
                                benchmarks
                            )
                    ) {
                        benchmark_name
                    } else {
                        NA_character_
                    },

                final_prepopulation_rank
            )
        ],
        by =
            "gene_id",
        all.x =
            TRUE,
        sort =
            FALSE
    )

} else {

    benchmark_local <- data.table()
}


# -----------------------------------------------------------------------------
# 18. Aggregate Ranger permutation importance
# -----------------------------------------------------------------------------

ranger_raw <- if (
    length(
        ranger_records
    ) >
        0L
) {

    rbindlist(
        ranger_records,
        fill = TRUE,
        use.names = TRUE
    )

} else {

    data.table()
}


if (
    nrow(
        ranger_raw
    ) >
        0L
) {

    ranger_importance <- ranger_raw[
        ,
        .(
            mean_permutation_importance =
                mean(
                    permutation_importance,
                    na.rm = TRUE
                ),

            sd_permutation_importance =
                stats::sd(
                    permutation_importance,
                    na.rm = TRUE
                ),

            n_bag_fold_models =
                .N
        ),
        by = .(
            model_set,
            algorithm,
            feature =
                source_feature
        )
    ]


    ranger_importance[
        ,
        importance_percentile :=
            rank01(
                mean_permutation_importance
            ),
        by = .(
            model_set,
            algorithm
        )
    ]


    ranger_importance[
        ,
        feature_rank :=
            frank(
                -mean_permutation_importance,
                ties.method =
                    "min"
            ),
        by = .(
            model_set,
            algorithm
        )
    ]


    data.table::setorder(
        ranger_importance,
        model_set,
        feature_rank,
        feature
    )

} else {

    ranger_importance <- data.table()
}


# -----------------------------------------------------------------------------
# 19. Cross-model consensus global importance
# -----------------------------------------------------------------------------

consensus_parts <- list()


part_counter <- 0L


if (
    nrow(
        xgb_global
    ) >
        0L
) {

    for (
        model_name in unique(
            xgb_global$
                model_set
        )
    ) {

        part_counter <-
            part_counter +
            1L


        z <- xgb_global[
            model_set ==
                model_name
        ]


        consensus_parts[[part_counter]] <- z[
            ,
            .(
                model_key =
                    paste(
                        model_set,
                        algorithm,
                        sep = "::"
                    ),

                feature,
                importance_percentile,
                raw_importance =
                    mean_abs_shap,

                importance_type =
                    "mean_abs_crossfitted_SHAP"
            )
        ]
    }
}


if (
    nrow(
        ranger_importance
    ) >
        0L
) {

    for (
        model_name in unique(
            ranger_importance$
                model_set
        )
    ) {

        part_counter <-
            part_counter +
            1L


        z <- ranger_importance[
            model_set ==
                model_name
        ]


        consensus_parts[[part_counter]] <- z[
            ,
            .(
                model_key =
                    paste(
                        model_set,
                        algorithm,
                        sep = "::"
                    ),

                feature,
                importance_percentile,
                raw_importance =
                    mean_permutation_importance,

                importance_type =
                    "crossfitted_bagged_permutation"
            )
        ]
    }
}


importance_long <- rbindlist(
    consensus_parts,
    fill = TRUE,
    use.names = TRUE
)


consensus_importance <- importance_long[
    ,
    .(
        consensus_importance_percentile =
            median(
                importance_percentile,
                na.rm = TRUE
            ),

        mean_importance_percentile =
            mean(
                importance_percentile,
                na.rm = TRUE
            ),

        min_importance_percentile =
            min(
                importance_percentile,
                na.rm = TRUE
            ),

        max_importance_percentile =
            max(
                importance_percentile,
                na.rm = TRUE
            ),

        n_selected_models_with_feature =
            uniqueN(
                model_key
            )
    ),
    by =
        feature
]


data.table::setorder(
    consensus_importance,
    -consensus_importance_percentile,
    -mean_importance_percentile,
    feature
)


consensus_importance[
    ,
    consensus_feature_rank :=
        seq_len(
            .N
        )
]


# -----------------------------------------------------------------------------
# 20. Refit summary
# -----------------------------------------------------------------------------

refit_summary <- rbindlist(
    refit_records,
    fill = TRUE,
    use.names = TRUE
)


# -----------------------------------------------------------------------------
# 21. Publication figures
# -----------------------------------------------------------------------------

# ---- Figure 11A: consensus global importance ----

plot_consensus <- consensus_importance[
    consensus_feature_rank <=
        min(
            15L,
            .N
        )
]


plot_consensus[
    ,
    feature :=
        factor(
            feature,
            levels =
                rev(
                    feature
                )
        )
]


p11a <- ggplot(
    plot_consensus,
    aes(
        x =
            consensus_importance_percentile,
        y =
            feature
    )
) +
    geom_col() +
    geom_point(
        aes(
            x =
                mean_importance_percentile
        )
    ) +
    labs(
        x =
            "Consensus importance percentile",

        y =
            NULL,

        title =
            "Cross-model consensus importance of leakage-safe PU predictors",

        subtitle =
            "Median rank-normalized importance across selected XGBoost and Ranger models"
    ) +
    theme_classic(
        base_size =
            11
    )


save_plot_pair(
    p11a,
    FIG11A_PNG,
    FIG11A_PDF,
    width =
        8.2,
    height =
        6.5
)


# ---- Figure 11B: M3-XGBoost global SHAP ----

m3_xgb <- xgb_global[
    model_set ==
        "M3_expression_orthology_annotation" &
        algorithm ==
        "xgboost"
]


if (
    nrow(
        m3_xgb
    ) >
        0L
) {

    m3_xgb_plot <- m3_xgb[
        feature_rank <=
            min(
                15L,
                .N
            )
    ]


    m3_xgb_plot[
        ,
        feature :=
            factor(
                feature,
                levels =
                    rev(
                        feature
                    )
            )
    ]


    p11b <- ggplot(
        m3_xgb_plot,
        aes(
            x =
                mean_abs_shap,
            y =
                feature
        )
    ) +
        geom_col() +
        labs(
            x =
                "Mean absolute cross-fitted SHAP",

            y =
                NULL,

            title =
                "Global feature importance in the selected M3 XGBoost PU model",

            subtitle =
                "SHAP values aggregated only from held-out genes"
        ) +
        theme_classic(
            base_size =
                11
        )


    save_plot_pair(
        p11b,
        FIG11B_PNG,
        FIG11B_PDF,
        width =
            8.2,
        height =
            6.5
    )

} else {

    warning(
        "Selected models do not contain M3 XGBoost; Figure 11B not written.",
        call. = FALSE
    )
}


# ---- Local heatmap feature set from top M3 XGBoost global features ----

local_features <- if (
    nrow(
        m3_xgb
    ) >
        0L
) {

    m3_xgb[
        order(
            feature_rank
        )
    ][
        1:min(
            8L,
            .N
        ),
        feature
    ]

} else {

    consensus_importance[
        1:min(
            8L,
            .N
        ),
        feature
    ]
}


# ---- Figure 11C: top-20 local SHAP heatmap ----

local_top20_plot <- top20_local[
    model_set ==
        "M3_expression_orthology_annotation" &
        algorithm ==
        "xgboost" &
        feature %in%
        local_features
]


if (
    nrow(
        local_top20_plot
    ) >
        0L
) {

    local_top20_plot <- merge(
        local_top20_plot,
        top20[
            ,
            .(
                gene_id,
                final_prepopulation_rank
            )
        ],
        by =
            "gene_id",
        all.x =
            TRUE,
        sort =
            FALSE,
        suffixes =
            c(
                "",
                ".top20"
            )
    )


    if (
        "final_prepopulation_rank.top20" %in%
            names(
                local_top20_plot
            )
    ) {

        local_top20_plot[
            ,
            final_prepopulation_rank :=
                final_prepopulation_rank.top20
        ]


        local_top20_plot[
            ,
            final_prepopulation_rank.top20 :=
                NULL
        ]
    }


    gene_order <- top20[
        order(
            final_prepopulation_rank
        ),
        gene_id
    ]


    local_top20_plot[
        ,
        gene_id :=
            factor(
                gene_id,
                levels =
                    rev(
                        gene_order
                    )
            )
    ]


    local_top20_plot[
        ,
        feature :=
            factor(
                feature,
                levels =
                    local_features
            )
    ]


    p11c <- ggplot(
        local_top20_plot,
        aes(
            x =
                feature,
            y =
                gene_id,
            fill =
                mean_shap
        )
    ) +
        geom_tile() +
        geom_text(
            aes(
                label =
                    sprintf(
                        "%.2f",
                        mean_shap
                    )
            ),
            size =
                2.8
        ) +
        scale_fill_gradient2(
            midpoint =
                0
        ) +
        labs(
            x =
                NULL,

            y =
                NULL,

            fill =
                "Mean SHAP",

            title =
                "Cross-fitted local SHAP profiles of the top 20 novel candidates",

            subtitle =
                "Selected M3 XGBoost PU model; positive values support the observed-positive class"
        ) +
        theme_minimal(
            base_size =
                10.5
        ) +
        theme(
            panel.grid =
                element_blank(),

            axis.text.x =
                element_text(
                    angle =
                        25,
                    hjust =
                        1
                )
        )


    save_plot_pair(
        p11c,
        FIG11C_PNG,
        FIG11C_PDF,
        width =
            9.5,
        height =
            8.5
    )

} else {

    warning(
        "No M3 XGBoost local SHAP records were available for the top-20 genes.",
        call. = FALSE
    )
}


# ---- Figure 11D: benchmark local SHAP heatmap ----

local_bench_plot <- benchmark_local[
    model_set ==
        "M3_expression_orthology_annotation" &
        algorithm ==
        "xgboost" &
        feature %in%
        local_features
]


if (
    nrow(
        local_bench_plot
    ) >
        0L
) {

    benchmark_labels <- benchmarks[
        ,
        .(
            gene_id,
            display_label =
                if (
                    "benchmark_name" %in%
                        names(
                            benchmarks
                        )
                ) {
                    paste0(
                        gene_id,
                        "\n",
                        benchmark_name
                    )
                } else {
                    gene_id
                }
        )
    ]


    local_bench_plot <- merge(
        local_bench_plot,
        benchmark_labels,
        by =
            "gene_id",
        all.x =
            TRUE,
        sort =
            FALSE
    )


    bench_order <- benchmark_labels$
        display_label


    local_bench_plot[
        ,
        display_label :=
            factor(
                display_label,
                levels =
                    rev(
                        bench_order
                    )
            )
    ]


    local_bench_plot[
        ,
        feature :=
            factor(
                feature,
                levels =
                    local_features
            )
    ]


    p11d <- ggplot(
        local_bench_plot,
        aes(
            x =
                feature,
            y =
                display_label,
            fill =
                mean_shap
        )
    ) +
        geom_tile() +
        geom_text(
            aes(
                label =
                    ifelse(
                        is.finite(
                            mean_shap
                        ),
                        sprintf(
                            "%.2f",
                            mean_shap
                        ),
                        "NA"
                    )
            ),
            size =
                2.8
        ) +
        scale_fill_gradient2(
            midpoint =
                0,
            na.value =
                "grey90"
        ) +
        labs(
            x =
                NULL,

            y =
                NULL,

            fill =
                "Mean SHAP",

            title =
                "External benchmark SHAP profiles",

            subtitle =
                "M3 XGBoost models never trained on benchmark genes"
        ) +
        theme_minimal(
            base_size =
                10.5
        ) +
        theme(
            panel.grid =
                element_blank(),

            axis.text.x =
                element_text(
                    angle =
                        25,
                    hjust =
                        1
                )
        )


    save_plot_pair(
        p11d,
        FIG11D_PNG,
        FIG11D_PDF,
        width =
            10,
        height =
            6.5
    )

} else {

    warning(
        "No M3 XGBoost benchmark SHAP records were available.",
        call. = FALSE
    )
}


# -----------------------------------------------------------------------------
# 22. QC
# -----------------------------------------------------------------------------

expected_xgb_fold_fits <-
    nrow(
        selected_xgb
    ) *
    N_REPEATS *
    N_FOLDS


expected_ranger_fold_fits <-
    nrow(
        selected_ranger
    ) *
    N_REPEATS *
    N_FOLDS


completed_xgb_fold_fits <-
    refit_summary[
        algorithm ==
            "xgboost",
        .N
    ]


completed_ranger_fold_fits <-
    refit_summary[
        algorithm ==
            "ranger",
        .N
    ]


if (
    completed_xgb_fold_fits !=
        expected_xgb_fold_fits
) {

    stop(
        paste0(
            "Expected ",
            expected_xgb_fold_fits,
            " XGBoost fold refits but completed ",
            completed_xgb_fold_fits,
            "."
        ),
        call. = FALSE
    )
}


if (
    completed_ranger_fold_fits !=
        expected_ranger_fold_fits
) {

    stop(
        paste0(
            "Expected ",
            expected_ranger_fold_fits,
            " Ranger fold refits but completed ",
            completed_ranger_fold_fits,
            "."
        ),
        call. = FALSE
    )
}


top20_explained_m3 <- if (
    nrow(
        top20_local
    ) >
        0L
) {

    uniqueN(
        top20_local[
            model_set ==
                "M3_expression_orthology_annotation" &
                algorithm ==
                "xgboost",
            gene_id
        ]
    )

} else {

    0L
}


bench_explained_m3 <- if (
    nrow(
        benchmark_local
    ) >
        0L
) {

    uniqueN(
        benchmark_local[
            model_set ==
                "M3_expression_orthology_annotation" &
                algorithm ==
                "xgboost" &
                n_external_model_explanations >
                0,
            gene_id
        ]
    )

} else {

    0L
}


qc <- data.table(
    metric = c(
        "Selected PU models",
        "Selected XGBoost models",
        "Selected Ranger models",
        "Expected XGBoost fold refits",
        "Completed XGBoost fold refits",
        "Expected Ranger fold refits",
        "Completed Ranger fold refits",
        "XGBoost global feature rows",
        "Ranger global feature rows",
        "Consensus feature rows",
        "Top-20 genes with M3 XGBoost cross-fitted SHAP",
        "Benchmarks with M3 XGBoost external SHAP",
        "Benchmark genes used in fitting",
        "FlyBase phenotype predictors added",
        "Step 04 ranking predictors added",
        "Step 05 CRISPR predictors added",
        "Population-genomics predictors added",
        "Figure 11A written",
        "Figure 11B written",
        "Figure 11C written",
        "Figure 11D written"
    ),

    value = c(
        nrow(
            selected_models
        ),

        nrow(
            selected_xgb
        ),

        nrow(
            selected_ranger
        ),

        expected_xgb_fold_fits,

        completed_xgb_fold_fits,

        expected_ranger_fold_fits,

        completed_ranger_fold_fits,

        nrow(
            xgb_global
        ),

        nrow(
            ranger_importance
        ),

        nrow(
            consensus_importance
        ),

        top20_explained_m3,

        bench_explained_m3,

        0,
        0,
        0,
        0,
        0,

        as.integer(
            file.exists(
                FIG11A_PNG
            ) &&
            file.exists(
                FIG11A_PDF
            )
        ),

        as.integer(
            file.exists(
                FIG11B_PNG
            ) &&
            file.exists(
                FIG11B_PDF
            )
        ),

        as.integer(
            file.exists(
                FIG11C_PNG
            ) &&
            file.exists(
                FIG11C_PDF
            )
        ),

        as.integer(
            file.exists(
                FIG11D_PNG
            ) &&
            file.exists(
                FIG11D_PDF
            )
        )
    )
)


# -----------------------------------------------------------------------------
# 23. Provenance
# -----------------------------------------------------------------------------

selected_model_text <- paste(
    selected_models$
        model_key,
    collapse =
        "; "
)


provenance <- data.table(
    component = c(
        "Selected models",
        "Refit design",
        "Fold assignment",
        "PU bagging",
        "XGBoost SHAP",
        "Global SHAP",
        "Local top-20 SHAP",
        "Benchmark SHAP",
        "Ranger permutation importance",
        "Engineered feature collapse",
        "Cross-model consensus importance",
        "Leakage control",
        "Population genomics"
    ),

    specification = c(
        selected_model_text,

        paste0(
            "Selected PU models are reconstructed using Step 07 modelling sets ",
            "and Step 08 hyperparameters."
        ),

        paste0(
            N_FOLDS,
            "-fold stratified folds repeated ",
            N_REPEATS,
            " times using the same Step 08 seed formula."
        ),

        paste0(
            N_BAGS,
            " bags per fold; each bag contains all observed positives plus ",
            "a fresh random unlabeled subset at ratio ",
            UNLABELED_TO_POSITIVE_RATIO,
            ":1."
        ),

        paste0(
            "XGBoost predcontrib values are computed from held-out fold genes ",
            "and external benchmark genes."
        ),

        paste0(
            "Mean absolute SHAP is aggregated across held-out genes, bags, folds, ",
            "and repeats."
        ),

        paste0(
            "Top-20 local SHAP values are cross-fitted: the explained gene is ",
            "never in the training subset of the model producing that SHAP value."
        ),

        paste0(
            "Six benchmark genes are excluded from all Step 07/08 fitting and are ",
            "explained externally across all refitted fold ensembles."
        ),

        paste0(
            "Ranger probability forests use permutation importance within each ",
            "PU bag; importances are averaged across bags, folds, and repeats."
        ),

        paste0(
            "Numeric values, missingness indicators, and one-hot factor levels are ",
            "collapsed back to their original biological predictor names."
        ),

        paste0(
            "Within-model global importance values are converted to percentiles; ",
            "consensus importance is the median percentile across selected models."
        ),

        paste0(
            "No FlyBase phenotype, Step 04 rank/score, Step 05 CRISPR, benchmark ",
            "identity, or Ag1000G variables are added to PU predictors."
        ),

        "Ag1000G remains pending and is not used in Step 11."
    )
)


# -----------------------------------------------------------------------------
# 24. Write outputs
# -----------------------------------------------------------------------------

fwrite(
    xgb_global,
    XGB_GLOBAL_FILE,
    na =
        "NA"
)


fwrite(
    top20_local,
    XGB_TOP20_LOCAL_FILE,
    na =
        "NA"
)


fwrite(
    benchmark_local,
    XGB_BENCH_LOCAL_FILE,
    na =
        "NA"
)


fwrite(
    ranger_importance,
    RANGER_IMPORTANCE_FILE,
    na =
        "NA"
)


fwrite(
    consensus_importance,
    CONSENSUS_IMPORTANCE_FILE,
    na =
        "NA"
)


fwrite(
    refit_summary,
    REFIT_SUMMARY_FILE,
    na =
        "NA"
)


fwrite(
    qc,
    QC_FILE,
    na =
        "NA"
)


fwrite(
    provenance,
    PROVENANCE_FILE,
    na =
        "NA"
)


capture.output(
    sessionInfo(),
    file =
        SESSION_INFO_FILE
)


checksum_files <- c(
    M2_FILE,
    M3_FILE,
    PREDICTOR_MANIFEST_FILE,
    MASTER_FILE,
    STEP08_CONFIG_FILE,
    SELECTED_MODELS_FILE,
    TOP20_FILE,
    BENCHMARK_FILE,
    XGB_GLOBAL_FILE,
    XGB_TOP20_LOCAL_FILE,
    XGB_BENCH_LOCAL_FILE,
    RANGER_IMPORTANCE_FILE,
    CONSENSUS_IMPORTANCE_FILE,
    REFIT_SUMMARY_FILE,
    QC_FILE,
    PROVENANCE_FILE,
    SESSION_INFO_FILE,
    FIG11A_PNG,
    FIG11A_PDF,
    FIG11B_PNG,
    FIG11B_PDF,
    FIG11C_PNG,
    FIG11C_PDF,
    FIG11D_PNG,
    FIG11D_PDF
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
# 25. Console summary
# -----------------------------------------------------------------------------

cat(
    "\n",
    "============================================================\n",
    "MOSQEDIT-R STEP 11 COMPLETED SUCCESSFULLY\n",
    "Cross-fitted model-specific explainability\n",
    "============================================================\n",
    "Selected PU models:                       ",
    nrow(
        selected_models
    ),
    "\n",
    "Selected XGBoost models:                  ",
    nrow(
        selected_xgb
    ),
    "\n",
    "Selected Ranger models:                   ",
    nrow(
        selected_ranger
    ),
    "\n",
    "XGBoost fold refits:                      ",
    completed_xgb_fold_fits,
    "\n",
    "Ranger fold refits:                       ",
    completed_ranger_fold_fits,
    "\n",
    "Top-20 genes with cross-fitted M3 SHAP:   ",
    top20_explained_m3,
    "\n",
    "Benchmarks with external M3 SHAP:          ",
    bench_explained_m3,
    "\n",
    "Consensus predictor features:              ",
    nrow(
        consensus_importance
    ),
    "\n",
    "Benchmark genes used for fitting:          NO\n",
    "FlyBase phenotype predictors added:        NO\n",
    "Step 04 score/rank predictors added:       NO\n",
    "Step 05 CRISPR predictors added:           NO\n",
    "Ag1000G predictors added:                  NO\n",
    "Publication figures written:               4 PNG + 4 PDF\n",
    "Python used:                               NO\n",
    "============================================================\n",
    sep =
        ""
)


cat(
    "\nTop global consensus features:\n"
)


print(
    consensus_importance[
        1:min(
            15L,
            .N
        )
    ]
)


cat(
    "\nM3 XGBoost global SHAP importance:\n"
)


if (
    nrow(
        m3_xgb
    ) >
        0L
) {

    print(
        m3_xgb[
            order(
                feature_rank
            )
        ][
            1:min(
                15L,
                .N
            )
        ]
    )
}


cat(
    "\nTop local M3 XGBoost SHAP features for the top-ranked novel gene:\n"
)


if (
    nrow(
        top20_local
    ) >
        0L
) {

    top_gene <- top20[
        order(
            final_prepopulation_rank
        ),
        gene_id
    ][[1]]


    print(
        top20_local[
            model_set ==
                "M3_expression_orthology_annotation" &
                algorithm ==
                "xgboost" &
                gene_id ==
                top_gene
        ][
            order(
                local_abs_importance_rank
            )
        ][
            1:min(
                10L,
                .N
            )
        ]
    )
}


cat(
    "\nQC summary:\n"
)


print(
    qc
)


log_step(
    "11",
    "Cross-fitted model-specific explainability completed successfully"
)

