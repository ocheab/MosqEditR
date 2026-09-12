# =============================================================================
# 08_train_bagged_pu_models.R
# MosqEdit-R Manuscript 1
# STEP 08 â€” repeated cross-fitted bagged positiveâ€“unlabeled learning
# =============================================================================

source("R/helpers.R")

required_packages <- c("data.table", "glmnet", "ranger", "xgboost", "e1071")
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

suppressPackageStartupMessages(library(data.table))

log_step("08", "Starting repeated cross-fitted bagged PU modelling")

# -----------------------------------------------------------------------------
# 1. Files and reproducibility
# -----------------------------------------------------------------------------

M1_FILE <- "data_processed/07_pu_M1_expression.csv"
M2_FILE <- "data_processed/07_pu_M2_expression_orthology.csv"
M3_FILE <- "data_processed/07_pu_M3_expression_orthology_annotation.csv"
PREDICTOR_MANIFEST_FILE <- "data_processed/07_predictor_manifest.csv"
MASTER_FILE <- "data_processed/master_gene_feature_matrix.csv"

CV_PREDICTIONS_FILE <- "data_processed/08_pu_cv_predictions.csv"
CV_METRICS_FILE <- "data_processed/08_pu_cv_metrics.csv"
METRICS_SUMMARY_FILE <- "data_processed/08_pu_metrics_summary.csv"
GENE_OOF_FILE <- "data_processed/08_pu_gene_oof_scores.csv"
BENCHMARK_SCORES_FILE <- "data_processed/08_external_benchmark_scores.csv"
CONFIG_FILE <- "data_processed/08_model_configuration.csv"
QC_FILE <- "data_processed/08_modelling_qc.csv"
PROVENANCE_FILE <- "data_processed/08_modelling_provenance.csv"
SESSION_INFO_FILE <- "logs/08_modelling_sessionInfo.txt"
CHECKSUM_FILE <- "logs/08_checksums.tsv"

dir.create("data_processed", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

MASTER_SEED <- 20260911L
N_FOLDS <- 5L
N_REPEATS <- 3L
N_BAGS <- 10L
UNLABELED_TO_POSITIVE_RATIO <- 1.0

GLMNET_ALPHA <- 0.50
GLMNET_NFOLDS <- 5L

RANGER_NUM_TREES <- 400L
RANGER_MIN_NODE_SIZE <- 5L

XGB_NROUNDS <- 150L
XGB_MAX_DEPTH <- 3L
XGB_ETA <- 0.05
XGB_SUBSAMPLE <- 0.80
XGB_COLSAMPLE <- 0.80
XGB_MIN_CHILD_WEIGHT <- 2

SVM_COST <- 1.0

MODEL_SETS <- c(
  "M1_expression",
  "M2_expression_orthology",
  "M3_expression_orthology_annotation"
)

ALGORITHMS <- c(
  "elastic_net",
  "ranger",
  "xgboost",
  "svm_radial"
)

# -----------------------------------------------------------------------------
# 2. Basic helpers
# -----------------------------------------------------------------------------

assert_file <- function(path) {
  if (!file.exists(path)) {
    stop(paste0("Required input file is missing:\n", path), call. = FALSE)
  }
  invisible(TRUE)
}

assert_unique_gene_ids <- function(dt, label) {
  if (!"gene_id" %in% names(dt)) {
    stop(paste0(label, " lacks gene_id."), call. = FALSE)
  }
  if (uniqueN(dt$gene_id) != nrow(dt)) {
    stop(paste0(label, " contains duplicate gene_id values."), call. = FALSE)
  }
  invisible(TRUE)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

# -----------------------------------------------------------------------------
# 3. PU proxy metrics
# -----------------------------------------------------------------------------
# IMPORTANT: labels are observed-positive versus unlabeled, not true negative.

roc_auc_rank <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])

  n_pos <- sum(y == 1L)
  n_neg <- sum(y == 0L)

  if (n_pos == 0L || n_neg == 0L) return(NA_real_)

  r <- rank(score, ties.method = "average")

  (
    sum(r[y == 1L]) -
      n_pos * (n_pos + 1) / 2
  ) / (n_pos * n_neg)
}

pr_auc_trapezoid <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])

  n_pos <- sum(y == 1L)
  if (n_pos == 0L) return(NA_real_)

  ord <- order(score, decreasing = TRUE)
  y <- y[ord]

  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)

  recall <- tp / n_pos
  precision <- tp / pmax(tp + fp, 1)

  recall <- c(0, recall)
  precision <- c(1, precision)

  sum(
    diff(recall) *
      (head(precision, -1L) + tail(precision, -1L)) / 2
  )
}

recall_at_k <- function(y, score, k) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])

  n_pos <- sum(y == 1L)
  if (n_pos == 0L) return(NA_real_)

  k <- min(as.integer(k), length(score))
  if (k <= 0L) return(NA_real_)

  ord <- order(score, decreasing = TRUE)
  sum(y[ord[seq_len(k)]] == 1L) / n_pos
}

mrr_score <- function(y, score) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])

  if (sum(y == 1L) == 0L) return(NA_real_)

  ord <- order(score, decreasing = TRUE)
  positive_ranks <- which(y[ord] == 1L)

  if (length(positive_ranks) == 0L) return(NA_real_)
  1 / positive_ranks[[1]]
}

ndcg_at_k <- function(y, score, k = 100L) {
  ok <- is.finite(score) & !is.na(y)
  y <- as.integer(y[ok])
  score <- as.numeric(score[ok])

  if (length(y) == 0L) return(NA_real_)

  k <- min(as.integer(k), length(y))
  ord <- order(score, decreasing = TRUE)
  rel <- y[ord[seq_len(k)]]
  discounts <- log2(seq_len(k) + 1)

  dcg <- sum(rel / discounts)

  ideal <- sort(y, decreasing = TRUE)[seq_len(k)]
  idcg <- sum(ideal / discounts)

  if (idcg <= 0) return(NA_real_)
  dcg / idcg
}

evaluate_scores <- function(y, score) {
  data.table(
    apparent_auprc = pr_auc_trapezoid(y, score),
    apparent_auroc = roc_auc_rank(y, score),
    recall_at_10 = recall_at_k(y, score, 10L),
    recall_at_50 = recall_at_k(y, score, 50L),
    recall_at_100 = recall_at_k(y, score, 100L),
    mrr = mrr_score(y, score),
    ndcg_at_100 = ndcg_at_k(y, score, 100L),
    observed_positive_prevalence = mean(y == 1L)
  )
}

# -----------------------------------------------------------------------------
# 4. Stratified repeated folds
# -----------------------------------------------------------------------------

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)

  fold <- integer(length(y))

  for (cls in c(0L, 1L)) {
    idx <- which(y == cls)
    idx <- sample(idx, length(idx), replace = FALSE)

    fold[idx] <- rep(
      seq_len(k),
      length.out = length(idx)
    )
  }

  fold
}

# -----------------------------------------------------------------------------
# 5. Fold-safe preprocessing
# -----------------------------------------------------------------------------
# Numeric/logical:
#   training median imputation + training mean/SD scaling + missing flag.
#
# Character/factor:
#   training levels, __MISSING__, __OTHER__, one-hot encoding.
#
# Zero-variance columns are removed using training data only.

fit_preprocessor <- function(train_dt, predictors) {
  specs <- list()
  pieces <- list()

  for (col in predictors) {
    x <- train_dt[[col]]

    if (is.numeric(x) || is.integer(x) || is.logical(x)) {
      raw <- if (is.logical(x)) as.numeric(x) else safe_numeric(x)

      missing <- is.na(raw) | !is.finite(raw)
      finite_values <- raw[!missing]

      med <- if (length(finite_values) > 0L) {
        median(finite_values, na.rm = TRUE)
      } else {
        0
      }

      imp <- raw
      imp[missing] <- med

      mu <- mean(imp)
      sdv <- stats::sd(imp)

      if (!is.finite(sdv) || sdv == 0) sdv <- 1

      mat <- cbind(
        (imp - mu) / sdv,
        as.numeric(missing)
      )

      colnames(mat) <- c(
        paste0("num__", col),
        paste0("missing__", col)
      )

      pieces[[length(pieces) + 1L]] <- mat

      specs[[col]] <- list(
        type = "numeric",
        median = med,
        mean = mu,
        sd = sdv
      )

    } else {
      raw <- as.character(x)
      raw[is.na(raw) | !nzchar(raw)] <- "__MISSING__"

      levels_train <- sort(unique(raw))
      levels_all <- unique(c(
        levels_train,
        "__MISSING__",
        "__OTHER__"
      ))

      mat <- matrix(
        0,
        nrow = length(raw),
        ncol = length(levels_all)
      )

      colnames(mat) <- paste0(
        "cat__",
        col,
        "__",
        make.names(levels_all, unique = TRUE)
      )

      for (j in seq_along(levels_all)) {
        mat[, j] <- as.numeric(raw == levels_all[[j]])
      }

      pieces[[length(pieces) + 1L]] <- mat

      specs[[col]] <- list(
        type = "categorical",
        levels = levels_all
      )
    }
  }

  X <- do.call(cbind, pieces)
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  colnames(X) <- make.names(
    colnames(X),
    unique = TRUE
  )

  keep <- vapply(
    seq_len(ncol(X)),
    function(j) {
      z <- X[, j]
      all(is.finite(z)) &&
        is.finite(stats::sd(z)) &&
        stats::sd(z) > 0
    },
    logical(1)
  )

  if (!any(keep)) {
    stop("Preprocessing produced no nonconstant predictors.", call. = FALSE)
  }

  X <- X[, keep, drop = FALSE]

  list(
    specs = specs,
    matrix_columns = colnames(X),
    X = X
  )
}

apply_preprocessor <- function(dt, prep) {
  pieces <- list()

  for (col in names(prep$specs)) {
    spec <- prep$specs[[col]]
    x <- dt[[col]]

    if (spec$type == "numeric") {
      raw <- if (is.logical(x)) as.numeric(x) else safe_numeric(x)

      missing <- is.na(raw) | !is.finite(raw)

      imp <- raw
      imp[missing] <- spec$median

      mat <- cbind(
        (imp - spec$mean) / spec$sd,
        as.numeric(missing)
      )

      colnames(mat) <- c(
        paste0("num__", col),
        paste0("missing__", col)
      )

      pieces[[length(pieces) + 1L]] <- mat

    } else {
      raw <- as.character(x)
      raw[is.na(raw) | !nzchar(raw)] <- "__MISSING__"
      raw[!raw %in% spec$levels] <- "__OTHER__"

      mat <- matrix(
        0,
        nrow = length(raw),
        ncol = length(spec$levels)
      )

      colnames(mat) <- paste0(
        "cat__",
        col,
        "__",
        make.names(spec$levels, unique = TRUE)
      )

      for (j in seq_along(spec$levels)) {
        mat[, j] <- as.numeric(raw == spec$levels[[j]])
      }

      pieces[[length(pieces) + 1L]] <- mat
    }
  }

  X <- do.call(cbind, pieces)
  X <- as.matrix(X)
  storage.mode(X) <- "double"

  colnames(X) <- make.names(
    colnames(X),
    unique = TRUE
  )

  absent <- setdiff(
    prep$matrix_columns,
    colnames(X)
  )

  if (length(absent) > 0L) {
    extra <- matrix(
      0,
      nrow = nrow(X),
      ncol = length(absent)
    )
    colnames(extra) <- absent
    X <- cbind(X, extra)
  }

  X[, prep$matrix_columns, drop = FALSE]
}

# -----------------------------------------------------------------------------
# 6. Base learners
# -----------------------------------------------------------------------------

fit_base_model <- function(algorithm, X, y, seed) {
  set.seed(seed)

  if (algorithm == "elastic_net") {
    nfolds_inner <- min(
      GLMNET_NFOLDS,
      sum(y == 1L),
      sum(y == 0L)
    )

    nfolds_inner <- max(3L, as.integer(nfolds_inner))

    fit <- glmnet::cv.glmnet(
      x = X,
      y = y,
      family = "binomial",
      alpha = GLMNET_ALPHA,
      nfolds = nfolds_inner,
      type.measure = "deviance",
      standardize = FALSE,
      intercept = TRUE
    )

    return(list(
      algorithm = algorithm,
      model = fit
    ))
  }

  if (algorithm == "ranger") {
    df <- as.data.frame(X)
    df$.y <- factor(y, levels = c(0, 1))

    fit <- ranger::ranger(
      dependent.variable.name = ".y",
      data = df,
      probability = TRUE,
      num.trees = RANGER_NUM_TREES,
      mtry = max(1L, floor(sqrt(ncol(X)))),
      min.node.size = RANGER_MIN_NODE_SIZE,
      importance = "none",
      seed = seed,
      num.threads = 1
    )

    return(list(
      algorithm = algorithm,
      model = fit
    ))
  }

  if (algorithm == "xgboost") {
    dtrain <- xgboost::xgb.DMatrix(
      data = X,
      label = y,
      nthread = 1
    )

    params <- list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      max_depth = XGB_MAX_DEPTH,
      eta = XGB_ETA,
      subsample = XGB_SUBSAMPLE,
      colsample_bytree = XGB_COLSAMPLE,
      min_child_weight = XGB_MIN_CHILD_WEIGHT,
      nthread = 1,
      seed = seed
    )

    fit <- xgboost::xgb.train(
      params = params,
      data = dtrain,
      nrounds = XGB_NROUNDS,
      verbose = 0
    )

    return(list(
      algorithm = algorithm,
      model = fit
    ))
  }

  if (algorithm == "svm_radial") {
    fit <- e1071::svm(
      x = X,
      y = factor(y, levels = c(0, 1)),
      kernel = "radial",
      cost = SVM_COST,
      gamma = 1 / max(1, ncol(X)),
      probability = TRUE,
      scale = FALSE
    )

    return(list(
      algorithm = algorithm,
      model = fit
    ))
  }

  stop(
    paste0("Unsupported algorithm: ", algorithm),
    call. = FALSE
  )
}

predict_base_model <- function(fit_object, X) {
  algorithm <- fit_object$algorithm
  model <- fit_object$model

  if (algorithm == "elastic_net") {
    return(
      as.numeric(
        predict(
          model,
          newx = X,
          s = "lambda.1se",
          type = "response"
        )
      )
    )
  }

  if (algorithm == "ranger") {
    pr <- predict(
      model,
      data = as.data.frame(X)
    )$predictions

    if (is.matrix(pr) || is.data.frame(pr)) {
      if (!is.null(colnames(pr)) && "1" %in% colnames(pr)) {
        return(as.numeric(pr[, "1"]))
      }
      return(as.numeric(pr[, ncol(pr)]))
    }

    return(as.numeric(pr))
  }

  if (algorithm == "xgboost") {
    dtest <- xgboost::xgb.DMatrix(
      data = X,
      nthread = 1
    )

    return(
      as.numeric(
        predict(model, dtest)
      )
    )
  }

  if (algorithm == "svm_radial") {
    pr <- predict(
      model,
      X,
      probability = TRUE
    )

    probs <- attr(pr, "probabilities")

    if (is.null(probs)) {
      stop("SVM probability predictions were not returned.", call. = FALSE)
    }

    if (!is.null(colnames(probs)) && "1" %in% colnames(probs)) {
      return(as.numeric(probs[, "1"]))
    }

    return(as.numeric(probs[, ncol(probs)]))
  }

  stop(
    paste0("Unsupported prediction algorithm: ", algorithm),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 7. Bagged PU learner
# -----------------------------------------------------------------------------

fit_bagged_pu <- function(
    X_train,
    y_train,
    algorithm,
    seed
) {
  positive_idx <- which(y_train == 1L)
  unlabeled_idx <- which(y_train == 0L)

  if (length(positive_idx) < 2L) {
    stop("Insufficient observed positives in training fold.", call. = FALSE)
  }

  if (length(unlabeled_idx) < 2L) {
    stop("Insufficient unlabeled genes in training fold.", call. = FALSE)
  }

  n_unlabeled_bag <- min(
    length(unlabeled_idx),
    max(
      length(positive_idx),
      as.integer(
        round(
          UNLABELED_TO_POSITIVE_RATIO *
            length(positive_idx)
        )
      )
    )
  )

  models <- vector("list", N_BAGS)

  for (b in seq_len(N_BAGS)) {
    bag_seed <- seed + b * 1009L
    set.seed(bag_seed)

    sampled_unlabeled <- sample(
      unlabeled_idx,
      size = n_unlabeled_bag,
      replace = FALSE
    )

    bag_idx <- c(
      positive_idx,
      sampled_unlabeled
    )

    bag_y <- c(
      rep(1L, length(positive_idx)),
      rep(0L, length(sampled_unlabeled))
    )

    models[[b]] <- fit_base_model(
      algorithm = algorithm,
      X = X_train[bag_idx, , drop = FALSE],
      y = bag_y,
      seed = bag_seed
    )
  }

  list(
    algorithm = algorithm,
    models = models,
    n_positive = length(positive_idx),
    n_unlabeled_per_bag = n_unlabeled_bag
  )
}

predict_bagged_pu <- function(bagged_fit, X) {
  pred_matrix <- vapply(
    bagged_fit$models,
    function(model_obj) {
      predict_base_model(model_obj, X)
    },
    numeric(nrow(X))
  )

  if (is.null(dim(pred_matrix))) {
    return(as.numeric(pred_matrix))
  }

  rowMeans(pred_matrix, na.rm = TRUE)
}

# -----------------------------------------------------------------------------
# 8. Load Step 07 datasets and manifest
# -----------------------------------------------------------------------------

for (path in c(
  M1_FILE,
  M2_FILE,
  M3_FILE,
  PREDICTOR_MANIFEST_FILE,
  MASTER_FILE
)) {
  assert_file(path)
}

datasets <- list(
  M1_expression = fread(M1_FILE),
  M2_expression_orthology = fread(M2_FILE),
  M3_expression_orthology_annotation = fread(M3_FILE)
)

predictor_manifest <- fread(PREDICTOR_MANIFEST_FILE)
master <- fread(MASTER_FILE)

for (nm in names(datasets)) {
  assert_unique_gene_ids(datasets[[nm]], nm)
}

assert_unique_gene_ids(master, "master_gene_feature_matrix")

# Recreate gene_length_bp because Step 07 created it in memory for M3.
if (
  !"gene_length_bp" %in% names(master) &&
  all(c("gene_start", "gene_end") %in% names(master))
) {
  master[
    ,
    gene_length_bp :=
      as.numeric(gene_end) -
      as.numeric(gene_start) +
      1
  ]
}

if ("is_external_landmark_validation" %in% names(master)) {
  master[
    ,
    external_validation_holdout :=
      !is.na(is_external_landmark_validation) &
      as.logical(is_external_landmark_validation)
  ]
} else if ("is_prespecified_benchmark" %in% names(master)) {
  master[
    ,
    external_validation_holdout :=
      !is.na(is_prespecified_benchmark) &
      as.logical(is_prespecified_benchmark)
  ]
} else {
  stop(
    "Could not identify benchmark holdout genes in master matrix.",
    call. = FALSE
  )
}

benchmark_master <- master[
  external_validation_holdout == TRUE
]

get_predictors <- function(model_set_name) {
  unique(
    predictor_manifest[
      modelling_set == model_set_name &
        allowed_for_pu_training == TRUE,
      predictor
    ]
  )
}

predictor_sets <- lapply(
  MODEL_SETS,
  get_predictors
)

names(predictor_sets) <- MODEL_SETS

for (nm in MODEL_SETS) {
  missing_predictors <- setdiff(
    predictor_sets[[nm]],
    names(datasets[[nm]])
  )

  if (length(missing_predictors) > 0L) {
    stop(
      paste0(
        nm,
        " is missing predictor(s): ",
        paste(missing_predictors, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

# -----------------------------------------------------------------------------
# 9. Repeated cross-fitted PU modelling
# -----------------------------------------------------------------------------

prediction_records <- list()
metric_records <- list()

prediction_counter <- 0L
metric_counter <- 0L

for (model_set_name in MODEL_SETS) {
  dt <- copy(datasets[[model_set_name]])
  predictors <- predictor_sets[[model_set_name]]
  y <- as.integer(dt$pu_label)

  log_step(
    "08",
    paste0(
      model_set_name,
      ": ",
      format(nrow(dt), big.mark = ","),
      " genes; ",
      sum(y == 1L),
      " observed positives"
    )
  )

  for (repeat_id in seq_len(N_REPEATS)) {
    repeat_seed <- MASTER_SEED +
      match(model_set_name, MODEL_SETS) * 100000L +
      repeat_id * 10000L

    fold_id <- make_stratified_folds(
      y = y,
      k = N_FOLDS,
      seed = repeat_seed
    )

    for (fold in seq_len(N_FOLDS)) {
      train_idx <- which(fold_id != fold)
      test_idx <- which(fold_id == fold)

      train_dt <- dt[train_idx]
      test_dt <- dt[test_idx]

      prep <- fit_preprocessor(
        train_dt = train_dt,
        predictors = predictors
      )

      X_train <- prep$X
      X_test <- apply_preprocessor(
        dt = test_dt,
        prep = prep
      )

      y_train <- as.integer(train_dt$pu_label)
      y_test <- as.integer(test_dt$pu_label)

      for (algorithm in ALGORITHMS) {
        fit_seed <- repeat_seed +
          fold * 1000L +
          match(algorithm, ALGORITHMS) * 100L

        log_step(
          "08",
          paste0(
            model_set_name,
            " | repeat ",
            repeat_id,
            "/",
            N_REPEATS,
            " | fold ",
            fold,
            "/",
            N_FOLDS,
            " | ",
            algorithm
          )
        )

        bagged_fit <- fit_bagged_pu(
          X_train = X_train,
          y_train = y_train,
          algorithm = algorithm,
          seed = fit_seed
        )

        score <- predict_bagged_pu(
          bagged_fit = bagged_fit,
          X = X_test
        )

        prediction_counter <- prediction_counter + 1L

        prediction_records[[prediction_counter]] <- data.table(
          model_set = model_set_name,
          algorithm = algorithm,
          repeat_id = repeat_id,
          fold = fold,
          gene_id = test_dt$gene_id,
          pu_label = y_test,
          pu_score = score
        )

        fold_metrics <- evaluate_scores(
          y = y_test,
          score = score
        )

        metric_counter <- metric_counter + 1L

        metric_records[[metric_counter]] <- cbind(
          data.table(
            model_set = model_set_name,
            algorithm = algorithm,
            repeat_id = repeat_id,
            fold = fold,
            n_test = length(y_test),
            n_observed_positive = sum(y_test == 1L),
            n_unlabeled = sum(y_test == 0L),
            n_bags = N_BAGS,
            n_unlabeled_per_bag =
              bagged_fit$n_unlabeled_per_bag
          ),
          fold_metrics
        )
      }
    }
  }
}

cv_predictions <- rbindlist(
  prediction_records,
  fill = TRUE,
  use.names = TRUE
)

cv_metrics <- rbindlist(
  metric_records,
  fill = TRUE,
  use.names = TRUE
)

# -----------------------------------------------------------------------------
# 10. Aggregate repeated out-of-fold scores
# -----------------------------------------------------------------------------

gene_oof <- cv_predictions[
  ,
  .(
    pu_label = unique(pu_label)[[1]],
    mean_oof_pu_score = mean(pu_score, na.rm = TRUE),
    sd_oof_pu_score = stats::sd(pu_score, na.rm = TRUE),
    min_oof_pu_score = min(pu_score, na.rm = TRUE),
    max_oof_pu_score = max(pu_score, na.rm = TRUE),
    n_oof_predictions = .N
  ),
  by = .(
    model_set,
    algorithm,
    gene_id
  )
]

gene_oof[
  !is.finite(sd_oof_pu_score),
  sd_oof_pu_score := 0
]

metrics_summary <- cv_metrics[
  ,
  .(
    apparent_auprc_mean =
      mean(apparent_auprc, na.rm = TRUE),
    apparent_auprc_sd =
      stats::sd(apparent_auprc, na.rm = TRUE),

    apparent_auroc_mean =
      mean(apparent_auroc, na.rm = TRUE),
    apparent_auroc_sd =
      stats::sd(apparent_auroc, na.rm = TRUE),

    recall_at_10_mean =
      mean(recall_at_10, na.rm = TRUE),
    recall_at_10_sd =
      stats::sd(recall_at_10, na.rm = TRUE),

    recall_at_50_mean =
      mean(recall_at_50, na.rm = TRUE),
    recall_at_50_sd =
      stats::sd(recall_at_50, na.rm = TRUE),

    recall_at_100_mean =
      mean(recall_at_100, na.rm = TRUE),
    recall_at_100_sd =
      stats::sd(recall_at_100, na.rm = TRUE),

    mrr_mean =
      mean(mrr, na.rm = TRUE),
    mrr_sd =
      stats::sd(mrr, na.rm = TRUE),

    ndcg_at_100_mean =
      mean(ndcg_at_100, na.rm = TRUE),
    ndcg_at_100_sd =
      stats::sd(ndcg_at_100, na.rm = TRUE),

    observed_positive_prevalence_mean =
      mean(observed_positive_prevalence, na.rm = TRUE)
  ),
  by = .(
    model_set,
    algorithm
  )
]

setorder(
  metrics_summary,
  -apparent_auprc_mean,
  -apparent_auroc_mean,
  model_set,
  algorithm
)

metrics_summary[
  ,
  apparent_auprc_rank := seq_len(.N)
]

# -----------------------------------------------------------------------------
# 11. Fit final ensembles and score external benchmark genes
# -----------------------------------------------------------------------------

benchmark_score_records <- list()
benchmark_counter <- 0L

for (model_set_name in MODEL_SETS) {
  train_dt <- copy(datasets[[model_set_name]])
  predictors <- predictor_sets[[model_set_name]]

  missing_benchmark_predictors <- setdiff(
    predictors,
    names(benchmark_master)
  )

  if (length(missing_benchmark_predictors) > 0L) {
    stop(
      paste0(
        "Benchmark master is missing predictor(s) for ",
        model_set_name,
        ": ",
        paste(missing_benchmark_predictors, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  prep <- fit_preprocessor(
    train_dt = train_dt,
    predictors = predictors
  )

  X_train <- prep$X
  y_train <- as.integer(train_dt$pu_label)

  X_benchmark <- apply_preprocessor(
    dt = benchmark_master,
    prep = prep
  )

  predictor_coverage <- rowSums(
    !is.na(
      benchmark_master[
        ,
        ..predictors
      ]
    )
  ) / length(predictors)

  for (algorithm in ALGORITHMS) {
    final_seed <- MASTER_SEED +
      900000L +
      match(model_set_name, MODEL_SETS) * 10000L +
      match(algorithm, ALGORITHMS) * 1000L

    log_step(
      "08",
      paste0(
        "External benchmark scoring | ",
        model_set_name,
        " | ",
        algorithm
      )
    )

    final_fit <- fit_bagged_pu(
      X_train = X_train,
      y_train = y_train,
      algorithm = algorithm,
      seed = final_seed
    )

    benchmark_score <- predict_bagged_pu(
      bagged_fit = final_fit,
      X = X_benchmark
    )

    # If a benchmark has no measured predictor at all for a given model set,
    # report NA rather than a score generated purely from imputation.
    benchmark_score[
      predictor_coverage == 0
    ] <- NA_real_

    # Use an explicit external loop variable name to avoid ambiguity with
    # the data.table column named algorithm.
    algorithm_name <- algorithm

    reference_scores <- gene_oof[
      model_set == model_set_name &
        algorithm == algorithm_name,
      mean_oof_pu_score
    ]

    percentile <- vapply(
      benchmark_score,
      function(s) {
        if (
          !is.finite(s) ||
          length(reference_scores) == 0L
        ) {
          return(NA_real_)
        }

        mean(
          reference_scores <= s,
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    benchmark_counter <- benchmark_counter + 1L

    benchmark_score_records[[benchmark_counter]] <- data.table(
      model_set = model_set_name,
      algorithm = algorithm_name,
      gene_id = benchmark_master$gene_id,

      benchmark_name =
        if ("benchmark_name" %in% names(benchmark_master)) {
          benchmark_master$benchmark_name
        } else {
          NA_character_
        },

      benchmark_class =
        if ("benchmark_class" %in% names(benchmark_master)) {
          benchmark_master$benchmark_class
        } else {
          NA_character_
        },

      predictor_coverage_fraction =
        predictor_coverage,

      external_pu_score =
        benchmark_score,

      percentile_vs_crossfitted_nonholdout_genes =
        percentile
    )
  }
}

benchmark_scores <- rbindlist(
  benchmark_score_records,
  fill = TRUE,
  use.names = TRUE
)

# -----------------------------------------------------------------------------
# 12. Configuration, QC, provenance
# -----------------------------------------------------------------------------

configuration <- data.table(
  parameter = c(
    "MASTER_SEED",
    "N_FOLDS",
    "N_REPEATS",
    "N_BAGS",
    "UNLABELED_TO_POSITIVE_RATIO",
    "GLMNET_ALPHA",
    "GLMNET_NFOLDS",
    "RANGER_NUM_TREES",
    "RANGER_MIN_NODE_SIZE",
    "XGB_NROUNDS",
    "XGB_MAX_DEPTH",
    "XGB_ETA",
    "XGB_SUBSAMPLE",
    "XGB_COLSAMPLE",
    "XGB_MIN_CHILD_WEIGHT",
    "SVM_COST"
  ),
  value = c(
    MASTER_SEED,
    N_FOLDS,
    N_REPEATS,
    N_BAGS,
    UNLABELED_TO_POSITIVE_RATIO,
    GLMNET_ALPHA,
    GLMNET_NFOLDS,
    RANGER_NUM_TREES,
    RANGER_MIN_NODE_SIZE,
    XGB_NROUNDS,
    XGB_MAX_DEPTH,
    XGB_ETA,
    XGB_SUBSAMPLE,
    XGB_COLSAMPLE,
    XGB_MIN_CHILD_WEIGHT,
    SVM_COST
  )
)

expected_fold_models <-
  length(MODEL_SETS) *
  length(ALGORITHMS) *
  N_REPEATS *
  N_FOLDS

if (nrow(cv_metrics) != expected_fold_models) {
  stop(
    paste0(
      "Expected ",
      expected_fold_models,
      " fold-level evaluations but obtained ",
      nrow(cv_metrics),
      "."
    ),
    call. = FALSE
  )
}

if (
  any(
    gene_oof$n_oof_predictions != N_REPEATS
  )
) {
  stop(
    paste0(
      "Cross-fitting integrity failure: every model-set/algorithm/gene ",
      "should have exactly ",
      N_REPEATS,
      " OOF predictions."
    ),
    call. = FALSE
  )
}

qc <- data.table(
  metric = c(
    "Model sets",
    "Algorithms",
    "CV folds",
    "CV repeats",
    "PU bags per ensemble",
    "Expected fold-level evaluations",
    "Completed fold-level evaluations",
    "CV prediction rows",
    "Aggregated OOF gene-score rows",
    "External benchmark genes",
    "External benchmark score rows",
    "Benchmark genes used in fitting",
    "FlyBase phenotype predictors used",
    "Step 04 rank/score predictors used",
    "Step 05 CRISPR genome-wide predictors used",
    "Population-genomics predictors used"
  ),
  value = c(
    length(MODEL_SETS),
    length(ALGORITHMS),
    N_FOLDS,
    N_REPEATS,
    N_BAGS,
    expected_fold_models,
    nrow(cv_metrics),
    nrow(cv_predictions),
    nrow(gene_oof),
    nrow(benchmark_master),
    nrow(benchmark_scores),
    0,
    0,
    0,
    0,
    0
  )
)

provenance <- data.table(
  component = c(
    "Learning paradigm",
    "Pseudo-negatives",
    "Cross-fitting",
    "Primary metric",
    "Metric interpretation",
    "Preprocessing",
    "Elastic net",
    "Ranger",
    "XGBoost",
    "SVM",
    "External benchmarks",
    "Population genomics",
    "CRISPR"
  ),
  specification = c(
    "Bagged positive-unlabeled learning.",
    paste0(
      "Each bag uses all observed positives plus a fresh random subset of ",
      "unlabeled genes at ratio ",
      UNLABELED_TO_POSITIVE_RATIO,
      ":1."
    ),
    paste0(
      N_FOLDS,
      "-fold stratified cross-fitting repeated ",
      N_REPEATS,
      " times; preprocessing is fit inside each training fold."
    ),
    "Apparent PR-AUC/AUPRC for observed-positive-vs-unlabeled discrimination.",
    paste0(
      "Metrics are PU proxy-discrimination metrics and are not conventional ",
      "positive-vs-true-negative diagnostic accuracy."
    ),
    paste0(
      "Training-fold median imputation and scaling for numeric/logical fields, ",
      "missingness indicators, training-derived categorical one-hot encoding, ",
      "and training-only zero-variance removal."
    ),
    paste0(
      "glmnet binomial elastic net alpha=",
      GLMNET_ALPHA,
      "; internal CV; lambda.1se prediction."
    ),
    paste0(
      "ranger probability forest; ",
      RANGER_NUM_TREES,
      " trees; min.node.size=",
      RANGER_MIN_NODE_SIZE,
      "."
    ),
    paste0(
      "xgboost binary:logistic; nrounds=",
      XGB_NROUNDS,
      "; max_depth=",
      XGB_MAX_DEPTH,
      "; eta=",
      XGB_ETA,
      "."
    ),
    paste0(
      "e1071 radial SVM; cost=",
      SVM_COST,
      "; gamma=1/p."
    ),
    "Six prespecified benchmark genes are external holdouts and never enter fitting.",
    "Ag1000G remains pending and unused.",
    "Step 05 CRISPR variables remain excluded from genome-wide PU fitting."
  )
)

# -----------------------------------------------------------------------------
# 13. Write outputs
# -----------------------------------------------------------------------------

fwrite(cv_predictions, CV_PREDICTIONS_FILE, na = "NA")
fwrite(cv_metrics, CV_METRICS_FILE, na = "NA")
fwrite(metrics_summary, METRICS_SUMMARY_FILE, na = "NA")
fwrite(gene_oof, GENE_OOF_FILE, na = "NA")
fwrite(benchmark_scores, BENCHMARK_SCORES_FILE, na = "NA")
fwrite(configuration, CONFIG_FILE, na = "NA")
fwrite(qc, QC_FILE, na = "NA")
fwrite(provenance, PROVENANCE_FILE, na = "NA")

capture.output(
  sessionInfo(),
  file = SESSION_INFO_FILE
)

checksum_files <- c(
  M1_FILE,
  M2_FILE,
  M3_FILE,
  PREDICTOR_MANIFEST_FILE,
  MASTER_FILE,
  CV_PREDICTIONS_FILE,
  CV_METRICS_FILE,
  METRICS_SUMMARY_FILE,
  GENE_OOF_FILE,
  BENCHMARK_SCORES_FILE,
  CONFIG_FILE,
  QC_FILE,
  PROVENANCE_FILE,
  SESSION_INFO_FILE
)

checksum_files <- checksum_files[
  file.exists(checksum_files)
]

write_checksum(
  checksum_files,
  CHECKSUM_FILE
)

# -----------------------------------------------------------------------------
# 14. Console summary
# -----------------------------------------------------------------------------

cat(
  "\n",
  "============================================================\n",
  "MOSQEDIT-R STEP 08 COMPLETED SUCCESSFULLY\n",
  "Repeated cross-fitted bagged PU model benchmarking\n",
  "============================================================\n",
  "Model sets:                              ",
  length(MODEL_SETS),
  "\n",
  "Algorithms:                              ",
  length(ALGORITHMS),
  "\n",
  "Cross-validation folds:                 ",
  N_FOLDS,
  "\n",
  "Cross-validation repeats:               ",
  N_REPEATS,
  "\n",
  "PU bags per ensemble:                   ",
  N_BAGS,
  "\n",
  "Fold-level evaluations:                 ",
  nrow(cv_metrics),
  "\n",
  "External benchmark genes:               ",
  nrow(benchmark_master),
  "\n",
  "Benchmark genes used for training:      NO\n",
  "FlyBase phenotype predictors used:      NO\n",
  "Step 04 score/rank predictors used:     NO\n",
  "Step 05 CRISPR predictors used:         NO\n",
  "Ag1000G predictors used:                NO\n",
  "Python used:                            NO\n",
  "============================================================\n",
  sep = ""
)

cat(
  "\nCross-validated metric summary, ordered by apparent AUPRC:\n"
)

print(
  metrics_summary[
    ,
    .(
      apparent_auprc_rank,
      model_set,
      algorithm,
      apparent_auprc_mean,
      apparent_auprc_sd,
      apparent_auroc_mean,
      apparent_auroc_sd,
      recall_at_50_mean,
      recall_at_100_mean,
      mrr_mean,
      ndcg_at_100_mean
    )
  ]
)

cat(
  "\nExternal benchmark PU scores:\n"
)

print(
  benchmark_scores[
    order(
      gene_id,
      model_set,
      algorithm
    )
  ]
)

cat(
  "\nQC summary:\n"
)

print(qc)

log_step(
  "08",
  "Repeated cross-fitted bagged PU model benchmarking completed successfully"
)

