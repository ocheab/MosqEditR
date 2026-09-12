source("R/helpers.R")
suppressPackageStartupMessages({library(data.table); library(ranger)})
log_step("07", "Training bagged positive-unlabelled target-prioritization model")
d <- fread("data_processed/master_gene_feature_matrix.csv")
# Training positives are transferred Drosophila female-sterility phenotype labels, following the biological logic
# of Hammond et al. but without treating every unstudied Anopheles gene as a confirmed negative.
if (!"flybase_female_sterile" %in% names(d)) stop("FlyBase sterility label missing.")
d[, y_pos := flybase_female_sterile %in% TRUE]
# Known Anopheles genetic-control landmarks are held out of label construction to preserve external validation.
d[known_vector_control_target %in% TRUE, y_pos := FALSE]
if (sum(d$y_pos, na.rm=TRUE) < 20) stop("Fewer than 20 transferred female-sterility positives are available. Supply/version the FlyBase phenotype export before modelling; the pipeline will not manufacture negative/positive labels.")
# Numeric predictors; explicitly exclude outcome-like or post-model fields.
exclude <- c("known_vector_control_target","y_pos","flybase_female_sterile")
num <- names(d)[vapply(d, is.numeric, logical(1))]
features <- setdiff(num, exclude)
features <- features[vapply(d[, ..features], function(z) mean(is.finite(z)) >= 0.70, logical(1))]
# Prevent phenotype-count variable from trivially leaking the label if it was constructed from the same FlyBase text.
features <- setdiff(features, "flybase_reproductive_phenotype_count")
if (length(features) < 3) stop("Too few complete non-leaking numeric features. Check source acquisition before modelling.")
B <- 500L
pred <- matrix(NA_real_, nrow=nrow(d), ncol=B)
set.seed(20260911)
for (b in seq_len(B)) {
  pidx <- which(d$y_pos)
  uidx <- which(!d$y_pos & !d$known_vector_control_target)
  nneg <- min(length(uidx), max(10L*length(pidx), 500L))
  neg <- sample(uidx, nneg, replace=FALSE)
  tr <- c(pidx, neg)
  X <- as.data.frame(d[, ..features])
  for (j in features) {
    med <- median(X[[j]][tr], na.rm=TRUE)
    if (!is.finite(med)) med <- 0
    X[[j]][!is.finite(X[[j]])] <- med
  }
  fit <- ranger(x=X[tr, , drop=FALSE], y=factor(d$y_pos[tr], levels=c(FALSE,TRUE)), probability=TRUE,
                num.trees=1000, importance="permutation", seed=20260911+b)
  pred[,b] <- predict(fit, data=X)$predictions[,"TRUE"]
}
d[, pu_probability := rowMeans(pred, na.rm=TRUE)]
d[, pu_probability_sd := apply(pred,1,sd,na.rm=TRUE)]
# Evolutionary/editability evidence is incorporated only when actually observed; missing layers are never replaced by hand scores.
robust_cols <- intersect(c("nucleotide_diversity_pi_proxy","high_freq_variant_density"), names(d))
if (length(robust_cols)) {
  # Convert lower-is-better risk features to percentile robustness and average them.
  rob <- lapply(robust_cols, function(j) 1 - frank(d[[j]], ties.method="average", na.last="keep")/sum(is.finite(d[[j]])))
  rob <- as.data.table(rob)
  d[, evolutionary_robustness := rowMeans(rob, na.rm=TRUE)]
  d[!is.finite(evolutionary_robustness), evolutionary_robustness := NA_real_]
} else d[, evolutionary_robustness := NA_real_]
# Prespecified ensemble: PU discovery evidence is primary. Evolutionary robustness contributes only when present.
d[, ensemble_score := ifelse(is.finite(evolutionary_robustness), 0.75*pu_probability + 0.25*evolutionary_robustness, pu_probability)]
d[, analysis_maturity := "MODEL"]
fwrite(d, "results/candidate_lists/model_scored_genes.csv")
saveRDS(list(features=features, pred=pred, seed=20260911, B=B), "results/model_bundle.rds")

