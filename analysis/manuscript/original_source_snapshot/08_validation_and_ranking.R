source("R/helpers.R")
suppressPackageStartupMessages({library(data.table)})
log_step("08", "Ranking candidates and landmark hold-out evaluation")
d <- fread("results/candidate_lists/model_scored_genes.csv")
d <- d[order(-ensemble_score)]
d[, rank := seq_len(.N)]
# Secondary non-coding landmarks are not forced into the protein-coding rank universe.
land <- fread("metadata/landmark_targets.csv")
val <- merge(land, d[, .(gene_id, ensemble_score, rank)], by="gene_id", all.x=TRUE)
fwrite(val, "results/tables/landmark_validation_ranks.csv")
# Bootstrap rank stability from saved PU prediction matrix.
bundle <- readRDS("results/model_bundle.rds")
pm <- bundle$pred
ranks <- apply(pm, 2, function(z) rank(-z, ties.method="average", na.last="keep"))
d[, bootstrap_rank_median := apply(ranks,1,median,na.rm=TRUE)]
d[, bootstrap_top50_prob := rowMeans(ranks <= 50, na.rm=TRUE)]
out <- d[, .(gene_id, ensemble_score, rank, bootstrap_rank_median, bootstrap_top50_prob, analysis_maturity)]
fwrite(out, "results/candidate_lists/ranked_targets.csv")

