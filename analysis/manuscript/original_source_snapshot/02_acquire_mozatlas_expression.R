source("R/helpers.R")
suppressPackageStartupMessages({library(GEOquery); library(Biobase); library(data.table); library(stringr)})
log_step("02", "Downloading GSE21689/MozAtlas processed expression")
gse <- getGEO("GSE21689", GSEMatrix=TRUE, getGPL=TRUE)
if (length(gse) < 1) stop("GSE21689 retrieval returned no ExpressionSet.")
es <- gse[[1]]
expr <- exprs(es); ph <- pData(es)
# Keep only Anopheles samples; ACPS controls are excluded from mosquito gene feature engineering.
titles <- as.character(ph$title)
keep <- !grepl("ACPS", titles, ignore.case=TRUE)
expr <- expr[, keep, drop=FALSE]; titles <- titles[keep]
# Determine signal scale without silently transforming.
q99 <- as.numeric(quantile(expr, 0.99, na.rm=TRUE))
expr_scale <- if (q99 > 100) "linear_to_log2p1" else "as_retrieved_GCRMA"
if (q99 > 100) expr <- log2(expr + 1)
# Platform mapping: search all feature annotation columns for AGAP identifiers.
fd <- fData(es)
probe_id <- rownames(fd)
blob <- apply(fd, 1, function(z) paste(z, collapse=" | "))
agap <- stringr::str_extract(blob, "AGAP\\d{6}")
map <- data.table(probe_id=probe_id, gene_id=agap)
map <- map[!is.na(gene_id)]
# Tissue labels from GEO titles.
lab <- tolower(titles)
lab <- gsub("anopheles[ _-]*", "", lab)
lab <- gsub("_[0-9]+$", "", lab)
lab <- gsub("[ .-]+", "_", lab)
# Aggregate replicate means at probe level then gene level.
probe_dt <- data.table(probe_id=rownames(expr))
for (t in unique(lab)) probe_dt[[paste0("expr_", t)]] <- rowMeans(expr[, lab==t, drop=FALSE], na.rm=TRUE)
probe_dt <- merge(probe_dt, map, by="probe_id")
numcols <- grep("^expr_", names(probe_dt), value=TRUE)
gene_expr <- probe_dt[, lapply(.SD, mean, na.rm=TRUE), by=gene_id, .SDcols=numcols]
# Harmonize expected names where present.
setnames(gene_expr, names(gene_expr), gsub("expr_female_ovary", "expr_ovary", names(gene_expr)))
setnames(gene_expr, names(gene_expr), gsub("expr_male_testis", "expr_testis", names(gene_expr)))
# Feature calculations.
expr_cols <- grep("^expr_", names(gene_expr), value=TRUE)
gene_expr[, tau_tissue_specificity := apply(.SD, 1, tau_index), .SDcols=expr_cols]
matched <- c("body","carcass","head","malpighian","midgut","salivary")
fcols <- paste0("expr_female_", matched); mcols <- paste0("expr_male_", matched)
fcols <- fcols[fcols %in% names(gene_expr)]; mcols <- mcols[mcols %in% names(gene_expr)]
if (length(fcols) && length(mcols)) {
  gene_expr[, female_mean_matched := rowMeans(.SD, na.rm=TRUE), .SDcols=fcols]
  gene_expr[, male_mean_matched := rowMeans(.SD, na.rm=TRUE), .SDcols=mcols]
  gene_expr[, female_male_delta := female_mean_matched - male_mean_matched]
}
if ("expr_ovary" %in% names(gene_expr)) {
  female_nonov <- setdiff(grep("^expr_female_", names(gene_expr), value=TRUE), "expr_ovary")
  if (length(female_nonov)) gene_expr[, ovary_specificity_delta := expr_ovary - rowMeans(.SD, na.rm=TRUE), .SDcols=female_nonov]
}
fwrite(gene_expr, "data_processed/02_mozatlas_gene_expression.csv")
writeLines(paste("expression_scale_decision", expr_scale, sep="\t"), "logs/02_expression_scale.tsv")
write_checksum("data_processed/02_mozatlas_gene_expression.csv", "logs/02_checksums.tsv")

