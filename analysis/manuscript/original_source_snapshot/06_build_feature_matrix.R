source("R/helpers.R")
suppressPackageStartupMessages({library(data.table)})
log_step("06", "Building auditable master gene-feature matrix")
files <- c(annotation="data_processed/01_gene_annotation.csv",
           expression="data_processed/02_mozatlas_gene_expression.csv",
           orthopheno="data_processed/03_gene_orthology_phenotypes.csv",
           variation="data_processed/04_ag1000g_gene_variation.csv",
           editability="data_processed/05_crispr_editability.csv")
miss <- files[!file.exists(files)]
if (length(miss)) stop("Missing upstream files: ", paste(miss, collapse=", "))
x <- fread(files[["annotation"]])
for (nm in names(files)[-1]) x <- merge(x, fread(files[[nm]]), by="gene_id", all.x=TRUE)
land <- fread("metadata/landmark_targets.csv")
x[, known_vector_control_target := gene_id %in% land$gene_id]
x <- merge(x, land[, .(gene_id, landmark_class=strategy)], by="gene_id", all.x=TRUE)
# Experimental landmarks are validation objects, not training labels.
x[, analysis_maturity := "RECOMPUTED"]
fwrite(x, "data_processed/master_gene_feature_matrix.csv")
write_checksum("data_processed/master_gene_feature_matrix.csv", "logs/06_checksums.tsv")

