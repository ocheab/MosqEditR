source("R/helpers.R")
suppressPackageStartupMessages({library(reticulate); library(data.table)})
log_step("04", "Computing resumable Ag1000G Phase 3 gene-level variation summaries")
# MalariaGEN's maintained Python API is called through reticulate. The companion routine
# uses Ag3.snp_allele_frequencies() per gene to avoid materialising chromosome-scale genotype arrays.
py <- Sys.getenv("MOSQEDIT_PYTHON", unset="")
if (nzchar(py)) use_python(py, required=TRUE)
if (!py_module_available("malariagen_data")) stop("Python package 'malariagen_data' is required. Create an environment per MalariaGEN documentation and set MOSQEDIT_PYTHON.")
source_python("R/ag1000g_gene_variation.py")
var <- compute_gene_variation_ag3(
  annotation_csv="data_processed/01_gene_annotation.csv",
  species=c("gambiae", "coluzzii"),
  min_minor_allele_frequency=0.01,
  cache_csv="data_raw/ag3_gene_variation_cache.csv"
)
var <- as.data.table(var)
if (!nrow(var)) stop("Ag1000G variation routine returned zero genes.")
fwrite(var, "data_processed/04_ag1000g_gene_variation.csv")
write_checksum(c("data_raw/ag3_gene_variation_cache.csv", "data_processed/04_ag1000g_gene_variation.csv"), "logs/04_checksums.tsv")

