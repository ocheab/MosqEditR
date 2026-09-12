# MosqEdit-R manuscript-pipeline compatibility helpers.
# Requires the installed MosqEditR package.
if (!requireNamespace("MosqEditR", quietly = TRUE)) {
  stop("Install MosqEditR before running the manuscript pipeline.", call. = FALSE)
}
log_step <- function(step, msg) MosqEditR:::.mosqedit_log_step(step, msg)
first_nonempty <- function(x) MosqEditR:::.mosqedit_first_nonempty(x)
tau_index <- function(x) MosqEditR::mosqedit_tau(x)
rank01 <- function(x, higher_is_better = TRUE) MosqEditR::mosqedit_percentile(x, higher_is_better)
write_checksum <- function(files, outfile) MosqEditR:::.mosqedit_write_checksum(files, outfile)

