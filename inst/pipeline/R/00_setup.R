source("R/helpers.R")
for (d in c("data_raw", "data_processed", "metadata", "figures", "logs", "results", "publication", "release")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
log_step("00", "MosqEdit-R analysis workspace initialized")
cat("MosqEdit-R project directories are ready.\n")

