suppressPackageStartupMessages(library(renv))
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
renv::snapshot(prompt=FALSE)
# Hash all final derived outputs for the audit trail.
fs <- list.files(c("data_processed","results"), recursive=TRUE, full.names=TRUE)
h <- data.frame(file=fs, sha256=vapply(fs, digest::digest, character(1), algo="sha256", file=TRUE))
write.table(h, "logs/final_output_sha256.tsv", sep="\t", row.names=FALSE, quote=FALSE)

