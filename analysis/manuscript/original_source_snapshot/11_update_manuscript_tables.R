# This script exports the exact tables consumed by the Word manuscript.
# It deliberately does not patch Word XML automatically; instead it writes machine-readable tables and a manuscript_value_manifest.csv
# so that every reported value can be traced to an output file and recomputed.
suppressPackageStartupMessages(library(data.table))
files <- c("results/tables/top50_candidates.csv", "results/tables/landmark_validation_ranks.csv")
rows <- list()
for (f in files[file.exists(files)]) {
  x <- fread(f)
  rows[[length(rows)+1]] <- data.table(source_file=f, n_rows=nrow(x), sha256=digest::digest(file=f, algo="sha256"))
}
if (length(rows)) fwrite(rbindlist(rows), "results/tables/manuscript_value_manifest.csv")

