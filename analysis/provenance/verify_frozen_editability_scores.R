# Read-only verification of the frozen 205-gene exon-aware reference editability score.
# Usage: Rscript verify_frozen_editability_scores.R /path/to/05_crispr_editability.csv
#        Rscript verify_frozen_editability_scores.R /path/to/data_processed
#        Rscript verify_frozen_editability_scores.R /path/to/MosqEditR_Phase2_Targeted_Extraction1.zip
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide a path to the processed CSV, its directory, or archive ZIP.")
source <- args[[1L]]
if (!file.exists(source) && !dir.exists(source)) stop("Input does not exist: ", source)
if (dir.exists(source)) source <- file.path(source, "05_crispr_editability.csv")
if (grepl("\\.zip$", source, ignore.case = TRUE)) {
  entries <- utils::unzip(source, list = TRUE)$Name
  matches <- entries[grepl("(^|/)05_crispr_editability\\.csv$", entries)]
  if (length(matches) != 1L) stop("Expected one 05_crispr_editability.csv inside ZIP; found ", length(matches))
  connection <- unz(source, matches)
  dat <- read.csv(connection, check.names = FALSE, stringsAsFactors = FALSE)
  close(connection)
} else {
  dat <- read.csv(source, check.names = FALSE, stringsAsFactors = FALSE)
}
required <- c(
  "gene_id", "editability_site_abundance_component",
  "editability_multiplex_component",
  "editability_transcript_coverage_component",
  "editability_gc_balance_component",
  "crispr_exon_aware_reference_editability_score"
)
stopifnot(all(required %in% names(dat)), nrow(dat) == 205L,
          !anyDuplicated(dat$gene_id))
components <- as.matrix(dat[,required[2:5]])
storage.mode(components) <- "double"
observed <- as.numeric(dat$crispr_exon_aware_reference_editability_score)
stopifnot(all(is.finite(components)), all(is.finite(observed)))
old_formula <- drop(components %*% c(.40, .25, .20, .15))
reconciled_formula <- drop(components %*% c(.35, .25, .25, .15))
d1 <- abs(observed - reconciled_formula)
d2 <- abs(observed - old_formula)
cat("Frozen genes assessed: ", length(observed), "\n", sep = "")
cat("35/25/25/15: matching = ", sum(d1 < 1e-12), "/205; max absolute diff = ",
    format(max(d1), digits = 8), "\n", sep = "")
cat("40/25/20/15: matching = ", sum(d2 < 1e-12), "/205; max absolute diff = ",
    format(max(d2), digits = 8), "\n", sep = "")
stopifnot(max(d1) < 1e-12, sum(d2 < 1e-12) < nrow(dat))
cat("PASS: Archived scores agree with 35/25/25/15; frozen CSV unchanged.\n")
