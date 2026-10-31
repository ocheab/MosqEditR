# Ag1000G Phase 2 AR1: descriptive QC of the frozen 530-site accessibility panel.
# Uses base R. Does not compute conservation, genotype missingness or allele frequencies.
# Execute from the MosqEditR_Phase2_Targeted_Extraction project directory:
#   Rscript 05_accessibility_summary.R
infile <- file.path("output", "AR1_frozen_site_accessibility_530.csv")
if (!file.exists(infile)) stop("Missing: ", infile, "; run 04_extract_accessibility.py first")
d <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(d) == 530L, length(unique(d$population_site_id)) == 530L)
stopifnot(sum(d$status == "EVALUATED_AG1000G_PHASE2_SITE_MASK") == 515L)
stopifnot(sum(d$status == "UNKN_COORDINATES_NOT_EVALUABLE") == 15L)
stopifnot(!anyNA(d$fraction_accessible[d$status == "EVALUATED_AG1000G_PHASE2_SITE_MASK"]))
arms <- c("2L", "2R", "3L", "3R", "X")
rows <- lapply(arms, function(a) {
  z <- d[d$genomic_seqid == a, , drop=FALSE]
  data.frame(chromosome = a,
             n_sites = nrow(z),
             n_fully_accessible_23bp = sum(z$fully_accessible_23bp == 1),
             pct_fully_accessible_23bp = round(100 * mean(z$fully_accessible_23bp == 1), 2),
             n_pam_fully_accessible = sum(z$pam_fully_accessible == 1),
             median_fraction_23bp_accessible = median(z$fraction_accessible),
             stringsAsFactors = FALSE)
})
sumtable <- do.call(rbind, rows)
outdir <- "output"
write.csv(sumtable, file.path(outdir, "AR1_accessibility_by_chromosome.csv"), row.names=FALSE)
print(sumtable, row.names=FALSE)
cat("\nAcross 515 evaluable frozen sites:\n")
cat("Fully accessible 23-bp targets:", sum(d$fully_accessible_23bp == 1, na.rm=TRUE), "\n")
cat("UNKN (unresolved):", sum(d$genomic_seqid == "UNKN"), "\n")
# Publication-sized 600 dpi plot, created in R (not Python).
png(filename=file.path(outdir, "AR1_frozen_targets_accessibility_by_arm_600dpi.png"),
    width=7.2, height=4.4, units="in", res=600)
par(mar=c(4.5, 4.8, 1.0, 0.8), family="sans", cex.axis=0.95)
bp <- barplot(sumtable$pct_fully_accessible_23bp,
              names.arg=sumtable$chromosome, ylim=c(0, 100),
              ylab="Fully accessible 23-bp targets (%)", xlab="Chromosome arm",
              col="grey55", border=NA)
text(bp, sumtable$pct_fully_accessible_23bp + 3,
     paste0(sumtable$n_fully_accessible_23bp, "/", sumtable$n_sites), cex=0.8)
box(bty="l")
dev.off()
cat("Wrote base-R summary table and 600-dpi figure to output/.\n")
