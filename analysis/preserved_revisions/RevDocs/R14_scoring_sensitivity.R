# MosqEdit-R R1.4 / R2.Q5
# CRISPR component-weight and integrated-ranking sensitivity
# Base R only; original source files remain unchanged.

root <- "../MosqEditR_Phase2_Targeted_Extraction/data_processed"
outdir <- "R14_results"
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

e <- read.csv(file.path(root,"05_crispr_editability.csv"))
d <- read.csv(file.path(root,"09_targeted_integrated_prioritization.csv"))

stopifnot(nrow(d)==205, nrow(e)==205,
          !anyDuplicated(d$gene_id),
          !anyDuplicated(e$gene_id))

e <- e[match(d$gene_id,e$gene_id),]
stopifnot(identical(d$gene_id,e$gene_id))

components <- c(
  "editability_site_abundance_component",
  "editability_multiplex_component",
  "editability_transcript_coverage_component",
  "editability_gc_balance_component"
)

x <- as.matrix(e[,components])
storage.mode(x) <- "double"
stopifnot(all(is.finite(x)))

# Archived output and selected alternative component weights
component_weights <- list(
  archived = c(.35,.25,.25,.15),
  supplied_code = c(.40,.25,.20,.15),
  equal_components = c(.25,.25,.25,.25),
  abundance_heavy = c(.60,.15,.15,.10),
  abundance_removed = c(0,.35,.35,.30)
)

domain_weights <- list(
  biology_led = c(.60,.25,.15),
  primary = c(.50,.30,.20),
  validation_led = c(.40,.35,.25)
)

rank_genes <- function(s) {
  stopifnot(all(is.finite(s)))
  order_idx <- order(-s,d$gene_id)
  r <- integer(length(s))
  r[order_idx] <- seq_along(s)
  r
}

domain_score <- function(crispr_pct,w) {
  m <- cbind(
    d$biological_discovery_percentile,
    d$pu_consensus_percentile,
    crispr_pct
  )
  valid <- is.finite(m)
  ww <- unname(w)
  num <- rowSums(sweep(ifelse(valid,m,0),2,ww,"*"))
  den <- rowSums(sweep(valid*1,2,ww,"*"))
  z <- num/den
  z[den==0] <- NA_real_
  z*(.8+.2*rowMeans(valid))
}

recalculate <- function(name,w) {
  score <- as.vector(x %*% w)

  # Preserve the exact archived score for baseline ranking;
  # formula reconstruction is checked independently below.
  if(name=="archived") {
    score <- e$crispr_exon_aware_reference_editability_score
  }

  crispr_pct <- (rank(score,ties.method="average")-1) /
                (length(score)-1)

  scenario_scores <- lapply(
    domain_weights,
    function(z) domain_score(crispr_pct,z)
  )

  scenario_ranks <- lapply(
    scenario_scores,
    rank_genes
  )

  median_rank <- apply(
    do.call(cbind,scenario_ranks),
    1,median
  )

  final_order <- order(
    median_rank,
    -scenario_scores$primary,
    d$gene_id
  )

  final_rank <- integer(nrow(d))
  final_rank[final_order] <- seq_len(nrow(d))

  list(
    score=score,
    crispr_pct=crispr_pct,
    scenario_scores=scenario_scores,
    scenario_ranks=scenario_ranks,
    final_rank=final_rank
  )
}

results <- lapply(
  names(component_weights),
  function(n) recalculate(n,component_weights[[n]])
)
names(results) <- names(component_weights)

original_formula <- as.vector(
  x %*% component_weights$archived
)

stopifnot(
  max(abs(original_formula -
    e$crispr_exon_aware_reference_editability_score)) < 1e-12,
  max(abs(results$archived$crispr_pct -
    d$crispr_tractability_percentile)) < 1e-12,
  max(abs(results$archived$scenario_scores$primary -
    d$final_integrated_score)) < 1e-12,
  identical(
    results$archived$scenario_ranks$primary,
    as.integer(d$final_integrated_rank)
  ),
  identical(
    results$archived$final_rank,
    as.integer(d$final_prepopulation_rank)
  )
)

cat("PASS: Archived CRISPR scores reproduced.\n")
cat("PASS: Original integrated scores and 205 consensus ranks reproduced.\n")

baseline <- results$archived$final_rank

summary <- do.call(rbind,lapply(names(results),function(n) {
  v <- results[[n]]
  data.frame(
    scenario=n,
    weights=paste(component_weights[[n]],collapse="/"),
    distinct_CRISPR_scores=length(unique(round(v$score,12))),
    genes_with_changed_final_rank=sum(v$final_rank != baseline),
    original_top20_retained=sum(
      v$final_rank <= 20 & baseline <= 20
    ),
    original_top50_retained=sum(
      v$final_rank <= 50 & baseline <= 50
    ),
    original_top100_retained=sum(
      v$final_rank <= 100 & baseline <= 100
    )
  )
}))

cat("\nSCORING SENSITIVITY SUMMARY:\n")
print(summary,row.names=FALSE)

details <- data.frame(
  gene_id=d$gene_id,
  original_rank=baseline,
  original_CRISPR_score=results$archived$score,
  original_CRISPR_percentile=results$archived$crispr_pct
)

for(n in names(results)) {
  details[[paste0(n,"_score")]] <- results[[n]]$score
  details[[paste0(n,"_final_rank")]] <- results[[n]]$final_rank
}

cat("\nGENES CHANGING FINAL RANK UNDER SUPPLIED CODE:\n")
changed <- subset(
  details,
  supplied_code_final_rank != original_rank
)

print(changed[,c(
  "gene_id","original_rank",
  "supplied_code_final_rank"
)],row.names=FALSE)

component_qc <- data.frame(
  component=components,
  min=apply(x,2,min),
  max=apply(x,2,max),
  distinct_values=apply(x,2,function(z)
    length(unique(z)))
)

cat("\nCOMPONENT VARIABILITY:\n")
print(component_qc,row.names=FALSE)

write.csv(
  summary,
  file.path(outdir,"R14_component_weight_sensitivity_summary.csv"),
  row.names=FALSE
)

write.csv(
  details,
  file.path(outdir,"R14_component_weight_sensitivity_205_genes.csv"),
  row.names=FALSE
)

write.csv(
  component_qc,
  file.path(outdir,"R14_component_variability_QC.csv"),
  row.names=FALSE
)

# R-generated publication-resolution diagnostic
png(
  file.path(outdir,"R14_component_weight_rank_sensitivity_600dpi.png"),
  width=3600,height=2700,res=600,type="cairo"
)

par(mar=c(5.5,5.5,3.5,1.5))

plot(
  c(1,205),c(1,205),
  type="n",
  xlim=c(1,205),ylim=c(205,1),
  xlab="Original consensus rank",
  ylab="Alternative component-weight consensus rank",
  main="CRISPR component-weight sensitivity"
)

abline(a=0,b=1,lty=2,col="grey60")

scenarios <- c(
  "supplied_code",
  "equal_components",
  "abundance_heavy",
  "abundance_removed"
)

cols <- c("#236A91","#C28A35","#397C5C","#9A4E6C")

for(i in seq_along(scenarios)) {
  points(
    baseline,
    results[[scenarios[i]]]$final_rank,
    pch=i,
    cex=.65,
    col=cols[i]
  )
}

legend(
  "bottomleft",
  legend=gsub("_"," ",scenarios),
  col=cols,pch=seq_along(scenarios),
  bty="n",cex=.75
)

dev.off()

cat("\nPASS: Three CSV tables and one 600-dpi R figure generated.\n")
cat("CAUTION: Abundance-removed scoring is a degenerate stress test.\n")
