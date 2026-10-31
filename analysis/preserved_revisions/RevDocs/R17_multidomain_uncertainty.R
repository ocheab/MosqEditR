# MosqEdit-R | Reviewer 1 comment 7
# Deterministic multi-domain scenario sensitivity + archived PU-resampling envelope.
# Base R only. No simulated confidence probabilities are assigned to scenarios.
# Run from RevDocs: Rscript R17_multidomain_uncertainty.R ../MosqEditR_Phase2_Targeted_Extraction/data_processed R17_results
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else '../MosqEditR_Phase2_Targeted_Extraction/data_processed'
out <- if(length(args)>=2) args[2] else 'R17_results'
dir.create(out,showWarnings=FALSE,recursive=TRUE)
read_in <- function(x) read.csv(file.path(root,x),stringsAsFactors=FALSE,check.names=FALSE)
d <- read_in('09_targeted_integrated_prioritization.csv')
e <- read_in('05_crispr_editability.csv')
g <- read_in('04_preliminary_genomewide_ranking.csv')
u <- read_in('12_bootstrap_rank_uncertainty.csv')
stopifnot(nrow(d)==205L,nrow(e)==205L,nrow(g)==13845L,nrow(u)==205L,
          !anyDuplicated(d$gene_id),!anyDuplicated(g$gene_id))
e <- e[match(d$gene_id,e$gene_id),]; u <- u[match(d$gene_id,u$gene_id),]
stopifnot(identical(d$gene_id,e$gene_id),identical(d$gene_id,u$gene_id))
perc <- function(z) {
 ans <- rep(NA_real_,length(z)); ok <- is.finite(z); n <- sum(ok)
 if(n==1L) ans[ok] <- .5 else if(n>1L) ans[ok] <- (rank(z[ok],ties.method='average')-1)/(n-1)
 ans
}
ranker <- function(z) {
 r <- rep(NA_integer_,length(z)); ok <- which(is.finite(z))
 ii <- ok[order(-z[ok],d$gene_id[ok])];r[ii] <- seq_along(ii);r
}
score3 <- function(b,p,c,w,floor=.8) {
 m <- cbind(b,p,c); obs <- is.finite(m); w <- unname(w)
 num <- rowSums(sweep(ifelse(obs,m,0),2,w,'*'))
 den <- rowSums(sweep(obs*1,2,w,'*'))
 s <- num/den;s[den==0] <- NA_real_
 s*(floor+(1-floor)*rowMeans(obs))
}
full_bio <- function(column) {
 z <- perc(g[[column]])
 z[match(d$gene_id,g$gene_id)]
}
bio <- list(primary=full_bio('score_primary'),
            expression_dominant=full_bio('score_expression_dominant'),
            functional_dominant=full_bio('score_functional_dominant'))
stopifnot(max(abs(bio$primary-d$biological_discovery_percentile),na.rm=TRUE)<1e-12,
          identical(is.na(bio$primary),is.na(d$biological_discovery_percentile)))
components <- c('editability_site_abundance_component','editability_multiplex_component',
                'editability_transcript_coverage_component','editability_gc_balance_component')
x <- as.matrix(e[,components]);storage.mode(x)<-'double'
stopifnot(all(is.finite(x)))
cw <- list(archived=c(.35,.25,.25,.15),supplied_code=c(.40,.25,.20,.15),
           equal=c(.25,.25,.25,.25),abundance_heavy=c(.60,.15,.15,.10))
stopifnot(max(abs(as.vector(x%*%cw$archived)-e$crispr_exon_aware_reference_editability_score))<1e-12)
cp <- lapply(cw,function(w) perc(as.vector(x%*%w)))
# The exact saved score is the baseline so its original ties and precision are maintained.
cp$archived <- perc(e$crispr_exon_aware_reference_editability_score)
stopifnot(max(abs(cp$archived-d$crispr_tractability_percentile))<1e-12)
dw <- list(biology_led=c(.60,.25,.15),primary=c(.50,.30,.20),validation_led=c(.40,.35,.25))
order_candidate <- function(b,p,c,view) {
 scores <- lapply(dw,function(w) score3(b,p,c,w))
 ranks <- lapply(scores,ranker)
 if(view=='consensus') {
   med <- apply(do.call(cbind,ranks),1,median,na.rm=TRUE)
   ix <- order(med,-ifelse(is.finite(scores$primary),scores$primary,-Inf),d$gene_id)
   r <- integer(nrow(d));r[ix] <- seq_len(nrow(d));return(r)
 }
 ranks[[view]]
}
base <- order_candidate(bio$primary,d$pu_consensus_percentile,cp$archived,'consensus')
stopifnot(identical(base,as.integer(d$final_prepopulation_rank)))
cat('PASS R17: Original Step 09 205-gene consensus ranking exactly reproduced.\n')
views <- c('consensus','primary','biology_led','validation_led')
settings <- expand.grid(biology=names(bio),crispr=names(cp),integration=views,stringsAsFactors=FALSE)
mat <- matrix(NA_integer_,nrow(d),nrow(settings))
for(i in seq_len(nrow(settings))) {
 z <- settings[i,]
 mat[,i] <- order_candidate(bio[[z$biology]],d$pu_consensus_percentile,
                             cp[[z$crispr]],z$integration)
}
stopifnot(ncol(mat)==48L,all(colSums(mat<=20)==20),all(colSums(mat<=100)==100))
scen_summary <- settings
scen_summary$scenario_id <- sprintf('S%02d',seq_len(nrow(settings)))
scen_summary$n_rank_changed <- colSums(mat!=base)
scen_summary$top20_retained <- colSums(mat<=20 & base<=20)
scen_summary$top50_retained <- colSums(mat<=50 & base<=50)
scen_summary$top100_retained <- colSums(mat<=100 & base<=100)
scen_summary$median_absolute_rank_shift <- apply(abs(mat-base),2,median)
scen_summary$max_absolute_rank_shift <- apply(abs(mat-base),2,max)
colnames(mat) <- scen_summary$scenario_id
write.csv(scen_summary,file.path(out,'R17_48_scenario_index_and_summary.csv'),row.names=FALSE)
write.csv(data.frame(gene_id=d$gene_id,original_rank=base,mat,check.names=FALSE),
          file.path(out,'R17_205_genes_x_48_ranks.csv'),row.names=FALSE)
lo <- apply(mat,1,min);hi <- apply(mat,1,max)
original_width <- u$bootstrap_rank_q975-u$bootstrap_rank_q025
combined_lo <- pmin(lo,u$bootstrap_rank_q025)
combined_hi <- pmax(hi,u$bootstrap_rank_q975)
comparison <- data.frame(gene_id=d$gene_id,is_external_benchmark=d$is_external_benchmark,
  original_rank=base,PU_resampling_q025=u$bootstrap_rank_q025,
  PU_resampling_q975=u$bootstrap_rank_q975,PU_resampling_width=original_width,
  deterministic_scenario_min_rank=lo,deterministic_scenario_max_rank=hi,
  deterministic_scenario_width=hi-lo,
  descriptive_combined_min_rank=combined_lo,descriptive_combined_max_rank=combined_hi,
  descriptive_combined_width=combined_hi-combined_lo,
  added_envelope_width=combined_hi-combined_lo-original_width,
  scenarios_remaining_top20=rowSums(mat<=20),scenarios_remaining_top50=rowSums(mat<=50),
  scenarios_remaining_top100=rowSums(mat<=100))
# Existing verified GC-window counts, if present, are included as a separate
# guide-availability sensitivity indicator; they are NOT treated as measured
# editing efficiencies or automatically converted to a modified CRISPR score.
gcfile <- 'R14_results/R14_gene_level_GC_sensitivity.csv'
if(file.exists(gcfile)) {
 gc <- read.csv(gcfile,stringsAsFactors=FALSE)
 stopifnot(!anyDuplicated(gc$gene_id),
           all(c('sites_GC_30_70','sites_GC_45_55') %in% names(gc)))
 ii <- match(comparison$gene_id,gc$gene_id)
 comparison$GC30_70_sites <- gc$sites_GC_30_70[ii]
 comparison$GC45_55_sites <- gc$sites_GC_45_55[ii]
 comparison$GC45_55_below_four_sites <- ifelse(
    is.na(comparison$GC45_55_sites),NA,comparison$GC45_55_sites<4)
 cat('GC threshold sensitivity linked for',sum(!is.na(ii)),'genes.\n')
} else {
 cat('GC threshold table not found: no guide-window sensitivity indicators added.\n')
}
write.csv(comparison,file.path(out,'R17_PU_interval_vs_multidomain_envelope.csv'),row.names=FALSE)
# Genome-wide selection stability under the THREE documented Step 04 rankings.
selected <- unique(g$gene_id[is.finite(g$preliminary_rank)&g$preliminary_rank<=200])
stopifnot(length(selected)==200L,all(selected %in% d$gene_id),sum(d$is_external_benchmark)==6)
genome_summary <- do.call(rbind,lapply(c('rank_primary','rank_expression_dominant',
                                          'rank_functional_dominant','preliminary_rank'),function(f) {
 ids <- unique(g$gene_id[is.finite(g[[f]]) & g[[f]]<=200]);
 data.frame(stage04_rank_method=f,genomewide_scored=sum(is.finite(g[[f]])),top200_n=length(ids),
  shared_with_frozen_stage04_top200=length(intersect(ids,selected)),
  changed_membership_vs_frozen_200=200-length(intersect(ids,selected)),
  included_in_targeted_205=length(intersect(ids,d$gene_id)))
}))
stopifnot(all(genome_summary$top200_n==200L))
write.csv(genome_summary,file.path(out,'R17_stage04_genomewide_shortlist_sensitivity.csv'),row.names=FALSE)
summary_stats <- data.frame(
 measure=c('genes','original_PU_resampling_median_interval_width',
           'deterministic_scenario_median_width','combined_descriptive_median_width',
           'median_added_width','max_added_width','original_top20_retained_in_all_48_scenarios',
           'original_top50_retained_in_all_48_scenarios','original_top100_retained_in_all_48_scenarios'),
 value=c(nrow(d),median(original_width),median(hi-lo),median(combined_hi-combined_lo),
         median(combined_hi-combined_lo-original_width),max(combined_hi-combined_lo-original_width),
         sum(base<=20 & rowSums(mat<=20)==48),
         sum(base<=50 & rowSums(mat<=50)==48),
         sum(base<=100 & rowSums(mat<=100)==48)))
write.csv(summary_stats,file.path(out,'R17_summary_metrics.csv'),row.names=FALSE)
cat('\nR17 UNCERTAINTY/SENSITIVITY METRICS\n');print(summary_stats,row.names=FALSE)
cat('\nR17 GENOME-WIDE STAGE04 TOP-200 MEMBERSHIP\n');print(genome_summary,row.names=FALSE)
cat('\nIMPORTANT: These 48 scenario settings are not random draws from a measured error distribution.\n')
cat('The combined envelope is descriptive and is NOT a 95% confidence interval.\n')
cat('No joint expression/orthology/guide/variant error model or population-level CI was inferred.\n')
# Scientific graphics are generated in R, never Python.
plot_rows <- order(base)[seq_len(35)]
png(file.path(out,'R17_uncertainty_envelopes_top35_600dpi.png'),width=4000,height=3200,res=600,type='cairo')
par(mar=c(5,11,4,1),cex.axis=.72)
y <- seq_along(plot_rows)
plot(NA,xlim=c(1,205),ylim=c(35.5,.5),yaxt='n',xlab='Rank among 205 genes (lower = higher)',ylab='',
     main='Archived PU interval versus expanded deterministic rank envelope')
segments(combined_lo[plot_rows],y,combined_hi[plot_rows],y,lwd=3,col='grey78')
segments(u$bootstrap_rank_q025[plot_rows],y,u$bootstrap_rank_q975[plot_rows],y,lwd=3,col='#3474A6')
points(base[plot_rows],y,pch=18,col='#A34C2F',cex=.8)
axis(2,at=y,labels=d$gene_id[plot_rows],las=1,cex.axis=.67)
legend('bottomright',legend=c('Descriptive combined range','Archived PU 2.5–97.5% interval','Frozen rank'),
       lty=c(1,1,NA),pch=c(NA,NA,18),col=c('grey65','#3474A6','#A34C2F'),bty='n',cex=.65)
dev.off()
png(file.path(out,'R17_stage04_shortlist_turnover_600dpi.png'),width=3200,height=2350,res=600,type='cairo')
par(mar=c(9,5,3,1))
barplot(genome_summary$changed_membership_vs_frozen_200,
        names.arg=c('Primary','Expression-led','Functional-led','Original consensus'),
        las=2,ylim=c(0,200),ylab='Genes outside original Stage 04 top 200',
        main='Pre-selection depends on biological evidence weighting')
dev.off()
cat('\nPASS R17: Four CSVs, metrics CSV and two 600-dpi R figures generated in ',out,'\n',sep='')
