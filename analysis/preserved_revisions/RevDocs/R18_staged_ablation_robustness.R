# MosqEdit-R | Reviewer 1 comment 8
# Reproduce original seven ablations; add absolute shifts and renormalization/floor sensitivity.
# Base R only. Analysis is CONDITIONAL on the pre-selected 205-gene panel.
# Run from RevDocs: Rscript R18_staged_ablation_robustness.R ../MosqEditR_Phase2_Targeted_Extraction/data_processed R18_results
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else '../MosqEditR_Phase2_Targeted_Extraction/data_processed'
out <- if(length(args)>=2) args[2] else 'R18_results'
dir.create(out,showWarnings=FALSE,recursive=TRUE)
read_in <- function(x) read.csv(file.path(root,x),stringsAsFactors=FALSE,check.names=FALSE)
d <- read_in('09_targeted_integrated_prioritization.csv')
a <- read_in('15_gene_level_ablation.csv')
g <- read_in('04_preliminary_genomewide_ranking.csv')
stopifnot(nrow(d)==205L,nrow(a)==1435L,nrow(g)==13845L,
          !anyDuplicated(d$gene_id),!anyDuplicated(g$gene_id))
components <- cbind(biology=d$biological_discovery_percentile,
                    pu=d$pu_consensus_percentile,crispr=d$crispr_tractability_percentile)
w <- c(biology=.5,pu=.3,crispr=.2)
ablations <- list(full=names(w),no_biology=c('pu','crispr'),
                  no_pu=c('biology','crispr'),no_crispr=c('biology','pu'),
                  biology_only='biology',pu_only='pu',crispr_only='crispr')
ranker <- function(z) {
 ans <- rep(NA_integer_,length(z));ok <- which(is.finite(z))
 oo <- ok[order(-z[ok],d$gene_id[ok])]
 ans[oo] <- seq_along(oo);ans
}
score <- function(keep,floor=.8,rule='available') {
 x <- components[,keep,drop=FALSE];ok <- is.finite(x)
 ww <- w[keep]
 numer <- rowSums(sweep(ifelse(ok,x,0),2,ww,'*'))
 den <- switch(rule,
   available=rowSums(sweep(ok*1,2,ww,'*')),
   zero_fill=rep(sum(ww),nrow(x)),
   original_denominator=rep(sum(w),nrow(x)),
   stop('Unknown rule'))
 v <- numer/den;v[den<=0] <- NA_real_
 v*(floor+(1-floor)*rowMeans(ok))
}
base <- as.integer(d$final_prepopulation_rank)
stopifnot(identical(sort(base),seq_len(nrow(d))))
# Exact reproduction of original Step 15 scores/ranks, including NA for missing biology.
original_rows <- list()
for(label in names(ablations)) {
 s <- score(ablations[[label]],floor=.8,rule='available')
 r <- ranker(s)
 archived <- a[a$ablation==label,]
 archived <- archived[match(d$gene_id,archived$gene_id),]
 stopifnot(nrow(archived)==205L,identical(archived$gene_id,d$gene_id),
           identical(is.na(s),is.na(archived$ablation_score)),
           max(abs(s-archived$ablation_score),na.rm=TRUE)<1e-12,
           identical(r,as.integer(archived$ablation_rank)))
 original_rows[[label]] <- data.frame(gene_id=d$gene_id,ablation=label,
   archived_full_consensus_rank=base,ablation_primary_rank=r,
   rank_shift=r-base,absolute_rank_shift=abs(r-base),
   is_external_benchmark=d$is_external_benchmark)
}
cat('PASS R18: All 1,435 original Step 15 ablation scores and ranks reproduced.\n')
original <- do.call(rbind,original_rows)
write.csv(original,file.path(out,'R18_original_ablation_1435_absolute_shifts.csv'),row.names=FALSE)
overlap <- function(r,k) sum(base<=k & r<=k,na.rm=TRUE)
metrics <- function(name,rule,floor,s,r) {
 delta <- r-base; z <- abs(delta)
 data.frame(ablation=name,renormalization=rule,completeness_floor=floor,
   genes_scoreable=sum(is.finite(s)),n_rank_changed=sum(delta!=0,na.rm=TRUE),
   n_shift_at_least_10=sum(z>=10,na.rm=TRUE),
   n_shift_at_least_25=sum(z>=25,na.rm=TRUE),
   median_absolute_rank_shift=median(z,na.rm=TRUE),
   maximum_absolute_rank_shift=max(z,na.rm=TRUE),
   n_in_original_top10_retained=overlap(r,10),
   n_in_original_top20_retained=overlap(r,20),
   n_in_original_top50_retained=overlap(r,50),
   n_in_original_top100_retained=overlap(r,100),
   spearman_vs_frozen=suppressWarnings(cor(base,r,method='spearman',use='complete.obs')),
   benchmark_median_rank=median(r[d$is_external_benchmark],na.rm=TRUE))
}
# Factorial stress test: 7 ablations x 3 floors x 3 renormalization rules = 63 conditions.
scenarios <- vector('list',63);counter <- 0L
for(floor in c(.6,.8,1.0)) for(rule in c('available','zero_fill','original_denominator')) {
 for(nm in names(ablations)) {
   counter <- counter+1L
   s <- score(ablations[[nm]],floor,rule)
   r <- ranker(s)
   scenarios[[counter]] <- metrics(nm,rule,floor,s,r)
 }
}
stopifnot(counter==63L)
tab <- do.call(rbind,scenarios)
write.csv(tab,file.path(out,'R18_63_ablation_floor_renormalization_scenarios.csv'),row.names=FALSE)
base7 <- tab[tab$renormalization=='available' & abs(tab$completeness_floor-.8)<1e-10,]
base7 <- base7[match(names(ablations),base7$ablation),]
write.csv(base7,file.path(out,'R18_seven_baseline_absolute_ablation_metrics.csv'),row.names=FALSE)
# Distinguish primary-score 'full' ablation from Step09 three-scenario CONSENSUS final rank.
full_shift <- subset(original,ablation=='full')
cat('Original Step15 full-primary versus Step09 consensus: ',
    sum(full_shift$absolute_rank_shift>0,na.rm=TRUE),' positions changed.\n',sep='')
# Stage 04 is genome-wide: the original 205 panel combines top 200 + five added landmarks.
selected <- unique(g$gene_id[is.finite(g$preliminary_rank) & g$preliminary_rank<=200])
stopifnot(length(selected)==200L,all(selected %in% d$gene_id),
          length(setdiff(d$gene_id,selected))==5L,
          sum(d$is_external_benchmark)==6L)
preselect <- data.frame(gene_id=d$gene_id,
  panel_final_rank=base,
  step04_genomewide_preliminary_rank=g$preliminary_rank[match(d$gene_id,g$gene_id)],
  entered_by_stage04_top200=d$gene_id %in% selected,
  is_external_benchmark=d$is_external_benchmark)
write.csv(preselect,file.path(out,'R18_stage04_preselection_audit_205.csv'),row.names=FALSE)
cat('\nR18 ORIGINAL SEVEN ABLATIONS — ABSOLUTE EFFECTS\n')
print(base7[,c('ablation','n_rank_changed','n_shift_at_least_10',
              'median_absolute_rank_shift','maximum_absolute_rank_shift',
              'n_in_original_top20_retained','n_in_original_top100_retained')],row.names=FALSE)
cat('\nR18 STAGE 04 PRESELECTION\n')
cat('Genome-wide annotated genes=',nrow(g),'\n',
    'Genome-wide scoreable preliminary genes=',sum(is.finite(g$preliminary_rank)),'\n',
    'Top-200 selected=',sum(preselect$entered_by_stage04_top200),'\n',
    'Added benchmarks=',sum(!preselect$entered_by_stage04_top200 & preselect$is_external_benchmark),'\n',sep='')
# Derive actual robustness to completeness factor & renormalization rather than assuming effects.
summary63 <- do.call(rbind,lapply(names(ablations),function(nm) {
 v <- tab[tab$ablation==nm,]
 data.frame(ablation=nm,min_top20_retained=min(v$n_in_original_top20_retained),
   max_top20_retained=max(v$n_in_original_top20_retained),
   min_top100_retained=min(v$n_in_original_top100_retained),
   max_top100_retained=max(v$n_in_original_top100_retained),
   min_n_rank_changed=min(v$n_rank_changed),max_n_rank_changed=max(v$n_rank_changed))
}))
write.csv(summary63,file.path(out,'R18_robustness_range_by_ablation.csv'),row.names=FALSE)
cat('\nR18 63-SCENARIO ROBUSTNESS RANGES\n');print(summary63,row.names=FALSE)
# Publication-resolution figures generated in R.
png(file.path(out,'R18_absolute_rank_shifts_600dpi.png'),width=3600,height=2500,res=600,type='cairo')
par(mar=c(9,5,3,1))
z <- base7[base7$ablation!='full',]
bp <- barplot(z$n_shift_at_least_10,names.arg=gsub('_',' ',z$ablation),las=2,
       ylab='Genes shifting by 10 or more ranks',ylim=c(0,210),
       main='Absolute rank changes after domain ablation')
text(bp,z$n_shift_at_least_10+7,labels=z$n_shift_at_least_10,cex=.8)
dev.off()
png(file.path(out,'R18_ablation_topk_robustness_600dpi.png'),width=3700,height=2500,res=600,type='cairo')
par(mar=c(9,5,3,1))
z <- summary63[summary63$ablation!='full',]
plot(seq_len(nrow(z)),z$min_top20_retained,ylim=c(0,20),xlim=c(.5,nrow(z)+.5),
     xaxt='n',xlab='',ylab='Original top-20 genes retained (out of 20)',
     pch=16,cex=1,main='Ablation robustness across floors and score normalization')
segments(seq_len(nrow(z)),z$min_top20_retained,seq_len(nrow(z)),z$max_top20_retained,lwd=5,col='grey60')
points(seq_len(nrow(z)),z$max_top20_retained,pch=16)
axis(1,at=seq_len(nrow(z)),labels=gsub('_',' ',z$ablation),las=2,cex.axis=.75)
dev.off()
cat('\nPASS R18: Five CSV outputs and two 600-dpi R figures generated in ',out,'\n',sep='')
cat('LIMIT: all Step15 ablations are conditional on the Step04-selected 205-gene panel.\n')
cat('Whole-genome CRISPR/PU domain ablation cannot be claimed without scores for unselected genes.\n')
