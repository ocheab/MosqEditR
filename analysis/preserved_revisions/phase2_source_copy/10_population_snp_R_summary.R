# MosqEdit-R Ag1000G Phase 2: descriptive summary of frozen-panel PASS SNP data.
# Base R, no packages; generates a 600-dpi PNG. This is NOT a guide haplotype
# conservation analysis. Execute from the kit root after 09_population_snp_metrics.py.
out <- 'output'
s <- read.csv(file.path(out,'AR1_SNP_site_summary_530.csv'), stringsAsFactors=FALSE)
p <- read.csv(file.path(out,'AR1_SNP_site_by_population.csv'), stringsAsFactors=FALSE)
expected <- c('2L','2R','3L','3R','X')
stopifnot(nrow(s)==530, nrow(p)>=530,
          all(expected %in% s$chromosome),
          sum(s$chromosome=='UNKN')==15,
          sum(s$fully_accessible_23bp==1,na.rm=TRUE)==491)
rows <- lapply(expected,function(arm) {
  x<-s[s$chromosome==arm,]
  data.frame(chromosome=arm, n_frozen_guide_sites=nrow(x),
    n_fully_accessible=sum(x$fully_accessible_23bp==1,na.rm=TRUE),
    n_sites_with_reported_PASS_SNP=sum(x$n_observed_PASS_SNP_records>0,na.rm=TRUE),
    n_sites_without_reported_PASS_SNP=sum(x$n_observed_PASS_SNP_records==0,na.rm=TRUE),
    n_sites_with_essential_PAM_GG_SNP_position=sum(x$n_essential_PAM_GG_SNP_positions>0,na.rm=TRUE),
    n_sites_with_protospacer_PAM_proximal_10_SNP_position=sum(x$n_protospacer_PAM_proximal10_SNP_positions>0,na.rm=TRUE))
})
tab <- do.call(rbind,rows)
write.csv(tab,file.path(out,'AR1_SNP_summary_by_chromosome.csv'),row.names=FALSE)
cat('Observed PASS SNP overlap with frozen 23-bp guides (not conservation):\n')
print(tab,row.names=FALSE)

# Aggregation denominator is *site-by-mosquito assessments at observed PASS SNP
# loci*, not the number of unique mosquitoes or full-23bp callable haplotypes.
# Guide sites without SNP rows and inaccessible sites are deliberately omitted.
use <- p[!is.na(p$n_observed_PASS_SNP_records) & p$n_observed_PASS_SNP_records>0 &
         !is.na(p$fully_accessible_23bp) & p$fully_accessible_23bp==1 &
         !is.na(p$n_complete_GT_at_reported_SNP_loci), ]
by_pop <- split(use,use$population)
ptab <- do.call(rbind,lapply(names(by_pop),function(pop) {
  d <- by_pop[[pop]]
  denominator <- sum(d$n_complete_GT_at_reported_SNP_loci,na.rm=TRUE)
  numerator <- sum(d$n_any_nonref_carriers,na.rm=TRUE)
  essential <- sum(d$n_essential_PAM_GG_nonref_carriers,na.rm=TRUE)
  denom_strict <- sum(d$n_complete_GT_strict_sensitivity,na.rm=TRUE)
  num_strict <- sum(d$n_any_nonref_strict_sensitivity,na.rm=TRUE)
  data.frame(population=pop,country=paste(unique(d$country),collapse=';'),
    n_snp_bearing_fully_accessible_guides=nrow(d),
    n_complete_site_sample_assessments=denominator,
    n_any_nonref_site_sample_assessments=numerator,
    fraction_any_nonref_site_sample_assessments=if(denominator) numerator/denominator else NA_real_,
    n_essential_PAM_GG_nonref_site_sample_assessments=essential,
    fraction_PAM_GG_nonref_site_sample_assessments=if(denominator) essential/denominator else NA_real_,
    n_complete_site_sample_assessments_GQ_DP_sensitivity=denom_strict,
    fraction_any_nonref_GQ_DP_sensitivity=if(denom_strict)num_strict/denom_strict else NA_real_)
}))
rownames(ptab)<-NULL
write.csv(ptab,file.path(out,'AR1_SNP_descriptive_by_population.csv'),row.names=FALSE)
cat('\nPopulation groups: ',nrow(ptab),'\n',sep='')
cat('Primary and strict sensitivity denominator definitions are in README.\n')

# Publication-resolution visual: only observed SNP positions, not inferred
# functional efficiency or callable whole-guide haplotypes.
png(file.path(out,'AR1_PASS_SNP_overlap_frozen_guides_600dpi.png'),width=8.5,height=5.4,
    units='in',res=600)
par(mar=c(5.0,4.7,3.4,1.2),family='sans')
mat <- rbind(tab$n_sites_with_reported_PASS_SNP,
             tab$n_sites_without_reported_PASS_SNP)
colnames(mat) <- tab$chromosome
barplot(mat,beside=FALSE,col=c('#4472C4','#B7C2D0'),border=NA,
        ylab='Number of frozen 23-bp guide sites',xlab='Chromosome arm',
        main='Observed Phase 2 PASS SNP overlap in frozen guide sites',
        ylim=c(0,max(colSums(mat))*1.17))
legend('topright',legend=c('At least one observed PASS SNP','No reported PASS SNP'),
       fill=c('#4472C4','#B7C2D0'),bty='n',cex=0.82)
mtext('Absence of a PASS SNP is NOT evidence of site invariance.',side=1,line=4,cex=0.73)
dev.off()
cat('Wrote output/AR1_SNP_summary_by_chromosome.csv, output/AR1_SNP_descriptive_by_population.csv,\n')
cat('and output/AR1_PASS_SNP_overlap_frozen_guides_600dpi.png\n')
