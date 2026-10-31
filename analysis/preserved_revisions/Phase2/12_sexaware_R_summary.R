#!/usr/bin/env Rscript
# MosqEdit-R Phase 2: descriptive, sex-aware PASS-SNP overlap statistics.
# Run from MosqEditR_Phase2_Targeted_Extraction, after Python step 11.
# Base R only; no claims of intact haplotypes, indel callability or editing efficacy.
options(stringsAsFactors=FALSE)
infile <- 'output/AR1_sexaware_site_by_population.csv'
stopifnot(file.exists(infile))
d <- read.csv(infile, check.names=FALSE)
required <- c('population','chromosome','n_samples','fully_accessible_23bp',
 'n_observed_PASS_SNP_records','n_complete_GT_at_reported_SNP_loci',
 'n_any_nonref_carriers','n_essential_PAM_GG_nonref_carriers',
 'n_complete_GT_strict_sensitivity','n_any_nonref_strict_sensitivity',
 'n_essential_PAM_GG_nonref_carriers_strict_sensitivity')
stopifnot(all(required %in% names(d)))
stopifnot(length(unique(d$population))==16L)
# Exclude UNKN and incompletely accessible guides in the primary summary.
# Guides with no reported SNPs have no observed-locus GT denominator: do not treat as invariant.
keep <- !is.na(d$fully_accessible_23bp) & d$fully_accessible_23bp==1 &
        !is.na(d$n_observed_PASS_SNP_records) & d$n_observed_PASS_SNP_records>0
p <- d[keep,,drop=FALSE]
stopifnot(nrow(p)==16L*478L)
num <- c('n_complete_GT_at_reported_SNP_loci','n_any_nonref_carriers',
 'n_essential_PAM_GG_nonref_carriers','n_complete_GT_strict_sensitivity',
 'n_any_nonref_strict_sensitivity','n_essential_PAM_GG_nonref_carriers_strict_sensitivity')
for (nm in num) p[[nm]] <- as.numeric(p[[nm]])
summ <- do.call(rbind,lapply(split(p,p$population), function(z) {
  n <- z$n_samples[1]; stopifnot(length(unique(z$n_samples))==1)
  denom <- sum(z$n_complete_GT_at_reported_SNP_loci,na.rm=TRUE)
  strict_denom <- sum(z$n_complete_GT_strict_sensitivity,na.rm=TRUE)
  safediv <- function(x,y) if(y>0) x/y else NA_real_
  data.frame(population=z$population[1],n_samples=n,n_guides= nrow(z),
    guide_sample_observations_called_at_reported_SNPs=denom,
    guide_sample_observations_strict_at_reported_SNPs=strict_denom,
    fraction_guide_sample_observations_any_nonref=safediv(sum(z$n_any_nonref_carriers,na.rm=TRUE),denom),
    fraction_guide_sample_observations_any_nonref_strict=safediv(sum(z$n_any_nonref_strict_sensitivity,na.rm=TRUE),strict_denom),
    fraction_guide_sample_observations_PAM_GG_nonref=safediv(sum(z$n_essential_PAM_GG_nonref_carriers,na.rm=TRUE),denom),
    fraction_guide_sample_observations_PAM_GG_nonref_strict=safediv(sum(z$n_essential_PAM_GG_nonref_carriers_strict_sensitivity,na.rm=TRUE),strict_denom),
    guides_with_any_nonref=sum(z$n_any_nonref_carriers>0,na.rm=TRUE),
    guides_with_PAM_GG_nonref=sum(z$n_essential_PAM_GG_nonref_carriers>0,na.rm=TRUE)
  )
}))
rownames(summ)<-NULL
summ <- summ[order(summ$population),,drop=FALSE]
out <- 'output/AR1_sexaware_descriptive_population_summary.csv'
write.csv(summ,out,row.names=FALSE,na='')
print(summ[,c('population','n_samples','guides_with_any_nonref',
    'guides_with_PAM_GG_nonref','fraction_guide_sample_observations_PAM_GG_nonref',
    'fraction_guide_sample_observations_PAM_GG_nonref_strict')],row.names=FALSE,digits=5)
# Focused, publication-resolution plot: no error bars; units correlated by guide.
figname <- 'output/AR1_sexaware_PAM_GG_observed_nonref_by_population_600dpi.png'
png(figname,width=9,height=5.6,units='in',res=600,pointsize=10)
op <- par(mar=c(6.6,4.8,2.5,1.2))
ord <- order(summ$fraction_guide_sample_observations_PAM_GG_nonref)
xx<-seq_along(ord)
y1 <- 100*summ$fraction_guide_sample_observations_PAM_GG_nonref[ord]
y2 <- 100*summ$fraction_guide_sample_observations_PAM_GG_nonref_strict[ord]
lim <- max(c(y1,y2),na.rm=TRUE)*1.17
plot(xx,y1,type='n',xaxt='n',xlab='',ylab='Guide-sample observations with PAM-GG non-reference (%)',
     main='Ag1000G Phase 2: observed PASS-SNP variation at essential PAM G positions',
     ylim=c(0,lim),xlim=c(0.5,length(xx)+0.5),las=1)
abline(h=pretty(c(0,lim)),col='grey92',lty=1)
segments(xx,y1,xx,y2,col='grey55',lwd=1)
points(xx,y1,pch=16,cex=1.15,col='#1f4e79')
points(xx,y2,pch=1,cex=1.15,col='#a15a3b')
axis(1,at=xx,labels=paste0(summ$population[ord],'\n(n=',summ$n_samples[ord],')'),
     las=2,cex.axis=0.70,tick=FALSE)
legend('topright',c('Primary: complete GT','Sensitivity: GQ >=20, DP >=10'),
       pch=c(16,1),col=c('#1f4e79','#a15a3b'),bty='n',cex=0.82)
mtext('Descriptive PASS-SNP results at reported loci only; site-sample observations are correlated.',
      side=1,line=5.4,cex=0.68)
par(op)
dev.off()
cat('\nWrote:',out,'\n',figname,'\n')
cat('Caution: No full 23-bp genotype callability, indel validation, haplotype conservation, or CRISPR editing efficacy inferred.\n')
