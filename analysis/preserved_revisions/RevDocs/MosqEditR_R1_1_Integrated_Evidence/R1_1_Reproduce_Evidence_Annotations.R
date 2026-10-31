#!/usr/bin/env Rscript
# MosqEdit-R R1.1 evidence annotations (base R). Run from the original extraction directory.
# INPUT: frozen 530-site manifest and previously generated sex-aware CSVs.
# Frozen gene order/ranks are not changed. NO interpretation as complete guide conservation.
options(stringsAsFactors = FALSE)
root <- normalizePath('.', mustWork=TRUE)
out <- file.path(root,'output')
req <- c('13_population_site_manifest_FROZEN.csv',
  'output/AR1_sexaware_site_summary_530.csv',
  'output/AR1_sexaware_site_by_population.csv',
  'output/AR1_sexaware_descriptive_population_summary.csv',
  'output/AR1_frozen_site_accessibility_530.csv',
  'output/AR1_frozen_site_FULL23bp_FASTA_QC_530.csv')
stopifnot(all(file.exists(file.path(root,req))))
manifest <- read.csv(file.path(root,req[1]),check.names=FALSE)
ss <- read.csv(file.path(root,req[2]),check.names=FALSE)
sp <- read.csv(file.path(root,req[3]),check.names=FALSE)
pp <- read.csv(file.path(root,req[4]),check.names=FALSE)
ac <- read.csv(file.path(root,req[5]),check.names=FALSE)
ref <- read.csv(file.path(root,req[6]),check.names=FALSE)
stopifnot(nrow(manifest)==530L,nrow(ss)==530L,length(unique(manifest$gene_id))==106L,
  all(table(manifest$gene_id)==5L),nrow(sp)==8480L,length(unique(sp$population))==16L,
  length(unique(manifest$population_site_id))==530L)
stopifnot(setequal(manifest$population_site_id,ss$population_site_id),
  setequal(manifest$population_site_id,ac$population_site_id),
  setequal(manifest$population_site_id,ref$population_site_id))
stopifnot(sum(ss$fully_accessible_23bp==1,na.rm=TRUE)==491L,
  sum(ss$n_observed_PASS_SNP_records>0,na.rm=TRUE)==479L,
  sum(ss$fully_accessible_23bp==1 & ss$n_observed_PASS_SNP_records>0,na.rm=TRUE)==478L,
  sum(ss$fully_accessible_23bp==1 & ss$n_essential_PAM_GG_SNP_positions>0,na.rm=TRUE)==220L)
# The sex-aware SNP-by-population file contains 16 rows per guide; undetected SNPs are NA.
ids <- unique(manifest$population_site_id)
sp_split <- split(sp,sp$population_site_id)
sums <- do.call(rbind,lapply(ids,function(id){
 z<-sp_split[[id]]
 cols <- c('n_any_nonref_carriers','n_essential_PAM_GG_nonref_carriers')
 if(is.null(z)) stop('Missing site population: ',id)
 if(all(is.na(z$n_observed_PASS_SNP_records))) {
  return(data.frame(population_site_id=id,
    total_nonref_carrier_guide_population_assessments=NA_real_,
    total_pam_GG_nonref_carrier_guide_population_assessments=NA_real_,
    n_populations_any_nonref=NA_real_, n_populations_PAM_GG_nonref=NA_real_,
    max_observed_population_PAM_GG_nonref_AF=NA_real_))
 }
 if(z$n_observed_PASS_SNP_records[1]==0) {
  return(data.frame(population_site_id=id,
    total_nonref_carrier_guide_population_assessments=NA_real_,
    total_pam_GG_nonref_carrier_guide_population_assessments=NA_real_,
    n_populations_any_nonref=NA_real_, n_populations_PAM_GG_nonref=NA_real_,
    max_observed_population_PAM_GG_nonref_AF=NA_real_))
 }
 freq <- z$max_PAM_essential_GG_nonref_AF_at_reported_SNP_loci
 data.frame(population_site_id=id,
  total_nonref_carrier_guide_population_assessments=sum(z$n_any_nonref_carriers,na.rm=TRUE),
  total_pam_GG_nonref_carrier_guide_population_assessments=sum(z$n_essential_PAM_GG_nonref_carriers,na.rm=TRUE),
  n_populations_any_nonref=sum(z$n_any_nonref_carriers>0,na.rm=TRUE),
  n_populations_PAM_GG_nonref=sum(z$n_essential_PAM_GG_nonref_carriers>0,na.rm=TRUE),
  max_observed_population_PAM_GG_nonref_AF=if(all(is.na(freq))) NA_real_ else max(freq,na.rm=TRUE))
}))
master<-merge(manifest,ss,by='population_site_id',all.x=TRUE,sort=FALSE,suffixes=c('_frozen','_qc'))
master<-merge(master,sums,by='population_site_id',all.x=TRUE,sort=FALSE)
master<-merge(master,ac[c('population_site_id','fraction_accessible','pam_fully_accessible')],by='population_site_id',all.x=TRUE,sort=FALSE)
master<-merge(master,ref[c('population_site_id','qc_status','n_matching_bases')],by='population_site_id',all.x=TRUE,sort=FALSE)
master<-master[match(ids,master$population_site_id),]
stopifnot(identical(as.character(master$population_site_id),as.character(ids)))
master$primary_guide_status <- ifelse(master$genomic_seqid=='UNKN','UNKN_COORDINATE_UNRESOLVED',
  ifelse(master$fully_accessible_23bp==0,'INCOMPLETELY_ACCESSIBLE',
  ifelse(master$n_observed_PASS_SNP_records==0,'NO_PASS_SNP_RECORDED_NOT_PROOF_OF_CONSERVATION',
  ifelse(master$total_pam_GG_nonref_carrier_guide_population_assessments>0,
  'PASS_SNP_PAM_GG_ALT_OBSERVED','PASS_SNP_NO_PAM_GG_ALT_OBSERVED_NOT_CONSERVATION'))))
write.csv(master,file.path(out,'R1_1_frozen_530_master_R_reproduction.csv'),row.names=FALSE,na='')
# R summary verifies the original population-level ratios; exact guide names remain fixed.
keep<-sp$fully_accessible_23bp==1 & !is.na(sp$n_observed_PASS_SNP_records) & sp$n_observed_PASS_SNP_records>0
keep[is.na(keep)]<-FALSE
x<-sp[keep,]
stopifnot(nrow(x)==16L*478L)
for (pop in unique(x$population)){
 z<-x[x$population==pop,]; s<-pp[pp$population==pop,]
 rate<-sum(z$n_essential_PAM_GG_nonref_carriers)/sum(z$n_complete_GT_at_reported_SNP_loci)
 stopifnot(length(s$population)==1L,abs(rate-s$fraction_guide_sample_observations_PAM_GG_nonref)<1e-7)
}
cat('R evidence annotations generated; 530 guides and 16 population summaries verified.\n')
cat('NOTE: old rankings retained. SNP overlap != conservation.\n')
