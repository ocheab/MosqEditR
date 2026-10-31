# MosqEdit-R frozen guide vs Phase 2 PASS VCF REF concordance QC.
# Base R only; deliberately *not* a complete genome reference/coordinate validation.
# Run in folder containing 13_population_site_manifest_FROZEN.csv and output/VCFs:
#   Rscript 07_vcf_ref_concordance.R

manifest_path <- '13_population_site_manifest_FROZEN.csv'
vcf_ref_path <- 'output/AR1_PASS_SNP_REF_5arms.tsv'
outdir <- 'output'
expected_sha <- '04687e2ea261caecee18065dcd2523313aca345ac865eb0483261976f703d85c'
if (!file.exists(manifest_path) || !file.exists(vcf_ref_path))
  stop('Missing manifest or VCF REF export; run 06_export_vcf_ref.sh first.')
sha <- system2('sha256sum', shQuote(manifest_path), stdout = TRUE)
sha <- strsplit(sha[1], '[[:space:]]+')[[1]][1]
if (!identical(sha, expected_sha)) stop('Frozen manifest SHA256 mismatch; STOP.')
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
v <- read.delim(vcf_ref_path, stringsAsFactors = FALSE, check.names = FALSE)
need_m <- c('population_site_id', 'genomic_seqid','genomic_start','genomic_end',
            'guide_genomic_strand','protospacer_20nt','pam')
need_v <- c('chromosome','pos','ref','alt')
stopifnot(all(need_m %in% names(manifest)), all(need_v %in% names(v)),
          nrow(manifest) == 530L, !anyDuplicated(manifest$population_site_id),
          nrow(v) == 2006L)
stopifnot(all(nchar(manifest$protospacer_20nt)==20),all(nchar(manifest$pam)==3),
          all(manifest$genomic_end-manifest$genomic_start+1==23),
          all(manifest$guide_genomic_strand %in% c('+','-')))
expected_per_arm <- c('2L'=727L,'2R'=395L,'3L'=239L,'3R'=576L,'X'=69L)
stopifnot(identical(as.integer(table(factor(v$chromosome, levels=names(expected_per_arm)))),
                    as.integer(expected_per_arm)))
comp <- function(s) chartr('ACGT','TGCA',s)
revcomp <- function(s) paste(rev(strsplit(comp(toupper(s)), '', fixed=TRUE)[[1]]),collapse='')
# 23-mer supplied in guide orientation. Genomic minus strand requires reverse complement.
guide <- paste0(toupper(manifest$protospacer_20nt),toupper(manifest$pam))
genomic_seq <- ifelse(manifest$guide_genomic_strand=='+', guide,
                       vapply(guide, revcomp, character(1)))
stopifnot(all(nchar(genomic_seq)==23))
idx_arms <- split(seq_len(nrow(v)), v$chromosome)
site_rows <- vector('list',nrow(manifest))
base_rows <- list(); b <- 0L
for (i in seq_len(nrow(manifest))) {
  r <- manifest[i, ]
  arm <- r$genomic_seqid
  if (arm == 'UNKN') {
    site_rows[[i]] <- data.frame(population_site_id=r$population_site_id,
        gene_id=r$gene_id, chromosome=arm, start=r$genomic_start,end=r$genomic_end,
        strand=r$guide_genomic_strand,n_snp_records_in_23bp=NA_integer_,
        n_ref_matches=NA_integer_,n_ref_mismatches=NA_integer_,
        qc_status='UNKN_NOT_EVALUABLE', stringsAsFactors=FALSE)
    next
  }
  ii <- idx_arms[[arm]]
  hits <- ii[v$pos[ii]>=r$genomic_start & v$pos[ii]<=r$genomic_end]
  # This script evaluates only single-base REF SNP records within the 23 bp interval.
  hits <- hits[nchar(v$ref[hits])==1L]
  matches <- logical(length(hits))
  for (j in seq_along(hits)) {
    k <- hits[j]; off <- v$pos[k]-r$genomic_start+1L
    expected_base <- substr(genomic_seq[i],off,off)
    observed_base <- toupper(v$ref[k]); matches[j] <- identical(expected_base,observed_base)
    b <- b+1L
    base_rows[[b]] <- data.frame(population_site_id=r$population_site_id,
      gene_id=r$gene_id, chromosome=arm, pos=v$pos[k],
      expected_manifest_ref=expected_base, phase2_vcf_ref=observed_base,
      match=matches[j],stringsAsFactors=FALSE)
  }
  status <- if (length(hits)==0) 'NO_SNP_REF_POSITION_TO_CHECK' else
    if (all(matches)) 'ALL_OBSERVED_SNP_REF_MATCH' else 'REF_MISMATCH_STOP_REVIEW'
  site_rows[[i]] <- data.frame(population_site_id=r$population_site_id,
      gene_id=r$gene_id, chromosome=arm, start=r$genomic_start,end=r$genomic_end,
      strand=r$guide_genomic_strand,n_snp_records_in_23bp=length(hits),
      n_ref_matches=sum(matches), n_ref_mismatches=sum(!matches),
      qc_status=status,stringsAsFactors=FALSE)
}
sites <- do.call(rbind,site_rows)
bases <- if(length(base_rows)) do.call(rbind,base_rows) else data.frame(
  population_site_id=character(),gene_id=character(),chromosome=character(),
  pos=integer(),expected_manifest_ref=character(),phase2_vcf_ref=character(),match=logical())
write.csv(sites,file.path(outdir,'AR1_frozen_site_REF_position_preflight_530.csv'),row.names=FALSE)
write.csv(bases,file.path(outdir,'AR1_manifest_vs_VCF_REF_base_comparisons.csv'),row.names=FALSE)
cat('Manifest SHA256:',sha,'\n')
cat('Frozen sites:',nrow(sites),'\n')
cat('Sites with >=1 VCF SNP REF comparison:',sum(sites$n_snp_records_in_23bp>0,na.rm=TRUE),'\n')
cat('Compared SNP REF positions:',nrow(bases),'\n')
cat('REF matches:',sum(bases$match),'mismatches:',sum(!bases$match),'\n')
cat('Site QC categories:\n');print(table(sites$qc_status))
cat('Output written to output/\n')
cat('\nIMPORTANT: Matching SNP REF alleles checks only loci reported as variant.\n')
cat('It does not establish concordance of all 23 guide bases or assembly compatibility.\n')
if (any(!bases$match)) {
  cat('\nSTOP: Reference mismatches found. Do not interpret Phase 2 genotypes as guide-level variation until resolved.\n')
  quit(save='no',status=3L)
}
