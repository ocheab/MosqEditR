#!/usr/bin/env Rscript
# External 2023 GSE239600 wild-type RNA-seq concordance check.
# Does not substitute the RNA-seq samples into the MozAtlas 15-tissue FRIS.
# Usage: Rscript 03_external_GSE239600_RNAseq_comparison.R <normalized_counts.txt.gz> <base_dir>
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2) stop("Usage: Rscript 03_external_GSE239600_RNAseq_comparison.R GSE239600_...counts.txt.gz .")
rnaseq_path <- args[1]; base <- args[2]
if(!file.exists(rnaseq_path)) stop("Missing RNA-seq file: ",rnaseq_path)
res <- file.path(base,"results");dir.create(res,recursive=TRUE,showWarnings=FALSE)
original <- read.csv(file.path(base,"inputs","02_mozatlas_gene_expression.csv"),stringsAsFactors=FALSE)
top <- read.csv(file.path(base,"inputs","09_novel_candidate_ranking.csv"),stringsAsFactors=FALSE)
stopifnot(anyDuplicated(original$gene_id)==0L,all(c("gene_id","expr_ovary","ovary_specificity_delta")%in%names(original)))
cat("Reading modern GEO table (GSE239600)\n")
r <- read.delim(gzfile(rnaseq_path),check.names=FALSE,stringsAsFactors=FALSE,quote="",comment.char="",fill=TRUE)
write.csv(data.frame(column_index=seq_along(names(r)),source_column=names(r)),file.path(res,"R2Q3_GSE239600_header_audit.csv"),row.names=FALSE)
cat("Modern GEO columns:\n");print(names(r));cat("Rows:",nrow(r),"\n")
# Corrected after source-schema audit (31 Oct 2026):
# Parent has exact and unique AGAP IDs in all 11,911 records.
# Gene Symbol also contains AGAP IDs, so automated identification is ambiguous.
gene_cand <- which(names(r) == "Parent")
if(length(gene_cand) != 1L) stop("Expected exactly one Parent column with exact AGAP IDs.")
ids <- trimws(as.character(r[[gene_cand]]))
if(anyNA(ids) || any(!grepl("^AGAP[0-9]{6}$",ids)))
  stop("Parent contains missing/non-exact AGAP IDs; manual annotation review required.")
if(anyDuplicated(ids)) stop("Duplicated AGAP IDs in Parent; resolve before merging.")
cat("Validated Parent gene IDs:", length(ids), "exact, unique AGAP identifiers\n")
if("Gene Symbol" %in% names(r)) {
  sym <- as.character(r[["Gene Symbol"]])
  prefix <- ifelse(grepl("^AGAP[0-9]{6}",sym), substr(sym,1L,10L), NA_character_)
  disagreement <- which(!is.na(prefix) & prefix != ids)
  write.csv(data.frame(Parent=ids, Gene_Symbol=sym, Symbol_agap_prefix=prefix,
                       prefix_matches_Parent=is.na(prefix) | prefix==ids),
            file.path(res,"R2Q3_GSE239600_gene_id_crosscheck.csv"),row.names=FALSE)
  if(length(disagreement)) warning("Gene Symbol/Parent prefix disagreement in ", length(disagreement),
                                   " rows. Check results/R2Q3_GSE239600_gene_id_crosscheck.csv")
}

norm<-tolower(gsub("[^[:alnum:]]","",names(r)))
# All known WT samples defined by GEO sample metadata, with negative-control KO samples excluded.
map <- list(
  ovary=c("GSM7668510","GSM7668511","GSM7668512","GSM7668513"),
  fatbody=c("GSM7668506","GSM7668507"),
  wholebody=c("GSM7668501","GSM7668502","GSM7668503")
)
map_titles<-list(ovary=c("wt1ovary","wt2ovary","wtovary4","wtovary5"),fatbody=c("wtfb1","wtfb2"),wholebody=c("wt4wb","wt5wb","wt6wb"))
chosen <- list()
for(group in names(map)) {
  indices<-integer(0)
  for(j in seq_along(map[[group]])){
    gsm <-tolower(map[[group]][j]);title<-map_titles[[group]][j]
    ii<-which(grepl(gsm,norm,fixed=TRUE) | grepl(title,norm,fixed=TRUE))
    if(length(ii)!=1L || ii==gene_cand) stop("Cannot uniquely match ",gsm," (",title,") to column; check header_audit.csv; no sample substitutions inferred.")
    indices<-c(indices,ii)
  }
  if(anyDuplicated(indices)) stop("Duplicate WT columns across ",group)
  chosen[[group]]<-indices
}
if(anyDuplicated(unlist(chosen)))stop("Sample groups overlap")
conv<-function(cols){
 x <- as.data.frame(lapply(r[cols],function(v)suppressWarnings(as.numeric(as.character(v)))))
 mat <- as.matrix(x)
 if(any(mat[is.finite(mat)]<0))stop("Negative normalized expression values; inspect data-processing scale before using the RNA-seq values.")
 if(any(colSums(is.finite(mat)) < 0.8*nrow(mat))) stop("Sample columns are not mostly numeric; inspect output header first.")
 rowMeans(mat,na.rm=FALSE)
}
rna <- data.frame(gene_id=ids,mean_wt_ovary=conv(chosen$ovary),mean_wt_fatbody=conv(chosen$fatbody),mean_wt_wholebody=conv(chosen$wholebody))
rna<-rna[grepl("^AGAP[0-9]{6}$",rna$gene_id),]
merged<-merge(original[,c("gene_id","expr_ovary","ovary_specificity_delta")],rna,by="gene_id",all=FALSE)
merged$original_ovary_log2p1 <-log2(pmax(0,merged$expr_ovary)+1)
merged$rna_ovary_log2p1<-log2(pmax(0,merged$mean_wt_ovary)+1)
merged$RNA_ovary_vs_fatbody_log2_ratio<-log2((merged$mean_wt_ovary+1)/(merged$mean_wt_fatbody+1))
merged$RNA_ovary_vs_wholebody_log2_ratio<-log2((merged$mean_wt_ovary+1)/(merged$mean_wt_wholebody+1))
matched <- with(merged,is.finite(expr_ovary)&is.finite(mean_wt_ovary))
if(sum(matched)<500L)stop("Fewer than 500 matched genes; inspect gene ID mapping and GEO table. No claimed cross-study comparison.")
rho<-cor(merged$expr_ovary[matched],merged$mean_wt_ovary[matched],method="spearman")
base_ids <- top$gene_id[top$novel_candidate_rank<=100 & !top$is_external_benchmark]
coverage<-sum(base_ids %in% merged$gene_id[matched])
# Top 10% rank concordance is calculated within the same *matched* gene universe.
ids_use<-merged$gene_id[matched]
old<-merged$expr_ovary[matched];new<-merged$mean_wt_ovary[matched]
k<-max(1,ceiling(.1*length(ids_use)))
oldtop<-ids_use[order(-old,ids_use)][1:k];newtop<-ids_use[order(-new,ids_use)][1:k]
summary <- data.frame(dataset="GSE239600_WT_ovary_vs_MozAtlas",n_matched_genes=sum(matched),spearman_rho_ovary=rho,n_original_top100_covered=coverage,top_decile_overlap_n=length(intersect(oldtop,newtop)),top_decile_denominator=k,n_WT_ovary=4,n_WT_fatbody=2,n_WT_wholebody=3)
write.csv(summary,file.path(res,"R2Q3_GSE239600_external_comparison_summary.csv"),row.names=FALSE)
write.csv(merged,file.path(res,"R2Q3_GSE239600_gene_level_matched.csv"),row.names=FALSE)
print(summary,row.names=FALSE)
cat("Scope: wild-type, blood-fed, 24 h post meal, limited three adult tissue sample groups; NOT full 15-tissue replacement; no PU refit.\n")
