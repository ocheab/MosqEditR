#!/usr/bin/env Rscript
# MosqEdit-R Reviewer R1.2 / R2.Q1: source-aware label provenance and OOF diagnostic
# Base R only. Requires original data_processed.zip or extracted data_processed directory.
# Usage: Rscript R12_Reproduce_Orthology_Audit.R /path/data_processed.zip /path/output_dir
args <- commandArgs(trailingOnly=TRUE)
input <- if(length(args)>=1) args[[1]] else "data_processed.zip"
out <- if(length(args)>=2) args[[2]] else "R12_outputs"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
read_src <- function(name) {
  if (grepl("[.]zip$", input, ignore.case=TRUE)) {
    read.csv(unz(input, paste0("data_processed/",name)), stringsAsFactors=FALSE, check.names=FALSE)
  } else {
    read.csv(file.path(input,name), stringsAsFactors=FALSE, check.names=FALSE)
  }
}
truth <- function(x) !is.na(x) & (x==TRUE | x=="TRUE" | x==1)
pr_auc_trapezoid <- function(y, score) {
  keep <- is.finite(score) & !is.na(y); y <- as.integer(y[keep]); score <- score[keep]
  np <- sum(y==1L); if(np==0L) return(NA_real_)
  ord <- order(score, decreasing=TRUE); y <- y[ord]
  tp <- cumsum(y==1L); fp <- cumsum(y==0L)
  recall <- c(0,tp/np); precision <- c(1,tp/pmax(tp+fp,1))
  sum(diff(recall)*(head(precision,-1L)+tail(precision,-1L))/2)
}
roc_auc_rank <- function(y,score) {
  keep <- is.finite(score) & !is.na(y); y <- as.integer(y[keep]); score <- score[keep]
  np <- sum(y==1L); nn <- sum(y==0L)
  if(np==0L || nn==0L) return(NA_real_)
  ranks <- rank(score,ties.method="average")
  (sum(ranks[y==1L])-np*(np+1)/2)/(np*nn)
}
master <- read_src("master_gene_feature_matrix.csv")
model <- read_src("07_pu_M3_expression_orthology_annotation.csv")
pred <- read_src("08_pu_cv_predictions.csv")
top100 <- read_src("09_top100_prepopulation_shortlist.csv")
stopifnot(nrow(master)==13845, !anyDuplicated(master$gene_id),nrow(model)==13839,
          length(unique(top100$gene_id))==100)
master$orth <- truth(master$has_dmel_ortholog)
master$positive <- truth(master$flybase_female_sterile)
stopifnot(sum(master$positive)==719, sum(master$orth)==8902,
          !any(master$positive & !master$orth))
master$label_status <- ifelse(master$positive,"comparative_FlyBase_positive",
                              ifelse(is.na(master$flybase_female_sterile),
                                     "no_mapped_FlyBase_phenotype_record",
                                     "FlyBase_record_without_strict_female_sterility"))
write.csv(master[,c("gene_id","orth","positive","label_status","is_prespecified_benchmark")],
          file.path(out,"R12_gene_label_provenance.csv"),row.names=FALSE)
partitions <- list(all=master, mapped=master[master$orth,], unmapped=master[!master$orth,])
count <- do.call(rbind,lapply(names(partitions),function(nm){
  d <- partitions[[nm]]
  data.frame(stratum=nm,n_genes=nrow(d),n_comparative_positive=sum(d$positive),
     n_with_mapped_phenotype_flag=sum(!is.na(d$flybase_female_sterile)),
     n_without_mapped_phenotype_flag=sum(is.na(d$flybase_female_sterile)),
     positive_fraction=sum(d$positive)/nrow(d))
}))
write.csv(count,file.path(out,"R12_label_provenance_summary.csv"),row.names=FALSE)
stopifnot(sum(top100$pu_label==1)==100)
write.csv(top100[,c("gene_id","novel_candidate_rank","pu_label","final_integrated_score")],
          file.path(out,"R12_top100_comparative_label_status.csv"),row.names=FALSE)
subset_pred <- pred[pred$model_set=="M3_expression_orthology_annotation" & pred$algorithm=="xgboost",]
idx <- match(subset_pred$gene_id,model$gene_id)
stopifnot(!anyNA(idx), nrow(subset_pred)==3*nrow(model))
subset_pred$orth <- truth(model$has_dmel_ortholog[idx])
keys <- unique(subset_pred[,c("repeat_id","fold")]); keys <- keys[order(keys$repeat_id,keys$fold),]
fold_metrics <- do.call(rbind,lapply(seq_len(nrow(keys)),function(i){
  g <- subset_pred[subset_pred$repeat_id==keys$repeat_id[i] & subset_pred$fold==keys$fold[i],]
  p <- g[g$orth,]; q <- g[!g$orth,]
  u <- g[order(-g$pu_score),]
  data.frame(repeat_id=keys$repeat_id[i],fold=keys$fold[i],n_genes=nrow(g),
   n_positives=sum(g$pu_label),n_orthology_mapped=nrow(p),n_orthology_unmapped=nrow(q),
   full_trapezoidal_AUPRC=pr_auc_trapezoid(g$pu_label,g$pu_score),
   orthology_mapped_trapezoidal_AUPRC=pr_auc_trapezoid(p$pu_label,p$pu_score),
   full_AUROC=roc_auc_rank(g$pu_label,g$pu_score),
   orthology_mapped_AUROC=roc_auc_rank(p$pu_label,p$pu_score),
   full_random_prevalence=mean(g$pu_label==1),
   mapped_random_prevalence=mean(p$pu_label==1),
   orthology_indicator_AUROC=roc_auc_rank(g$pu_label,as.integer(g$orth)),
   orthology_unmapped_in_top1000=sum(!head(u$orth,1000)),
   min_mapped_score=min(p$pu_score),max_unmapped_score=max(q$pu_score))
}))
stopifnot(nrow(fold_metrics)==15L)
write.csv(fold_metrics,file.path(out,"R12_crossfit_folds.csv"),row.names=FALSE)
means <- sapply(fold_metrics[,vapply(fold_metrics,is.numeric,logical(1)),drop=FALSE],mean)
mean_table <- data.frame(metric=names(means),mean_15_folds=as.numeric(means))
write.csv(mean_table,file.path(out,"R12_crossfit_means.csv"),row.names=FALSE)
# Manual, DOI-linked evidence: separate biological endpoints, never joined to training labels.
native <- data.frame(gene_id=c("AGAP005958","AGAP011377","AGAP007280","AGAP004050","AGAP013051","AGAP006241"),
  phenotype=c(rep("experimental female sterility after disruption",3),
   "female-specific sterility after dsx disruption","female lethality after fle disruption",
   "germ-cell development/egg production impairment"),
  doi=c(rep("10.1038/nbt.3439",3),"10.1038/nbt.4245",
    "10.1126/sciadv.ade8903","10.1186/s12864-023-09887-0"),
  endpoint_group=c(rep("strict_female_sterility",4),"female_lethality_not_identical_to_sterility",
    "reproductive_function_not_uniform_sterility"),stringsAsFactors=FALSE)
m <- match(native$gene_id,master$gene_id);stopifnot(!anyNA(m))
native$mapped_Dmel_orthology <- master$orth[m];native$FlyBase_positive <- master$positive[m]
native$preexisting_landmark <- truth(master$is_prespecified_benchmark[m])
write.csv(native,file.path(out,"R12_Anopheles_native_evidence.csv"),row.names=FALSE)
# Figure, base R only. The two panels compare within the same prediction files;
# no cross-validation retraining and no claim of Anopheles-specific accuracy.
png(file.path(out,"R12_orthology_bias_diagnostic_600dpi.png"),
    width=8.0,height=3.6,units="in",res=600)
par(mfrow=c(1,2),mar=c(6,4.3,2.5,1.2),oma=c(0,0,1,0),cex=0.86)
barplot(c(mean(fold_metrics$full_AUROC),mean(fold_metrics$orthology_mapped_AUROC)),
  names.arg=c("Full held-out","Mapped-orthology\nonly"),ylim=c(0,1),
  ylab="AUROC: comparative PU labels",col=c("grey38","grey72"),
  main="Held-out discrimination",las=1)
barplot(c(mean(fold_metrics$full_trapezoidal_AUPRC),mean(fold_metrics$orthology_mapped_trapezoidal_AUPRC),
          mean(fold_metrics$mapped_random_prevalence)),
 names.arg=c("Full held-out","Mapped-orthology\nonly","Mapped\nprevalence"),
 ylim=c(0,0.20),ylab="PR AUC / conditional baseline",col=c("grey38","grey72","white"),
 main="Comparative-label retrieval",las=1)
mtext("Source: original MosqEdit-R fivefold x three repeat cross-fitted M3-XGBoost predictions",outer=TRUE,cex=0.7)
dev.off()
cat("Source audit complete. Full positive count:",sum(master$positive),"\n")
cat("Orthology-mapped positives:",sum(master$positive & master$orth),"/",sum(master$positive),"\n")
cat("Original vs orthology-only AUROC:",round(mean(fold_metrics$full_AUROC),6),"vs",round(mean(fold_metrics$orthology_mapped_AUROC),6),"\n")
cat("Original vs orthology-only trapezoidal AUPRC:",round(mean(fold_metrics$full_trapezoidal_AUPRC),6),"vs",round(mean(fold_metrics$orthology_mapped_trapezoidal_AUPRC),6),"\n")
cat("Important: All outcomes are comparative-label proxy measures, NOT Anopheles sterility accuracy.\n")
