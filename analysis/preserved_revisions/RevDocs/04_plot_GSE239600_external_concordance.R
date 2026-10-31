#!/usr/bin/env Rscript
# Base-R, source-derived publication figure. Execute inside RevDocs after script 03.
# Rscript 04_plot_GSE239600_external_concordance.R .
a <- commandArgs(trailingOnly = TRUE)
base <- if(length(a)) a[1] else "."
indir <- file.path(base, "results")
f <- file.path(indir, "R2Q3_GSE239600_gene_level_matched.csv")
if(!file.exists(f)) stop("Missing results/R2Q3_GSE239600_gene_level_matched.csv; run corrected script 03 first")
d <- read.csv(f, stringsAsFactors = FALSE)
needed <- c("gene_id", "expr_ovary", "mean_wt_ovary")
if(!all(needed %in% names(d))) stop("Missing required columns: ", paste(setdiff(needed,names(d)), collapse=", "))
stopifnot(!anyDuplicated(d$gene_id))
ok <- is.finite(d$expr_ovary) & is.finite(d$mean_wt_ovary)
d <- d[ok, ,drop=FALSE]
if(nrow(d)<500 || any(d$expr_ovary < 0) || any(d$mean_wt_ovary < 0)) stop("Unexpected comparison dataset; inspect source")
n <- nrow(d); k <- ceiling(.10*n)
old <- d$gene_id[order(-d$expr_ovary, d$gene_id)][seq_len(k)]
new <- d$gene_id[order(-d$mean_wt_ovary, d$gene_id)][seq_len(k)]
both <- d$gene_id %in% intersect(old,new)
oldonly <- d$gene_id %in% setdiff(old,new)
newonly <- d$gene_id %in% setdiff(new,old)
rho <- cor(d$expr_ovary,d$mean_wt_ovary,method="spearman")
ninter <- sum(both)
cat(sprintf("Matched N=%d; Spearman rho=%.7f; top decile overlap=%d/%d (%.2f%%)\n", n, rho, ninter,k,100*ninter/k))
result <- data.frame(n_matched=n, spearman_rho=rho, k_top_decile=k,n_overlap=ninter, overlap_fraction=ninter/k, independent_random_expected_overlap=k*k/n)
write.csv(result,file.path(indir,"R2Q3_GSE239600_external_concordance_figure_stats.csv"),row.names=FALSE)
fig <- file.path(indir, "R2Q3_GSE239600_ovarian_expression_concordance_600dpi.png")
png(fig,width=3300,height=2700,res=600,type="cairo")
par(mar=c(5.1,5.3,2.7,1),mgp=c(3.1,0.9,0),cex=0.82)
x<-log2(d$expr_ovary + 1);y<-log2(d$mean_wt_ovary + 1)
plot(x,y,type="n",xlab="MozAtlas ovarian expression (log2[value + 1])",ylab="GSE239600 WT ovarian expression (log2[normalized count + 1])",main="Cross-study ovarian-expression concordance",las=1)
points(x[!both & !oldonly & !newonly],y[!both & !oldonly & !newonly],pch=16,cex=.26,col=grDevices::adjustcolor("#6b7280",alpha.f=.22))
points(x[oldonly],y[oldonly],pch=16,cex=.45,col=grDevices::adjustcolor("#2563eb",alpha.f=.45))
points(x[newonly],y[newonly],pch=16,cex=.45,col=grDevices::adjustcolor("#c2410c",alpha.f=.45))
points(x[both],y[both],pch=16,cex=.6,col=grDevices::adjustcolor("#157347",alpha.f=.70))
legend("topleft",bty="n",cex=.74,pch=16,col=c("#157347","#2563eb","#c2410c","#6b7280"),legend=c(sprintf("Both top deciles: %d genes",ninter),"MozAtlas top-decile only","RNA-seq top-decile only","Neither top decile"))
mtext(sprintf("N = %s matched genes; Spearman rho = %.3f; top-decile overlap = %d/%d (%.1f%%)",format(n,big.mark=","),rho,ninter,k,100*ninter/k),side=3,cex=.79,line=.2)
dev.off()
cat("Wrote:",fig,"\n")
