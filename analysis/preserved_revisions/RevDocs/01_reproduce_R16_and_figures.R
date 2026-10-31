# R1.6 / R2.Q4: reproducible landmark rank and weight-sensitivity analysis
# Base R only. Input: original data_processed directory or original data_processed.zip
# Usage: Rscript 01_reproduce_R16_and_figures.R <data_processed_folder_or_zip> <output_folder>
a <- commandArgs(trailingOnly=TRUE)
if(length(a)<1L) stop('Usage: Rscript 01_reproduce_R16_and_figures.R <folder_or_zip> [output_folder]')
source_path <- a[1]; dest <- if(length(a)>=2L) a[2] else 'R16_results'
dir.create(dest,showWarnings=FALSE,recursive=TRUE)
read_original <- function(name) {
 if(dir.exists(source_path)) {
   f <- file.path(source_path,name)
   if(!file.exists(f)) stop('Missing required CSV: ',f)
   return(read.csv(f,check.names=FALSE,stringsAsFactors=FALSE))
 }
 if(!file.exists(source_path)) stop('Input does not exist: ',source_path)
 n <- paste0('data_processed/',name)
 if(!(n %in% unzip(source_path,list=TRUE)$Name)) stop('Missing archive entry: ',n)
 read.csv(unz(source_path,n),check.names=FALSE,stringsAsFactors=FALSE)
}
d <- read_original('09_targeted_integrated_prioritization.csv')
ab <- read_original('15_benchmark_ablation.csv')
stopifnot(nrow(d)==205L,!anyDuplicated(d$gene_id),sum(d$is_external_benchmark)==6L)
mat <- as.matrix(d[,c('biological_discovery_percentile','pu_consensus_percentile','crispr_tractability_percentile')]); storage.mode(mat)<-'double'
valid <- is.finite(mat)
calculated <- function(w) {
 stopifnot(length(w)==3L,all(w>=0),abs(sum(w)-1)<1e-8)
 numer <- rowSums(sweep(ifelse(valid,mat,0),2,w,'*'))
 denom <- rowSums(sweep(valid*1,2,w,'*'))
 raw <- numer / denom;raw[denom<1e-12]<-NA_real_
 raw*(0.8+0.2*rowMeans(valid))
}
drank <- function(s) {
 r <- rep(NA_integer_,length(s));good <- which(is.finite(s))
 ord <- good[order(-s[good],d$gene_id[good])]
 r[ord] <- seq_along(ord);r
}
w <- list(biology_led=c(.60,.25,.15),primary=c(.50,.30,.20),validation_led=c(.40,.35,.25))
r <- lapply(w,function(ww) drank(calculated(ww)))
stopifnot(max(abs(calculated(w$primary)-d$final_integrated_score),na.rm=TRUE)<1e-12,
 identical(as.integer(r$primary),as.integer(d$final_integrated_rank)),
 identical(as.integer(r$biology_led),as.integer(d$rank_biology_led)),
 identical(as.integer(r$validation_led),as.integer(d$rank_validation_led)))
medianr <- apply(do.call(cbind,r),1,median)
ordered <- order(medianr,-calculated(w$primary),d$gene_id)
finalr <- integer(nrow(d));finalr[ordered] <- seq_len(nrow(d))
stopifnot(identical(finalr,as.integer(d$final_prepopulation_rank)))
b <- d[d$is_external_benchmark,]
write.csv(b[,c('gene_id','benchmark_name','benchmark_class','pu_label','biological_discovery_percentile','pu_consensus_percentile','crispr_tractability_percentile','final_evidence_domain_coverage','final_integrated_score','final_integrated_rank','final_prepopulation_rank','rank_biology_led','rank_validation_led')],file.path(dest,'R16_original_six_landmarks.csv'),row.names=FALSE,na='')
sc <- expand.grid(biology_units=0:20,pu_units=0:20)
sc <- sc[sc$biology_units+sc$pu_units<=20,]
sc$crispr_units <- 20-sc$biology_units-sc$pu_units
stopifnot(nrow(sc)==231L)
res <- vector('list',nrow(sc))
for(i in seq_len(nrow(sc))) {
 ww <- unlist(sc[i,c('biology_units','pu_units','crispr_units')])/20
 rr <- drank(calculated(ww))
 res[[i]] <- data.frame(scenario_id=sprintf('B%02d_P%02d_C%02d',sc$biology_units[i],sc$pu_units[i],sc$crispr_units[i]),gene_id=b$gene_id,benchmark_name=b$benchmark_name,benchmark_class=b$benchmark_class,biology_weight=ww[1],pu_weight=ww[2],crispr_weight=ww[3],rank=rr[match(b$gene_id,d$gene_id)])
}
long <- do.call(rbind,res)
write.csv(long,file.path(dest,'R16_231_weight_grid_landmark_ranks.csv'),row.names=FALSE,na='')
summ <- do.call(rbind,lapply(split(long,long$gene_id),function(x) {
 good <- x$rank[!is.na(x$rank)]
 data.frame(gene_id=x$gene_id[1],n_scenarios=nrow(x),n_ranks=length(good),min_rank=min(good),median_rank=median(good),max_rank=max(good),fraction_top50=mean(good<=50),fraction_top100=mean(good<=100))
}))
summ$original_rank <- d$final_prepopulation_rank[match(summ$gene_id,d$gene_id)]
summ <- summ[order(summ$original_rank),]
write.csv(summ,file.path(dest,'R16_231_weight_grid_landmark_summary.csv'),row.names=FALSE,na='')
write.csv(ab,file.path(dest,'R16_original_step15_benchmark_ablations.csv'),row.names=FALSE,na='')
labels <- c('AGAP007280 (fertility)','AGAP005958 (fertility)','AGAP007031 (FREP1)','AGAP011377 (fertility)','AGAP004050 (doublesex)','AGAP028779 (mir-184)')
ordgenes <- b$gene_id[order(b$final_prepopulation_rank)]
# Visual 1: weight-grid range (descriptive; do not optimize)
png(file.path(dest,'R16_landmark_weight_sensitivity_600dpi.png'),width=3600,height=2700,res=600,type='cairo')
par(mar=c(4.7,13,3.5,2),cex.axis=.85)
plot(NA,xlim=c(1,205),ylim=c(0.5,6.5),yaxt='n',xlab='Rank among 205 candidates (smaller = higher priority)',ylab='',main='Landmark ranks across 231 exploratory domain-weight settings')
abline(v=c(20,50,100),lty=c(3,2,2),col='grey75')
for(i in seq_along(ordgenes)) {
 g<-ordgenes[i];y<-7-i
 vals<-long$rank[long$gene_id==g];vals<-vals[is.finite(vals)]
 points(vals,rep(y,length(vals)),pch=16,cex=.12,col=rgb(.18,.35,.65,.20))
 points(d$final_prepopulation_rank[d$gene_id==g],y,pch=18,cex=1.25,col='#AD4935')
}
axis(2,at=6:1,labels=labels,las=1,cex.axis=.76)
legend('bottomright',legend=c('Exploratory scenario ranks','Original frozen rank'),pch=c(16,18),col=c('#6484AA','#AD4935'),bty='n',cex=.77)
mtext('Original rank is not an independently validated efficacy measure.',side=1,line=3.5,cex=.75)
dev.off()
# Visual 2: domain-score structural difference
png(file.path(dest,'R16_landmark_domain_percentiles_600dpi.png'),width=3500,height=2600,res=600,type='cairo')
par(mar=c(5,13,3.5,2),cex.axis=.85)
plot(NA,xlim=c(0,1.03),ylim=c(.5,6.5),yaxt='n',xlab='Original domain percentile (0–1; missing not plotted)',ylab='',main='Original evidence-domain support for six strategy-diverse landmarks')
cl<-c('#2F6F9E','#BE7032','#45815C');pchs<-c(16,17,15)
for(i in seq_along(ordgenes)) {
 ix <- match(ordgenes[i],d$gene_id);y<-7-i
 for(j in 1:3) if(is.finite(mat[ix,j])) points(mat[ix,j],y+(.16*(j-2)),pch=pchs[j],col=cl[j],cex=1.1)
}
axis(2,at=6:1,labels=labels,las=1,cex.axis=.76)
legend('bottomleft',legend=c('Biology','PU comparative-label score','Reference CRISPR tractability'),pch=pchs,col=cl,bty='n',cex=.77)
dev.off()
cat('PASS: exact reconstruction of 205 original Step09 primary + scenario + consensus ranks.\n')
print(summ,row.names=FALSE,digits=3)
cat('Created: 4 CSVs and two 600-dpi R figures in',dest,'\n')
cat('Interpretation: descriptive sensitivity only; no independent efficacy validation, no weight calibration.\n')
