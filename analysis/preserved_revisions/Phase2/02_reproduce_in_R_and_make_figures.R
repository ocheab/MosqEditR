# MosqEdit-R R1.3 / R2.Q2 -- reproducible base-R benchmarking and figures
# Usage: Rscript 02_reproduce_in_R_and_make_figures.R /path/to/data_processed.zip [output_folder]
# Requires only base R. No FlyBase-outcome variables enter the simple comparators.
# Reuses the ORIGINAL 3 x 5 held-out M3/XGBoost predictions, not a fresh CV split.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1L) stop('Usage: Rscript 02_reproduce_in_R_and_make_figures.R data_processed.zip [output_folder]')
archive <- args[1]
outdir <- if(length(args)>=2) args[2] else '.'
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
read_source <- function(n) read.csv(unz(archive, paste0('data_processed/', n)), check.names=FALSE, stringsAsFactors=FALSE)
P <- read_source('08_pu_cv_predictions.csv')
P <- P[P$model_set=='M3_expression_orthology_annotation' & P$algorithm=='xgboost',]
F <- read_source('07_pu_M3_expression_orthology_annotation.csv')
stopifnot(nrow(P)==41517L, nrow(F)==13839L, sum(F$pu_label)==718L)
stopifnot(!anyDuplicated(F$gene_id), !anyDuplicated(P[c('repeat_id','fold','gene_id')]))

score_pct <- function(x) {
 x <- suppressWarnings(as.numeric(x))
 v <- rank(x, ties.method='average', na.last='keep') / sum(!is.na(x))
 v[is.na(v)] <- 0.5
 v
}
feature_names <- c('expr_ovary','ovary_specificity_delta','female_male_delta','tau_tissue_specificity','dmel_ortholog_count')
for(n in feature_names) F[[paste0('rank_',n)]] <- score_pct(F[[n]])
F$rank_has_dmel_ortholog <- as.numeric(F$has_dmel_ortholog)
F$rank_has_dmel_ortholog[is.na(F$rank_has_dmel_ortholog)] <- .5
F$simple_expression_mean4 <- rowMeans(F[,paste0('rank_',feature_names[1:4])])
F$simple_expression_orthology_mean5 <- (4*F$simple_expression_mean4+F$rank_dmel_ortholog_count)/5

feature_cols <- c(paste0('rank_',feature_names),'rank_has_dmel_ortholog','simple_expression_mean4','simple_expression_orthology_mean5')
P <- merge(P,F[,c('gene_id','pu_label',feature_cols)],by='gene_id',all.x=TRUE,sort=FALSE,suffixes=c('','_source'))
stopifnot(all(P$pu_label==P$pu_label_source),nrow(P)==41517L)
P$pu_label_source <- NULL

methods <- c('M3_XGBoost_PU_out_of_fold'='pu_score',
 'Single_ovary_expression'='rank_expr_ovary',
 'Single_ovary_specificity'='rank_ovary_specificity_delta',
 'Single_female_male_delta'='rank_female_male_delta',
 'Single_Dmel_ortholog_count'='rank_dmel_ortholog_count',
 'Simple_4_expression_mean'='simple_expression_mean4',
 'Simple_4_expression_plus_orthology_mean5'='simple_expression_orthology_mean5',
 'Ortholog_presence_diagnostic'='rank_has_dmel_ortholog')

measure <- function(y,s) {
 y <- as.integer(y);s<-as.numeric(s)
 stopifnot(all(is.finite(s)), all(y %in% c(0L,1L)),sum(y)>0,sum(y==0)>0)
 n <- length(y); npos<-sum(y); prev<-npos/n
 # AUROC with average rank at equal scores
 rk <- rank(s, ties.method='average')
 auroc <- (sum(rk[y==1]) - npos*(npos+1)/2) /(npos*(n-npos))
 o <- order(-s,seq_along(s));so <- s[o];yo<-y[o]
 cumpos<-cumsum(yo)
 # Tie-aware AP: one threshold per distinct score, no original-order artefact.
 thr <- c(which(diff(so)!=0),n)
 prec <- cumpos[thr]/thr
 rec <- cumpos[thr]/npos
 ap <- sum(diff(c(0,rec))*prec)
 trap_tie <- sum(diff(c(0,rec))*(head(c(1,prec),-1)+prec)/2)
 # Original legacy PR trapezoid -- ONLY for verifying the M3 source result.
 source_rec<-c(0,cumpos/npos)
 source_prec<-c(1,cumpos/seq_len(n))
 trap_source<-sum(diff(source_rec)*(head(source_prec,-1)+tail(source_prec,-1))/2)
 ans <- c(auprc_trapezoidal=trap_source, auprc_tieaware_trapezoidal=trap_tie,
   average_precision_tieaware=ap,auroc=auroc,n=n,n_positive=npos,prevalence=prev)
 for(k in c(10,20,50,100,200)){
  kk<-min(k,n);cutoff<-so[kk]
  higher<-s>cutoff; tied<-s==cutoff
  # Mean over all equally plausible selections at the cutoff tie.
  hits<-sum(y[higher])+(kk-sum(higher))*mean(y[tied])
  ans[paste0('positive_at_',k)]<-hits
  ans[paste0('precision_at_',k)]<-hits/kk
  ans[paste0('recall_at_',k)]<-hits/npos
  ans[paste0('enrichment_at_',k)]<-(hits/kk)/prev
 }
 ans
}

set.seed(20261031L)
groups <- split(seq_len(nrow(P)), paste(P$repeat_id,P$fold,sep=':'))
all_res <- list();j<-0L
for(ix in groups){
 df<-P[ix,];y<-df$pu_label
 for(m in names(methods)){
  j<-j+1L
  all_res[[j]]<-data.frame(repeat_id=df$repeat_id[1],fold=df$fold[1],method=m,
      null_draw=NA_integer_,as.list(measure(y,df[[methods[m]]])),check.names=FALSE)
 }
 for(q in seq_len(300L)){
  j<-j+1L
  all_res[[j]]<-data.frame(repeat_id=df$repeat_id[1],fold=df$fold[1],method='Permuted_scores_null',
      null_draw=q,as.list(measure(y,sample(df$pu_score))),check.names=FALSE)
 }
}
D <- do.call(rbind,all_res)
write.csv(D,file.path(outdir,'R13_R_fold_baseline_detail.csv'),row.names=FALSE)
numvars<-names(D)[vapply(D,is.numeric,logical(1))]
numvars<-setdiff(numvars,c('repeat_id','fold','null_draw'))
method_order<-c(names(methods),'Permuted_scores_null')
S <- do.call(rbind,lapply(method_order,function(m){
 a<-D[D$method==m,]
 ans<-data.frame(method=m)
 for(v in numvars){ans[[paste0(v,'_mean')]]<-mean(a[[v]],na.rm=TRUE)
  ans[[paste0(v,'_sd')]]<-sd(a[[v]],na.rm=TRUE)}
 ans
}))
write.csv(S,file.path(outdir,'R13_R_fold_matched_summary.csv'),row.names=FALSE)

# Full gene-universe: three cross-fitted M3 scores per gene averaged, not true validation.
M <- aggregate(pu_score~gene_id+pu_label,data=P,FUN=mean)
nms <- names(M);nms[nms=='pu_score'] <- 'mean_M3_crossfit_PU_score'; names(M)<-nms
G <- merge(M,F[,c('gene_id','has_dmel_ortholog',feature_cols)],by='gene_id',all.x=TRUE,sort=FALSE)
stopifnot(nrow(G)==13839L,sum(G$pu_label)==718L)
write.csv(G,file.path(outdir,'R13_R_genomewide_OOF_scores.csv'),row.names=FALSE)
topmethods<-c('M3_XGBoost_mean_OOF'='mean_M3_crossfit_PU_score',
 'Single_ovary_expression'='rank_expr_ovary',
 'Simple_4_expression_mean'='simple_expression_mean4',
 'Simple_4_expression_plus_orthology_mean5'='simple_expression_orthology_mean5')
toprecords <- list();j<-0L
for(m in names(topmethods)){
 score<-G[[topmethods[m]]]
 o<-order(-score,G$gene_id)
 for(k in c(20,50,100,200,500,1000)){
  ix<-o[seq_len(k)]
  hit<-sum(G$pu_label[ix]);prev<-mean(G$pu_label)
  j<-j+1L;toprecords[[j]]<-data.frame(method=m,k=k,n_eligible=nrow(G),
   n_observed_positive=sum(G$pu_label),expected_random_positive_count=k*prev,
   observed_positive_count=hit,unlabeled_count=k-hit,
   observed_positive_precision=hit/k,observed_positive_recall=hit/sum(G$pu_label),
   enrichment_over_prevalence=(hit/k)/prev)
 }
}
T <- do.call(rbind,toprecords)
write.csv(T,file.path(outdir,'R13_R_full_universe_topk_utility.csv'),row.names=FALSE)

# Publication figure 1: tie-aware AP comparison. 15 folds are correlated.
labels <- c('M3 XGBoost (out-of-fold)','Ovary expression','Ovary specificity',
 'Female-male contrast','Ortholog count','4-expression composite',
 'Expression + orthology composite','Ortholog-present indicator','Permuted-score null')
vals<-S$average_precision_tieaware_mean
ord <- order(vals)
png(file.path(outdir,'R13_R_AveragePrecision_Baselines_600dpi.png'),width=7200,height=4700,res=600)
par(mar=c(5.4,16.8,3.5,2.1),family='sans')
cols<-ifelse(S$method[ord]=='M3_XGBoost_PU_out_of_fold','#315A78','#9BAAB7')
mids<-barplot(vals[ord],horiz=TRUE,names.arg=labels[ord],las=1,col=cols,border=NA,
 xlim=c(0,max(vals)*1.19),xlab='Tie-aware average precision (observed-positive vs unlabeled)',
 main='Descriptive PU-label recovery: original M3 versus simple baselines',cex.names=.75)
axis(1,at=pretty(c(0,max(vals))))
text(vals[ord]+0.003,mids,labels=sprintf('%.3f',vals[ord]),cex=.73,pos=4)
mtext('Original 3 x 5 out-of-fold M3 comparisons; 300 null permutations/fold. Not biological validation.',side=1,line=4,cex=.65)
dev.off()

# Publication figure 2: label vs unlabeled composition in genome-wide OOF top-k.
tm<-T[T$method=='M3_XGBoost_mean_OOF',];xs<-seq_len(nrow(tm))
png(file.path(outdir,'R13_R_Genomewide_OOF_TopK_600dpi.png'),width=7200,height=4600,res=600)
par(mar=c(5,5,3.5,1.5),family='sans')
barplot(rbind(tm$observed_positive_count,tm$unlabeled_count),names.arg=tm$k,
 col=c('#315A78','#CFD6D8'),border=NA,legend.text=c('Observed FlyBase positive','Unlabeled (not negative)'),
 args.legend=list(x='topleft',bty='n',cex=.76),xlab='Top-k genes by mean M3 out-of-fold PU score',
 ylab='Gene count',main='Genome-wide out-of-fold ranking is distinct from the preselected CRISPR shortlist')
dev.off()
cat('Done. Main reproduced metrics:\n')
print(S[,c('method','average_precision_tieaware_mean','auroc_mean','enrichment_at_50_mean')],row.names=FALSE)
cat('\nImportant: no experimental Anopheles sterility validation or true recall in unlabeled genes inferred.\n')
