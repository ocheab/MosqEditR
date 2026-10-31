# MosqEdit-R: R1.6 and R2.Q4 landmark analysis
# Requires base R only.

a <- commandArgs(trailingOnly=TRUE)

if(length(a)<1)
  stop("Provide data_processed folder and optional output folder")

source_path <- a[1]
dest <- if(length(a)>=2) a[2] else "R16_results"
dir.create(dest,showWarnings=FALSE,recursive=TRUE)

read_original <- function(name) {
  if(dir.exists(source_path)) {
    f <- file.path(source_path,name)
    if(!file.exists(f)) stop("Missing file: ",f)
    return(read.csv(f,check.names=FALSE,
                    stringsAsFactors=FALSE))
  }
  if(!file.exists(source_path))
    stop("Input not found: ",source_path)
  n <- paste0("data_processed/",name)
  if(!(n %in% unzip(source_path,list=TRUE)$Name))
    stop("Missing archive entry: ",n)
  read.csv(unz(source_path,n),check.names=FALSE,
           stringsAsFactors=FALSE)
}

d <- read_original(
  "09_targeted_integrated_prioritization.csv"
)
ab <- read_original("15_benchmark_ablation.csv")

stopifnot(
  nrow(d)==205,
  !anyDuplicated(d$gene_id),
  sum(d$is_external_benchmark)==6
)

domains <- c(
  "biological_discovery_percentile",
  "pu_consensus_percentile",
  "crispr_tractability_percentile"
)

mat <- as.matrix(d[,domains])
storage.mode(mat) <- "double"
valid <- is.finite(mat)

calculate <- function(w) {
  stopifnot(length(w)==3,all(w>=0),
            abs(sum(w)-1)<1e-8)
  numer <- rowSums(sweep(
    ifelse(valid,mat,0),2,w,"*"
  ))
  denom <- rowSums(sweep(valid*1,2,w,"*"))
  raw <- numer/denom
  raw[denom<1e-12] <- NA_real_
  raw*(0.8+0.2*rowMeans(valid))
}

rank_genes <- function(s) {
  r <- rep(NA_integer_,length(s))
  good <- which(is.finite(s))
  ord <- good[order(-s[good],d$gene_id[good])]
  r[ord] <- seq_along(ord)
  r
}

weights <- list(
  biology_led=c(.60,.25,.15),
  primary=c(.50,.30,.20),
  validation_led=c(.40,.35,.25)
)

rr <- lapply(weights,function(w)
  rank_genes(calculate(w))
)

stopifnot(
  max(abs(calculate(weights$primary)-
          d$final_integrated_score),
      na.rm=TRUE)<1e-12,
  identical(as.integer(rr$primary),
            as.integer(d$final_integrated_rank)),
  identical(as.integer(rr$biology_led),
            as.integer(d$rank_biology_led)),
  identical(as.integer(rr$validation_led),
            as.integer(d$rank_validation_led))
)

med <- apply(do.call(cbind,rr),1,median)
ord <- order(
  med,-calculate(weights$primary),d$gene_id
)
final_rank <- integer(nrow(d))
final_rank[ord] <- seq_len(nrow(d))

stopifnot(identical(
  final_rank,
  as.integer(d$final_prepopulation_rank)
))

cat("PASS: All 205 original gene ranks reproduced.\n")

b <- d[d$is_external_benchmark,]

write.csv(
  b,
  file.path(dest,"R16_original_six_landmarks.csv"),
  row.names=FALSE
)

sc <- expand.grid(
  biology_units=0:20,
  pu_units=0:20
)
sc <- sc[
  sc$biology_units+sc$pu_units<=20,
]
sc$crispr_units <- 20-
  sc$biology_units-sc$pu_units

stopifnot(nrow(sc)==231)

res <- vector("list",nrow(sc))

for(i in seq_len(nrow(sc))) {
  w <- unlist(sc[
    i,c("biology_units","pu_units","crispr_units")
  ])/20
  ranks <- rank_genes(calculate(w))

  res[[i]] <- data.frame(
    scenario_id=sprintf(
      "B%02d_P%02d_C%02d",
      sc$biology_units[i],
      sc$pu_units[i],
      sc$crispr_units[i]
    ),
    gene_id=b$gene_id,
    benchmark_name=b$benchmark_name,
    benchmark_class=b$benchmark_class,
    biology_weight=w[1],
    pu_weight=w[2],
    crispr_weight=w[3],
    rank=ranks[match(b$gene_id,d$gene_id)]
  )
}

long <- do.call(rbind,res)

write.csv(
  long,
  file.path(
    dest,"R16_231_weight_grid_landmark_ranks.csv"
  ),
  row.names=FALSE
)

summ <- do.call(rbind,
  lapply(split(long,long$gene_id),function(x) {
    z <- x$rank[!is.na(x$rank)]
    data.frame(
      gene_id=x$gene_id[1],
      n_scenarios=nrow(x),
      min_rank=min(z),
      median_rank=median(z),
      max_rank=max(z),
      fraction_top50=mean(z<=50),
      fraction_top100=mean(z<=100)
    )
  })
)

summ$original_rank <-
  d$final_prepopulation_rank[
    match(summ$gene_id,d$gene_id)
  ]

summ <- summ[order(summ$original_rank),]

write.csv(
  summ,
  file.path(
    dest,"R16_231_weight_grid_landmark_summary.csv"
  ),
  row.names=FALSE
)

write.csv(
  ab,
  file.path(
    dest,"R16_original_step15_benchmark_ablations.csv"
  ),
  row.names=FALSE
)

labels <- c(
  "AGAP007280 (fertility)",
  "AGAP005958 (fertility)",
  "AGAP007031 (FREP1)",
  "AGAP011377 (fertility)",
  "AGAP004050 (doublesex)",
  "AGAP028779 (mir-184)"
)

genes <- b$gene_id[
  order(b$final_prepopulation_rank)
]

# Figure 1: landmark rank sensitivity
png(
  file.path(
    dest,"R16_landmark_weight_sensitivity_600dpi.png"
  ),
  width=3600,height=2700,res=600,type="cairo"
)

par(mar=c(4.7,13,3.5,2),cex.axis=.85)

plot(
  NA,xlim=c(1,205),ylim=c(.5,6.5),
  yaxt="n",
  xlab="Rank among 205 candidates (lower is better)",
  ylab="",
  main="Landmark sensitivity across 231 weight settings"
)

abline(v=c(20,50,100),
       lty=c(3,2,2),col="grey75")

for(i in seq_along(genes)) {
  g <- genes[i]
  y <- 7-i
  vals <- long$rank[long$gene_id==g]
  vals <- vals[is.finite(vals)]

  points(
    vals,rep(y,length(vals)),
    pch=16,cex=.12,
    col=rgb(.18,.35,.65,.20)
  )

  points(
    d$final_prepopulation_rank[d$gene_id==g],
    y,pch=18,cex=1.25,col="#AD4935"
  )
}

axis(2,at=6:1,labels=labels,
     las=1,cex.axis=.76)

legend(
  "bottomright",
  legend=c("Alternative weights","Original rank"),
  pch=c(16,18),
  col=c("#6484AA","#AD4935"),
  bty="n",cex=.77
)

dev.off()

# Figure 2: original evidence-domain scores
png(
  file.path(
    dest,"R16_landmark_domain_percentiles_600dpi.png"
  ),
  width=3500,height=2600,res=600,type="cairo"
)

par(mar=c(5,13,3.5,2),cex.axis=.85)

plot(
  NA,xlim=c(0,1.03),ylim=c(.5,6.5),
  yaxt="n",
  xlab="Original evidence-domain percentile",
  ylab="",
  main="Evidence-domain support for six landmarks"
)

cl <- c("#2F6F9E","#BE7032","#45815C")
pchs <- c(16,17,15)

for(i in seq_along(genes)) {
  ix <- match(genes[i],d$gene_id)
  y <- 7-i

  for(j in 1:3) {
    if(is.finite(mat[ix,j])) {
      points(
        mat[ix,j],y+(.16*(j-2)),
        pch=pchs[j],col=cl[j],cex=1.1
      )
    }
  }
}

axis(2,at=6:1,labels=labels,
     las=1,cex.axis=.76)

legend(
  "bottomleft",
  legend=c(
    "Biology",
    "PU comparative-label score",
    "Reference CRISPR tractability"
  ),
  pch=pchs,col=cl,bty="n",cex=.77
)

dev.off()

cat("\nLANDMARK WEIGHT-SENSITIVITY SUMMARY\n")
print(summ,row.names=FALSE,digits=3)

cat("\nSUCCESS: Created four CSV tables and",
    "two 600-dpi R figures in",dest,"\n")
cat("Interpretation: descriptive sensitivity,",
    "not independent experimental validation.\n")
