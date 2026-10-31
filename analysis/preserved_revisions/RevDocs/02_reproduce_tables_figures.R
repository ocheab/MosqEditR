#!/usr/bin/env Rscript
# Uses audited, executed result CSVs and renders publication-quality 600 dpi figures.
# Numerical sensitivity is computed in 01_audit_sensitivity.py. No Python figures.
args <- commandArgs(trailingOnly=TRUE)
base <- if (length(args)>=1) args[1] else "."
res <- file.path(base,"results")
stopifnot(file.exists(file.path(res,"R2Q3_step04_domain_weight_sensitivity.csv")))
w <- read.csv(file.path(res,"R2Q3_step04_domain_weight_sensitivity.csv"),check.names=FALSE)
d <- read.csv(file.path(res,"R2Q3_expression_masking_replicates.csv"),check.names=FALSE)
c <- read.csv(file.path(res,"R2Q3_expression_coverage_by_cohort.csv"),check.names=FALSE)
stopifnot(all(c("expression_weight","top100_overlap_n","top20_overlap_n") %in% names(w)), nrow(w)==8)
stopifnot(nrow(d)==450, all(c("fraction_of_available_expression_masked","top100_overlap_n") %in% names(d)))
stopifnot(any(c$cohort=="Whole annotated gene universe"))

png(file.path(res,"Figure_R2Q3_Weight_Sensitivity_600dpi.png"),width=3900,height=2400,res=600,pointsize=10,type="cairo")
par(mar=c(4.2,4.2,2.3,0.6),family="sans")
x <- w$expression_weight * 100
plot(x,w$top100_overlap_n,type="b",pch=19,col="#185A9D",lwd=2,ylim=c(0,105),xlim=c(0,100),
     xlab="Weight assigned to the MozAtlas expression domain (%)",
     ylab="Original preliminary top 100 retained (n)", main="Step 04 biological-ranking sensitivity")
abline(v=60,col="gray40",lty=2)
points(60,100,pch=21,bg="#B0502B",cex=1.4)
text(60,90,"Original: 60%",pos=4,cex=0.83)
mtext("Fixed original FRIS and functional scores; not a PU model refit or new RNA-seq validation",side=1,line=3.1,cex=0.73)
dev.off()

png(file.path(res,"Figure_R2Q3_Expression_Masking_600dpi.png"),width=3900,height=2400,res=600,pointsize=10,type="cairo")
par(mar=c(4.4,4.2,2.3,0.6),family="sans")
boxplot(top100_overlap_n~fraction_of_available_expression_masked,data=d,
        names=c("10%","25%","50%"),ylim=c(0,105),col=c("#ACD5E7","#78B2D2","#4A8EB8"),
        outline=FALSE,xlab="Fraction of originally available expression profiles masked",
        ylab="Original preliminary top 100 retained (n)",
        main="Illustrative expression-availability stress test")
mtext("150 random masks per fraction; ranks recalculated with original completeness adjustment",side=1,line=3.1,cex=0.73)
dev.off()

write.csv(w,file.path(res,"R2Q3_R_reproduced_weight_table.csv"),row.names=FALSE)
write.csv(c,file.path(res,"R2Q3_R_reproduced_coverage_table.csv"),row.names=FALSE)
cat("DONE: two 600 dpi figures; eight weighting scenarios; 450 masking replicates\n")
print(w[,c("expression_weight","spearman_common_with_original","top20_overlap_n","top100_overlap_n")],row.names=FALSE)
