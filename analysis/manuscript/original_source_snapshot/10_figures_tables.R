source("R/helpers.R")
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
log_step("10", "Generating publication tables and figures from computed outputs")
r <- fread("results/candidate_lists/ranked_targets.csv")
top <- r[1:min(50,.N)]
fwrite(top, "results/tables/top50_candidates.csv")
p <- ggplot(top[1:min(20,.N)], aes(x=reorder(gene_id, ensemble_score), y=ensemble_score)) +
  geom_col() + coord_flip() + labs(x=NULL, y="MosqEdit-R score", title="Top MosqEdit-R candidate genes") + theme_bw(base_size=11)
ggsave("figures/Figure4_top_candidates.pdf", p, width=7, height=6)
ggsave("figures/Figure4_top_candidates.png", p, width=7, height=6, dpi=600)
ggsave("figures/Figure4_top_candidates.tiff", p, width=7, height=6, dpi=600, compression="lzw")

