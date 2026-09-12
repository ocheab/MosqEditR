source("R/helpers.R")
suppressPackageStartupMessages({library(biomaRt); library(data.table); library(stringr); library(tidyr); library(curl)})
log_step("03", "Acquiring Anopheles-Drosophila orthology and pinned FlyBase phenotype labels")

# ---- Anopheles -> D. melanogaster orthology ----
mart <- useEnsembl(biomart="metazoa_mart", host="https://metazoa.ensembl.org")
ds <- listDatasets(mart)
ag_ds <- ds$dataset[grepl("Anopheles gambiae", ds$description, ignore.case=TRUE)][1]
if (is.na(ag_ds)) stop("Anopheles gambiae dataset not found.")
ag <- useDataset(ag_ds, mart=mart)
a <- listAttributes(ag)
hom_gene <- a$name[grepl("dmel.*homolog.*ensembl_gene", a$name, ignore.case=TRUE)][1]
hom_name <- a$name[grepl("dmel.*homolog.*associated_gene_name", a$name, ignore.case=TRUE)][1]
if (is.na(hom_gene)) stop("Could not identify D. melanogaster homolog attributes; inspect listAttributes().")
attrs <- unique(c("ensembl_gene_id", hom_gene, hom_name[!is.na(hom_name)]))
orth <- getBM(attrs, mart=ag)
setDT(orth)
setnames(orth, 1, "gene_id")
setnames(orth, 2, "dmel_ortholog")
if (ncol(orth) >= 3) setnames(orth, 3, "dmel_gene_name") else orth[, dmel_gene_name := NA_character_]
orth <- orth[grepl("^AGAP", gene_id)]
fwrite(orth, "data_processed/03_anopheles_dmel_orthology.csv")

# ---- FlyBase phenotype evidence ----
# Pin the release used for this analysis. Override only deliberately, and record the replacement in provenance.
fb_release <- Sys.getenv("MOSQEDIT_FLYBASE_RELEASE", unset="FB2026_02")
fb_tag <- tolower(sub("^FB", "", fb_release))
fb_tag <- gsub("_", "_", fb_tag)
phen_file <- file.path("data_raw", paste0("genotype_phenotype_data_fb_", fb_tag, ".tsv.gz"))
map_file  <- file.path("data_raw", paste0("fbal_to_fbgn_fb_", fb_tag, ".tsv.gz"))
base_url <- paste0("https://s3ftp.flybase.org/releases/", fb_release, "/precomputed_files/alleles/")
phen_url <- paste0(base_url, basename(phen_file))
map_url  <- paste0(base_url, basename(map_file))

for (z in list(c(phen_url, phen_file), c(map_url, map_file))) {
  if (!file.exists(z[2])) {
    message("Downloading pinned FlyBase file: ", z[1])
    curl_download(z[1], destfile=z[2], quiet=FALSE)
  }
}
writeLines(c(paste0("FlyBase_release\t", fb_release), paste0("phenotype_url\t", phen_url), paste0("allele_gene_url\t", map_url)),
           "logs/03_flybase_release.tsv")

# FlyBase precomputed files are commented TSVs and schemas can evolve. Parse defensively from row text.
read_fb <- function(path) fread(path, sep="\t", header=TRUE, comment.char="#", fill=TRUE, quote="", data.table=TRUE)
ph <- read_fb(phen_file)
am <- read_fb(map_file)
if (!nrow(ph) || !nrow(am)) stop("Pinned FlyBase files were downloaded but contained no data rows.")

ph_txt <- names(ph)[vapply(ph, function(x) is.character(x) || is.factor(x), logical(1))]
am_txt <- names(am)[vapply(am, function(x) is.character(x) || is.factor(x), logical(1))]
ph[, phenotype_blob := do.call(paste, c(lapply(.SD, as.character), sep=" | ")), .SDcols=ph_txt]
am[, map_blob := do.call(paste, c(lapply(.SD, as.character), sep=" | ")), .SDcols=am_txt]

# Each genotype may contain multiple alleles. Expand all FBal identifiers, then map alleles to FBgn genes.
ph_long <- ph[, .(phenotype_blob, fbal=str_extract_all(phenotype_blob, "FBal\\d+")[[1]]), by=seq_len(nrow(ph))]
ph_long <- ph_long[!is.na(fbal) & nzchar(fbal)]
am[, fbal := str_extract(map_blob, "FBal\\d+")]
am[, fbgn := str_extract(map_blob, "FBgn\\d+")]
am <- unique(am[!is.na(fbal) & !is.na(fbgn), .(fbal, fbgn)])
ph_long <- merge(ph_long, am, by="fbal", all.x=FALSE, all.y=FALSE)
if (!nrow(ph_long)) stop("Could not map FlyBase genotype phenotype rows to FBgn identifiers; inspect pinned file schemas.")

sterile_re <- paste(c("female sterile", "female sterility", "grandchildless", "maternal[- ]effect",
                      "oogenesis defective", "egg laying defective", "reduced female fertility"), collapse="|")
repro_re <- paste(c(sterile_re, "fertility", "fecundity", "ovary", "oogenesis", "egg", "germline"), collapse="|")
ph_long[, female_sterile_hit := grepl(sterile_re, phenotype_blob, ignore.case=TRUE)]
ph_long[, reproductive_hit := grepl(repro_re, phenotype_blob, ignore.case=TRUE)]
fb_ev <- ph_long[, .(flybase_female_sterile=any(female_sterile_hit, na.rm=TRUE),
                     flybase_reproductive_phenotype_count=sum(reproductive_hit, na.rm=TRUE)), by=fbgn]

gene_ev <- merge(orth, fb_ev, by.x="dmel_ortholog", by.y="fbgn", all.x=TRUE)
gene_ev <- gene_ev[, .(
  dmel_ortholog=first_nonempty(dmel_ortholog),
  dmel_gene_name=first_nonempty(dmel_gene_name),
  flybase_female_sterile=if (all(is.na(flybase_female_sterile))) NA else any(flybase_female_sterile, na.rm=TRUE),
  flybase_reproductive_phenotype_count=if (all(is.na(flybase_reproductive_phenotype_count))) NA_integer_ else sum(flybase_reproductive_phenotype_count, na.rm=TRUE)
), by=gene_id]
fwrite(gene_ev, "data_processed/03_gene_orthology_phenotypes.csv")
write_checksum(c("data_processed/03_anopheles_dmel_orthology.csv", "data_processed/03_gene_orthology_phenotypes.csv", phen_file, map_file), "logs/03_checksums.tsv")

