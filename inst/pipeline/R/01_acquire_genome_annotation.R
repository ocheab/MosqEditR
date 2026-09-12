source("R/helpers.R")

required <- c("data.table", "curl")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required package(s): ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages({library(data.table); library(curl)})

log_step("01", "Acquiring release-pinned AgamP4 annotation from Ensembl Metazoa release 63")

dir.create("data_raw/ensembl_release63_agamp4", recursive = TRUE, showWarnings = FALSE)
dir.create("data_processed", recursive = TRUE, showWarnings = FALSE)

url <- paste0(
  "https://ftp.ebi.ac.uk/ensemblgenomes/pub/metazoa/release-63/gff3/",
  "anopheles_gambiae/Anopheles_gambiae.AgamP4.63.gff3.gz"
)
dest <- "data_raw/ensembl_release63_agamp4/Anopheles_gambiae.AgamP4.63.gff3.gz"
if (!file.exists(dest)) curl::curl_download(url, dest, quiet = FALSE)

gff <- data.table::fread(
  cmd = paste("gzip -dc", shQuote(dest)),
  sep = "\t", header = FALSE, comment.char = "#", fill = TRUE,
  col.names = c("seqid","source","type","start","end","score","strand","phase","attributes")
)

attr_value <- function(x, key) {
  pat <- paste0("(?:^|;)", key, "=([^;]+)")
  m <- regexec(pat, x, perl = TRUE)
  z <- regmatches(x, m)
  vapply(z, function(y) if (length(y) >= 2L) URLdecode(y[[2]]) else NA_character_, character(1))
}

# All gene-like rows are retained, matching the frozen 13,845-locus universe.
gene_rows <- gff[grepl("gene$", type, ignore.case = TRUE)]
gene_rows[, gene_id := attr_value(attributes, "ID")]
gene_rows[, gene_id := sub("^gene:", "", gene_id)]
gene_rows[, gene_symbol := attr_value(attributes, "Name")]
gene_rows[, biotype := attr_value(attributes, "biotype")]
gene_rows[is.na(biotype), biotype := attr_value(attributes, "gene_biotype")]
gene_rows <- gene_rows[grepl("^AGAP", gene_id)]

ann <- gene_rows[, .(
  gene_id,
  gene_symbol,
  chromosome = seqid,
  gene_start = as.integer(start),
  gene_end = as.integer(end),
  strand,
  feature_type = type,
  biotype
)]
setorder(ann, chromosome, gene_start, gene_id)

# Transcript audit table.
trans <- gff[grepl("transcript$|mRNA$|RNA$", type, ignore.case = TRUE)]
trans[, transcript_id := sub("^transcript:", "", attr_value(attributes, "ID"))]
trans[, parent_gene := sub("^gene:", "", attr_value(attributes, "Parent"))]
trans <- trans[grepl("^AGAP", parent_gene)]
trans_out <- trans[, .(
  gene_id = parent_gene,
  transcript_id,
  transcript_type = type,
  chromosome = seqid,
  transcript_start = as.integer(start),
  transcript_end = as.integer(end),
  strand
)]

counts <- trans_out[, .(transcript_count = uniqueN(transcript_id)), by = gene_id]
ann <- merge(ann, counts, by = "gene_id", all.x = TRUE, sort = FALSE)
ann[is.na(transcript_count), transcript_count := 0L]

fwrite(trans_out, "data_raw/ensembl_anopheles_annotation_transcripts.csv")
fwrite(ann, "data_processed/01_gene_annotation.csv")

prov <- data.table(
  item = c("source", "release", "genome_build", "universe_rule", "gene_count"),
  value = c(url, "Ensembl Metazoa 63", "AgamP4", "GFF3 rows whose feature type ends in gene and ID begins AGAP", as.character(nrow(ann)))
)
fwrite(prov, "data_processed/01_gene_annotation_provenance.csv")

qc <- data.table(metric = c("AGAP loci", "Duplicate gene IDs"), value = c(nrow(ann), anyDuplicated(ann$gene_id)))
fwrite(qc, "data_processed/01_gene_annotation_qc.csv")
write_checksum(c(dest, "data_processed/01_gene_annotation.csv", "data_processed/01_gene_annotation_provenance.csv", "data_processed/01_gene_annotation_qc.csv"), "logs/01_checksums.tsv")

cat("\nMOSQEDIT-R STEP 01 COMPLETED\nAGAP loci: ", nrow(ann), "\n", sep = "")
log_step("01", "Annotation acquisition completed")

