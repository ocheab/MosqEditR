source("R/helpers.R")
suppressPackageStartupMessages({library(biomaRt); library(data.table)})
log_step("01", "Querying Ensembl Metazoa/Anopheles annotation")
# Dataset names may change between Ensembl releases; discover rather than hard-code silently.
mart <- useEnsembl(biomart="metazoa_mart", host="https://metazoa.ensembl.org")
ds <- listDatasets(mart)
hit <- ds[grepl("Anopheles gambiae", ds$description, ignore.case=TRUE), , drop=FALSE]
if (!nrow(hit)) stop("Could not locate an Anopheles gambiae dataset in Ensembl Metazoa.")
mart_ag <- useDataset(hit$dataset[1], mart=mart)
attrs <- c("ensembl_gene_id","external_gene_name","chromosome_name","gene_biotype","ensembl_transcript_id")
ann <- getBM(attributes=attrs, mart=mart_ag)
setDT(ann)
setnames(ann, c("ensembl_gene_id","external_gene_name","chromosome_name","gene_biotype","ensembl_transcript_id"),
         c("gene_id","gene_symbol","chromosome","biotype","transcript_id"))
# Retain AGAP IDs and collapse transcripts.
ann <- ann[grepl("^AGAP", gene_id)]
genes <- ann[, .(gene_symbol=first_nonempty(gene_symbol), chromosome=first_nonempty(chromosome),
                 biotype=first_nonempty(biotype), transcript_count=uniqueN(transcript_id)), by=gene_id]
fwrite(ann, "data_raw/ensembl_anopheles_annotation_transcripts.csv")
fwrite(genes, "data_processed/01_gene_annotation.csv")
write_checksum(c("data_raw/ensembl_anopheles_annotation_transcripts.csv","data_processed/01_gene_annotation.csv"), "logs/01_checksums.tsv")

