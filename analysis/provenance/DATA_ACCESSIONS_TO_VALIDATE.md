# Upstream source-data accession and license verification

- GEO GSE21689 (MozAtlas), GPL1321 — archive raw GEO SOFT, cite accession, and record retrieval date/checksum.
- GEO GSE239600 (RNA-seq sensitivity) — record exact series matrix/count files, raw sequence coverage and source license before redistribution.
- Ensembl Metazoa release 63, *Anopheles gambiae* PEST AgamP4 GFF3/cDNA/ncRNA, orthology TSV — record exact URLs, release and SHA-256.
- FlyBase FB2026_02 phenotype/mapping files — verify redistribution and credit/citation terms; don't assume MIT package licence applies to these datasets.
- MalariaGEN Ag1000G **Phase 2 AR1** (1,142 wild mosquitoes analysed) — cite and pin phase/release and SNP subset; note no whole-genome dataset copy and no per-haplotype 23-bp validation. Do not silently replace with Ag3.0 Phase 3.
- AgamP4 genome toplevel FASTA for Bowtie specificity — original reference available from Ensembl (exact URI/version to be recorded); genome index is regenerated using `bowtie-build` rather than deposited as evidence. Reference SHA-256 is required for a genuinely identical recreation.

Before public posting, inspect provider reuse conditions for copied third-party source files. Derivative tables and original source code have different reuse rights.
