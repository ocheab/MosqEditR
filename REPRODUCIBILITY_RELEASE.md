# Reproducibility and archival status — MosqEditR

**Prepared: 2026-10-31. Status: LOCAL STAGING; NOT YET PUBLIC DEPOSITED.**

The `R/` directory contains installable package helpers; `inst/pipeline/R/` contains project workflow templates. Neither by itself constitutes the full executed manuscript environment.

The `analysis/preserved_revisions/` directory contains byte-preserved R/Python/shell scripts collected from the two user-supplied study archives. Their original source archive, member path and SHA-256 are recorded in `PRESERVED_CODE_SHA256.csv`. This is preservation, **not a claim that every script can reproduce every historical result from an empty directory without reconstruction**.

## Important historical discrepancies

1. The frozen 205 gene-level Step 05 scores and the original `05_crispr_editability_provenance.csv` document component weights of **0.35 / 0.25 / 0.25 / 0.15**, while `inst/pipeline/R/05b_crispr_exon_aware_validation.R` has a **0.40 / 0.25 / 0.20 / 0.15** formula. Preserve both; the independently validated `R14_scoring_sensitivity.R` reconstructs archived outputs exactly using documented original weights. Do not present the supplied 05b implementation as the exact executed source of the archived scores.
2. `analysis/manuscript/original_source_snapshot/05_crispr_editability.R` explicitly generates **placeholder NA** editability values. The other original-source-snapshot scripts describe an earlier unexecuted prototype, not the final research output.
3. `analysis/renv.lock` in the uploaded package lists **zero package entries** and cannot reproduce the full dependency state. Replace with an `renv::snapshot()` from the actually executed R environment and include session info, system versions (R, Bowtie1, samtools, Python), and required Python packages.
4. The uploaded `analysis/analysis_config.yml` says Ag1000G Phase 3 and 500 PU resamples, whereas the frozen processed data and extended analysis used Phase 2 AR1, three 5-fold repeats, and 2000 PU-resampling replicates. A release-pinned, *executed* configuration must be reconstructed, tested and supplied; do not silently retcon the older config.
5. The archives contain selected interval-level Ag1000G SNP results and PAM/protospacer alignment summaries, not an experimental cutting-efficiency dataset or phase-resolved 23-bp validation per mosquito.

## Separate data records (Zenodo recommended)

- `MosqEditR_Phase2_Targeted_Extraction1.zip` — exact raw acquisition snapshots, processed analysis matrices, model outputs, targeted variant calls and supporting outputs as uploaded.
- `MosqEditR_Revision_Outputs_No_Regenerable_Indexes.zip` — analysis evidence and derived outputs from `RevDocs1.zip`, minus regenerable Bowtie indexes and large AgamP4 FASTA. The excluded reference binary requires upstream accession/URL, version and reference SHA-256. See the exclusion manifest.

The original full archives' SHA-256 values are in the attached top-level `source_archive_sha256.json`; the member-level inventory is in `source_archives_file_manifest.csv`.

## Requirements for a claim of full reproducibility

- Source-data URLs/accessions, original release labels, licensing and checksums verified.
- Executed source script and parameter provenance reconciled for all stages; document reconstruction when exact script not preserved.
- Exact R package/environment lock and OS/system binary versions captured and independently tested.
- A clean machine should be able to install dependencies and replay full pipeline or explicitly documented frozen-result verification with a source-data retrieval map.
- Byte-identical frozen outputs, figure-specific source data, script-to-output table, tests and hashes linked to the manuscript.
- The GitHub commit/tag and both Zenodo record DOIs recorded **after** public publication and verification.

Nothing in this staging bundle is itself an issued DOI, immutable Git tag, or verified public GitHub commit.
