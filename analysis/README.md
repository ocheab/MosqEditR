# Manuscript-analysis area

The installable package API lives in the repository-level `R/` directory.
The manuscript workflow is intentionally kept separate so package installation
does not execute downloads or long-running analyses.

- `manuscript/original_source_snapshot/` preserves the earlier uploaded analysis scripts for audit.
- `MosqEditR::mosqedit_init()` copies the curated current pipeline shipped in `inst/pipeline/R/` into a fresh analysis workspace.
- Large frozen result matrices are not committed as package data. Publish them as GitHub Release assets and/or Zenodo records.
- Population-genomic Step 14 must only be run after real site-level population metrics have been populated.

