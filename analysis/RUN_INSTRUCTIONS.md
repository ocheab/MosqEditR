# MosqEdit-R execution guide

## 1. Prerequisites
- R 4.4 or newer recommended.
- Python 3.10+ with packages in `requirements-python.txt`.
- Network access to Ensembl Metazoa, NCBI GEO, FlyBase, and MalariaGEN cloud data.

## 2. Create the Python environment
Example with a virtual environment:

```bash
python -m venv .venv
. .venv/bin/activate          # Windows PowerShell: .venv\\Scripts\\Activate.ps1
pip install -r requirements-python.txt
```

Set `MOSQEDIT_PYTHON` to the absolute path of that environment's Python executable before running R.

## 3. Smoke-test Ag1000G before the full genome-wide run
The Ag1000G routine is resumable. Start with a small validation set:

```bash
export MOSQEDIT_MAX_GENES=20
```

Run steps 00-04, inspect `data_raw/ag3_gene_variation_cache.csv`, then unset the variable for the full run.

## 4. Run the R pipeline
From the project root:

```r
source("R/00_setup.R")
source("R/run_all.R")
```

By default the FlyBase phenotype layer is pinned to `FB2026_02`. A different release can be used only deliberately via `MOSQEDIT_FLYBASE_RELEASE`, and the release is written to `logs/03_flybase_release.tsv`.

## 5. Pre-submission gates
Do not treat the current PENDING candidate file as a result. Before submission, confirm:
- all acquisition steps completed without silent fallback;
- the exact releases/API versions are recorded;
- known Anopheles landmarks were withheld from label construction;
- `results/candidate_lists/ranked_targets.csv` contains only MODEL-generated ranks;
- model diagnostics, top-k validation, and bootstrap stability were generated;
- `R/99_capture_session.R` has replaced the bootstrap `renv.lock` and `sessionInfo.txt`;
- the manuscript's locked Results section was updated only from code-generated outputs.

