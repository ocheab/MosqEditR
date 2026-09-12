# Runtime requirements and pre-submission gates

The scripts in this directory are designed for an external R-enabled environment with network access. They have not been executed in the current build environment because R is unavailable here.

Before manuscript submission:
1. Run the full pipeline from `R/run_all.R` in a clean environment.
2. Pin the exact FlyBase release used for phenotype labels and record it in provenance.
3. Confirm the current `malariagen_data` API used by `R/ag1000g_gene_variation.py` and validate gene-level population-genetic calculations against a small known region.
4. Recompute or import gene-level CRISPR editability metrics rather than relying only on the Schmidt et al. historical genome-wide benchmark.
5. Confirm that all landmark targets are excluded from label construction before external hold-out validation.
6. Replace PENDING rows in `results/candidate_lists/ranked_targets.csv` only with code-generated values.
7. Recreate all final tables/figures and update the manuscript from code-generated outputs.
8. Run `R/99_capture_session.R` so `sessionInfo.txt` and `renv.lock` represent the executed environment.

