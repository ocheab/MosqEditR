# MosqEditR

**MosqEditR** is an installable R package plus reproducible analysis templates for explainable prioritization of candidate gene-editing targets in *Anopheles gambiae*.

The framework keeps distinct evidence layers separate:

1. biological discovery evidence (adult tissue expression and comparative functional evidence),
2. leakage-safe positive-unlabeled (PU) model support,
3. exon-aware reference-sequence CRISPR tractability, and
4. optional targeted population-genomic robustness when real population data are available.

> Population-genomic metrics are never fabricated or silently imputed. The public manuscript workflow treats them as a second-stage validation layer.

## Installation from GitHub

```r
install.packages("pak")
pak::pak("https://github.com/ocheab/MosqEditR")
```

or:

```r
install.packages("remotes")
remotes::install_github("https://github.com/ocheab/MosqEditR")
```

Then:

```r
library(MosqEditR)
mosqedit_version()
```

## Quick example

```r
score <- mosqedit_integrated_score(
  biology = c(0.99, 0.82, 0.61),
  pu      = c(0.98, 0.75, 0.80),
  crispr  = c(0.95, 0.91, 0.40)
)
score

mosqedit_rank(
  score$adjusted_score,
  gene_id = c("AGAP_A", "AGAP_B", "AGAP_C")
)
```

A small frozen manuscript example is bundled:

```r
head(mosqedit_example_candidates())
```

## Create a complete analysis workspace

```r
mosqedit_init("MosqEditR-analysis")
setwd("MosqEditR-analysis")
mosqedit_pipeline_steps()

# Install the optional dependencies needed by the full workflow:
mosqedit_install_pipeline_deps()
```

Run one step:

```r
mosqedit_run_step("04")
```

Or a range:

```r
mosqedit_run(from = "00", to = "04")
```

Later stages require their upstream files. Several acquisition steps require internet access and optional Bioconductor/CRAN packages. Step 14 intentionally stops unless real targeted population-genomic metrics exist.

## Repository layout

```text
MosqEditR/
DESCRIPTION
NAMESPACE
R/                         # installable package functions
man/                       # help pages
tests/testthat/            # unit tests
vignettes/                 # tutorials
inst/
extdata/                   # small examples/landmarks
pipeline/R/                # manuscript workflow templates
templates/                 # config + metadata templates
analysis/                  # manuscript reproducibility material
data-raw/                  # package-data preparation scripts
.github/workflows/         # R CMD check + pkgdown
CITATION.cff
CONTRIBUTING.md
â””â”€â”€ _pkgdown.yml
```

## Scientific interpretation

PU scores are ranking scores, not calibrated probabilities of biological success. Unlabeled genes are not confirmed negatives. The Drosophila/FlyBase label system introduces ascertainment structure that should be discussed when orthology predictors are used. CRISPR scores quantify reference-sequence tractability and do not substitute for population conservation, genome-wide off-target analysis, or experimental validation.

## Manuscript reproducibility

The manuscript analysis used an AgamP4 gene-like universe, MozAtlas/GSE21689 expression, Ensembl/FlyBase comparative evidence, bagged PU learning, exon-aware SpCas9 reference target assessment, explainability, and PU-resampling rank uncertainty. Large frozen result files should be deposited in a release/archival repository rather than installed as package data.

## Citation

After installation:

```r
citation("MosqEditR")
```

## License

MIT by default in this repository template. Change it before release if your institution or collaborators require another license.

