# Publishing MosqEditR on GitHub

## 1. Edit repository metadata first

Replace these placeholders:

- `ocheab` in `DESCRIPTION`, `README.md`, `_pkgdown.yml`, and `CITATION.cff`.
- `ocheab1@gmail.com` in `DESCRIPTION`.
- authors/contributors in `DESCRIPTION` and `CITATION.cff`.
- DOI fields after Zenodo/manuscript DOI assignment.
- license if MIT is not the agreed project license.

## 2. Create and push the repository

From a terminal in this folder:

```bash
git init
git add .
git commit -m "Initial public MosqEditR package"
git branch -M main
git remote add origin https://github.com/ocheab/MosqEditR.git
git push -u origin main
```

Or with GitHub CLI:

```bash
gh repo create MosqEditR --public --source=. --remote=origin --push
```

## 3. Users install it in R

```r
install.packages("pak")
pak::pak("https://github.com/ocheab/MosqEditR")
```

or:

```r
install.packages("remotes")
remotes::install_github("https://github.com/ocheab/MosqEditR")
```

## 4. Recommended release workflow

- Tag the first public package as `v0.1.0`.
- Connect the repository to Zenodo if you want a citable software DOI.
- Upload the large frozen manuscript outputs as a GitHub Release asset and/or Zenodo dataset rather than package data.
- Enable GitHub Pages from the `gh-pages` branch if using the included pkgdown action.

## 5. Local package checks

The repository includes GitHub Actions for R CMD check. Locally, after installing development tools:

```r
install.packages(c("devtools", "testthat", "roxygen2"))
devtools::document()
devtools::test()
devtools::check()
```

This package bundle was assembled in an environment without an R runtime, so the first GitHub Actions run should be treated as the authoritative package-build check.

