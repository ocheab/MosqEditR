#' List manuscript-pipeline steps shipped with MosqEditR
#'
#' @return A data.frame describing the available analysis scripts.
#' @export
mosqedit_pipeline_steps <- function() {
  data.frame(
    step = c("00", "01", "02", "03", "04", "05a", "05b", "06", "07", "08",
             "09", "10", "11", "12", "13", "14", "15", "16", "17"),
    file = c(
      "00_setup.R",
      "01_acquire_genome_annotation.R",
      "02_acquire_mozatlas_expression.R",
      "03_acquire_orthology_phenotypes.R",
      "04_preliminary_genomewide_ranking.R",
      "05a_crispr_reference_scan.R",
      "05b_crispr_exon_aware_validation.R",
      "06_build_feature_matrix.R",
      "07_build_leakage_safe_modelling_sets.R",
      "08_train_bagged_pu_models.R",
      "09_integrate_pu_crispr_final_prioritization.R",
      "10_explain_final_ranking_and_figures.R",
      "11_model_specific_explainability.R",
      "12_bootstrap_rank_uncertainty.R",
      "13_freeze_population_validation_panel.R",
      "14_integrate_population_genomics.R",
      "15_ablation_and_validation.R",
      "16_build_publication_tables_and_figures.R",
      "17_freeze_reproducible_release_bundle.R"
    ),
    description = c(
      "Create analysis directories and verify package availability",
      "Acquire release-pinned AgamP4 annotation",
      "Acquire and summarize MozAtlas expression",
      "Acquire orthology and FlyBase phenotype evidence",
      "Preliminary genome-wide biological ranking",
      "Targeted transcript-level SpCas9 reference scan",
      "Exon-aware CRISPR target refinement",
      "Build auditable master feature matrix",
      "Build leakage-safe PU modelling sets",
      "Repeated cross-fitted bagged PU modelling",
      "Integrate biology, PU support, and CRISPR tractability",
      "Explain integrated ranking and produce figures",
      "Cross-fitted model-specific explainability",
      "PU-resampling rank uncertainty",
      "Freeze targeted population-validation panel",
      "Integrate real population-genomic target-site metrics",
      "Evidence-domain ablation and validation",
      "Assemble publication tables and figures",
      "Freeze reproducible analysis release"
    ),
    stringsAsFactors = FALSE
  )
}

#' Initialize a MosqEdit-R analysis project
#'
#' Creates a manuscript-analysis workspace and copies the pipeline templates,
#' configuration, and small metadata files bundled with the installed package.
#'
#' @param path Destination directory.
#' @param overwrite Logical; allow replacing existing template files.
#' @return Invisibly returns the normalized project path.
#' @export
mosqedit_init <- function(path = "MosqEditR-analysis", overwrite = FALSE) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE)) && !overwrite) {
    stop("Destination is not empty. Use overwrite = TRUE or choose a new directory.", call. = FALSE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  for (d in c("R", "metadata", "data_raw", "data_processed", "figures", "logs", "results", "publication", "release")) {
    dir.create(file.path(path, d), recursive = TRUE, showWarnings = FALSE)
  }

  pipeline <- system.file("pipeline", "R", package = "MosqEditR")
  templates <- system.file("templates", package = "MosqEditR")
  if (!nzchar(pipeline)) stop("Installed package does not contain pipeline templates.", call. = FALSE)

  pipe_files <- list.files(pipeline, full.names = TRUE)
  file.copy(pipe_files, file.path(path, "R", basename(pipe_files)), overwrite = overwrite)

  if (nzchar(templates)) {
    cfg <- file.path(templates, "analysis_config.yml")
    if (file.exists(cfg)) file.copy(cfg, file.path(path, "analysis_config.yml"), overwrite = overwrite)
    mdir <- file.path(templates, "metadata")
    if (dir.exists(mdir)) {
      mfiles <- list.files(mdir, full.names = TRUE)
      file.copy(mfiles, file.path(path, "metadata", basename(mfiles)), overwrite = overwrite)
    }
  }

  readme <- c(
    "# MosqEdit-R analysis workspace",
    "",
    "Created by `MosqEditR::mosqedit_init()`.",
    "",
    "List steps with `MosqEditR::mosqedit_pipeline_steps()`.",
    "Run one step with `MosqEditR::mosqedit_run_step(\"04\", project = \".\")`.",
    "",
    "Population-genomic integration (Step 14) intentionally stops unless real site-level metrics are supplied."
  )
  writeLines(readme, file.path(path, "README.md"))
  invisible(path)
}

#' Validate a MosqEdit-R project skeleton
#'
#' @param project Project root.
#' @return A data.frame of required-path checks.
#' @export
mosqedit_validate_project <- function(project = ".") {
  required <- c("R", "metadata", "data_raw", "data_processed", "figures", "logs")
  paths <- file.path(project, required)
  data.frame(
    component = required,
    path = paths,
    exists = file.exists(paths) | dir.exists(paths),
    stringsAsFactors = FALSE
  )
}

.resolve_step <- function(step) {
  tab <- mosqedit_pipeline_steps()
  step <- as.character(step)
  hit <- which(tab$step == step | tab$file == step)
  if (length(hit) != 1L) {
    stop(
      "Unknown or ambiguous step. Use mosqedit_pipeline_steps() to list valid steps.",
      call. = FALSE
    )
  }
  tab[hit, , drop = FALSE]
}

#' Run one MosqEdit-R manuscript-analysis step
#'
#' @param step Step identifier (for example, "04" or "12") or exact script name.
#' @param project Project root created by `mosqedit_init()`.
#' @param echo Passed to `source()`.
#' @return Invisibly returns the sourced script result.
#' @export
mosqedit_run_step <- function(step, project = ".", echo = FALSE) {
  spec <- .resolve_step(step)
  script <- file.path(project, "R", spec$file)
  if (!file.exists(script)) {
    stop(
      paste0("Step script is missing: ", script, "\nRun mosqedit_init() or copy the pipeline templates."),
      call. = FALSE
    )
  }
  old <- getwd()
  on.exit(setwd(old), add = TRUE)
  setwd(project)
  source(file.path("R", spec$file), local = new.env(parent = globalenv()), echo = echo)
}

#' Run a contiguous range of MosqEdit-R steps
#'
#' @param from First step identifier.
#' @param to Last step identifier.
#' @param project Project root.
#' @param echo Passed to `source()`.
#' @return Invisibly returns TRUE after successful completion.
#' @export
mosqedit_run <- function(from = "00", to = "13", project = ".", echo = FALSE) {
  tab <- mosqedit_pipeline_steps()
  i <- match(as.character(from), tab$step)
  j <- match(as.character(to), tab$step)
  if (is.na(i) || is.na(j) || i > j) {
    stop("from/to must identify an ordered range in mosqedit_pipeline_steps().", call. = FALSE)
  }
  for (k in seq.int(i, j)) {
    message("\n===== MosqEdit-R step ", tab$step[k], ": ", tab$description[k], " =====")
    mosqedit_run_step(tab$step[k], project = project, echo = echo)
  }
  invisible(TRUE)
}


#' Install optional dependencies for the full MosqEdit-R analysis pipeline
#'
#' The installed package itself is deliberately lightweight. This helper
#' installs the CRAN and Bioconductor packages required by the manuscript
#' acquisition, modelling, explainability, and figure-generation templates.
#'
#' @param ask Logical; passed to Bioconductor installation where supported.
#' @param update Logical; whether Bioconductor may update old packages.
#' @return Invisibly returns a list of requested CRAN and Bioconductor packages.
#' @export
mosqedit_install_pipeline_deps <- function(ask = FALSE, update = FALSE) {
  cran <- c(
    "data.table", "digest", "yaml", "curl", "stringr",
    "glmnet", "ranger", "xgboost", "e1071", "ggplot2"
  )
  bioc <- c("Biostrings", "GEOquery", "Biobase")

  cran_missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(cran_missing)) {
    utils::install.packages(cran_missing)
  }

  bioc_missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(bioc_missing)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      utils::install.packages("BiocManager")
    }
    BiocManager::install(bioc_missing, ask = ask, update = update)
  }

  invisible(list(cran = cran, bioconductor = bioc))
}

#' Locate a bundled example or metadata file
#'
#' @param name File name in `inst/extdata`.
#' @return Full installed path.
#' @export
mosqedit_extdata <- function(name) {
  path <- system.file("extdata", name, package = "MosqEditR")
  if (!nzchar(path)) stop("Unknown bundled extdata file: ", name, call. = FALSE)
  path
}

#' Load the bundled top-candidate example table
#'
#' @return A data.frame with a small frozen example from the manuscript analysis.
#' @export
mosqedit_example_candidates <- function() {
  utils::read.csv(mosqedit_extdata("example_top20_candidates.csv"), stringsAsFactors = FALSE)
}

#' Return the installed MosqEditR version
#'
#' @return A package_version object.
#' @export
mosqedit_version <- function() {
  utils::packageVersion("MosqEditR")
}

