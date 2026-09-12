# =============================================================================
# 17_freeze_reproducible_release_bundle.R
# MosqEdit-R Manuscript 1
# STEP 17 â€” Freeze reproducible analysis release bundle
# PURE R
# =============================================================================
#
# Creates a final release directory and ZIP containing:
#   - R scripts
#   - metadata
#   - processed result tables
#   - publication tables/figures
#   - logs/session information/checksums
#
# It does NOT silently include large raw genomic files.
# =============================================================================

source("R/helpers.R")

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
}
suppressPackageStartupMessages(library(data.table))

log_step("17", "Freezing reproducible MosqEdit-R analysis release")

RELEASE_DIR <- "release/MosqEditR_Manuscript1_Analysis"
ZIP_FILE <- "release/MosqEditR_Manuscript1_Analysis.zip"

dir.create(RELEASE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("release", recursive = TRUE, showWarnings = FALSE)

copy_tree_files <- function(source_dir, destination_dir, pattern = NULL) {
    if (!dir.exists(source_dir)) {
        return(character())
    }

    files <- list.files(
        source_dir,
        pattern = pattern,
        recursive = TRUE,
        full.names = TRUE
    )

    if (length(files) == 0L) {
        return(character())
    }

    copied <- character()

    for (src in files) {
        rel <- substring(
            src,
            nchar(source_dir) + 2L
        )

        dest <- file.path(
            destination_dir,
            rel
        )

        dir.create(
            dirname(dest),
            recursive = TRUE,
            showWarnings = FALSE
        )

        ok <- file.copy(
            src,
            dest,
            overwrite = TRUE
        )

        if (isTRUE(ok)) {
            copied <- c(copied, dest)
        }
    }

    copied
}

# -------------------------------------------------------------------------
# 1. R scripts
# -------------------------------------------------------------------------

r_dest <- file.path(RELEASE_DIR, "R")
dir.create(r_dest, recursive = TRUE, showWarnings = FALSE)

r_files <- list.files(
    "R",
    pattern = "\\.[Rr]$",
    full.names = TRUE
)

for (src in r_files) {
    file.copy(
        src,
        file.path(r_dest, basename(src)),
        overwrite = TRUE
    )
}

# -------------------------------------------------------------------------
# 2. Metadata
# -------------------------------------------------------------------------

metadata_dest <- file.path(RELEASE_DIR, "metadata")
dir.create(metadata_dest, recursive = TRUE, showWarnings = FALSE)

if (dir.exists("metadata")) {
    metadata_files <- list.files(
        "metadata",
        full.names = TRUE
    )

    for (src in metadata_files[file.info(metadata_files)$isdir == FALSE]) {
        file.copy(
            src,
            file.path(metadata_dest, basename(src)),
            overwrite = TRUE
        )
    }
}

# -------------------------------------------------------------------------
# 3. Processed analysis outputs
# -------------------------------------------------------------------------

processed_dest <- file.path(RELEASE_DIR, "data_processed")
dir.create(processed_dest, recursive = TRUE, showWarnings = FALSE)

processed_patterns <- c(
    "^0[1-9]_",
    "^1[0-7]_",
    "^master_gene_feature_matrix\\.csv$"
)

processed_files <- list.files(
    "data_processed",
    full.names = TRUE
)

keep_processed <- processed_files[
    vapply(
        basename(processed_files),
        function(x) {
            any(
                vapply(
                    processed_patterns,
                    grepl,
                    logical(1),
                    x = x
                )
            )
        },
        logical(1)
    )
]

for (src in keep_processed[file.info(keep_processed)$isdir == FALSE]) {
    file.copy(
        src,
        file.path(processed_dest, basename(src)),
        overwrite = TRUE
    )
}

# -------------------------------------------------------------------------
# 4. Publication outputs
# -------------------------------------------------------------------------

publication_dest <- file.path(RELEASE_DIR, "publication")
copy_tree_files(
    "publication",
    publication_dest
)

# -------------------------------------------------------------------------
# 5. Logs
# -------------------------------------------------------------------------

logs_dest <- file.path(RELEASE_DIR, "logs")
dir.create(logs_dest, recursive = TRUE, showWarnings = FALSE)

log_files <- list.files(
    "logs",
    full.names = TRUE
)

for (src in log_files[file.info(log_files)$isdir == FALSE]) {
    file.copy(
        src,
        file.path(logs_dest, basename(src)),
        overwrite = TRUE
    )
}

# -------------------------------------------------------------------------
# 6. README
# -------------------------------------------------------------------------

population_available <- file.exists(
    "data_processed/14_final_population_validated_prioritization.csv"
)

readme <- c(
    "# MosqEdit-R Manuscript 1 reproducible analysis release",
    "",
    paste0("Created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "",
    "## Scope",
    "",
    "This release freezes the analysis outputs, scripts, metadata, publication tables, figures, and audit logs used for Manuscript 1.",
    "",
    "## Analysis architecture",
    "",
    "1. Step 01 â€” AgamP4 annotation universe",
    "2. Step 02 â€” MozAtlas expression",
    "3. Step 03 â€” Drosophila orthology/FlyBase phenotype evidence",
    "4. Step 04 â€” preliminary biological ranking",
    "5. Step 05 â€” exon-aware SpCas9 reference editability",
    "6. Step 06 â€” auditable master feature matrix",
    "7. Step 07 â€” leakage-safe PU modelling sets",
    "8. Step 08 â€” repeated cross-fitted bagged PU learning",
    "9. Step 09 â€” robust PU consensus + CRISPR-aware ranking",
    "10. Step 10 â€” integrated ranking explainability",
    "11. Step 11 â€” model-specific SHAP/permutation explainability",
    "12. Step 12 â€” PU-resampling rank uncertainty",
    "13. Step 13 â€” targeted population validation panel",
    "14. Step 14 â€” targeted population-genomic robustness integration",
    "15. Step 15 â€” ablation and validation analysis",
    "16. Step 16 â€” publication tables and figures",
    "17. Step 17 â€” frozen reproducible release",
    "",
    "## Population genomics",
    "",
    paste0(
        "Population-genomic integration status: ",
        if (population_available) "COMPLETED" else "NOT COMPLETED / PENDING REAL DATA"
    ),
    "",
    "Population metrics are never fabricated or silently inferred. Step 14 only runs after real site-level population-genomic metrics are supplied.",
    "",
    "## Raw data",
    "",
    "Large raw population-genomic files are intentionally not bundled. Refer to provenance files and source manifests for accession/version details."
)

writeLines(
    readme,
    file.path(RELEASE_DIR, "README.md"),
    useBytes = TRUE
)

# -------------------------------------------------------------------------
# 7. Final release manifest + checksums
# -------------------------------------------------------------------------

release_files <- list.files(
    RELEASE_DIR,
    recursive = TRUE,
    full.names = TRUE
)

manifest <- data.table(
    relative_path = substring(
        release_files,
        nchar(RELEASE_DIR) + 2L
    ),
    size_bytes = file.info(release_files)$size
)

manifest <- manifest[
    file.info(release_files)$isdir == FALSE
]

fwrite(
    manifest,
    file.path(RELEASE_DIR, "release_manifest.csv"),
    na = "NA"
)

release_files <- list.files(
    RELEASE_DIR,
    recursive = TRUE,
    full.names = TRUE
)

release_files <- release_files[
    file.info(release_files)$isdir == FALSE
]

write_checksum(
    release_files,
    file.path(RELEASE_DIR, "release_checksums.tsv")
)

# -------------------------------------------------------------------------
# 8. ZIP
# -------------------------------------------------------------------------

if (file.exists(ZIP_FILE)) {
    file.remove(ZIP_FILE)
}

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)

setwd("release")

utils::zip(
    zipfile = basename(ZIP_FILE),
    files = basename(RELEASE_DIR)
)

setwd(old_wd)

if (!file.exists(ZIP_FILE)) {
    stop("Release ZIP was not created.", call. = FALSE)
}

cat(
    "\n============================================================\n",
    "MOSQEDIT-R STEP 17 COMPLETED SUCCESSFULLY\n",
    "Reproducible analysis release frozen\n",
    "============================================================\n",
    "Release directory:                        ", RELEASE_DIR, "\n",
    "Release ZIP:                              ", ZIP_FILE, "\n",
    "Population genomics integrated:           ",
    if (population_available) "YES" else "NO",
    "\n",
    "============================================================\n",
    sep = ""
)

log_step("17", "Reproducible analysis release completed successfully")

