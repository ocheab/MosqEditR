# =============================================================================
# 06_build_feature_matrix_FULL_CORRECTED.R
#
# MosqEdit-R Manuscript 1
# STEP 06 â€” AUTHORITATIVE MASTER FEATURE MATRIX
#
# Builds one auditable gene-level matrix on the COMPLETE Step 01 AgamP4
# universe. Step 05 CRISPR variables are targeted-shortlist-only and remain NA
# outside the shortlist. Ag1000G population genomics is intentionally deferred.
# =============================================================================

source("R/helpers.R")

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop(
        "Package 'data.table' is required. Install with install.packages('data.table').",
        call. = FALSE
    )
}

suppressPackageStartupMessages({
    library(data.table)
})

log_step(
    "06",
    "Building auditable master gene-feature matrix from Steps 01-05"
)

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------

ANNOTATION_FILE  <- "data_processed/01_gene_annotation.csv"
EXPRESSION_FILE  <- "data_processed/02_mozatlas_gene_expression.csv"
ORTHOPHENO_FILE  <- "data_processed/03_gene_orthology_phenotypes.csv"
RANKING_FILE     <- "data_processed/04_preliminary_genomewide_ranking.csv"
EDITABILITY_FILE <- "data_processed/05_crispr_editability.csv"

LANDMARK_FILE <- "metadata/landmark_targets.csv"

MASTER_FILE           <- "data_processed/master_gene_feature_matrix.csv"
FEATURE_MANIFEST_FILE <- "data_processed/06_feature_manifest.csv"
QC_FILE               <- "data_processed/06_master_matrix_qc.csv"
PROVENANCE_FILE       <- "data_processed/06_master_matrix_provenance.csv"
SESSION_INFO_FILE     <- "logs/06_master_matrix_sessionInfo.txt"
CHECKSUM_FILE         <- "logs/06_checksums.tsv"

dir.create("data_processed", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. Utility functions
# -----------------------------------------------------------------------------

assert_file <- function(path) {
    if (!file.exists(path)) {
        stop(
            paste0("Required upstream file is missing:\n", path),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

assert_gene_table <- function(dt, label) {
    if (!"gene_id" %in% names(dt)) {
        stop(
            paste0(label, " does not contain gene_id."),
            call. = FALSE
        )
    }

    if (any(is.na(dt$gene_id) | !nzchar(dt$gene_id))) {
        stop(
            paste0(label, " contains missing or blank gene_id values."),
            call. = FALSE
        )
    }

    dups <- dt[
        duplicated(gene_id),
        unique(gene_id)
    ]

    if (length(dups) > 0L) {
        stop(
            paste0(
                label,
                " contains duplicate gene_id values. Example(s): ",
                paste(head(dups, 10L), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}

# Add ONLY columns that do not already exist in the master table.
# Later-stage files re-carry earlier-stage columns; this prevents .x/.y
# duplication and preserves the earliest authoritative version of each field.
merge_new_columns <- function(master, incoming, source_label) {
    assert_gene_table(incoming, source_label)

    new_columns <- setdiff(
        names(incoming),
        names(master)
    )

    if (length(new_columns) == 0L) {
        message(source_label, ": no new columns to add.")
        return(
            list(
                data = master,
                added = character()
            )
        )
    }

    incoming_small <- incoming[
        ,
        c("gene_id", new_columns),
        with = FALSE
    ]

    n_before <- nrow(master)

    master <- merge(
        master,
        incoming_small,
        by = "gene_id",
        all.x = TRUE,
        sort = FALSE
    )

    if (nrow(master) != n_before) {
        stop(
            paste0(
                "Merge with ",
                source_label,
                " changed master row count from ",
                n_before,
                " to ",
                nrow(master),
                "."
            ),
            call. = FALSE
        )
    }

    message(
        source_label,
        ": added ",
        length(new_columns),
        " new column(s)."
    )

    list(
        data = master,
        added = new_columns
    )
}

# -----------------------------------------------------------------------------
# 3. Validate and load upstream files
# -----------------------------------------------------------------------------

required_files <- c(
    ANNOTATION_FILE,
    EXPRESSION_FILE,
    ORTHOPHENO_FILE,
    RANKING_FILE,
    EDITABILITY_FILE
)

for (path in required_files) {
    assert_file(path)
}

annotation  <- fread(ANNOTATION_FILE)
expression  <- fread(EXPRESSION_FILE)
orthopheno  <- fread(ORTHOPHENO_FILE)
ranking     <- fread(RANKING_FILE)
editability <- fread(EDITABILITY_FILE)

assert_gene_table(annotation,  "Step 01 annotation")
assert_gene_table(expression,  "Step 02 expression")
assert_gene_table(orthopheno,  "Step 03 orthology/phenotype")
assert_gene_table(ranking,     "Step 04 preliminary ranking")
assert_gene_table(editability, "Step 05 CRISPR editability")

# -----------------------------------------------------------------------------
# 4. Freeze Step 01 as the master universe
# -----------------------------------------------------------------------------

master <- copy(annotation)

n_annotation <- nrow(master)
step01_gene_ids <- master$gene_id

if (uniqueN(step01_gene_ids) != n_annotation) {
    stop("Step 01 universe is not gene_id-unique.", call. = FALSE)
}

check_outside_universe <- function(dt, label) {
    outside <- setdiff(
        dt$gene_id,
        step01_gene_ids
    )

    if (length(outside) > 0L) {
        stop(
            paste0(
                label,
                " contains ",
                length(outside),
                " gene(s) outside Step 01. Example(s): ",
                paste(head(outside, 10L), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}

check_outside_universe(expression,  "Step 02")
check_outside_universe(orthopheno,  "Step 03")
check_outside_universe(ranking,     "Step 04")
check_outside_universe(editability, "Step 05")

# -----------------------------------------------------------------------------
# 5. Explicit source-availability flags
# -----------------------------------------------------------------------------

master[
    ,
    has_step02_expression_record :=
        gene_id %in% expression$gene_id
]

master[
    ,
    has_step03_orthopheno_record :=
        gene_id %in% orthopheno$gene_id
]

master[
    ,
    has_step04_ranking_record :=
        gene_id %in% ranking$gene_id
]

master[
    ,
    selected_for_step05_crispr :=
        gene_id %in% editability$gene_id
]

# -----------------------------------------------------------------------------
# 6. Merge Steps 02â€“05, adding only genuinely new fields
# -----------------------------------------------------------------------------

res02 <- merge_new_columns(
    master,
    expression,
    "Step 02 MozAtlas expression"
)
master <- res02$data
step02_added <- res02$added

res03 <- merge_new_columns(
    master,
    orthopheno,
    "Step 03 orthology/FlyBase phenotype"
)
master <- res03$data
step03_added <- res03$added

res04 <- merge_new_columns(
    master,
    ranking,
    "Step 04 preliminary ranking"
)
master <- res04$data
step04_added <- res04$added

res05 <- merge_new_columns(
    master,
    editability,
    "Step 05 targeted CRISPR editability"
)
master <- res05$data
step05_added <- res05$added

# -----------------------------------------------------------------------------
# 7. Benchmark / landmark handling
# -----------------------------------------------------------------------------
# Benchmarks are validation objects, NOT training labels.

landmark_columns_added <- character()

if (file.exists(LANDMARK_FILE)) {
    landmarks <- fread(LANDMARK_FILE)
    assert_gene_table(landmarks, "Landmark target metadata")

    landmark_ids <- intersect(
        landmarks$gene_id,
        step01_gene_ids
    )

    master[
        ,
        is_external_landmark_validation :=
            gene_id %in% landmark_ids
    ]

    landmark_nonkey <- setdiff(
        names(landmarks),
        "gene_id"
    )

    if (length(landmark_nonkey) > 0L) {
        landmark_small <- copy(
            landmarks[
                gene_id %in% step01_gene_ids
            ]
        )

        setnames(
            landmark_small,
            old = landmark_nonkey,
            new = paste0(
                "landmark_meta_",
                landmark_nonkey
            )
        )

        res_landmark <- merge_new_columns(
            master,
            landmark_small,
            "Optional landmark metadata"
        )

        master <- res_landmark$data
        landmark_columns_added <- res_landmark$added
    }

    landmark_metadata_status <-
        "metadata/landmark_targets.csv merged as validation metadata"

} else {
    if ("is_prespecified_benchmark" %in% names(master)) {
        master[
            ,
            is_external_landmark_validation :=
                as.logical(is_prespecified_benchmark)
        ]
    } else if ("benchmark_name" %in% names(master)) {
        master[
            ,
            is_external_landmark_validation :=
                !is.na(benchmark_name)
        ]
    } else {
        master[
            ,
            is_external_landmark_validation :=
                FALSE
        ]
    }

    landmark_metadata_status <-
        "No optional landmark metadata file; Step 04 benchmark annotation used"
}

master[
    ,
    landmark_validation_only :=
        is_external_landmark_validation
]

# -----------------------------------------------------------------------------
# 8. Evidence and stage-maturity flags
# -----------------------------------------------------------------------------
#
# IMPORTANT:
# Step 03 is a FULL-UNIVERSE table, so membership in the Step 03 CSV does NOT
# mean that a gene has orthology or FlyBase phenotype evidence. Distinguish
# structural row presence from actual usable biological evidence.

if ("preliminary_rank" %in% names(master)) {
    master[
        ,
        has_preliminary_rank :=
            !is.na(preliminary_rank)
    ]
} else {
    master[
        ,
        has_preliminary_rank :=
            FALSE
    ]
}

# Step 02 expression evidence.
if ("expression_domain_available" %in% names(master)) {
    master[
        ,
        has_step02_expression_evidence :=
            !is.na(expression_domain_available) &
            as.logical(expression_domain_available)
    ]
} else {
    master[
        ,
        has_step02_expression_evidence :=
            has_step02_expression_record
    ]
}

# Step 03 D. melanogaster orthology evidence.
if ("has_dmel_ortholog" %in% names(master)) {
    master[
        ,
        has_step03_dmel_orthology_evidence :=
            !is.na(has_dmel_ortholog) &
            as.logical(has_dmel_ortholog)
    ]
} else {
    master[
        ,
        has_step03_dmel_orthology_evidence :=
            FALSE
    ]
}

# Step 03 FlyBase phenotype-record evidence.
if ("has_flybase_phenotype_records" %in% names(master)) {
    master[
        ,
        has_step03_flybase_phenotype_evidence :=
            !is.na(has_flybase_phenotype_records) &
            as.logical(has_flybase_phenotype_records)
    ]
} else {
    master[
        ,
        has_step03_flybase_phenotype_evidence :=
            FALSE
    ]
}

# Functional domain actually used in Step 04 ranking.
if ("functional_domain_available" %in% names(master)) {
    master[
        ,
        has_step03_functional_ranking_evidence :=
            !is.na(functional_domain_available) &
            as.logical(functional_domain_available)
    ]
} else {
    master[
        ,
        has_step03_functional_ranking_evidence :=
            has_step03_flybase_phenotype_evidence
    ]
}

# At least one evidence domain capable of contributing to Step 04 ranking.
master[
    ,
    has_genomewide_ranking_evidence :=
        has_step02_expression_evidence |
        has_step03_functional_ranking_evidence
]

if ("reference_sequence_available" %in% names(master)) {
    master[
        ,
        has_step05_reference_sequence :=
            fifelse(
                selected_for_step05_crispr,
                as.logical(reference_sequence_available),
                NA
            )
    ]
} else {
    master[
        ,
        has_step05_reference_sequence :=
            fifelse(
                selected_for_step05_crispr,
                TRUE,
                NA
            )
    ]
}

if (
    "crispr_exon_aware_reference_editability_score" %in%
        names(master)
) {
    master[
        ,
        has_step05_editability_score :=
            selected_for_step05_crispr &
            !is.na(
                crispr_exon_aware_reference_editability_score
            )
    ]
} else {
    master[
        ,
        has_step05_editability_score :=
            FALSE
    ]
}

# Population genomics is intentionally still pending.
master[
    ,
    population_genomics_assessed :=
        FALSE
]

master[
    ,
    population_genomics_scope :=
        "PENDING_TARGETED_VALIDATION"
]

master[
    ,
    analysis_maturity :=
        fcase(
            selected_for_step05_crispr &
                has_step05_editability_score,
            "STEP05_TARGETED_CRISPR_ASSESSED",

            has_preliminary_rank,
            "STEP04_PRELIMINARY_RANKED",

            has_genomewide_ranking_evidence,
            "UPSTREAM_EVIDENCE_ONLY",

            default =
                "ANNOTATION_ONLY"
        )
]

# -----------------------------------------------------------------------------
# 9. Evidence-domain coverage
# -----------------------------------------------------------------------------

master[
    ,
    n_genomewide_evidence_domains :=
        as.integer(has_step02_expression_evidence) +
        as.integer(has_step03_functional_ranking_evidence)
]

master[
    ,
    genomewide_evidence_domain_coverage :=
        n_genomewide_evidence_domains / 2
]

master[
    ,
    targeted_validation_domains_assessed :=
        as.integer(selected_for_step05_crispr) +
        as.integer(population_genomics_assessed)
]

# -----------------------------------------------------------------------------
# 10. Stable ordering
# -----------------------------------------------------------------------------

if ("preliminary_rank" %in% names(master)) {
    master[
        ,
        preliminary_rank_missing :=
            is.na(preliminary_rank)
    ]

    setorder(
        master,
        preliminary_rank_missing,
        preliminary_rank,
        gene_id
    )

    master[
        ,
        preliminary_rank_missing :=
            NULL
    ]
} else {
    setorder(
        master,
        gene_id
    )
}

# -----------------------------------------------------------------------------
# 11. Feature manifest
# -----------------------------------------------------------------------------

annotation_columns <- names(annotation)

source_for_column <- function(col) {
    if (col %in% annotation_columns) {
        return("Step 01 annotation")
    }

    if (col %in% step02_added) {
        return("Step 02 MozAtlas expression")
    }

    if (col %in% step03_added) {
        return("Step 03 orthology/FlyBase phenotype")
    }

    if (col %in% step04_added) {
        return("Step 04 preliminary ranking")
    }

    if (col %in% step05_added) {
        return("Step 05 targeted CRISPR editability")
    }

    if (col %in% landmark_columns_added) {
        return("Optional landmark metadata")
    }

    "Step 06 audit/derived field"
}

scope_for_column <- function(col) {
    src <- source_for_column(col)

    if (src == "Step 05 targeted CRISPR editability") {
        return("targeted_shortlist_only")
    }

    if (src == "Optional landmark metadata") {
        return("validation_metadata_only")
    }

    if (src == "Step 06 audit/derived field") {
        return("derived_audit_field")
    }

    "genomewide_or_partial_genomewide"
}

class_for_column <- function(x) {
    paste(
        class(x),
        collapse = "/"
    )
}

feature_manifest <- data.table(
    column_name =
        names(master),

    source_stage =
        vapply(
            names(master),
            source_for_column,
            character(1)
        ),

    analysis_scope =
        vapply(
            names(master),
            scope_for_column,
            character(1)
        ),

    r_class =
        vapply(
            master,
            class_for_column,
            character(1)
        )
)

feature_manifest[
    ,
    non_missing_n :=
        vapply(
            names(master),
            function(col) {
                sum(
                    !is.na(
                        master[[col]]
                    )
                )
            },
            integer(1)
        )
]

feature_manifest[
    ,
    missing_n :=
        n_annotation - non_missing_n
]

feature_manifest[
    ,
    non_missing_percent :=
        round(
            100 *
                non_missing_n /
                n_annotation,
            3
        )
]

feature_manifest[
    ,
    modelling_note :=
        fcase(
            source_stage ==
                "Step 05 targeted CRISPR editability",
            "NA outside Step 05 shortlist means not assessed, not biological absence.",

            source_stage ==
                "Optional landmark metadata",
            "Validation metadata only; not a training label.",

            column_name %in%
                c(
                    "benchmark_name",
                    "benchmark_class",
                    "benchmark_reference",
                    "is_prespecified_benchmark",
                    "is_external_landmark_validation",
                    "landmark_validation_only"
                ),
            "Benchmark/validation information; exclude from predictor training.",

            column_name %in%
                c(
                    "preliminary_rank",
                    "rank_primary",
                    "rank_expression_dominant",
                    "rank_functional_dominant",
                    "rank_consensus_median"
                ),
            "Derived ranking output; avoid leakage when modelling component evidence.",

            default =
                "Review before modelling."
        )
]

# -----------------------------------------------------------------------------
# 12. QC summary
# -----------------------------------------------------------------------------

n_expression_rows <- master[
    has_step02_expression_record == TRUE,
    .N
]

n_expression_evidence <- master[
    has_step02_expression_evidence == TRUE,
    .N
]

n_step03_rows <- master[
    has_step03_orthopheno_record == TRUE,
    .N
]

n_dmel_orthology <- master[
    has_step03_dmel_orthology_evidence == TRUE,
    .N
]

n_flybase_phenotype <- master[
    has_step03_flybase_phenotype_evidence == TRUE,
    .N
]

n_functional_evidence <- master[
    has_step03_functional_ranking_evidence == TRUE,
    .N
]

n_any_ranking_evidence <- master[
    has_genomewide_ranking_evidence == TRUE,
    .N
]

n_ranked <- master[
    has_preliminary_rank == TRUE,
    .N
]

n_step05_selected <- master[
    selected_for_step05_crispr == TRUE,
    .N
]

n_step05_scored <- master[
    has_step05_editability_score == TRUE,
    .N
]

n_landmarks <- master[
    is_external_landmark_validation == TRUE,
    .N
]

n_annotation_only <- master[
    analysis_maturity == "ANNOTATION_ONLY",
    .N
]

n_upstream_only <- master[
    analysis_maturity == "UPSTREAM_EVIDENCE_ONLY",
    .N
]

n_step04_ranked <- master[
    analysis_maturity == "STEP04_PRELIMINARY_RANKED",
    .N
]

n_step05_assessed <- master[
    analysis_maturity == "STEP05_TARGETED_CRISPR_ASSESSED",
    .N
]

qc <- data.table(
    metric = c(
        "Step 01 master gene universe",
        "Genes represented in Step 02 expression table",
        "Genes with usable Step 02 expression evidence",
        "Genes represented in Step 03 full-universe table",
        "Genes with D. melanogaster orthology evidence",
        "Genes with FlyBase phenotype-record evidence",
        "Genes with Step 03 functional ranking evidence",
        "Genes with at least one Step 04 ranking evidence domain",
        "Genes with preliminary rank",
        "Genes selected for Step 05 CRISPR assessment",
        "Genes with Step 05 exon-aware editability score",
        "External landmark/benchmark validation genes",
        "Annotation-only genes",
        "Upstream-evidence-only genes",
        "Step 04 preliminary-ranked genes not Step 05 assessed",
        "Step 05 targeted CRISPR-assessed genes",
        "Population-genomics assessed genes",
        "Master matrix columns",
        "Feature manifest rows"
    ),

    value = c(
        n_annotation,
        n_expression_rows,
        n_expression_evidence,
        n_step03_rows,
        n_dmel_orthology,
        n_flybase_phenotype,
        n_functional_evidence,
        n_any_ranking_evidence,
        n_ranked,
        n_step05_selected,
        n_step05_scored,
        n_landmarks,
        n_annotation_only,
        n_upstream_only,
        n_step04_ranked,
        n_step05_assessed,
        0,
        ncol(master),
        nrow(feature_manifest)
    )
)

# -----------------------------------------------------------------------------
# 13. Provenance
# -----------------------------------------------------------------------------

provenance <- data.table(
    component = c(
        "Universe",
        "Step 02 integration",
        "Step 03 integration",
        "Step 04 integration",
        "Step 05 integration",
        "Duplicate-column handling",
        "Missing Step 05 values",
        "Benchmarks/landmarks",
        "Population genomics",
        "Master-matrix interpretation",
        "Optional landmark metadata"
    ),

    specification = c(
        paste0(
            "Master matrix is anchored to all ",
            format(n_annotation, big.mark = ","),
            " unique Step 01 AGAP loci."
        ),

        "Step 02 contributes only columns not already present in the Step 01 base table.",

        paste0(
            "Step 03 contributes only new orthology/phenotype fields. The Step 03 ",
            "table contains all Step 01 genes, so table membership is structural ",
            "and is NOT interpreted as evidence. Actual orthology and FlyBase ",
            "phenotype evidence are tracked with dedicated logical flags."
        ),

        paste0(
            "Step 04 contributes preliminary scores/ranks and benchmark annotation ",
            "without replacing upstream Step 01-03 fields."
        ),

        paste0(
            "Step 05 contributes exon-aware CRISPR reference-editability variables ",
            "only for the targeted shortlist."
        ),

        paste0(
            "Later-stage copies of earlier columns are not re-merged, preventing ",
            "ambiguous .x/.y duplicates."
        ),

        paste0(
            "CRISPR values outside the Step 05 shortlist remain NA and ",
            "selected_for_step05_crispr=FALSE. NA therefore means not assessed."
        ),

        paste0(
            "Benchmarks/landmarks are retained as external validation objects. ",
            "Step 06 creates no supervised training label."
        ),

        paste0(
            "Ag1000G population genomics is intentionally not required or merged ",
            "at Step 06; population_genomics_assessed remains FALSE."
        ),

        paste0(
            "The matrix contains both genome-wide discovery evidence and targeted ",
            "validation evidence. Use the feature manifest before modelling."
        ),

        landmark_metadata_status
    )
)

# -----------------------------------------------------------------------------
# 14. Final integrity checks
# -----------------------------------------------------------------------------

if (nrow(master) != n_annotation) {
    stop(
        paste0(
            "Final master row count is ",
            nrow(master),
            " but Step 01 contains ",
            n_annotation,
            "."
        ),
        call. = FALSE
    )
}

if (uniqueN(master$gene_id) != nrow(master)) {
    stop(
        "Final master matrix contains duplicate gene_id values.",
        call. = FALSE
    )
}

if (any(!master$gene_id %in% step01_gene_ids)) {
    stop(
        "Final master matrix contains genes outside Step 01.",
        call. = FALSE
    )
}

if (n_step05_selected != nrow(editability)) {
    stop(
        paste0(
            "Step 05 selection flag count (",
            n_step05_selected,
            ") does not equal Step 05 rows (",
            nrow(editability),
            ")."
        ),
        call. = FALSE
    )
}

# Step 05 values must never appear on genes outside the Step 05 shortlist.
step05_metric_candidates <- intersect(
    c(
        "crispr_exon_aware_reference_editability_score",
        "crispr_exon_aware_reference_editability_rank",
        "n_unique_functional_region_reference_quality_sites"
    ),
    names(master)
)

if (length(step05_metric_candidates) > 0L) {
    leak_check <- master[
        selected_for_step05_crispr == FALSE,
        rowSums(
            !is.na(.SD)
        ),
        .SDcols = step05_metric_candidates
    ]

    if (any(leak_check > 0L)) {
        stop(
            paste0(
                "Step 05 values were found on genes outside the Step 05 shortlist. ",
                "This indicates a merge-scope error."
            ),
            call. = FALSE
        )
    }
}

# No duplicate output column names.
if (anyDuplicated(names(master)) > 0L) {
    stop(
        "Duplicate column names detected in master matrix.",
        call. = FALSE
    )
}

# -----------------------------------------------------------------------------
# 15. Write outputs
# -----------------------------------------------------------------------------

fwrite(
    master,
    MASTER_FILE,
    na = "NA"
)

fwrite(
    feature_manifest,
    FEATURE_MANIFEST_FILE,
    na = "NA"
)

fwrite(
    qc,
    QC_FILE,
    na = "NA"
)

fwrite(
    provenance,
    PROVENANCE_FILE,
    na = "NA"
)

capture.output(
    sessionInfo(),
    file = SESSION_INFO_FILE
)

checksum_files <- c(
    ANNOTATION_FILE,
    EXPRESSION_FILE,
    ORTHOPHENO_FILE,
    RANKING_FILE,
    EDITABILITY_FILE,
    MASTER_FILE,
    FEATURE_MANIFEST_FILE,
    QC_FILE,
    PROVENANCE_FILE,
    SESSION_INFO_FILE
)

if (file.exists(LANDMARK_FILE)) {
    checksum_files <- c(
        checksum_files,
        LANDMARK_FILE
    )
}

checksum_files <- unique(
    checksum_files[
        file.exists(checksum_files)
    ]
)

write_checksum(
    checksum_files,
    CHECKSUM_FILE
)

# -----------------------------------------------------------------------------
# 16. Console summary
# -----------------------------------------------------------------------------

cat(
    "\n",
    "============================================================\n",
    "MOSQEDIT-R STEP 06 COMPLETED SUCCESSFULLY\n",
    "Auditable master gene-feature matrix\n",
    "============================================================\n",
    "Step 01 master universe:                  ",
    format(n_annotation, big.mark = ","),
    "\n",
    "Genes with Step 02 expression evidence:   ",
    format(n_expression_evidence, big.mark = ","),
    "\n",
    "Genes with Dmel orthology evidence:        ",
    format(n_dmel_orthology, big.mark = ","),
    "\n",
    "Genes with FlyBase phenotype evidence:     ",
    format(n_flybase_phenotype, big.mark = ","),
    "\n",
    "Genes with functional ranking evidence:    ",
    format(n_functional_evidence, big.mark = ","),
    "\n",
    "Genes with any Step 04 ranking evidence:   ",
    format(n_any_ranking_evidence, big.mark = ","),
    "\n",
    "Genes with preliminary rank:              ",
    format(n_ranked, big.mark = ","),
    "\n",
    "Genes selected for Step 05 CRISPR:        ",
    format(n_step05_selected, big.mark = ","),
    "\n",
    "Genes with Step 05 editability score:     ",
    format(n_step05_scored, big.mark = ","),
    "\n",
    "External landmark/benchmark genes:        ",
    format(n_landmarks, big.mark = ","),
    "\n",
    "Population genomics merged:               NO\n",
    "Master matrix rows:                       ",
    format(nrow(master), big.mark = ","),
    "\n",
    "Master matrix columns:                    ",
    format(ncol(master), big.mark = ","),
    "\n",
    "Duplicate gene IDs:                      ",
    nrow(master) - uniqueN(master$gene_id),
    "\n",
    "Python used:                              NO\n",
    "============================================================\n",
    sep = ""
)

cat(
    "\nQC summary:\n"
)

print(qc)

cat(
    "\nAnalysis maturity distribution:\n"
)

print(
    master[
        ,
        .N,
        by = analysis_maturity
    ][
        order(analysis_maturity)
    ]
)

cat(
    "\nTop 20 preliminary-ranked genes in master matrix:\n"
)

display_candidates <- c(
    "preliminary_rank",
    "gene_id",
    "benchmark_name",
    "fris",
    "functional_evidence_score",
    "score_primary",
    "selected_for_step05_crispr",
    "crispr_exon_aware_reference_editability_score",
    "crispr_exon_aware_reference_editability_rank",
    "population_genomics_assessed",
    "analysis_maturity"
)

display_columns <- display_candidates[
    display_candidates %in%
        names(master)
]

if ("preliminary_rank" %in% names(master)) {
    print(
        master[
            !is.na(preliminary_rank)
        ][
            1:min(20L, .N),
            ..display_columns
        ]
    )
}

cat(
    "\nMaster feature matrix written to:\n",
    MASTER_FILE,
    "\n",
    sep = ""
)

cat(
    "\nFeature manifest written to:\n",
    FEATURE_MANIFEST_FILE,
    "\n",
    sep = ""
)

log_step(
    "06",
    "Auditable master gene-feature matrix completed successfully"
)

