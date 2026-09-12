# =============================================================================
# 07_build_leakage_safe_modelling_sets.R
#
# MosqEdit-R Manuscript 1
#
# STEP 07
# Leakage-safe modelling datasets and benchmark holdout construction
#
# PURE R
#
# PURPOSE
# -------
# Convert the Step 06 master matrix into clearly separated modelling/validation
# datasets without leaking benchmark status, FlyBase-derived outcome variables,
# Step 04 ranking outputs, or targeted-only Step 05 CRISPR information into the
# genome-wide discovery model.
#
# IMPORTANT CONCEPT
# -----------------
# The supervised signal available at this stage is POSITIVEâ€“UNLABELED (PU):
#
#   positive = strict mapped FlyBase female-sterility evidence
#   unlabeled = every other eligible gene
#
# "Unlabeled" MUST NOT be interpreted as a confirmed negative.
#
# Because FlyBase female-sterility evidence defines the PU outcome, ALL
# FlyBase phenotype-derived variables are excluded from PU predictors.
#
# Prespecified benchmark/landmark genes are held out completely from model
# fitting and are reserved for external recovery validation.
#
# Step 04 ranking outputs are also excluded from PU predictors because they are
# derived from upstream evidence and would cause circular evaluation.
#
# Step 05 CRISPR variables are restricted to 205 shortlisted genes and are NOT
# valid genome-wide predictors. They form a separate targeted-validation table.
#
# POPULATION GENOMICS
# -------------------
# Ag1000G remains pending and is not imputed as zero.
#
# INPUTS
# ------
# data_processed/master_gene_feature_matrix.csv
# data_processed/06_feature_manifest.csv
#
# OUTPUTS
# -------
# data_processed/07_pu_M1_expression.csv
# data_processed/07_pu_M2_expression_orthology.csv
# data_processed/07_pu_M3_expression_orthology_annotation.csv
# data_processed/07_targeted_crispr_validation_panel.csv
# data_processed/07_external_benchmark_holdout.csv
# data_processed/07_predictor_manifest.csv
# data_processed/07_modelling_qc.csv
# data_processed/07_modelling_provenance.csv
# logs/07_modelling_sessionInfo.txt
# logs/07_checksums.tsv
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Helpers and package checks
# -----------------------------------------------------------------------------

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
    "07",
    "Building leakage-safe PU modelling and targeted-validation datasets"
)


# -----------------------------------------------------------------------------
# 2. Configuration
# -----------------------------------------------------------------------------

MASTER_FILE <-
    "data_processed/master_gene_feature_matrix.csv"


STEP06_MANIFEST_FILE <-
    "data_processed/06_feature_manifest.csv"


M1_FILE <-
    "data_processed/07_pu_M1_expression.csv"


M2_FILE <-
    "data_processed/07_pu_M2_expression_orthology.csv"


M3_FILE <-
    "data_processed/07_pu_M3_expression_orthology_annotation.csv"


TARGETED_CRISPR_FILE <-
    "data_processed/07_targeted_crispr_validation_panel.csv"


BENCHMARK_HOLDOUT_FILE <-
    "data_processed/07_external_benchmark_holdout.csv"


PREDICTOR_MANIFEST_FILE <-
    "data_processed/07_predictor_manifest.csv"


QC_FILE <-
    "data_processed/07_modelling_qc.csv"


PROVENANCE_FILE <-
    "data_processed/07_modelling_provenance.csv"


SESSION_INFO_FILE <-
    "logs/07_modelling_sessionInfo.txt"


CHECKSUM_FILE <-
    "logs/07_checksums.tsv"


dir.create(
    "data_processed",
    recursive = TRUE,
    showWarnings = FALSE
)


dir.create(
    "logs",
    recursive = TRUE,
    showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# 3. Utility functions
# -----------------------------------------------------------------------------

assert_file <- function(path) {

    if (!file.exists(path)) {

        stop(
            paste0(
                "Required input file is missing:\n",
                path
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


assert_unique_gene_ids <- function(
        dt,
        label
) {

    if (!"gene_id" %in% names(dt)) {

        stop(
            paste0(
                label,
                " lacks gene_id."
            ),
            call. = FALSE
        )
    }


    if (
        any(
            is.na(dt$gene_id) |
                !nzchar(dt$gene_id)
        )
    ) {

        stop(
            paste0(
                label,
                " contains missing/blank gene_id."
            ),
            call. = FALSE
        )
    }


    if (
        uniqueN(dt$gene_id) !=
            nrow(dt)
    ) {

        stop(
            paste0(
                label,
                " contains duplicate gene_id values."
            ),
            call. = FALSE
        )
    }


    invisible(TRUE)
}


existing_columns <- function(
        preferred,
        dt
) {

    preferred[
        preferred %in%
            names(dt)
    ]
}


all_missing_rows <- function(
        dt,
        cols
) {

    if (length(cols) == 0L) {

        return(
            rep(
                TRUE,
                nrow(dt)
            )
        )
    }


    rowSums(
        !is.na(
            dt[
                ,
                ..cols
            ]
        )
    ) ==
        0L
}


# Add explicit missingness indicators without imputing values.
add_missingness_flags <- function(
        dt,
        predictor_cols
) {

    out <- copy(
        dt
    )


    for (
        col in predictor_cols
    ) {

        missing_col <- paste0(
            "missing__",
            col
        )


        out[
            ,
            (missing_col) :=
                is.na(
                    get(
                        col
                    )
                )
        ]
    }


    out
}


# -----------------------------------------------------------------------------
# 4. Load Step 06 master matrix
# -----------------------------------------------------------------------------

assert_file(
    MASTER_FILE
)


assert_file(
    STEP06_MANIFEST_FILE
)


master <- fread(
    MASTER_FILE
)


step06_manifest <- fread(
    STEP06_MANIFEST_FILE
)


assert_unique_gene_ids(
    master,
    "Step 06 master matrix"
)


n_master <- nrow(
    master
)


# -----------------------------------------------------------------------------
# 5. Define external benchmark/landmark holdout
# -----------------------------------------------------------------------------

if (
    "is_external_landmark_validation" %in%
    names(master)
) {

    master[
        ,
        external_validation_holdout :=
            !is.na(
                is_external_landmark_validation
            ) &
            as.logical(
                is_external_landmark_validation
            )
    ]

} else if (
    "is_prespecified_benchmark" %in%
    names(master)
) {

    master[
        ,
        external_validation_holdout :=
            !is.na(
                is_prespecified_benchmark
            ) &
            as.logical(
                is_prespecified_benchmark
            )
    ]

} else if (
    "benchmark_name" %in%
    names(master)
) {

    master[
        ,
        external_validation_holdout :=
            !is.na(
                benchmark_name
            ) &
            nzchar(
                benchmark_name
            )
    ]

} else {

    master[
        ,
        external_validation_holdout :=
            FALSE
    ]
}


n_external_holdout <- master[
    external_validation_holdout ==
        TRUE,
    .N
]


if (n_external_holdout == 0L) {

    warning(
        paste0(
            "No external benchmark/landmark genes were identified. ",
            "The expected workflow contains six prespecified validation genes."
        ),
        call. = FALSE
    )
}


# -----------------------------------------------------------------------------
# 6. Construct PU outcome
# -----------------------------------------------------------------------------
#
# Preferred source: strict Step 04 score_strict_female_sterility.
# Fallback: Step 03 flybase_female_sterile.
#
# Label semantics:
#   1 = observed positive
#   0 = unlabeled, NOT confirmed negative.

if (
    "score_strict_female_sterility" %in%
    names(master)
) {

    master[
        ,
        pu_positive :=
            !is.na(
                score_strict_female_sterility
            ) &
            as.numeric(
                score_strict_female_sterility
            ) ==
                1
    ]


    pu_label_source <-
        "score_strict_female_sterility"

} else if (
    "flybase_female_sterile" %in%
    names(master)
) {

    master[
        ,
        pu_positive :=
            !is.na(
                flybase_female_sterile
            ) &
            as.logical(
                flybase_female_sterile
            )
    ]


    pu_label_source <-
        "flybase_female_sterile"

} else {

    stop(
        paste0(
            "Cannot construct PU outcome: neither score_strict_female_sterility ",
            "nor flybase_female_sterile is present."
        ),
        call. = FALSE
    )
}


master[
    ,
    pu_label :=
        as.integer(
            pu_positive
        )
]


master[
    ,
    pu_label_semantics :=
        data.table::fifelse(
            pu_positive,
            "observed_positive",
            "unlabeled_not_negative"
        )
]


# -----------------------------------------------------------------------------
# 7. Explicit leakage exclusions
# -----------------------------------------------------------------------------

phenotype_leakage_patterns <- c(
    "flybase",
    "female_steril",
    "reproductive_phenotype",
    "functional_evidence",
    "score_strict",
    "score_broader"
)


ranking_leakage_patterns <- c(
    "^preliminary_rank$",
    "^rank_",
    "^score_primary",
    "^score_expression",
    "^score_functional",
    "^rank_consensus",
    "^rank_sensitivity",
    "^selection_reason$"
)


benchmark_leakage_patterns <- c(
    "benchmark",
    "landmark",
    "external_validation"
)


targeted_leakage_patterns <- c(
    "crispr",
    "step05",
    "selected_for_step05"
)


all_master_columns <- names(
    master
)


matches_any_pattern <- function(
        x,
        patterns
) {

    Reduce(
        `|`,
        lapply(
            patterns,
            function(p) {

                grepl(
                    p,
                    x,
                    ignore.case = TRUE,
                    perl = TRUE
                )
            }
        )
    )
}


phenotype_leakage_columns <- all_master_columns[
    matches_any_pattern(
        all_master_columns,
        phenotype_leakage_patterns
    )
]


ranking_leakage_columns <- all_master_columns[
    matches_any_pattern(
        all_master_columns,
        ranking_leakage_patterns
    )
]


benchmark_leakage_columns <- all_master_columns[
    matches_any_pattern(
        all_master_columns,
        benchmark_leakage_patterns
    )
]


targeted_leakage_columns <- all_master_columns[
    matches_any_pattern(
        all_master_columns,
        targeted_leakage_patterns
    )
]


# -----------------------------------------------------------------------------
# 8. Curated leakage-safe predictor families
# -----------------------------------------------------------------------------

expression_predictor_candidates <- c(
    "expr_ovary",
    "expr_testis",
    "expr_female_body",
    "expr_male_body",
    "expr_female_midgut",
    "expr_female_salivary",
    "tau_tissue_specificity",
    "female_male_delta",
    "ovary_specificity_delta"
)


expression_predictors <- existing_columns(
    expression_predictor_candidates,
    master
)


if (length(expression_predictors) == 0L) {

    stop(
        "No expected Step 02 expression predictors were found.",
        call. = FALSE
    )
}


orthology_predictor_candidates <- c(
    "has_dmel_ortholog",
    "dmel_ortholog_count"
)


orthology_predictors <- existing_columns(
    orthology_predictor_candidates,
    master
)


# Annotation-derived features.
if (
    all(
        c(
            "gene_start",
            "gene_end"
        ) %in%
        names(master)
    )
) {

    master[
        ,
        gene_length_bp :=
            as.numeric(
                gene_end
            ) -
            as.numeric(
                gene_start
            ) +
            1
    ]

} else {

    master[
        ,
        gene_length_bp :=
            NA_real_
    ]
}


annotation_predictor_candidates <- c(
    "gene_length_bp",
    "strand",
    "gene_biotype",
    "biotype"
)


annotation_predictors <- existing_columns(
    annotation_predictor_candidates,
    master
)


# Prefer one biotype field if both aliases exist.
if (
    all(
        c(
            "gene_biotype",
            "biotype"
        ) %in%
        annotation_predictors
    )
) {

    annotation_predictors <- setdiff(
        annotation_predictors,
        "biotype"
    )
}


# Safety check: curated predictor families must not contain leakage fields.
all_safe_candidates <- unique(
    c(
        expression_predictors,
        orthology_predictors,
        annotation_predictors
    )
)


forbidden_candidates <- unique(
    c(
        phenotype_leakage_columns,
        ranking_leakage_columns,
        benchmark_leakage_columns,
        targeted_leakage_columns
    )
)


safe_predictors <- setdiff(
    all_safe_candidates,
    forbidden_candidates
)


if (
    !setequal(
        safe_predictors,
        all_safe_candidates
    )
) {

    removed <- setdiff(
        all_safe_candidates,
        safe_predictors
    )


    warning(
        paste0(
            "The following curated predictors were removed by leakage rules: ",
            paste(
                removed,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}


expression_predictors <- intersect(
    expression_predictors,
    safe_predictors
)


orthology_predictors <- intersect(
    orthology_predictors,
    safe_predictors
)


annotation_predictors <- intersect(
    annotation_predictors,
    safe_predictors
)


# -----------------------------------------------------------------------------
# 9. Shared modelling metadata
# -----------------------------------------------------------------------------

shared_columns <- existing_columns(
    c(
        "gene_id",
        "gene_name",
        "chromosome",
        "gene_start",
        "gene_end",
        "pu_label",
        "pu_label_semantics",
        "external_validation_holdout",
        "has_step02_expression_evidence",
        "has_step03_dmel_orthology_evidence",
        "has_genomewide_ranking_evidence"
    ),
    master
)


# -----------------------------------------------------------------------------
# 10. M1 â€” expression-only PU matrix
# -----------------------------------------------------------------------------

m1_predictors <- unique(
    expression_predictors
)


m1_rows <-
    !master$external_validation_holdout &
    (
        if (
            "has_step02_expression_evidence" %in%
            names(master)
        ) {

            master$has_step02_expression_evidence

        } else {

            !all_missing_rows(
                master,
                m1_predictors
            )
        }
    )


M1 <- master[
    m1_rows,
    c(
        shared_columns,
        m1_predictors
    ),
    with = FALSE
]


M1 <- add_missingness_flags(
    M1,
    m1_predictors
)


M1[
    ,
    modelling_set :=
        "M1_expression"
]


M1[
    ,
    training_role :=
        "PU_training_pool"
]


# -----------------------------------------------------------------------------
# 11. M2 â€” expression + orthology structure
# -----------------------------------------------------------------------------

m2_predictors <- unique(
    c(
        expression_predictors,
        orthology_predictors
    )
)


m2_available <-
    !all_missing_rows(
        master,
        m2_predictors
    )


M2 <- master[
    !external_validation_holdout &
        m2_available,
    c(
        shared_columns,
        m2_predictors
    ),
    with = FALSE
]


M2 <- add_missingness_flags(
    M2,
    m2_predictors
)


M2[
    ,
    modelling_set :=
        "M2_expression_orthology"
]


M2[
    ,
    training_role :=
        "PU_training_pool"
]


# -----------------------------------------------------------------------------
# 12. M3 â€” expression + orthology + neutral annotation
# -----------------------------------------------------------------------------

m3_predictors <- unique(
    c(
        expression_predictors,
        orthology_predictors,
        annotation_predictors
    )
)


m3_biological_predictors <- unique(
    c(
        expression_predictors,
        orthology_predictors
    )
)


# Requiring some expression/orthology evidence prevents annotation-only genes
# from entering merely because coordinates/length are universally available.
m3_available <-
    !all_missing_rows(
        master,
        m3_biological_predictors
    )


M3 <- master[
    !external_validation_holdout &
        m3_available,
    c(
        shared_columns,
        m3_predictors
    ),
    with = FALSE
]


M3 <- add_missingness_flags(
    M3,
    m3_predictors
)


M3[
    ,
    modelling_set :=
        "M3_expression_orthology_annotation"
]


M3[
    ,
    training_role :=
        "PU_training_pool"
]


# -----------------------------------------------------------------------------
# 13. External benchmark holdout table
# -----------------------------------------------------------------------------

benchmark_display_candidates <- c(
    "gene_id",
    "gene_name",
    "benchmark_name",
    "benchmark_class",
    "benchmark_reference",
    "preliminary_rank",
    "fris",
    "functional_evidence_score",
    "score_primary",
    "crispr_exon_aware_reference_editability_score",
    "crispr_exon_aware_reference_editability_rank",
    "pu_label",
    "pu_label_semantics",
    "population_genomics_assessed"
)


benchmark_columns <- existing_columns(
    benchmark_display_candidates,
    master
)


benchmark_holdout <- master[
    external_validation_holdout ==
        TRUE,
    ..benchmark_columns
]


benchmark_holdout[
    ,
    evaluation_role :=
        "EXTERNAL_VALIDATION_ONLY"
]


# -----------------------------------------------------------------------------
# 14. Targeted Step 05 CRISPR validation panel
# -----------------------------------------------------------------------------
#
# This table may contain Step 04 ranking outputs because it is for validation
# and re-ranking, NOT genome-wide PU model fitting.

targeted_columns_candidates <- c(
    "gene_id",
    "gene_name",
    "chromosome",
    "gene_start",
    "gene_end",

    "benchmark_name",
    "benchmark_class",
    "benchmark_reference",
    "is_external_landmark_validation",

    "preliminary_rank",
    "rank_primary",
    "rank_expression_dominant",
    "rank_functional_dominant",
    "rank_consensus_median",
    "rank_sensitivity_range",

    "fris",
    "functional_evidence_score",
    "score_primary",
    "evidence_domain_coverage",

    "gene_has_cds",
    "n_transcripts_scanned",
    "n_unique_exon_contiguous_ngg_sites",
    "n_unique_exon_contiguous_reference_quality_sites",
    "n_unique_functional_region_reference_quality_sites",
    "n_unique_CDS_reference_quality_sites",
    "n_unique_UTR_reference_quality_sites",
    "n_unique_noncoding_exon_reference_quality_sites",
    "n_unique_functional_reference_guide_sequences",
    "fraction_raw_guide_instances_exon_contiguous",
    "fraction_raw_guide_instances_junction_or_unmapped",
    "fraction_transcripts_with_functional_region_reference_quality_guide",
    "best_functional_guide_gc_percent",
    "median_functional_guide_gc_percent",
    "crispr_exon_aware_reference_editability_score",
    "crispr_exon_aware_reference_editability_rank",
    "step05_status",

    "population_conservation_status",
    "population_genomics_assessed",
    "genomewide_offtarget_status",

    "pu_label",
    "pu_label_semantics"
)


targeted_columns <- existing_columns(
    targeted_columns_candidates,
    master
)


targeted_crispr <- master[
    selected_for_step05_crispr ==
        TRUE,
    ..targeted_columns
]


targeted_crispr[
    ,
    analysis_role :=
        "TARGETED_CRISPR_VALIDATION_AND_RERANKING"
]


if (
    "preliminary_rank" %in%
    names(targeted_crispr)
) {

    targeted_crispr[
        ,
        preliminary_rank_missing :=
            is.na(
                preliminary_rank
            )
    ]


    data.table::setorder(
        targeted_crispr,
        preliminary_rank_missing,
        preliminary_rank,
        gene_id
    )


    targeted_crispr[
        ,
        preliminary_rank_missing :=
            NULL
    ]
}


# -----------------------------------------------------------------------------
# 15. Predictor manifest
# -----------------------------------------------------------------------------

manifest_rows <- list()


add_manifest_rows <- function(
        model_set,
        predictors,
        family
) {

    if (length(predictors) == 0L) {
        return(NULL)
    }


    data.table(
        modelling_set =
            model_set,

        predictor =
            predictors,

        predictor_family =
            family,

        allowed_for_pu_training =
            TRUE,

        leakage_status =
            "SAFE",

        missingness_handling =
            "Preserve NA + explicit missingness flag; imputation deferred to modelling step"
    )
}


manifest_rows[[1]] <- add_manifest_rows(
    "M1_expression",
    expression_predictors,
    "expression"
)


manifest_rows[[2]] <- add_manifest_rows(
    "M2_expression_orthology",
    expression_predictors,
    "expression"
)


manifest_rows[[3]] <- add_manifest_rows(
    "M2_expression_orthology",
    orthology_predictors,
    "orthology_structure"
)


manifest_rows[[4]] <- add_manifest_rows(
    "M3_expression_orthology_annotation",
    expression_predictors,
    "expression"
)


manifest_rows[[5]] <- add_manifest_rows(
    "M3_expression_orthology_annotation",
    orthology_predictors,
    "orthology_structure"
)


manifest_rows[[6]] <- add_manifest_rows(
    "M3_expression_orthology_annotation",
    annotation_predictors,
    "annotation"
)


predictor_manifest <- rbindlist(
    manifest_rows,
    fill = TRUE,
    use.names = TRUE
)


excluded_manifest <- data.table(
    modelling_set =
        "ALL_PU_MODELS",

    predictor =
        unique(
            c(
                phenotype_leakage_columns,
                ranking_leakage_columns,
                benchmark_leakage_columns,
                targeted_leakage_columns
            )
        )
)


excluded_manifest[
    ,
    predictor_family :=
        fcase(
            predictor %in%
                phenotype_leakage_columns,
            "phenotype_label_or_derivative",

            predictor %in%
                ranking_leakage_columns,
            "derived_ranking_output",

            predictor %in%
                benchmark_leakage_columns,
            "benchmark_validation_metadata",

            predictor %in%
                targeted_leakage_columns,
            "targeted_only_crispr",

            default =
                "other"
        )
]


excluded_manifest[
    ,
    `:=`(
        allowed_for_pu_training =
            FALSE,

        leakage_status =
            "EXCLUDED",

        missingness_handling =
            "Not applicable"
    )
]


predictor_manifest <- rbindlist(
    list(
        predictor_manifest,
        excluded_manifest
    ),
    fill = TRUE,
    use.names = TRUE
)


predictor_manifest <- unique(
    predictor_manifest,
    by = c(
        "modelling_set",
        "predictor"
    )
)


data.table::setorder(
    predictor_manifest,
    modelling_set,
    -allowed_for_pu_training,
    predictor_family,
    predictor
)


# -----------------------------------------------------------------------------
# 16. Integrity checks â€” no benchmark leakage
# -----------------------------------------------------------------------------

if (
    M1[
        external_validation_holdout ==
            TRUE,
        .N
    ] >
        0L ||
    M2[
        external_validation_holdout ==
            TRUE,
        .N
    ] >
        0L ||
    M3[
        external_validation_holdout ==
            TRUE,
        .N
    ] >
        0L
) {

    stop(
        "Benchmark leakage detected: external validation genes entered a PU training pool.",
        call. = FALSE
    )
}


# No forbidden predictor may be in the actual model predictor sets.
used_predictors <- unique(
    c(
        m1_predictors,
        m2_predictors,
        m3_predictors
    )
)


forbidden_used <- intersect(
    used_predictors,
    forbidden_candidates
)


if (length(forbidden_used) > 0L) {

    stop(
        paste0(
            "Leakage failure: forbidden predictors entered PU sets: ",
            paste(
                forbidden_used,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}


# Positive counts should be nonzero.
n_positive_master <- master[
    pu_positive ==
        TRUE,
    .N
]


if (n_positive_master == 0L) {

    stop(
        "PU outcome contains zero observed positives.",
        call. = FALSE
    )
}


# PU labels must be only 0/1.
for (
    z in list(
        M1,
        M2,
        M3
    )
) {

    if (
        !all(
            unique(
                z$pu_label
            ) %in%
                c(
                    0L,
                    1L
                )
        )
    ) {

        stop(
            "PU label integrity failure: values outside {0,1}.",
            call. = FALSE
        )
    }
}


# -----------------------------------------------------------------------------
# 17. QC
# -----------------------------------------------------------------------------

n_positive_holdout <- master[
    external_validation_holdout ==
        TRUE &
        pu_positive ==
        TRUE,
    .N
]


n_positive_available_training <- master[
    external_validation_holdout ==
        FALSE &
        pu_positive ==
        TRUE,
    .N
]


qc <- data.table(
    metric = c(
        "Step 06 master genes",
        "Observed PU-positive genes in master",
        "External benchmark/landmark holdout genes",
        "Observed PU positives among external holdouts",
        "Observed PU positives available outside holdout",
        "M1 expression training-pool genes",
        "M1 observed positives",
        "M1 unlabeled genes",
        "M1 predictors before missingness flags",
        "M2 expression+orthology training-pool genes",
        "M2 observed positives",
        "M2 unlabeled genes",
        "M2 predictors before missingness flags",
        "M3 expression+orthology+annotation training-pool genes",
        "M3 observed positives",
        "M3 unlabeled genes",
        "M3 predictors before missingness flags",
        "Targeted CRISPR validation genes",
        "Population-genomics predictor used",
        "Benchmark metadata used as predictor",
        "FlyBase phenotype variable used as predictor",
        "Step 04 rank/score used as PU predictor",
        "Step 05 CRISPR used as genome-wide PU predictor"
    ),

    value = c(
        n_master,
        n_positive_master,
        n_external_holdout,
        n_positive_holdout,
        n_positive_available_training,

        nrow(
            M1
        ),
        M1[
            pu_label ==
                1L,
            .N
        ],
        M1[
            pu_label ==
                0L,
            .N
        ],
        length(
            m1_predictors
        ),

        nrow(
            M2
        ),
        M2[
            pu_label ==
                1L,
            .N
        ],
        M2[
            pu_label ==
                0L,
            .N
        ],
        length(
            m2_predictors
        ),

        nrow(
            M3
        ),
        M3[
            pu_label ==
                1L,
            .N
        ],
        M3[
            pu_label ==
                0L,
            .N
        ],
        length(
            m3_predictors
        ),

        nrow(
            targeted_crispr
        ),
        0,
        0,
        0,
        0,
        0
    )
)


# -----------------------------------------------------------------------------
# 18. Provenance
# -----------------------------------------------------------------------------

provenance <- data.table(
    component = c(
        "PU outcome",
        "Unlabeled interpretation",
        "External holdout",
        "Phenotype leakage control",
        "Ranking leakage control",
        "CRISPR scope control",
        "M1",
        "M2",
        "M3",
        "Missing predictors",
        "Population genomics",
        "Targeted validation panel"
    ),

    specification = c(
        paste0(
            "Observed-positive label derived from ",
            pu_label_source,
            "."
        ),

        paste0(
            "pu_label=0 denotes unlabeled, not experimentally or genetically ",
            "confirmed negative."
        ),

        paste0(
            "All prespecified benchmark/landmark genes are removed from PU ",
            "training pools and written separately for external validation."
        ),

        paste0(
            "All FlyBase phenotype/female-sterility/reproductive-phenotype ",
            "variables and functional_evidence derivatives are excluded from PU ",
            "predictors because they define or directly encode the outcome."
        ),

        paste0(
            "Step 04 preliminary ranks and integrated scores are excluded from ",
            "PU predictors to avoid circular prediction/evaluation."
        ),

        paste0(
            "Step 05 CRISPR variables are excluded from genome-wide PU predictors ",
            "because they were measured only for the targeted shortlist."
        ),

        paste0(
            "M1 uses curated MozAtlas expression/tissue-specificity predictors."
        ),

        paste0(
            "M2 adds D. melanogaster orthology structure but excludes FlyBase ",
            "phenotype variables."
        ),

        paste0(
            "M3 adds neutral annotation variables such as gene length/strand/",
            "biotype where available."
        ),

        paste0(
            "No biological predictor is zero-imputed in Step 07. Missing values ",
            "remain NA and receive explicit missingness flags. Model-specific ",
            "imputation is deferred to the fitting step."
        ),

        paste0(
            "Ag1000G remains pending and is not used as a predictor."
        ),

        paste0(
            "The 205-gene CRISPR table is a targeted validation/re-ranking panel, ",
            "not a genome-wide supervised training matrix."
        )
    )
)


# -----------------------------------------------------------------------------
# 19. Write outputs
# -----------------------------------------------------------------------------

fwrite(
    M1,
    M1_FILE,
    na = "NA"
)


fwrite(
    M2,
    M2_FILE,
    na = "NA"
)


fwrite(
    M3,
    M3_FILE,
    na = "NA"
)


fwrite(
    targeted_crispr,
    TARGETED_CRISPR_FILE,
    na = "NA"
)


fwrite(
    benchmark_holdout,
    BENCHMARK_HOLDOUT_FILE,
    na = "NA"
)


fwrite(
    predictor_manifest,
    PREDICTOR_MANIFEST_FILE,
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
    MASTER_FILE,
    STEP06_MANIFEST_FILE,
    M1_FILE,
    M2_FILE,
    M3_FILE,
    TARGETED_CRISPR_FILE,
    BENCHMARK_HOLDOUT_FILE,
    PREDICTOR_MANIFEST_FILE,
    QC_FILE,
    PROVENANCE_FILE,
    SESSION_INFO_FILE
)


checksum_files <- checksum_files[
    file.exists(
        checksum_files
    )
]


write_checksum(
    checksum_files,
    CHECKSUM_FILE
)


# -----------------------------------------------------------------------------
# 20. Console summary
# -----------------------------------------------------------------------------

cat(
    "\n",
    "============================================================\n",
    "MOSQEDIT-R STEP 07 COMPLETED SUCCESSFULLY\n",
    "Leakage-safe PU modelling datasets\n",
    "============================================================\n",
    "Master genes:                              ",
    format(
        n_master,
        big.mark = ","
    ),
    "\n",
    "Observed PU-positive genes:                ",
    format(
        n_positive_master,
        big.mark = ","
    ),
    "\n",
    "External benchmark holdout genes:          ",
    format(
        n_external_holdout,
        big.mark = ","
    ),
    "\n",
    "Observed positives outside holdout:        ",
    format(
        n_positive_available_training,
        big.mark = ","
    ),
    "\n",
    "M1 expression training pool:               ",
    format(
        nrow(
            M1
        ),
        big.mark = ","
    ),
    "\n",
    "M2 expression+orthology training pool:     ",
    format(
        nrow(
            M2
        ),
        big.mark = ","
    ),
    "\n",
    "M3 +annotation training pool:              ",
    format(
        nrow(
            M3
        ),
        big.mark = ","
    ),
    "\n",
    "Targeted CRISPR validation panel:          ",
    format(
        nrow(
            targeted_crispr
        ),
        big.mark = ","
    ),
    "\n",
    "FlyBase phenotype predictors used:         NO\n",
    "Step 04 score/rank predictors used:        NO\n",
    "Step 05 CRISPR genome-wide predictors:     NO\n",
    "Benchmark genes used for model fitting:    NO\n",
    "Ag1000G predictor used:                    NO\n",
    "Python used:                               NO\n",
    "============================================================\n",
    sep = ""
)


cat(
    "\nQC summary:\n"
)


print(
    qc
)


cat(
    "\nLeakage-safe predictor sets:\n"
)


cat(
    "\nM1 expression predictors:\n",
    paste(
        m1_predictors,
        collapse = "\n"
    ),
    "\n",
    sep = ""
)


cat(
    "\nM2 additional orthology predictors:\n",
    paste(
        setdiff(
            m2_predictors,
            m1_predictors
        ),
        collapse = "\n"
    ),
    "\n",
    sep = ""
)


cat(
    "\nM3 additional annotation predictors:\n",
    paste(
        setdiff(
            m3_predictors,
            m2_predictors
        ),
        collapse = "\n"
    ),
    "\n",
    sep = ""
)


cat(
    "\nExternal benchmark holdout:\n"
)


print(
    benchmark_holdout
)


log_step(
    "07",
    "Leakage-safe PU modelling datasets completed successfully"
)

