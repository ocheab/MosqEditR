# =============================================================================
# 04_preliminary_genomewide_ranking.R
#
# MosqEdit-R Manuscript 1
#
# STEP 04 (REVISED LOW-DATA DESIGN)
# Genome-wide preliminary prioritisation using ONLY Steps 01-03.
#
# PURPOSE
# -------
# Build an explainable, genome-wide preliminary ranking from:
#   Step 01: current AgamP4 gene-like annotation universe
#   Step 02: MozAtlas adult tissue/sex expression
#   Step 03: D. melanogaster orthology + FlyBase phenotype evidence
#
# Then produce a small candidate set for later TARGETED Ag1000G validation:
#   - top N preliminary candidates
#   - all prespecified benchmark/control loci
#
# IMPORTANT INTERPRETATION
# ------------------------
# This is an EVIDENCE PRIORITISATION SCORE, not a probability of gene-drive
# success and not a fitted machine-learning model.
#
# Missing evidence is never silently converted to biological absence.
# Scores are calculated from available evidence and are accompanied by explicit
# coverage/completeness variables. Final ranking uses a modest completeness
# adjustment so a gene supported by one isolated feature does not automatically
# outrank a similarly scoring gene supported by both evidence domains.
#
# Ag1000G population-genomic evidence is intentionally deferred to a targeted
# second-stage validation because downloading chromosome-scale phased VCFs is
# unnecessarily expensive for the current discovery phase.
#
# INPUTS
# ------
# data_processed/01_gene_annotation.csv
# data_processed/02_mozatlas_gene_expression.csv
# data_processed/03_gene_orthology_phenotypes.csv
#
# OUTPUTS
# -------
# data_processed/04_preliminary_genomewide_ranking.csv
# data_processed/04_targeted_ag3_shortlist.csv
# data_processed/04_benchmark_recovery.csv
# data_processed/04_preliminary_ranking_qc.csv
# data_processed/04_preliminary_ranking_provenance.csv
# logs/04_preliminary_ranking_sessionInfo.txt
# logs/04_preliminary_ranking_checksums.tsv
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
    "04",
    "Starting genome-wide preliminary prioritisation from Steps 01-03"
)


# -----------------------------------------------------------------------------
# 2. Configuration
# -----------------------------------------------------------------------------

ANNOTATION_FILE <-
    "data_processed/01_gene_annotation.csv"

EXPRESSION_FILE <-
    "data_processed/02_mozatlas_gene_expression.csv"

ORTHOLOGY_FILE <-
    "data_processed/03_gene_orthology_phenotypes.csv"


RANKING_FILE <-
    "data_processed/04_preliminary_genomewide_ranking.csv"

SHORTLIST_FILE <-
    "data_processed/04_targeted_ag3_shortlist.csv"

BENCHMARK_FILE <-
    "data_processed/04_benchmark_recovery.csv"

QC_FILE <-
    "data_processed/04_preliminary_ranking_qc.csv"

PROVENANCE_FILE <-
    "data_processed/04_preliminary_ranking_provenance.csv"

SESSION_INFO_FILE <-
    "logs/04_preliminary_ranking_sessionInfo.txt"

CHECKSUM_FILE <-
    "logs/04_preliminary_ranking_checksums.tsv"


TOP_N_FOR_TARGETED_AG3 <- 200L


# Primary integrated weighting.
#
# Expression/reproductive relevance is deliberately the larger domain because
# the objective is discovery of A. gambiae female-reproduction targets.
#
# FlyBase is used as cross-species functional evidence, not as a required label.
PRIMARY_EXPRESSION_WEIGHT <- 0.60
PRIMARY_FUNCTIONAL_WEIGHT <- 0.40


# Female Reproductive Importance Score (FRIS) weights.
#
# Components:
#   1. ovary abundance
#   2. ovary specificity
#   3. female-vs-male bias
#   4. tissue specificity (Tau)
#
# Weights are renormalised per gene over the components actually available.
FRIS_WEIGHTS <- c(
    ovary_abundance = 0.35,
    ovary_specificity = 0.30,
    female_bias = 0.20,
    tissue_specificity = 0.15
)


# Functional-evidence weights.
#
# Explicit FlyBase female-sterility evidence is given greatest weight.
FUNCTIONAL_WEIGHTS <- c(
    strict_female_sterility = 0.75,
    broader_reproductive_evidence = 0.25
)


# Completeness adjustment:
#
# final scenario score =
#   observed-domain weighted score *
#   (COMPLETENESS_FLOOR + (1-COMPLETENESS_FLOOR)*domain_coverage)
#
# This does NOT classify missing evidence as negative; it only prevents a
# candidate supported by one isolated domain from automatically outranking an
# equally scoring candidate supported by both expression and functional data.
COMPLETENESS_FLOOR <- 0.75


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
# 3. Prespecified benchmark / external-control loci
# -----------------------------------------------------------------------------

# These are not used to train the score.
#
# They are carried through independently to assess whether the preliminary
# ranking recovers biologically established targets and to ensure that known
# noncoding / modification targets are not lost from the targeted validation
# stage.

benchmarks <- data.table(
    gene_id = c(
        "AGAP005958",
        "AGAP011377",
        "AGAP007280",
        "AGAP004050",
        "AGAP028779",
        "AGAP007031"
    ),

    benchmark_name = c(
        "Hammond female-fertility target 1",
        "Hammond female-fertility target 2",
        "Hammond female-fertility target 3",
        "doublesex",
        "mir-184",
        "FREP1"
    ),

    benchmark_class = c(
        "population-suppression female-fertility",
        "population-suppression female-fertility",
        "population-suppression female-fertility",
        "population-suppression sex-determination",
        "suppression-modification noncoding RNA",
        "population-modification / parasite-interaction"
    ),

    benchmark_reference = c(
        "Hammond et al. Nature Biotechnology 2016",
        "Hammond et al. Nature Biotechnology 2016",
        "Hammond et al. Nature Biotechnology 2016",
        "Kyrou et al. Nature Biotechnology 2018",
        "Verkuijl et al. Nature Communications 2025",
        "malaria-transmission literature; FREP1 = AGAP007031"
    )
)


# -----------------------------------------------------------------------------
# 4. Utility functions
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


assert_unique_gene_ids <- function(dt, label) {

    if (!"gene_id" %in% names(dt)) {
        stop(
            paste0(
                label,
                " does not contain a 'gene_id' column."
            ),
            call. = FALSE
        )
    }

    bad <- dt[
        !is.na(gene_id) &
            duplicated(gene_id),
        unique(gene_id)
    ]

    if (length(bad) > 0L) {
        stop(
            paste0(
                label,
                " contains duplicate gene_id values. Example(s): ",
                paste(
                    head(bad, 10L),
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
    }

    invisible(TRUE)
}


as_numeric_safe <- function(x) {

    suppressWarnings(
        as.numeric(x)
    )
}


# Percentile score among observed values only.
#
# Output is 0-1.
# NA stays NA.
rank01 <- function(x) {

    x <- as_numeric_safe(x)

    out <- rep(
        NA_real_,
        length(x)
    )

    ok <- is.finite(x)

    n_ok <- sum(ok)

    if (n_ok == 0L) {
        return(out)
    }

    if (n_ok == 1L) {
        out[ok] <- 0.5
        return(out)
    }

    r <- rank(
        x[ok],
        ties.method = "average",
        na.last = "keep"
    )

    out[ok] <- (r - 1) / (n_ok - 1)

    out
}


clip01 <- function(x) {

    x <- as_numeric_safe(x)

    x[
        is.finite(x)
    ] <- pmin(
        1,
        pmax(
            0,
            x[
                is.finite(x)
            ]
        )
    )

    x
}


row_weighted_available <- function(
        dt,
        columns,
        weights
) {

    stopifnot(
        length(columns) ==
            length(weights)
    )

    weights <- as.numeric(
        weights
    )

    m <- as.matrix(
        dt[
            ,
            ..columns
        ]
    )

    storage.mode(m) <- "double"

    observed <- is.finite(
        m
    )

    weighted_sum <- rowSums(
        sweep(
            ifelse(
                observed,
                m,
                0
            ),
            2,
            weights,
            `*`
        )
    )

    observed_weight <- rowSums(
        sweep(
            observed * 1,
            2,
            weights,
            `*`
        )
    )

    score <- weighted_sum /
        observed_weight

    score[
        observed_weight <= 0
    ] <- NA_real_

    list(
        score = score,
        observed_weight = observed_weight,
        fraction = observed_weight /
            sum(weights)
    )
}


scenario_score <- function(
        expression_score,
        functional_score,
        expression_weight,
        functional_weight,
        completeness_floor = COMPLETENESS_FLOOR
) {

    m <- cbind(
        expression_score,
        functional_score
    )

    w <- c(
        expression_weight,
        functional_weight
    )

    observed <- is.finite(
        m
    )

    numerator <- rowSums(
        sweep(
            ifelse(
                observed,
                m,
                0
            ),
            2,
            w,
            `*`
        )
    )

    denominator <- rowSums(
        sweep(
            observed * 1,
            2,
            w,
            `*`
        )
    )

    raw <- numerator /
        denominator

    raw[
        denominator <= 0
    ] <- NA_real_

    coverage <- denominator /
        sum(w)

    adjustment <-
        completeness_floor +
        (
            1 -
            completeness_floor
        ) *
        coverage

    adjusted <- raw *
        adjustment

    list(
        raw = raw,
        domain_coverage = coverage,
        completeness_adjustment = adjustment,
        adjusted = adjusted
    )
}


rank_candidates <- function(
        score,
        evidence_completeness,
        gene_id
) {

    # Lower rank is better.
    #
    # Score is primary.
    # Evidence completeness is a tie-breaker only.
    # Stable gene_id provides deterministic ordering.
    order_dt <- data.table(
        row_index = seq_along(score),
        score = score,
        evidence_completeness = evidence_completeness,
        gene_id = gene_id
    )

    order_dt[
        ,
        score_for_sort :=
            fifelse(
                is.finite(score),
                score,
                -Inf
            )
    ]

    order_dt[
        ,
        completeness_for_sort :=
            fifelse(
                is.finite(evidence_completeness),
                evidence_completeness,
                -Inf
            )
    ]

    setorder(
        order_dt,
        -score_for_sort,
        -completeness_for_sort,
        gene_id
    )

    order_dt[
        ,
        rank := seq_len(.N)
    ]

    result <- integer(
        nrow(order_dt)
    )

    result[
        order_dt$row_index
    ] <- order_dt$rank

    result[
        !is.finite(score)
    ] <- NA_integer_

    result
}


# -----------------------------------------------------------------------------
# 5. Load inputs
# -----------------------------------------------------------------------------

assert_file(
    ANNOTATION_FILE
)

assert_file(
    EXPRESSION_FILE
)

assert_file(
    ORTHOLOGY_FILE
)


annotation <- fread(
    ANNOTATION_FILE
)

expression <- fread(
    EXPRESSION_FILE
)

orthology <- fread(
    ORTHOLOGY_FILE
)


assert_unique_gene_ids(
    annotation,
    "Step 01 annotation"
)

assert_unique_gene_ids(
    expression,
    "Step 02 MozAtlas gene expression"
)

assert_unique_gene_ids(
    orthology,
    "Step 03 orthology/phenotype table"
)


annotation <- annotation[
    !is.na(gene_id) &
        grepl(
            "^AGAP",
            gene_id
        )
]


n_annotation <- nrow(
    annotation
)


log_step(
    "04",
    paste(
        "Loaded",
        format(
            n_annotation,
            big.mark = ","
        ),
        "Step 01 AGAP loci"
    )
)


# -----------------------------------------------------------------------------
# 6. Retain the Step 02 fields needed for preliminary ranking
# -----------------------------------------------------------------------------

step02_preferred <- c(
    "gene_id",
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


step02_keep <- intersect(
    step02_preferred,
    names(expression)
)


if (!"gene_id" %in% step02_keep) {
    stop(
        "Step 02 table lacks gene_id.",
        call. = FALSE
    )
}


expr_small <- expression[
    ,
    ..step02_keep
]


# -----------------------------------------------------------------------------
# 7. Retain the Step 03 fields needed for preliminary ranking
# -----------------------------------------------------------------------------

step03_preferred <- c(
    "gene_id",
    "dmel_ortholog_count",
    "dmel_ortholog",
    "dmel_gene_name",
    "has_dmel_ortholog",
    "has_flybase_phenotype_records",
    "flybase_female_sterile",
    "flybase_female_sterile_record_count",
    "flybase_reproductive_phenotype_count"
)


step03_keep <- intersect(
    step03_preferred,
    names(orthology)
)


if (!"gene_id" %in% step03_keep) {
    stop(
        "Step 03 table lacks gene_id.",
        call. = FALSE
    )
}


orth_small <- orthology[
    ,
    ..step03_keep
]


# -----------------------------------------------------------------------------
# 8. Merge to the complete Step 01 universe
# -----------------------------------------------------------------------------

x <- merge(
    annotation,
    expr_small,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


x <- merge(
    x,
    orth_small,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


if (nrow(x) != n_annotation) {
    stop(
        paste0(
            "Merge changed the Step 01 universe size.\n",
            "Before merge: ",
            n_annotation,
            "\nAfter merge: ",
            nrow(x)
        ),
        call. = FALSE
    )
}


if (uniqueN(x$gene_id) != nrow(x)) {
    stop(
        "Merged Step 04 table contains duplicate gene_id values.",
        call. = FALSE
    )
}


# -----------------------------------------------------------------------------
# 9. Expression-derived evidence variables
# -----------------------------------------------------------------------------

# Use log2(x+1) for intensity contrasts where raw intensity fields are present.
#
# Ranking ovary abundance itself is monotonic under log transformation, but
# transformed values are retained because they are easier to interpret and
# safer for later modelling.

if ("expr_ovary" %in% names(x)) {

    x[
        ,
        log2_expr_ovary :=
            log2(
                pmax(
                    as_numeric_safe(
                        expr_ovary
                    ),
                    0
                ) +
                1
            )
    ]

} else {

    x[
        ,
        log2_expr_ovary :=
            NA_real_
    ]
}


if (
    all(
        c(
            "expr_female_body",
            "expr_male_body"
        ) %in%
        names(x)
    )
) {

    x[
        ,
        log2_female_male_contrast :=
            log2(
                pmax(
                    as_numeric_safe(
                        expr_female_body
                    ),
                    0
                ) +
                1
            ) -
            log2(
                pmax(
                    as_numeric_safe(
                        expr_male_body
                    ),
                    0
                ) +
                1
            )
    ]

} else if (
    "female_male_delta" %in%
    names(x)
) {

    # Fall back to the Step 02 contrast if body intensities are not present.
    x[
        ,
        log2_female_male_contrast :=
            as_numeric_safe(
                female_male_delta
            )
    ]

} else {

    x[
        ,
        log2_female_male_contrast :=
            NA_real_
    ]
}


if (
    all(
        c(
            "expr_ovary",
            "expr_testis",
            "expr_female_body",
            "expr_male_body",
            "expr_female_midgut",
            "expr_female_salivary"
        ) %in%
        names(x)
    )
) {

    non_ovary_matrix <- cbind(
        log2(
            pmax(
                as_numeric_safe(
                    x$expr_testis
                ),
                0
            ) +
            1
        ),
        log2(
            pmax(
                as_numeric_safe(
                    x$expr_female_body
                ),
                0
            ) +
            1
        ),
        log2(
            pmax(
                as_numeric_safe(
                    x$expr_male_body
                ),
                0
            ) +
            1
        ),
        log2(
            pmax(
                as_numeric_safe(
                    x$expr_female_midgut
                ),
                0
            ) +
            1
        ),
        log2(
            pmax(
                as_numeric_safe(
                    x$expr_female_salivary
                ),
                0
            ) +
            1
        )
    )


    max_non_ovary <- apply(
        non_ovary_matrix,
        1,
        function(z) {

            if (all(!is.finite(z))) {
                return(NA_real_)
            }

            max(
                z[
                    is.finite(z)
                ]
            )
        }
    )


    x[
        ,
        log2_ovary_specificity_contrast :=
            log2_expr_ovary -
            max_non_ovary
    ]

} else if (
    "ovary_specificity_delta" %in%
    names(x)
) {

    x[
        ,
        log2_ovary_specificity_contrast :=
            as_numeric_safe(
                ovary_specificity_delta
            )
    ]

} else {

    x[
        ,
        log2_ovary_specificity_contrast :=
            NA_real_
    ]
}


if (
    "tau_tissue_specificity" %in%
    names(x)
) {

    x[
        ,
        tau_for_scoring :=
            clip01(
                tau_tissue_specificity
            )
    ]

} else {

    x[
        ,
        tau_for_scoring :=
            NA_real_
    ]
}


# Percentile components.
x[
    ,
    score_ovary_abundance :=
        rank01(
            log2_expr_ovary
        )
]


x[
    ,
    score_ovary_specificity :=
        rank01(
            log2_ovary_specificity_contrast
        )
]


x[
    ,
    score_female_bias :=
        rank01(
            log2_female_male_contrast
        )
]


x[
    ,
    score_tissue_specificity :=
        tau_for_scoring
]


fris_result <- row_weighted_available(
    dt = x,

    columns = c(
        "score_ovary_abundance",
        "score_ovary_specificity",
        "score_female_bias",
        "score_tissue_specificity"
    ),

    weights = FRIS_WEIGHTS
)


x[
    ,
    fris :=
        fris_result$score
]


x[
    ,
    fris_feature_coverage :=
        fris_result$fraction
]


# -----------------------------------------------------------------------------
# 10. Functional evidence from Drosophila/FlyBase
# -----------------------------------------------------------------------------

if (
    "has_flybase_phenotype_records" %in%
    names(x)
) {

    x[
        ,
        has_flybase_records_for_scoring :=
            fifelse(
                is.na(
                    has_flybase_phenotype_records
                ),
                NA,
                as.logical(
                    has_flybase_phenotype_records
                )
            )
    ]

} else {

    x[
        ,
        has_flybase_records_for_scoring :=
            NA
    ]
}


if (
    "flybase_female_sterile" %in%
    names(x)
) {

    # IMPORTANT:
    # FALSE is meaningful only if phenotype records were actually mapped.
    x[
        ,
        score_strict_female_sterility :=
            fifelse(
                has_flybase_records_for_scoring == TRUE,
                fifelse(
                    as.logical(
                        flybase_female_sterile
                    ) == TRUE,
                    1,
                    0
                ),
                NA_real_
            )
    ]

} else {

    x[
        ,
        score_strict_female_sterility :=
            NA_real_
    ]
}


if (
    "flybase_reproductive_phenotype_count" %in%
    names(x)
) {

    repro_count <- as_numeric_safe(
        x$flybase_reproductive_phenotype_count
    )


    repro_count[
        x$has_flybase_records_for_scoring !=
            TRUE
    ] <- NA_real_


    x[
        ,
        score_broader_reproductive_evidence :=
            rank01(
                log1p(
                    pmax(
                        repro_count,
                        0
                    )
                )
            )
    ]

} else {

    x[
        ,
        score_broader_reproductive_evidence :=
            NA_real_
    ]
}


functional_result <- row_weighted_available(
    dt = x,

    columns = c(
        "score_strict_female_sterility",
        "score_broader_reproductive_evidence"
    ),

    weights = FUNCTIONAL_WEIGHTS
)


x[
    ,
    functional_evidence_score :=
        functional_result$score
]


x[
    ,
    functional_feature_coverage :=
        functional_result$fraction
]


# -----------------------------------------------------------------------------
# 11. Evidence-domain availability and completeness
# -----------------------------------------------------------------------------

x[
    ,
    expression_domain_available :=
        is.finite(
            fris
        )
]


x[
    ,
    functional_domain_available :=
        is.finite(
            functional_evidence_score
        )
]


x[
    ,
    n_evidence_domains_available :=
        as.integer(
            expression_domain_available
        ) +
        as.integer(
            functional_domain_available
        )
]


x[
    ,
    evidence_domain_coverage :=
        n_evidence_domains_available /
        2
]


# A finer-grained completeness indicator for descriptive QC.
x[
    ,
    evidence_feature_completeness :=
        rowMeans(
            cbind(
                fris_feature_coverage,
                functional_feature_coverage
            ),
            na.rm = TRUE
        )
]


x[
    !is.finite(
        evidence_feature_completeness
    ),
    evidence_feature_completeness :=
        NA_real_
]


# -----------------------------------------------------------------------------
# 12. Three transparent weighting scenarios
# -----------------------------------------------------------------------------

# Scenario A: primary/balanced discovery model
sc_primary <- scenario_score(
    expression_score =
        x$fris,

    functional_score =
        x$functional_evidence_score,

    expression_weight =
        PRIMARY_EXPRESSION_WEIGHT,

    functional_weight =
        PRIMARY_FUNCTIONAL_WEIGHT
)


x[
    ,
    score_primary_raw :=
        sc_primary$raw
]


x[
    ,
    score_primary :=
        sc_primary$adjusted
]


# Scenario B: expression-dominant
sc_expr <- scenario_score(
    expression_score =
        x$fris,

    functional_score =
        x$functional_evidence_score,

    expression_weight = 0.75,
    functional_weight = 0.25
)


x[
    ,
    score_expression_dominant :=
        sc_expr$adjusted
]


# Scenario C: phenotype-dominant
sc_func <- scenario_score(
    expression_score =
        x$fris,

    functional_score =
        x$functional_evidence_score,

    expression_weight = 0.40,
    functional_weight = 0.60
)


x[
    ,
    score_functional_dominant :=
        sc_func$adjusted
]


# -----------------------------------------------------------------------------
# 13. Ranks and weight-sensitivity consensus
# -----------------------------------------------------------------------------

x[
    ,
    rank_primary :=
        rank_candidates(
            score_primary,
            evidence_feature_completeness,
            gene_id
        )
]


x[
    ,
    rank_expression_dominant :=
        rank_candidates(
            score_expression_dominant,
            evidence_feature_completeness,
            gene_id
        )
]


x[
    ,
    rank_functional_dominant :=
        rank_candidates(
            score_functional_dominant,
            evidence_feature_completeness,
            gene_id
        )
]


rank_matrix <- cbind(
    x$rank_primary,
    x$rank_expression_dominant,
    x$rank_functional_dominant
)


x[
    ,
    rank_consensus_median :=
        apply(
            rank_matrix,
            1,
            function(z) {

                z <- z[
                    is.finite(z)
                ]

                if (length(z) == 0L) {
                    return(
                        NA_real_
                    )
                }

                median(z)
            }
        )
]


x[
    ,
    rank_sensitivity_range :=
        apply(
            rank_matrix,
            1,
            function(z) {

                z <- z[
                    is.finite(z)
                ]

                if (length(z) == 0L) {
                    return(
                        NA_real_
                    )
                }

                max(z) -
                    min(z)
            }
        )
]


# Deterministic final preliminary ordering:
#   1. median rank across weighting scenarios
#   2. primary score
#   3. evidence completeness
#   4. gene_id
#
# IMPORTANT:
# Assign preliminary ranks back to x by the stable key gene_id, NOT by row
# position. Using .I from a filtered data.table can mis-map ranks to unrelated
# genes.

ranking_order <- x[
    !is.na(
        rank_consensus_median
    ),
    .(
        gene_id,
        rank_consensus_median,
        score_primary,
        evidence_feature_completeness
    )
]


if (
    uniqueN(ranking_order$gene_id) !=
        nrow(ranking_order)
) {

    stop(
        "Duplicate gene_id values detected while constructing preliminary ranking.",
        call. = FALSE
    )
}


ranking_order[
    ,
    score_primary_for_sort :=
        data.table::fifelse(
            is.finite(score_primary),
            score_primary,
            -Inf
        )
]


ranking_order[
    ,
    evidence_completeness_for_sort :=
        data.table::fifelse(
            is.finite(evidence_feature_completeness),
            evidence_feature_completeness,
            -Inf
        )
]


data.table::setorder(
    ranking_order,
    rank_consensus_median,
    -score_primary_for_sort,
    -evidence_completeness_for_sort,
    gene_id
)


ranking_order[
    ,
    preliminary_rank :=
        seq_len(.N)
]


x[
    ,
    preliminary_rank :=
        NA_integer_
]


# Stable key-based assignment.
x[
    ranking_order,
    on = "gene_id",
    preliminary_rank :=
        i.preliminary_rank
]


# Ranking integrity checks.
if (
    any(
        is.na(x$rank_consensus_median) &
            !is.na(x$preliminary_rank)
    )
) {

    stop(
        paste0(
            "Ranking integrity failure: at least one gene without a consensus ",
            "rank received a preliminary rank."
        ),
        call. = FALSE
    )
}


if (
    any(
        !is.na(x$rank_consensus_median) &
            is.na(x$preliminary_rank)
    )
) {

    stop(
        paste0(
            "Ranking integrity failure: at least one gene with a consensus ",
            "rank did not receive a preliminary rank."
        ),
        call. = FALSE
    )
}


if (
    uniqueN(
        x[
            !is.na(preliminary_rank),
            preliminary_rank
        ]
    ) !=
    x[
        !is.na(preliminary_rank),
        .N
    ]
) {

    stop(
        "Ranking integrity failure: preliminary ranks are not unique.",
        call. = FALSE
    )
}


ranking_order[
    ,
    c(
        "score_primary_for_sort",
        "evidence_completeness_for_sort"
    ) := NULL
]


# -----------------------------------------------------------------------------
# 14. Add benchmark annotations
# -----------------------------------------------------------------------------

x <- merge(
    x,
    benchmarks,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


x[
    ,
    is_prespecified_benchmark :=
        !is.na(
            benchmark_name
        )
]


# -----------------------------------------------------------------------------
# 15. Create targeted Ag1000G validation shortlist
# -----------------------------------------------------------------------------

top_candidate_ids <- x[
    !is.na(
        preliminary_rank
    ) &
        preliminary_rank <=
        TOP_N_FOR_TARGETED_AG3,
    gene_id
]


benchmark_ids <- benchmarks$gene_id


shortlist_ids <- unique(
    c(
        top_candidate_ids,
        benchmark_ids
    )
)


shortlist <- x[
    gene_id %in%
        shortlist_ids
]


benchmark_ids_present <- intersect(
    benchmark_ids,
    x$gene_id
)


benchmark_ids_missing <- setdiff(
    benchmark_ids,
    x$gene_id
)


shortlist[
    ,
    selection_reason :=
        data.table::fifelse(
            (gene_id %in% top_candidate_ids) &
                (gene_id %in% benchmark_ids_present),
            paste0(
                "top_",
                TOP_N_FOR_TARGETED_AG3,
                "_preliminary_and_benchmark"
            ),
            data.table::fifelse(
                gene_id %in%
                    top_candidate_ids,
                paste0(
                    "top_",
                    TOP_N_FOR_TARGETED_AG3,
                    "_preliminary"
                ),
                "benchmark_added_independently"
            )
        )
]


shortlist[
    ,
    selected_as_top_candidate :=
        gene_id %in%
            top_candidate_ids
]


shortlist[
    ,
    selected_as_benchmark :=
        gene_id %in%
            benchmark_ids_present
]


# Create chromosome region strings for later targeted regional queries.
coord_fields <- c(
    "chromosome",
    "gene_start",
    "gene_end"
)


if (
    all(
        coord_fields %in%
        names(shortlist)
    )
) {

    shortlist[
        ,
        ag3_query_region :=
            fifelse(
                !is.na(
                    chromosome
                ) &
                    !is.na(
                        gene_start
                    ) &
                    !is.na(
                        gene_end
                    ),
                paste0(
                    chromosome,
                    ":",
                    gene_start,
                    "-",
                    gene_end
                ),
                NA_character_
            )
    ]

} else {

    shortlist[
        ,
        ag3_query_region :=
            NA_character_
    ]
}


shortlist[
    ,
    preliminary_rank_missing :=
        is.na(preliminary_rank)
]


data.table::setorder(
    shortlist,
    preliminary_rank_missing,
    preliminary_rank,
    gene_id
)


shortlist[
    ,
    preliminary_rank_missing :=
        NULL
]


# -----------------------------------------------------------------------------
# 16. Benchmark recovery table
# -----------------------------------------------------------------------------

benchmark_recovery <- merge(
    benchmarks,
    x[
        ,
        .(
            gene_id,
            gene_present_in_step01 = TRUE,
            preliminary_rank,
            rank_primary,
            rank_expression_dominant,
            rank_functional_dominant,
            rank_consensus_median,
            rank_sensitivity_range,
            fris,
            functional_evidence_score,
            score_primary,
            expression_domain_available,
            functional_domain_available,
            evidence_domain_coverage
        )
    ],
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


benchmark_recovery[
    is.na(
        gene_present_in_step01
    ),
    gene_present_in_step01 :=
        FALSE
]


# -----------------------------------------------------------------------------
# 17. Stable output ordering
# -----------------------------------------------------------------------------

x[
    ,
    preliminary_rank_missing :=
        is.na(preliminary_rank)
]


data.table::setorder(
    x,
    preliminary_rank_missing,
    preliminary_rank,
    gene_id
)


x[
    ,
    preliminary_rank_missing :=
        NULL
]


# -----------------------------------------------------------------------------
# 18. QC tables
# -----------------------------------------------------------------------------

n_with_expression <- x[
    expression_domain_available == TRUE,
    .N
]


n_with_functional <- x[
    functional_domain_available == TRUE,
    .N
]


n_with_both <- x[
    expression_domain_available == TRUE &
        functional_domain_available == TRUE,
    .N
]


n_with_any <- x[
    n_evidence_domains_available > 0,
    .N
]


n_without_any <- x[
    n_evidence_domains_available == 0,
    .N
]


n_flybase_female_sterile <- x[
    score_strict_female_sterility == 1,
    .N
]


n_benchmark_present <- benchmark_recovery[
    gene_present_in_step01 == TRUE,
    .N
]


n_benchmark_top_n <- benchmark_recovery[
    !is.na(
        preliminary_rank
    ) &
        preliminary_rank <=
        TOP_N_FOR_TARGETED_AG3,
    .N
]


qc <- data.table(
    metric = c(
        "Step 01 AGAP loci",
        "Genes with MozAtlas-derived expression-domain score",
        "Genes with FlyBase functional-domain score",
        "Genes with both expression and functional domains",
        "Genes with at least one evidence domain",
        "Genes with neither evidence domain",
        "Genes with strict FlyBase female-sterility evidence",
        "Prespecified benchmark loci",
        "Benchmark loci present in current Step 01 universe",
        paste0(
            "Benchmark loci already in top ",
            TOP_N_FOR_TARGETED_AG3
        ),
        paste0(
            "Top preliminary genes requested for targeted Ag1000G"
        ),
        "Total targeted Ag1000G shortlist after adding benchmarks",
        "Primary expression-domain weight",
        "Primary functional-domain weight",
        "Completeness adjustment floor"
    ),

    value = c(
        n_annotation,
        n_with_expression,
        n_with_functional,
        n_with_both,
        n_with_any,
        n_without_any,
        n_flybase_female_sterile,
        nrow(benchmarks),
        n_benchmark_present,
        n_benchmark_top_n,
        TOP_N_FOR_TARGETED_AG3,
        nrow(shortlist),
        PRIMARY_EXPRESSION_WEIGHT,
        PRIMARY_FUNCTIONAL_WEIGHT,
        COMPLETENESS_FLOOR
    )
)


# -----------------------------------------------------------------------------
# 19. Provenance / scoring specification
# -----------------------------------------------------------------------------

provenance <- data.table(
    component = c(
        "Genome universe",
        "Expression evidence",
        "Functional evidence",
        "FRIS",
        "Functional score",
        "Primary integrated score",
        "Weight sensitivity",
        "Missing-evidence handling",
        "Targeted Ag1000G shortlist",
        "Population genomics"
    ),

    specification = c(
        "Current AGAP gene-like loci from Step 01 AgamP4 annotation.",
        paste0(
            "MozAtlas Step 02 adult tissue/sex evidence. ",
            "Raw intensities are log2(x+1)-transformed only for derived contrasts; ",
            "percentile components are calculated among observed genes."
        ),
        paste0(
            "Step 03 D. melanogaster orthology + FlyBase phenotype evidence. ",
            "No phenotype record is represented as missing evidence, not FALSE."
        ),
        paste0(
            "Female Reproductive Importance Score = available-weight-normalised ",
            "combination of ovary abundance (0.35), ovary specificity (0.30), ",
            "female-vs-male bias (0.20), and Tau tissue specificity (0.15)."
        ),
        paste0(
            "Available-weight-normalised functional score = strict mapped ",
            "FlyBase female-sterility evidence (0.75) + broader reproductive ",
            "phenotype evidence (0.25)."
        ),
        paste0(
            "Primary preliminary evidence score uses 0.60 FRIS + 0.40 functional ",
            "evidence over observed domains, then applies a modest evidence-domain ",
            "completeness adjustment with floor 0.75."
        ),
        paste0(
            "Sensitivity scenarios: expression-dominant 0.75/0.25 and ",
            "functional-dominant 0.40/0.60. Preliminary ordering is based on ",
            "median scenario rank, then primary score and evidence completeness."
        ),
        paste0(
            "Missing evidence is not imputed as biological absence. All domain ",
            "and feature coverage variables are retained explicitly."
        ),
        paste0(
            "Top ",
            TOP_N_FOR_TARGETED_AG3,
            " preliminary candidates plus all prespecified benchmark loci."
        ),
        paste0(
            "Deferred to targeted second-stage validation. Population-genomic ",
            "variables are therefore NOT genome-wide predictors in this version ",
            "of MosqEdit-R."
        )
    )
)


# -----------------------------------------------------------------------------
# 20. Write outputs
# -----------------------------------------------------------------------------

fwrite(
    x,
    RANKING_FILE,
    na = "NA"
)


fwrite(
    shortlist,
    SHORTLIST_FILE,
    na = "NA"
)


fwrite(
    benchmark_recovery,
    BENCHMARK_FILE,
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
    ORTHOLOGY_FILE,
    RANKING_FILE,
    SHORTLIST_FILE,
    BENCHMARK_FILE,
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
# 21. Console summary
# -----------------------------------------------------------------------------

cat(
    "\n",
    "============================================================\n",
    "MOSQEDIT-R STEP 04 COMPLETED SUCCESSFULLY\n",
    "Genome-wide preliminary ranking + targeted-validation shortlist\n",
    "============================================================\n",
    "Step 01 AGAP loci:                      ",
    format(
        n_annotation,
        big.mark = ","
    ),
    "\n",
    "Genes with expression evidence:         ",
    format(
        n_with_expression,
        big.mark = ","
    ),
    "\n",
    "Genes with functional evidence:         ",
    format(
        n_with_functional,
        big.mark = ","
    ),
    "\n",
    "Genes with both evidence domains:        ",
    format(
        n_with_both,
        big.mark = ","
    ),
    "\n",
    "Genes with any ranking evidence:         ",
    format(
        n_with_any,
        big.mark = ","
    ),
    "\n",
    "Strict FlyBase female-sterility loci:    ",
    format(
        n_flybase_female_sterile,
        big.mark = ","
    ),
    "\n",
    "Requested top candidates for Ag1000G:    ",
    TOP_N_FOR_TARGETED_AG3,
    "\n",
    "Final shortlist after benchmark union:   ",
    nrow(shortlist),
    "\n",
    "Population genomics used genome-wide:    NO\n",
    "Python used:                             NO\n",
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
    "\nTop 20 preliminary candidates:\n"
)

print(
    x[
        !is.na(
            preliminary_rank
        )
    ][
        1:min(
            20L,
            .N
        ),
        .(
            preliminary_rank,
            gene_id,
            fris,
            functional_evidence_score,
            score_primary,
            rank_consensus_median,
            rank_sensitivity_range,
            evidence_domain_coverage,
            benchmark_name
        )
    ]
)


cat(
    "\nBenchmark recovery:\n"
)

print(
    benchmark_recovery[
        ,
        .(
            gene_id,
            benchmark_name,
            gene_present_in_step01,
            preliminary_rank,
            rank_consensus_median,
            score_primary,
            expression_domain_available,
            functional_domain_available
        )
    ]
)


cat(
    "\nTargeted Ag1000G shortlist written to:\n",
    SHORTLIST_FILE,
    "\n",
    sep = ""
)


log_step(
    "04",
    "Preliminary genome-wide prioritisation completed successfully"
)

