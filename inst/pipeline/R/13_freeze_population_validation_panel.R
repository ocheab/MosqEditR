# =============================================================================
# 13_freeze_population_validation_panel.R
# MosqEdit-R Manuscript 1
# STEP 13 â€” Freeze robust targeted population-genomic validation panel
# PURE R
# =============================================================================

source("R/helpers.R")

if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.", call. = FALSE)
}
suppressPackageStartupMessages(library(data.table))

log_step("13", "Freezing targeted population-genomic validation panel")

UNCERTAINTY_FILE <- "data_processed/12_bootstrap_rank_uncertainty.csv"
CRISPR_SITES_FILE <- "data_processed/05_crispr_candidate_sites_exon_aware.csv"

PANEL_FILE <- "data_processed/13_population_validation_panel.csv"
SITE_MANIFEST_FILE <- "data_processed/13_population_site_manifest.csv"
POP_TEMPLATE_FILE <- "data_processed/13_population_site_metrics_TEMPLATE.csv"
QC_FILE <- "data_processed/13_population_panel_qc.csv"
PROVENANCE_FILE <- "data_processed/13_population_panel_provenance.csv"
SESSION_INFO_FILE <- "logs/13_population_panel_sessionInfo.txt"
CHECKSUM_FILE <- "logs/13_checksums.tsv"

NOVEL_N <- 100L
MAX_SITES_PER_GENE <- 5L

dir.create("data_processed", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

assert_file <- function(path) {
    if (!file.exists(path)) {
        stop(paste0("Required input file is missing:\n", path), call. = FALSE)
    }
}

assert_file(UNCERTAINTY_FILE)
assert_file(CRISPR_SITES_FILE)

u <- fread(UNCERTAINTY_FILE)
sites <- fread(CRISPR_SITES_FILE)

if (uniqueN(u$gene_id) != nrow(u)) {
    stop("Step 12 uncertainty table contains duplicate gene IDs.", call. = FALSE)
}
if (!"gene_id" %in% names(sites)) {
    stop("Step 05 site table lacks gene_id.", call. = FALSE)
}

required_u <- c(
    "gene_id",
    "final_prepopulation_rank",
    "is_external_benchmark",
    "probability_top20",
    "probability_top50",
    "probability_top100",
    "bootstrap_rank_median",
    "bootstrap_rank_q025",
    "bootstrap_rank_q975"
)

missing_u <- setdiff(required_u, names(u))
if (length(missing_u) > 0L) {
    stop(
        paste0(
            "Step 12 table is missing: ",
            paste(missing_u, collapse = ", ")
        ),
        call. = FALSE
    )
}

novel <- u[
    is_external_benchmark == FALSE &
    !is.na(final_prepopulation_rank)
][
    order(final_prepopulation_rank, gene_id)
][
    1:min(NOVEL_N, .N)
]

bench <- u[
    is_external_benchmark == TRUE
]

panel <- rbindlist(
    list(novel, bench),
    fill = TRUE,
    use.names = TRUE
)

panel[
    ,
    population_validation_role :=
        fifelse(
            is_external_benchmark,
            "EXTERNAL_BENCHMARK",
            fifelse(
                probability_top50 >= 0.75,
                "ROBUST_CORE_NOVEL",
                "EXTENDED_TOP100_NOVEL"
            )
        )
]

panel[
    ,
    population_validation_priority :=
        fcase(
            is_external_benchmark, 3L,
            probability_top20 >= 0.75, 1L,
            probability_top50 >= 0.75, 1L,
            default = 2L
        )
]

setorder(
    panel,
    population_validation_priority,
    final_prepopulation_rank,
    gene_id
)

# -------------------------------------------------------------------------
# Select representative CRISPR target sites for population validation.
# Population validation should be site-specific, not merely gene-wide.
# -------------------------------------------------------------------------

site_pool <- sites[
    gene_id %in% panel$gene_id
]

# Prefer the authoritative Step 05 functional-quality flag when present.
if ("passes_functional_region_reference_quality" %in% names(site_pool)) {
    site_pool <- site_pool[
        passes_functional_region_reference_quality == TRUE
    ]
} else if ("passes_reference_sequence_quality" %in% names(site_pool)) {
    site_pool <- site_pool[
        passes_reference_sequence_quality == TRUE
    ]
}

# Derive a deterministic site priority using only reference properties.
site_pool[
    ,
    site_gc_distance :=
        if ("gc_percent" %in% names(site_pool)) {
            abs(as.numeric(gc_percent) - 50)
        } else if ("guide_gc_percent" %in% names(site_pool)) {
            abs(as.numeric(guide_gc_percent) - 50)
        } else {
            0
        }
]

site_pool[
    ,
    site_unique_priority :=
        if ("sequence_unique_within_gene" %in% names(site_pool)) {
            fifelse(
                is.na(sequence_unique_within_gene),
                1L,
                fifelse(as.logical(sequence_unique_within_gene), 0L, 1L)
            )
        } else {
            0L
        }
]

site_pool[
    ,
    site_cds_priority :=
        if ("target_region" %in% names(site_pool)) {
            fifelse(target_region == "CDS", 0L, 1L)
        } else if ("region_class" %in% names(site_pool)) {
            fifelse(region_class == "CDS", 0L, 1L)
        } else {
            0L
        }
]

coord_start <- intersect(
    c("genomic_start", "target_genomic_start", "start"),
    names(site_pool)
)
coord_end <- intersect(
    c("genomic_end", "target_genomic_end", "end"),
    names(site_pool)
)

if (length(coord_start) == 0L || length(coord_end) == 0L) {
    stop(
        paste0(
            "Could not identify genomic target coordinates in ",
            CRISPR_SITES_FILE,
            "."
        ),
        call. = FALSE
    )
}

setorderv(
    site_pool,
    cols = c(
        "gene_id",
        "site_unique_priority",
        "site_cds_priority",
        "site_gc_distance",
        coord_start[[1]],
        coord_end[[1]]
    ),
    order = c(1L, 1L, 1L, 1L, 1L, 1L),
    na.last = TRUE
)

site_pool[
    ,
    population_site_rank_within_gene :=
        seq_len(.N),
    by = gene_id
]

site_manifest <- site_pool[
    population_site_rank_within_gene <= MAX_SITES_PER_GENE
]

# Join panel-level context.
panel_context_cols <- intersect(
    c(
        "gene_id",
        "final_prepopulation_rank",
        "probability_top20",
        "probability_top50",
        "probability_top100",
        "bootstrap_rank_median",
        "bootstrap_rank_q025",
        "bootstrap_rank_q975",
        "population_validation_role",
        "population_validation_priority",
        "benchmark_name",
        "benchmark_class"
    ),
    names(panel)
)

site_manifest <- merge(
    site_manifest,
    panel[, ..panel_context_cols],
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)

setorder(
    site_manifest,
    population_validation_priority,
    final_prepopulation_rank,
    gene_id,
    population_site_rank_within_gene
)

site_manifest[
    ,
    population_site_id :=
        paste0(
            gene_id,
            "__SITE",
            sprintf("%02d", population_site_rank_within_gene)
        )
]

# -------------------------------------------------------------------------
# Population metrics template.
# Fill ONLY from real targeted population-genomic analysis.
# -------------------------------------------------------------------------

template_id_cols <- intersect(
    c(
        "population_site_id",
        "gene_id",
        "genomic_seqid",
        "chromosome",
        "genomic_start",
        "target_genomic_start",
        "genomic_end",
        "target_genomic_end",
        "guide_genomic_strand",
        "protospacer",
        "pam",
        "protospacer_pam",
        "guide_sequence",
        "target_region",
        "region_class",
        "population_site_rank_within_gene",
        "final_prepopulation_rank",
        "population_validation_role"
    ),
    names(site_manifest)
)

population_template <- site_manifest[
    ,
    ..template_id_cols
]

population_template[
    ,
    `:=`(
        population_source_release = NA_character_,
        population_sample_set = NA_character_,
        population_taxa = NA_character_,
        population_accessibility_mask = NA_character_,
        population_n_samples = NA_integer_,
        population_n_populations = NA_integer_,
        population_callable_fraction = NA_real_,
        pam_intact_fraction = NA_real_,
        protospacer_exact_match_fraction = NA_real_,
        target_23bp_exact_match_fraction = NA_real_,
        max_target_variant_alt_af = NA_real_,
        mean_target_heterozygosity_proxy = NA_real_,
        populations_with_target_23bp_exact_match_ge_0_95 = NA_integer_,
        populations_assessed = NA_integer_,
        population_target_site_robustness_score = NA_real_,
        population_metric_status = "PENDING_REAL_DATA"
    )
]

fwrite(panel, PANEL_FILE, na = "NA")
fwrite(site_manifest, SITE_MANIFEST_FILE, na = "NA")
fwrite(population_template, POP_TEMPLATE_FILE, na = "NA")

qc <- data.table(
    metric = c(
        "Novel genes frozen for population validation",
        "External benchmarks included",
        "Total genes in population panel",
        "Selected CRISPR sites",
        "Maximum sites per gene",
        "Genes with at least one selected site",
        "Genes without selected reference-quality site",
        "Population metrics populated"
    ),
    value = c(
        nrow(novel),
        nrow(bench),
        nrow(panel),
        nrow(site_manifest),
        MAX_SITES_PER_GENE,
        uniqueN(site_manifest$gene_id),
        nrow(panel) - uniqueN(site_manifest$gene_id),
        0
    )
)

provenance <- data.table(
    component = c(
        "Population panel",
        "Novel candidate rule",
        "Benchmark rule",
        "Site specificity",
        "Reference-site prioritization",
        "Population data status",
        "Interpretation"
    ),
    specification = c(
        "Targeted validation only; not genome-wide population modelling.",
        paste0(
            "Top ",
            NOVEL_N,
            " novel genes by frozen Step 09 final_prepopulation_rank."
        ),
        "All six prespecified external benchmarks are retained regardless of rank.",
        paste0(
            "Up to ",
            MAX_SITES_PER_GENE,
            " reference-quality exon-contiguous CRISPR sites per gene are exported."
        ),
        paste0(
            "Within gene, reference sites are prioritized by sequence uniqueness, ",
            "CDS status where available, GC proximity to 50%, then coordinates."
        ),
        "Population metrics remain NA until produced from real Ag1000G or equivalent population-genomic data.",
        "Step 13 freezes the targeted panel; it does not infer population conservation."
    )
)

fwrite(qc, QC_FILE, na = "NA")
fwrite(provenance, PROVENANCE_FILE, na = "NA")
capture.output(sessionInfo(), file = SESSION_INFO_FILE)

checksum_files <- c(
    UNCERTAINTY_FILE,
    CRISPR_SITES_FILE,
    PANEL_FILE,
    SITE_MANIFEST_FILE,
    POP_TEMPLATE_FILE,
    QC_FILE,
    PROVENANCE_FILE,
    SESSION_INFO_FILE
)
checksum_files <- checksum_files[file.exists(checksum_files)]
write_checksum(checksum_files, CHECKSUM_FILE)

cat(
    "\n============================================================\n",
    "MOSQEDIT-R STEP 13 COMPLETED SUCCESSFULLY\n",
    "Targeted population-genomic validation panel frozen\n",
    "============================================================\n",
    "Novel genes:                              ", nrow(novel), "\n",
    "External benchmarks:                      ", nrow(bench), "\n",
    "Total genes:                              ", nrow(panel), "\n",
    "Selected CRISPR sites:                    ", nrow(site_manifest), "\n",
    "Genes with selected sites:                ", uniqueN(site_manifest$gene_id), "\n",
    "Population metrics populated:             NO\n",
    "============================================================\n",
    sep = ""
)

print(qc)

log_step("13", "Targeted population-genomic validation panel completed")

