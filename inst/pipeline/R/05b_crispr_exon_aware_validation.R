# =============================================================================
# 05b_crispr_exon_aware_validation.R
#
# MosqEdit-R Manuscript 1
#
# STEP 05b
# Exon-aware refinement of the Step 05 transcript-based SpCas9 scan.
#
# PURPOSE
# -------
# Remove transcript-derived candidate sites that cross exon-exon junctions and
# distinguish:
#   - CDS-contained targets in protein-coding transcripts,
#   - UTR/other exonic targets,
#   - noncoding-exon targets.
#
# This uses the same release-pinned AgamP4.63 GFF3 annotation used in Step 01.
# It DOES NOT require a whole-genome FASTA or Ag1000G VCF download.
#
# INPUTS
# ------
# data_processed/05_crispr_candidate_guides_transcript_level.csv
# data_processed/05_crispr_editability.csv
# release-63 AgamP4 GFF3 (auto-discovered or downloaded if absent)
#
# OUTPUTS
# -------
# data_processed/05b_crispr_guides_exon_aware.csv
# data_processed/05b_crispr_editability_exon_aware.csv
# data_processed/05b_crispr_editability_benchmarks.csv
# data_processed/05b_crispr_exon_aware_qc.csv
# data_processed/05b_crispr_exon_aware_provenance.csv
# logs/05b_crispr_exon_aware_sessionInfo.txt
# logs/05b_crispr_exon_aware_checksums.tsv
#
# INTERPRETATION
# --------------
# "functional-region guide" means:
#   - full 23-nt target lies within CDS for genes with CDS annotation;
#   - full 23-nt target lies within an exon for genes without CDS annotation.
#
# Population conservation and genome-wide off-target specificity remain pending.
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Helpers and packages
# -----------------------------------------------------------------------------

source("R/helpers.R")


required_cran <- c(
    "data.table",
    "curl"
)


missing_cran <- required_cran[
    !vapply(
        required_cran,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]


if (length(missing_cran) > 0L) {

    stop(
        paste0(
            "Missing CRAN package(s): ",
            paste(missing_cran, collapse = ", "),
            "\nInstall with:\n",
            "install.packages(c(",
            paste0('"', missing_cran, '"', collapse = ", "),
            "))"
        ),
        call. = FALSE
    )
}


suppressPackageStartupMessages({

    library(data.table)
    library(curl)

})


log_step(
    "05b",
    "Starting exon-aware CRISPR target refinement"
)


# -----------------------------------------------------------------------------
# 2. Configuration
# -----------------------------------------------------------------------------

GUIDE_INSTANCE_FILE <-
    "data_processed/05_crispr_candidate_guides_transcript_level.csv"


STEP05_EDITABILITY_FILE <-
    "data_processed/05_crispr_editability.csv"


GUIDE_OUTPUT_FILE <-
    "data_processed/05b_crispr_guides_exon_aware.csv"


EDITABILITY_OUTPUT_FILE <-
    "data_processed/05b_crispr_editability_exon_aware.csv"


BENCHMARK_OUTPUT_FILE <-
    "data_processed/05b_crispr_editability_benchmarks.csv"


QC_FILE <-
    "data_processed/05b_crispr_exon_aware_qc.csv"


PROVENANCE_FILE <-
    "data_processed/05b_crispr_exon_aware_provenance.csv"


SESSION_INFO_FILE <-
    "logs/05b_crispr_exon_aware_sessionInfo.txt"


CHECKSUM_FILE <-
    "logs/05b_crispr_exon_aware_checksums.tsv"


GFF_URL <-
    paste0(
        "https://ftp.ebi.ac.uk/ensemblgenomes/pub/metazoa/",
        "release-63/gff3/anopheles_gambiae/",
        "Anopheles_gambiae.AgamP4.63.gff3.gz"
    )


GFF_FALLBACK_FILE <-
    "data_raw/ensembl_release63_agamp4/Anopheles_gambiae.AgamP4.63.gff3.gz"


MULTIPLEX_TARGET_2 <- 2L
MULTIPLEX_TARGET_4 <- 4L


dir.create(
    dirname(GFF_FALLBACK_FILE),
    recursive = TRUE,
    showWarnings = FALSE
)


# -----------------------------------------------------------------------------
# 3. Basic checks
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


assert_file(
    GUIDE_INSTANCE_FILE
)


assert_file(
    STEP05_EDITABILITY_FILE
)


# -----------------------------------------------------------------------------
# 4. Locate or obtain the release-63 GFF3
# -----------------------------------------------------------------------------

gff_candidates <- list.files(
    "data_raw",
    pattern = "Anopheles_gambiae\\.AgamP4\\.63\\.gff3(\\.gz)?$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = FALSE
)


if (length(gff_candidates) > 0L) {

    gff_candidates <- gff_candidates[
        order(
            file.info(gff_candidates)$size,
            decreasing = TRUE
        )
    ]


    GFF_FILE <- gff_candidates[[1]]


    message(
        "Using existing release-63 GFF3:\n",
        GFF_FILE
    )

} else {

    GFF_FILE <- GFF_FALLBACK_FILE


    message(
        "Release-63 GFF3 not found locally; downloading only the small annotation file:\n",
        GFF_URL
    )


    tmp <- paste0(
        GFF_FILE,
        ".part"
    )


    if (file.exists(tmp)) {
        unlink(tmp)
    }


    tryCatch(

        curl::curl_download(
            url = GFF_URL,
            destfile = tmp,
            quiet = FALSE,
            mode = "wb"
        ),

        error = function(e) {

            if (file.exists(tmp)) {
                unlink(tmp)
            }


            stop(
                paste0(
                    "Could not download release-63 GFF3.\n\n",
                    conditionMessage(e)
                ),
                call. = FALSE
            )
        }
    )


    if (
        !file.exists(tmp) ||
        file.info(tmp)$size < 100000L
    ) {

        stop(
            "Downloaded GFF3 is unexpectedly small.",
            call. = FALSE
        )
    }


    if (file.exists(GFF_FILE)) {
        unlink(GFF_FILE)
    }


    if (!file.rename(tmp, GFF_FILE)) {

        stop(
            "Could not move downloaded GFF3 into place.",
            call. = FALSE
        )
    }
}


# -----------------------------------------------------------------------------
# 5. Read GFF3
# -----------------------------------------------------------------------------

read_gff_lines <- function(path) {

    con <- if (
        grepl(
            "\\.gz$",
            path,
            ignore.case = TRUE
        )
    ) {

        gzfile(
            path,
            open = "rt"
        )

    } else {

        file(
            path,
            open = "rt"
        )
    }


    on.exit(
        close(con),
        add = TRUE
    )


    lines <- readLines(
        con,
        warn = FALSE
    )


    lines[
        nzchar(lines) &
            !startsWith(
                lines,
                "#"
            )
    ]
}


log_step(
    "05b",
    "Reading release-63 GFF3 exon/CDS annotation"
)


gff_lines <- read_gff_lines(
    GFF_FILE
)


gff <- fread(
    text = paste(
        gff_lines,
        collapse = "\n"
    ),
    sep = "\t",
    header = FALSE,
    quote = "",
    fill = TRUE,
    showProgress = FALSE
)


if (ncol(gff) < 9L) {

    stop(
        "GFF3 did not parse into the expected 9 columns.",
        call. = FALSE
    )
}


setnames(
    gff,
    names(gff)[1:9],
    c(
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes"
    )
)


gff <- gff[
    type %in%
        c(
            "exon",
            "CDS"
        )
]


if (nrow(gff) == 0L) {

    stop(
        "No exon/CDS records were found in the release-63 GFF3.",
        call. = FALSE
    )
}


# -----------------------------------------------------------------------------
# 6. Parse transcript Parent identifiers
# -----------------------------------------------------------------------------

extract_attribute <- function(
        x,
        key
) {

    pattern <- paste0(
        "(?:^|;)",
        key,
        "=([^;]+)"
    )


    m <- regexec(
        pattern,
        x,
        perl = TRUE
    )


    hits <- regmatches(
        x,
        m
    )


    out <- vapply(
        hits,
        function(z) {

            if (length(z) >= 2L) {
                z[[2]]
            } else {
                NA_character_
            }
        },
        character(1)
    )


    out
}


clean_parent_id <- function(x) {

    x <- as.character(
        x
    )


    # If multiple parents are ever present, use the first here. Ensembl normally
    # emits transcript-specific exon/CDS records for this annotation.
    x <- sub(
        ",.*$",
        "",
        x
    )


    x <- sub(
        "^transcript:",
        "",
        x
    )


    x
}


gff[
    ,
    transcript_id :=
        clean_parent_id(
            extract_attribute(
                attributes,
                "Parent"
            )
        )
]


gff <- gff[
    !is.na(transcript_id) &
        nzchar(transcript_id)
]


gff[
    ,
    `:=`(
        start = as.integer(start),
        end = as.integer(end)
    )
]


# -----------------------------------------------------------------------------
# 7. Load Step 05 guide instances and final Step 05 gene table
# -----------------------------------------------------------------------------

guides <- fread(
    GUIDE_INSTANCE_FILE
)


step05 <- fread(
    STEP05_EDITABILITY_FILE
)


required_guide_fields <- c(
    "gene_id",
    "transcript_id",
    "strand",
    "transcript_start",
    "transcript_end",
    "protospacer_20nt",
    "pam",
    "guide_gc_percent",
    "gc_balance_score",
    "passes_reference_sequence_quality"
)


missing_guide_fields <- setdiff(
    required_guide_fields,
    names(guides)
)


if (length(missing_guide_fields) > 0L) {

    stop(
        paste0(
            "Step 05 guide file is missing required fields: ",
            paste(
                missing_guide_fields,
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}


target_transcripts <- unique(
    guides$transcript_id
)


target_genes <- unique(
    step05$gene_id
)


# Keep only annotation relevant to this Step 05 shortlist.
gff <- gff[
    transcript_id %in%
        target_transcripts
]


# -----------------------------------------------------------------------------
# 8. Build transcript-coordinate exon map
# -----------------------------------------------------------------------------

exons <- gff[
    type == "exon"
]


if (nrow(exons) == 0L) {

    stop(
        "No exon annotations matched the Step 05 transcript IDs.",
        call. = FALSE
    )
}


# Transcript orientation:
# + strand exons progress from low to high genomic coordinate.
# - strand exons progress from high to low genomic coordinate.
exons[
    ,
    transcript_order_key :=
        data.table::fifelse(
            strand == "-",
            -end,
            start
        )
]


data.table::setorder(
    exons,
    transcript_id,
    transcript_order_key,
    start,
    end
)


exons[
    ,
    exon_length_bp :=
        end -
        start +
        1L
]


exons[
    ,
    exon_rank :=
        seq_len(.N),
    by = transcript_id
]


exons[
    ,
    exon_cdna_end :=
        cumsum(
            exon_length_bp
        ),
    by = transcript_id
]


exons[
    ,
    exon_cdna_start :=
        exon_cdna_end -
        exon_length_bp +
        1L
]


# One transcript should not switch chromosome or gene strand.
exon_structure_qc <- exons[
    ,
    .(
        n_seqids = uniqueN(seqid),
        n_strands = uniqueN(strand),
        n_exons = .N,
        reconstructed_cdna_length =
            max(exon_cdna_end)
    ),
    by = transcript_id
]


bad_structure <- exon_structure_qc[
    n_seqids != 1L |
        n_strands != 1L,
    transcript_id
]


if (length(bad_structure) > 0L) {

    stop(
        paste0(
            "Unexpected multi-seqid/multi-strand exon structure for transcript(s): ",
            paste(
                head(
                    bad_structure,
                    10L
                ),
                collapse = ", "
            )
        ),
        call. = FALSE
    )
}


# -----------------------------------------------------------------------------
# 9. Build CDS intervals in transcript coordinates
# -----------------------------------------------------------------------------

cds <- gff[
    type == "CDS"
]


cds_map_list <- list()


if (nrow(cds) > 0L) {

    cds_split <- split(
        cds,
        by = "transcript_id",
        keep.by = TRUE
    )


    for (
        tx in names(cds_split)
    ) {

        tx_cds <- cds_split[[tx]]


        tx_exons <- exons[
            transcript_id == tx
        ]


        if (nrow(tx_exons) == 0L) {
            next
        }


        mapped <- vector(
            "list",
            nrow(tx_cds)
        )


        for (
            i in seq_len(
                nrow(tx_cds)
            )
        ) {

            c0 <- tx_cds[i]


            ex <- tx_exons[
                start <= c0$start &
                    end >= c0$end
            ]


            if (nrow(ex) == 0L) {
                next
            }


            ex <- ex[1L]


            if (ex$strand == "+") {

                cdna_start <-
                    ex$exon_cdna_start +
                    (
                        c0$start -
                        ex$start
                    )


                cdna_end <-
                    ex$exon_cdna_start +
                    (
                        c0$end -
                        ex$start
                    )

            } else {

                cdna_start <-
                    ex$exon_cdna_start +
                    (
                        ex$end -
                        c0$end
                    )


                cdna_end <-
                    ex$exon_cdna_start +
                    (
                        ex$end -
                        c0$start
                    )
            }


            mapped[[i]] <- data.table(
                transcript_id = tx,
                cds_cdna_start =
                    as.integer(
                        cdna_start
                    ),
                cds_cdna_end =
                    as.integer(
                        cdna_end
                    )
            )
        }


        cds_map_list[[tx]] <- rbindlist(
            mapped,
            fill = TRUE
        )
    }
}


cds_map <- rbindlist(
    cds_map_list,
    fill = TRUE,
    use.names = TRUE
)


transcripts_with_cds <- unique(
    cds_map$transcript_id
)


# -----------------------------------------------------------------------------
# 10. Exon-aware annotation of every transcript-level guide instance
# -----------------------------------------------------------------------------

guides[
    ,
    `:=`(
        is_genomically_contiguous_exonic =
            FALSE,

        exon_rank =
            NA_integer_,

        genomic_seqid =
            NA_character_,

        genomic_start =
            NA_integer_,

        genomic_end =
            NA_integer_,

        transcript_gene_strand =
            NA_character_,

        guide_genomic_strand =
            NA_character_,

        target_region_class =
            "exon_junction_or_unmapped"
    )
]


guide_split_indices <- split(
    seq_len(
        nrow(guides)
    ),
    guides$transcript_id
)


for (
    tx in names(guide_split_indices)
) {

    idx <- guide_split_indices[[tx]]


    tx_exons <- exons[
        transcript_id == tx
    ]


    if (nrow(tx_exons) == 0L) {
        next
    }


    data.table::setorder(
        tx_exons,
        exon_cdna_start
    )


    g_start <- as.integer(
        guides$transcript_start[idx]
    )


    g_end <- as.integer(
        guides$transcript_end[idx]
    )


    # Find the exon whose cDNA start is the last <= guide start.
    exon_idx <- findInterval(
        g_start,
        tx_exons$exon_cdna_start
    )


    valid_idx <-
        exon_idx >= 1L &
        exon_idx <= nrow(tx_exons)


    contained <- rep(
        FALSE,
        length(idx)
    )


    contained[valid_idx] <-
        g_end[valid_idx] <=
        tx_exons$exon_cdna_end[
            exon_idx[valid_idx]
        ]


    if (!any(contained)) {
        next
    }


    local_rows <- which(
        contained
    )


    global_rows <- idx[
        local_rows
    ]


    ex_selected <- tx_exons[
        exon_idx[
            local_rows
        ]
    ]


    guides[
        global_rows,
        is_genomically_contiguous_exonic :=
            TRUE
    ]


    guides[
        global_rows,
        exon_rank :=
            ex_selected$exon_rank
    ]


    guides[
        global_rows,
        genomic_seqid :=
            ex_selected$seqid
    ]


    guides[
        global_rows,
        transcript_gene_strand :=
            ex_selected$strand
    ]


    # Map the 23-nt transcript interval back to an ascending genomic interval.
    plus_gene <- ex_selected$strand == "+"


    mapped_start <- integer(
        length(global_rows)
    )


    mapped_end <- integer(
        length(global_rows)
    )


    mapped_start[plus_gene] <-
        ex_selected$start[plus_gene] +
        (
            g_start[local_rows][plus_gene] -
            ex_selected$exon_cdna_start[plus_gene]
        )


    mapped_end[plus_gene] <-
        ex_selected$start[plus_gene] +
        (
            g_end[local_rows][plus_gene] -
            ex_selected$exon_cdna_start[plus_gene]
        )


    minus_gene <- !plus_gene


    mapped_start[minus_gene] <-
        ex_selected$end[minus_gene] -
        (
            g_end[local_rows][minus_gene] -
            ex_selected$exon_cdna_start[minus_gene]
        )


    mapped_end[minus_gene] <-
        ex_selected$end[minus_gene] -
        (
            g_start[local_rows][minus_gene] -
            ex_selected$exon_cdna_start[minus_gene]
        )


    guides[
        global_rows,
        genomic_start :=
            as.integer(
                mapped_start
            )
    ]


    guides[
        global_rows,
        genomic_end :=
            as.integer(
                mapped_end
            )
    ]


    # The '+'/'-' from Step 05 is relative to transcript sequence orientation.
    # Convert to reference-genome strand orientation.
    guides[
        global_rows,
        guide_genomic_strand :=
            data.table::fifelse(
                ex_selected$strand == "+",
                strand,
                data.table::fifelse(
                    strand == "+",
                    "-",
                    "+"
                )
            )
    ]


    # Determine whether the full 23-nt target is contained in CDS.
    tx_cds <- cds_map[
        transcript_id == tx
    ]


    if (nrow(tx_cds) > 0L) {

        data.table::setorder(
            tx_cds,
            cds_cdna_start
        )


        cds_idx <- findInterval(
            g_start[local_rows],
            tx_cds$cds_cdna_start
        )


        cds_valid <-
            cds_idx >= 1L &
            cds_idx <= nrow(tx_cds)


        in_cds <- rep(
            FALSE,
            length(local_rows)
        )


        in_cds[cds_valid] <-
            g_end[local_rows][cds_valid] <=
            tx_cds$cds_cdna_end[
                cds_idx[cds_valid]
            ]


        guides[
            global_rows,
            target_region_class :=
                data.table::fifelse(
                    in_cds,
                    "CDS",
                    "UTR_or_other_exonic"
                )
        ]

    } else {

        guides[
            global_rows,
            target_region_class :=
                "noncoding_exon"
        ]
    }
}


# -----------------------------------------------------------------------------
# 11. Gene coding status and functional-region definition
# -----------------------------------------------------------------------------

gene_coding_status <- guides[
    ,
    .(
        gene_has_cds =
            any(
                transcript_id %in%
                    transcripts_with_cds
            )
    ),
    by = gene_id
]


guides <- merge(
    guides,
    gene_coding_status,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


guides[
    ,
    is_functional_region_target :=
        data.table::fifelse(
            gene_has_cds == TRUE,
            target_region_class == "CDS",
            target_region_class == "noncoding_exon"
        )
]


guides[
    ,
    passes_exon_aware_reference_quality :=
        is_genomically_contiguous_exonic &
        passes_reference_sequence_quality
]


guides[
    ,
    passes_functional_region_reference_quality :=
        passes_exon_aware_reference_quality &
        is_functional_region_target
]


# -----------------------------------------------------------------------------
# 12. Unique exon-aware guide table
# -----------------------------------------------------------------------------

valid_instances <- guides[
    is_genomically_contiguous_exonic == TRUE
]


if (nrow(valid_instances) > 0L) {

    unique_exon_guides <- valid_instances[
        ,
        .(
            n_transcripts_supporting_exon_contiguous_guide =
                uniqueN(transcript_id),

            any_reference_quality =
                any(
                    passes_reference_sequence_quality,
                    na.rm = TRUE
                ),

            any_functional_region =
                any(
                    is_functional_region_target,
                    na.rm = TRUE
                ),

            any_functional_region_reference_quality =
                any(
                    passes_functional_region_reference_quality,
                    na.rm = TRUE
                ),

            any_CDS =
                any(
                    target_region_class == "CDS",
                    na.rm = TRUE
                ),

            any_UTR_or_other_exonic =
                any(
                    target_region_class == "UTR_or_other_exonic",
                    na.rm = TRUE
                ),

            any_noncoding_exon =
                any(
                    target_region_class == "noncoding_exon",
                    na.rm = TRUE
                ),

            best_gc_balance_score =
                max(
                    gc_balance_score,
                    na.rm = TRUE
                ),

            best_gc_percent =
                guide_gc_percent[
                    which.min(
                        abs(
                            guide_gc_percent -
                            50
                        )
                    )
                ],

            representative_seqid =
                genomic_seqid[[1]],

            representative_start =
                genomic_start[[1]],

            representative_end =
                genomic_end[[1]],

            representative_genomic_strand =
                guide_genomic_strand[[1]]
        ),
        by = .(
            gene_id,
            protospacer_20nt,
            pam
        )
    ]


    unique_exon_guides[
        !is.finite(
            best_gc_balance_score
        ),
        best_gc_balance_score :=
            NA_real_
    ]

} else {

    unique_exon_guides <- data.table(
        gene_id = character(),
        protospacer_20nt = character(),
        pam = character(),
        n_transcripts_supporting_exon_contiguous_guide = integer(),
        any_reference_quality = logical(),
        any_functional_region = logical(),
        any_functional_region_reference_quality = logical(),
        any_CDS = logical(),
        any_UTR_or_other_exonic = logical(),
        any_noncoding_exon = logical(),
        best_gc_balance_score = numeric(),
        best_gc_percent = numeric(),
        representative_seqid = character(),
        representative_start = integer(),
        representative_end = integer(),
        representative_genomic_strand = character()
    )
}


# -----------------------------------------------------------------------------
# 13. Gene-level exon-aware metrics
# -----------------------------------------------------------------------------

guide_gene_summary <- guides[
    ,
    .(
        n_transcript_guide_instances =
            .N,

        n_exon_contiguous_guide_instances =
            sum(
                is_genomically_contiguous_exonic,
                na.rm = TRUE
            ),

        n_exon_junction_or_unmapped_guide_instances =
            sum(
                !is_genomically_contiguous_exonic,
                na.rm = TRUE
            ),

        n_reference_quality_exon_contiguous_instances =
            sum(
                passes_exon_aware_reference_quality,
                na.rm = TRUE
            ),

        n_functional_region_reference_quality_instances =
            sum(
                passes_functional_region_reference_quality,
                na.rm = TRUE
            ),

        n_transcripts_with_functional_region_reference_quality_guide =
            uniqueN(
                transcript_id[
                    passes_functional_region_reference_quality ==
                        TRUE
                ]
            )
    ),
    by = gene_id
]


unique_gene_summary <- unique_exon_guides[
    ,
    .(
        n_unique_exon_contiguous_ngg_guides =
            .N,

        n_unique_exon_contiguous_reference_quality_guides =
            sum(
                any_reference_quality,
                na.rm = TRUE
            ),

        n_unique_functional_region_reference_quality_guides =
            sum(
                any_functional_region_reference_quality,
                na.rm = TRUE
            ),

        n_unique_CDS_reference_quality_guides =
            sum(
                any_CDS &
                    any_reference_quality,
                na.rm = TRUE
            ),

        n_unique_UTR_reference_quality_guides =
            sum(
                any_UTR_or_other_exonic &
                    any_reference_quality,
                na.rm = TRUE
            ),

        n_unique_noncoding_exon_reference_quality_guides =
            sum(
                any_noncoding_exon &
                    any_reference_quality,
                na.rm = TRUE
            ),

        best_functional_region_gc_balance_score =
            if (
                any(
                    any_functional_region_reference_quality,
                    na.rm = TRUE
                )
            ) {

                max(
                    best_gc_balance_score[
                        any_functional_region_reference_quality ==
                            TRUE
                    ],
                    na.rm = TRUE
                )

            } else {

                NA_real_
            },

        best_functional_region_gc_percent =
            if (
                any(
                    any_functional_region_reference_quality,
                    na.rm = TRUE
                )
            ) {

                z <- best_gc_percent[
                    any_functional_region_reference_quality ==
                        TRUE
                ]

                z[
                    which.min(
                        abs(
                            z -
                            50
                        )
                    )
                ]

            } else {

                NA_real_
            }
    ),
    by = gene_id
]


# -----------------------------------------------------------------------------
# 14. Merge into Step 05 gene table
# -----------------------------------------------------------------------------

editability <- merge(
    step05,
    gene_coding_status,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


editability <- merge(
    editability,
    guide_gene_summary,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


editability <- merge(
    editability,
    unique_gene_summary,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


count_fields <- c(
    "n_transcript_guide_instances",
    "n_exon_contiguous_guide_instances",
    "n_exon_junction_or_unmapped_guide_instances",
    "n_reference_quality_exon_contiguous_instances",
    "n_functional_region_reference_quality_instances",
    "n_transcripts_with_functional_region_reference_quality_guide",
    "n_unique_exon_contiguous_ngg_guides",
    "n_unique_exon_contiguous_reference_quality_guides",
    "n_unique_functional_region_reference_quality_guides",
    "n_unique_CDS_reference_quality_guides",
    "n_unique_UTR_reference_quality_guides",
    "n_unique_noncoding_exon_reference_quality_guides"
)


for (
    col in count_fields
) {

    if (col %in% names(editability)) {

        editability[
            is.na(
                get(col)
            ),
            (col) :=
                0L
        ]
    }
}


editability[
    ,
    fraction_raw_guide_instances_exon_contiguous :=
        data.table::fifelse(
            n_transcript_guide_instances > 0L,
            n_exon_contiguous_guide_instances /
                n_transcript_guide_instances,
            NA_real_
        )
]


editability[
    ,
    fraction_raw_guide_instances_junction_or_unmapped :=
        data.table::fifelse(
            n_transcript_guide_instances > 0L,
            n_exon_junction_or_unmapped_guide_instances /
                n_transcript_guide_instances,
            NA_real_
        )
]


editability[
    ,
    fraction_transcripts_with_functional_region_reference_quality_guide :=
        data.table::fifelse(
            n_transcripts_scanned > 0L,
            n_transcripts_with_functional_region_reference_quality_guide /
                n_transcripts_scanned,
            NA_real_
        )
]


editability[
    ,
    has_functional_region_reference_quality_guide :=
        n_unique_functional_region_reference_quality_guides >=
        1L
]


editability[
    ,
    has_at_least_2_functional_region_reference_quality_guides :=
        n_unique_functional_region_reference_quality_guides >=
        MULTIPLEX_TARGET_2
]


editability[
    ,
    has_at_least_4_functional_region_reference_quality_guides :=
        n_unique_functional_region_reference_quality_guides >=
        MULTIPLEX_TARGET_4
]


# -----------------------------------------------------------------------------
# 15. Exon-aware reference editability score
# -----------------------------------------------------------------------------
#
# More conservative than Step 05:
#   40% functional-region guide abundance
#   25% multiplexability (saturates at 4)
#   20% transcript coverage by functional-region guides
#   15% best GC balance
#
# For coding genes, "functional region" = CDS.
# For genes without CDS, "functional region" = noncoding exon.

rank01 <- function(x) {

    x <- suppressWarnings(
        as.numeric(x)
    )


    out <- rep(
        NA_real_,
        length(x)
    )


    ok <- is.finite(
        x
    )


    n_ok <- sum(
        ok
    )


    if (n_ok == 0L) {
        return(out)
    }


    if (n_ok == 1L) {

        out[ok] <- 0.5

        return(out)
    }


    r <- rank(
        x[ok],
        ties.method = "average"
    )


    out[ok] <-
        (
            r -
            1
        ) /
        (
            n_ok -
            1
        )


    out
}


editability[
    ,
    exon_aware_abundance_component :=
        rank01(
            log1p(
                n_unique_functional_region_reference_quality_guides
            )
        )
]


editability[
    ,
    exon_aware_multiplex_component :=
        pmin(
            1,
            n_unique_functional_region_reference_quality_guides /
                MULTIPLEX_TARGET_4
        )
]


editability[
    ,
    exon_aware_transcript_coverage_component :=
        pmin(
            1,
            pmax(
                0,
                fraction_transcripts_with_functional_region_reference_quality_guide
            )
        )
]


editability[
    ,
    exon_aware_gc_balance_component :=
        pmin(
            1,
            pmax(
                0,
                best_functional_region_gc_balance_score
            )
        )
]


editability[
    ,
    crispr_exon_aware_reference_editability_score :=
        0.40 *
        exon_aware_abundance_component +
        0.25 *
        exon_aware_multiplex_component +
        0.20 *
        exon_aware_transcript_coverage_component +
        0.15 *
        exon_aware_gc_balance_component
]


editability[
    !has_functional_region_reference_quality_guide,
    crispr_exon_aware_reference_editability_score :=
        0
]


editability[
    ,
    crispr_exon_aware_reference_editability_rank :=
        data.table::frank(
            -crispr_exon_aware_reference_editability_score,
            ties.method = "first"
        )
]


editability[
    ,
    step05b_status :=
        data.table::fcase(
            n_unique_functional_region_reference_quality_guides >= 4L,
            "functional_region_multiplex_4plus",

            n_unique_functional_region_reference_quality_guides >= 2L,
            "functional_region_multiplex_2plus",

            n_unique_functional_region_reference_quality_guides >= 1L,
            "functional_region_single_plus",

            n_unique_exon_contiguous_reference_quality_guides >= 1L,
            "exonic_guides_only_no_functional_region_guide",

            default =
                "no_exon_contiguous_reference_quality_guide"
        )
]


# -----------------------------------------------------------------------------
# 16. Stable ordering
# -----------------------------------------------------------------------------

editability[
    ,
    preliminary_rank_missing :=
        is.na(preliminary_rank)
]


data.table::setorder(
    editability,
    preliminary_rank_missing,
    preliminary_rank,
    gene_id
)


editability[
    ,
    preliminary_rank_missing :=
        NULL
]


data.table::setorder(
    guides,
    gene_id,
    transcript_id,
    -passes_functional_region_reference_quality,
    -passes_exon_aware_reference_quality,
    transcript_start,
    strand
)


# -----------------------------------------------------------------------------
# 17. Benchmark table
# -----------------------------------------------------------------------------

if (
    "is_prespecified_benchmark" %in%
    names(editability)
) {

    benchmark_editability <- editability[
        is_prespecified_benchmark == TRUE
    ]

} else if (
    "selected_as_benchmark" %in%
    names(editability)
) {

    benchmark_editability <- editability[
        selected_as_benchmark == TRUE
    ]

} else {

    benchmark_editability <- editability[
        0
    ]
}


# -----------------------------------------------------------------------------
# 18. QC
# -----------------------------------------------------------------------------

n_genes <- nrow(
    editability
)


n_with_functional_guide <- editability[
    has_functional_region_reference_quality_guide == TRUE,
    .N
]


n_with_2 <- editability[
    has_at_least_2_functional_region_reference_quality_guides == TRUE,
    .N
]


n_with_4 <- editability[
    has_at_least_4_functional_region_reference_quality_guides == TRUE,
    .N
]


n_coding <- editability[
    gene_has_cds == TRUE,
    .N
]


n_noncoding <- editability[
    gene_has_cds == FALSE,
    .N
]


n_raw_instances <- nrow(
    guides
)


n_contiguous_instances <- guides[
    is_genomically_contiguous_exonic == TRUE,
    .N
]


n_junction_instances <- guides[
    is_genomically_contiguous_exonic == FALSE,
    .N
]


qc <- data.table(
    metric = c(
        "Step 05 shortlist genes",
        "Genes classified as having CDS",
        "Genes classified as noncoding/no CDS",
        "Raw transcript-level guide instances",
        "Exon-contiguous guide instances",
        "Junction-spanning or unmapped guide instances",
        "Percent raw guide instances exon-contiguous",
        "Genes with >=1 functional-region reference-quality guide",
        "Genes with >=2 functional-region reference-quality guides",
        "Genes with >=4 functional-region reference-quality guides",
        "Population conservation evaluated",
        "Genome-wide off-target specificity evaluated"
    ),

    value = c(
        n_genes,
        n_coding,
        n_noncoding,
        n_raw_instances,
        n_contiguous_instances,
        n_junction_instances,
        round(
            100 *
            n_contiguous_instances /
            n_raw_instances,
            3
        ),
        n_with_functional_guide,
        n_with_2,
        n_with_4,
        0,
        0
    )
)


# -----------------------------------------------------------------------------
# 19. Provenance
# -----------------------------------------------------------------------------

provenance <- data.table(
    component = c(
        "Guide instances",
        "Annotation",
        "Exon-junction control",
        "Functional-region definition",
        "Exon-aware score",
        "Population conservation",
        "Off-target specificity"
    ),

    specification = c(
        "Step 05 transcript-level 20-nt+NGG candidate instances.",

        paste0(
            "Ensembl Metazoa release 63 AgamP4 GFF3: ",
            basename(GFF_FILE),
            "."
        ),

        paste0(
            "A 23-nt protospacer+PAM is retained as genomically contiguous only ",
            "when its complete transcript-coordinate interval lies inside one ",
            "annotated exon."
        ),

        paste0(
            "For genes with CDS annotation, functional-region guides must lie ",
            "fully within an annotated CDS segment. For genes without CDS, ",
            "functional-region guides must lie fully within a noncoding exon."
        ),

        paste0(
            "40% functional-region guide abundance, 25% multiplexability, ",
            "20% transcript coverage, 15% best GC balance. This is an ",
            "explainable reference-editability score, not a probability of ",
            "successful editing."
        ),

        "PENDING_TARGETED_AG1000G.",

        paste0(
            "Not evaluated. The present filter does not establish genome-wide ",
            "off-target safety."
        )
    )
)


# -----------------------------------------------------------------------------
# 20. Write outputs
# -----------------------------------------------------------------------------

fwrite(
    guides,
    GUIDE_OUTPUT_FILE,
    na = "NA"
)


fwrite(
    editability,
    EDITABILITY_OUTPUT_FILE,
    na = "NA"
)


fwrite(
    benchmark_editability,
    BENCHMARK_OUTPUT_FILE,
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
    GUIDE_INSTANCE_FILE,
    STEP05_EDITABILITY_FILE,
    GFF_FILE,
    GUIDE_OUTPUT_FILE,
    EDITABILITY_OUTPUT_FILE,
    BENCHMARK_OUTPUT_FILE,
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
# 21. Integrity checks
# -----------------------------------------------------------------------------

stopifnot(
    uniqueN(editability$gene_id) ==
        nrow(editability)
)


stopifnot(
    nrow(editability) ==
        nrow(step05)
)


stopifnot(
    all(
        editability$
            n_unique_functional_region_reference_quality_guides >=
            0L
    )
)


stopifnot(
    all(
        guides[
            is_genomically_contiguous_exonic == TRUE,
            genomic_end -
                genomic_start +
                1L
        ] ==
            23L
    )
)


# -----------------------------------------------------------------------------
# 22. Console summary
# -----------------------------------------------------------------------------

cat(
    "\n",
    "============================================================\n",
    "MOSQEDIT-R STEP 05b COMPLETED SUCCESSFULLY\n",
    "Exon-aware SpCas9 reference-editability refinement\n",
    "============================================================\n",
    "Shortlisted genes:                             ",
    format(n_genes, big.mark = ","),
    "\n",
    "Genes with CDS annotation:                    ",
    format(n_coding, big.mark = ","),
    "\n",
    "Genes without CDS annotation:                 ",
    format(n_noncoding, big.mark = ","),
    "\n",
    "Raw transcript-level guide instances:         ",
    format(n_raw_instances, big.mark = ","),
    "\n",
    "Exon-contiguous guide instances:              ",
    format(n_contiguous_instances, big.mark = ","),
    "\n",
    "Junction-spanning/unmapped guide instances:    ",
    format(n_junction_instances, big.mark = ","),
    "\n",
    "Genes with >=1 functional-region guide:       ",
    format(n_with_functional_guide, big.mark = ","),
    "\n",
    "Genes with >=2 functional-region guides:      ",
    format(n_with_2, big.mark = ","),
    "\n",
    "Genes with >=4 functional-region guides:      ",
    format(n_with_4, big.mark = ","),
    "\n",
    "Population conservation evaluated:            NO\n",
    "Genome-wide off-target specificity evaluated: NO\n",
    "Python used:                                  NO\n",
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
    "\nTop 20 preliminary-ranked genes after exon-aware refinement:\n"
)


display_columns <- c(
    "preliminary_rank",
    "gene_id",
    "benchmark_name",
    "gene_has_cds",
    "n_unique_reference_quality_guides",
    "n_unique_exon_contiguous_reference_quality_guides",
    "n_unique_functional_region_reference_quality_guides",
    "fraction_raw_guide_instances_junction_or_unmapped",
    "crispr_reference_editability_score",
    "crispr_exon_aware_reference_editability_score",
    "crispr_exon_aware_reference_editability_rank",
    "step05b_status"
)


display_columns <- display_columns[
    display_columns %in%
        names(editability)
]


print(
    editability[
        !is.na(preliminary_rank)
    ][
        1:min(
            20L,
            .N
        ),
        ..display_columns
    ]
)


if (nrow(benchmark_editability) > 0L) {

    cat(
        "\nBenchmark exon-aware editability:\n"
    )


    benchmark_columns <- c(
        "gene_id",
        "benchmark_name",
        "preliminary_rank",
        "gene_has_cds",
        "n_unique_reference_quality_guides",
        "n_unique_exon_contiguous_reference_quality_guides",
        "n_unique_functional_region_reference_quality_guides",
        "best_functional_region_gc_percent",
        "crispr_exon_aware_reference_editability_score",
        "crispr_exon_aware_reference_editability_rank",
        "step05b_status",
        "population_conservation_status"
    )


    benchmark_columns <- benchmark_columns[
        benchmark_columns %in%
            names(benchmark_editability)
    ]


    print(
        benchmark_editability[
            ,
            ..benchmark_columns
        ]
    )
}


log_step(
    "05b",
    "Exon-aware CRISPR target refinement completed successfully"
)


# -----------------------------------------------------------------------------
# MosqEditR package compatibility aliases
# -----------------------------------------------------------------------------
# The frozen manuscript pipeline subsequently expects the authoritative
# exon-aware outputs under the Step-05 canonical filenames.
file.copy(
    GUIDE_OUTPUT_FILE,
    "data_processed/05_crispr_candidate_sites_exon_aware.csv",
    overwrite = TRUE
)
file.copy(
    EDITABILITY_OUTPUT_FILE,
    "data_processed/05_crispr_editability.csv",
    overwrite = TRUE
)


