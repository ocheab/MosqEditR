# =============================================================================
# 05_crispr_editability_TARGETED_R.R
#
# MosqEdit-R Manuscript 1
#
# STEP 05
# Data-light, targeted SpCas9 reference-sequence editability analysis for the
# Step 04 shortlist.
#
# PURE R.
#
# PURPOSE
# -------
# Scan the transcript sequences of the Step 04 shortlisted genes for canonical
# Streptococcus pyogenes Cas9 (SpCas9) 20-nt protospacers followed by NGG PAMs.
#
# This step intentionally evaluates REFERENCE-SEQUENCE EDITABILITY only.
# It does NOT claim:
#   - population conservation,
#   - Ag1000G DRA frequency,
#   - genome-wide off-target specificity,
#   - experimentally measured cutting efficiency.
#
# Those layers remain independent validation tasks.
#
# The design mirrors the reference-target component of Schmidt et al. (2020):
# potential SpCas9 targets are 20-nt protospacers adjacent to an NGG PAM and
# guide GC content of 30-70% is retained as a Schmidt-compatible sequence
# filter. We additionally flag 4+ nt homopolymers and TTTT motifs because these
# are common guide-design liabilities. Genome-wide off-target filtering is NOT
# reproduced here because that requires a full-genome specificity search.
#
# DATA-LIGHT STRATEGY
# -------------------
# Instead of downloading Ag1000G chromosome VCFs or the complete genomic FASTA,
# download only release-pinned AgamP4 transcript FASTA resources:
#
#   Ensembl Metazoa release 63
#   Anopheles gambiae PEST / AgamP4
#   cDNA FASTA
#   ncRNA FASTA
#
# Only transcripts belonging to the Step 04 shortlist are retained in memory.
#
# INPUT
# -----
# data_processed/04_targeted_ag3_shortlist.csv
#
# OUTPUTS
# -------
# data_raw/ensembl_release63_agamp4/
#   release-63 cDNA and ncRNA FASTA files
#
# data_processed/
#   05_crispr_candidate_guides_transcript_level.csv
#   05_crispr_candidate_guides_unique.csv
#   05_crispr_editability.csv
#   05_crispr_editability_benchmarks.csv
#   05_crispr_editability_qc.csv
#   05_crispr_editability_provenance.csv
#
# logs/
#   05_crispr_editability_sessionInfo.txt
#   05_crispr_editability_checksums.tsv
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Helpers and packages
# -----------------------------------------------------------------------------

source("R/helpers.R")


required_cran <- c(
    "data.table",
    "curl"
)

required_bioc <- c(
    "Biostrings"
)


missing_cran <- required_cran[
    !vapply(
        required_cran,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]


missing_bioc <- required_bioc[
    !vapply(
        required_bioc,
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


if (length(missing_bioc) > 0L) {

    stop(
        paste0(
            "Missing Bioconductor package(s): ",
            paste(missing_bioc, collapse = ", "),
            "\nInstall with:\n",
            "if (!requireNamespace(\"BiocManager\", quietly=TRUE)) ",
            "install.packages(\"BiocManager\")\n",
            "BiocManager::install(c(",
            paste0('"', missing_bioc, '"', collapse = ", "),
            "))"
        ),
        call. = FALSE
    )
}


suppressPackageStartupMessages({

    library(data.table)
    library(curl)
    library(Biostrings)

})


log_step(
    "05",
    "Starting targeted reference-sequence SpCas9 editability analysis"
)


# -----------------------------------------------------------------------------
# 2. Configuration
# -----------------------------------------------------------------------------

SHORTLIST_FILE <-
    "data_processed/04_targeted_ag3_shortlist.csv"


RAW_DIR <-
    "data_raw/ensembl_release63_agamp4"


TRANSCRIPT_GUIDE_FILE <-
    "data_processed/05_crispr_candidate_guides_transcript_level.csv"


UNIQUE_GUIDE_FILE <-
    "data_processed/05_crispr_candidate_guides_unique.csv"


EDITABILITY_FILE <-
    "data_processed/05_crispr_editability.csv"


BENCHMARK_FILE <-
    "data_processed/05_crispr_editability_benchmarks.csv"


QC_FILE <-
    "data_processed/05_crispr_editability_qc.csv"


PROVENANCE_FILE <-
    "data_processed/05_crispr_editability_provenance.csv"


SESSION_INFO_FILE <-
    "logs/05_crispr_editability_sessionInfo.txt"


CHECKSUM_FILE <-
    "logs/05_crispr_editability_checksums.tsv"


ENSEMBL_GENOMES_RELEASE <- "63"

ENSEMBL_ASSEMBLY <- "AgamP4"

ENSEMBL_SPECIES <- "anopheles_gambiae"


FASTA_ROOT <- paste0(
    "https://ftp.ebi.ac.uk/ensemblgenomes/pub/metazoa/release-",
    ENSEMBL_GENOMES_RELEASE,
    "/fasta/",
    ENSEMBL_SPECIES
)


CDNA_DIR_URL <- paste0(
    FASTA_ROOT,
    "/cdna/"
)


NCRNA_DIR_URL <- paste0(
    FASTA_ROOT,
    "/ncrna/"
)


# SpCas9 reference-sequence rules.
GUIDE_LENGTH <- 20L
PAM_LENGTH <- 3L

MIN_GUIDE_GC <- 30
MAX_GUIDE_GC <- 70

HOMOPOLYMER_RUN <- 4L

# Reference multiplexability thresholds.
MULTIPLEX_TARGET_2 <- 2L
MULTIPLEX_TARGET_4 <- 4L


dir.create(
    RAW_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

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


rank01 <- function(x) {

    x <- suppressWarnings(
        as.numeric(x)
    )

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

    x <- suppressWarnings(
        as.numeric(x)
    )

    ok <- is.finite(x)

    x[ok] <- pmin(
        1,
        pmax(
            0,
            x[ok]
        )
    )

    x
}


# -----------------------------------------------------------------------------
# 4. Discover release-pinned Ensembl FASTA files dynamically
# -----------------------------------------------------------------------------

discover_fasta_file <- function(
        directory_url,
        preferred_regex,
        label
) {

    message(
        "Discovering ",
        label,
        " FASTA from:\n",
        directory_url
    )


    response <- tryCatch(

        curl::curl_fetch_memory(
            directory_url
        ),

        error = function(e) {

            stop(
                paste0(
                    "Could not read Ensembl FASTA directory for ",
                    label,
                    ".\n\nURL:\n",
                    directory_url,
                    "\n\nOriginal error:\n",
                    conditionMessage(e)
                ),
                call. = FALSE
            )
        }
    )


    if (
        is.null(response$status_code) ||
        response$status_code < 200L ||
        response$status_code >= 300L
    ) {

        stop(
            paste0(
                "Ensembl directory request returned HTTP ",
                response$status_code,
                " for ",
                directory_url
            ),
            call. = FALSE
        )
    }


    html <- rawToChar(
        response$content
    )


    m <- gregexpr(
        'href="[^"]+\\.fa\\.gz"',
        html,
        perl = TRUE
    )[[1]]


    if (
        length(m) == 1L &&
        m[[1]] == -1L
    ) {

        stop(
            paste0(
                "No .fa.gz links were found in Ensembl directory:\n",
                directory_url
            ),
            call. = FALSE
        )
    }


    hits <- regmatches(
        html,
        list(m)
    )[[1]]


    files <- sub(
        '^href="',
        "",
        hits
    )


    files <- sub(
        '"$',
        "",
        files
    )


    files <- basename(
        files
    )


    files <- unique(
        files
    )


    preferred <- files[
        grepl(
            preferred_regex,
            files,
            ignore.case = TRUE,
            perl = TRUE
        )
    ]


    if (length(preferred) == 0L) {

        stop(
            paste0(
                "Could not identify the expected ",
                label,
                " FASTA file.\n\nAvailable .fa.gz files:\n",
                paste(
                    files,
                    collapse = "\n"
                )
            ),
            call. = FALSE
        )
    }


    # Prefer the shortest matching filename to avoid ancillary subsets.
    preferred <- preferred[
        order(
            nchar(preferred),
            preferred
        )
    ]


    selected <- preferred[[1]]


    list(
        filename = selected,
        url = paste0(
            directory_url,
            selected
        )
    )
}


# Expected Ensembl naming conventions:
#   *.cdna.all.fa.gz
#   *.ncrna.fa.gz
#
# Dynamic discovery prevents a silent failure if minor filename details differ.
cdna_resource <- discover_fasta_file(
    directory_url = CDNA_DIR_URL,
    preferred_regex = "cdna\\.all\\.fa\\.gz$",
    label = "cDNA"
)


ncrna_resource <- discover_fasta_file(
    directory_url = NCRNA_DIR_URL,
    preferred_regex = "ncrna\\.fa\\.gz$",
    label = "ncRNA"
)


cat(
    "\nSelected Ensembl sequence resources:\n",
    "cDNA:  ",
    cdna_resource$url,
    "\n",
    "ncRNA: ",
    ncrna_resource$url,
    "\n",
    sep = ""
)


# -----------------------------------------------------------------------------
# 5. Download small release-pinned transcript FASTA resources
# -----------------------------------------------------------------------------

download_if_needed <- function(
        url,
        destfile,
        minimum_bytes = 1000L
) {

    if (
        file.exists(destfile) &&
        !is.na(file.info(destfile)$size) &&
        file.info(destfile)$size >= minimum_bytes
    ) {

        message(
            "Using existing local file: ",
            destfile
        )

        return(
            invisible(destfile)
        )
    }


    tmp <- paste0(
        destfile,
        ".part"
    )


    if (file.exists(tmp)) {
        unlink(tmp)
    }


    message(
        "Downloading: ",
        url
    )


    tryCatch(

        curl::curl_download(
            url = url,
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
                    "Unable to download Ensembl transcript FASTA.\n\n",
                    "URL:\n",
                    url,
                    "\n\nOriginal error:\n",
                    conditionMessage(e)
                ),
                call. = FALSE
            )
        }
    )


    if (
        !file.exists(tmp) ||
        is.na(file.info(tmp)$size) ||
        file.info(tmp)$size < minimum_bytes
    ) {

        if (file.exists(tmp)) {
            unlink(tmp)
        }

        stop(
            paste0(
                "Downloaded FASTA is unexpectedly small:\n",
                url
            ),
            call. = FALSE
        )
    }


    if (file.exists(destfile)) {
        unlink(destfile)
    }


    if (!file.rename(tmp, destfile)) {

        stop(
            paste0(
                "Could not move downloaded FASTA into place:\n",
                destfile
            ),
            call. = FALSE
        )
    }


    invisible(destfile)
}


CDNA_FASTA <- file.path(
    RAW_DIR,
    cdna_resource$filename
)


NCRNA_FASTA <- file.path(
    RAW_DIR,
    ncrna_resource$filename
)


download_if_needed(
    cdna_resource$url,
    CDNA_FASTA
)


download_if_needed(
    ncrna_resource$url,
    NCRNA_FASTA
)


cat(
    "\nDownloaded/available FASTA sizes:\n",
    "cDNA:  ",
    round(
        file.info(CDNA_FASTA)$size / 1024^2,
        2
    ),
    " MB\n",
    "ncRNA: ",
    round(
        file.info(NCRNA_FASTA)$size / 1024^2,
        2
    ),
    " MB\n",
    sep = ""
)


# -----------------------------------------------------------------------------
# 6. Load Step 04 shortlist
# -----------------------------------------------------------------------------

assert_file(
    SHORTLIST_FILE
)


shortlist <- fread(
    SHORTLIST_FILE
)


if (
    !"gene_id" %in%
    names(shortlist)
) {

    stop(
        "Step 04 shortlist lacks gene_id.",
        call. = FALSE
    )
}


if (
    uniqueN(shortlist$gene_id) !=
    nrow(shortlist)
) {

    stop(
        "Step 04 shortlist contains duplicate gene_id values.",
        call. = FALSE
    )
}


shortlist_gene_ids <- shortlist[
    !is.na(gene_id) &
        nzchar(gene_id),
    gene_id
]


log_step(
    "05",
    paste(
        "Loaded",
        format(length(shortlist_gene_ids), big.mark = ","),
        "shortlisted genes from Step 04"
    )
)


# -----------------------------------------------------------------------------
# 7. FASTA header parser
# -----------------------------------------------------------------------------

extract_gene_id <- function(header) {

    header <- as.character(
        header
    )


    out <- rep(
        NA_character_,
        length(header)
    )


    has_gene_tag <- grepl(
        "gene:AGAP[0-9]+",
        header,
        perl = TRUE
    )


    out[has_gene_tag] <- sub(
        ".*gene:(AGAP[0-9]+).*",
        "\\1",
        header[has_gene_tag],
        perl = TRUE
    )


    # Fallback: transcript IDs normally begin with the stable AGAP gene ID.
    missing <- is.na(out)


    if (any(missing)) {

        first_token <- sub(
            "\\s.*$",
            "",
            header[missing]
        )


        has_agap <- grepl(
            "^AGAP[0-9]+",
            first_token
        )


        temp <- rep(
            NA_character_,
            length(first_token)
        )


        temp[has_agap] <- sub(
            "^(AGAP[0-9]+).*$",
            "\\1",
            first_token[has_agap]
        )


        out[missing] <- temp
    }


    out
}


extract_transcript_id <- function(header) {

    sub(
        "\\s.*$",
        "",
        as.character(header)
    )
}


# -----------------------------------------------------------------------------
# 8. Read FASTA and retain only shortlist transcripts
# -----------------------------------------------------------------------------

read_target_transcripts <- function(
        fasta_file,
        source_label,
        target_gene_ids
) {

    message(
        "Reading ",
        source_label,
        " FASTA..."
    )


    dna <- Biostrings::readDNAStringSet(
        fasta_file,
        format = "fasta",
        use.names = TRUE
    )


    hdr <- names(
        dna
    )


    gene_id <- extract_gene_id(
        hdr
    )


    transcript_id <- extract_transcript_id(
        hdr
    )


    keep <- !is.na(gene_id) &
        gene_id %in%
        target_gene_ids


    if (!any(keep)) {

        return(
            data.table(
                gene_id = character(),
                transcript_id = character(),
                sequence_source = character(),
                sequence = character()
            )
        )
    }


    selected <- dna[
        keep
    ]


    dt <- data.table(
        gene_id = gene_id[keep],
        transcript_id = transcript_id[keep],
        sequence_source = source_label,
        sequence = toupper(
            as.character(
                selected
            )
        )
    )


    dt[
        ,
        transcript_length_bp :=
            nchar(sequence)
    ]


    dt
}


cdna_transcripts <- read_target_transcripts(
    fasta_file = CDNA_FASTA,
    source_label = "Ensembl_cDNA",
    target_gene_ids = shortlist_gene_ids
)


ncrna_transcripts <- read_target_transcripts(
    fasta_file = NCRNA_FASTA,
    source_label = "Ensembl_ncRNA",
    target_gene_ids = shortlist_gene_ids
)


transcripts <- rbindlist(
    list(
        cdna_transcripts,
        ncrna_transcripts
    ),
    fill = TRUE,
    use.names = TRUE
)


if (nrow(transcripts) == 0L) {

    stop(
        "No Step 04 shortlist transcripts were recovered from the Ensembl FASTA files.",
        call. = FALSE
    )
}


# If the same transcript is present in both sequence resources, keep one copy.
# Prefer cDNA simply for deterministic ordering; sequence equality is checked.
setorder(
    transcripts,
    transcript_id,
    sequence_source
)


dup_ids <- transcripts[
    duplicated(transcript_id) |
        duplicated(
            transcript_id,
            fromLast = TRUE
        ),
    unique(transcript_id)
]


if (length(dup_ids) > 0L) {

    dup_check <- transcripts[
        transcript_id %in%
            dup_ids,
        .(
            n_unique_sequences =
                uniqueN(sequence)
        ),
        by = transcript_id
    ]


    bad_dup <- dup_check[
        n_unique_sequences > 1L,
        transcript_id
    ]


    if (length(bad_dup) > 0L) {

        stop(
            paste0(
                "Duplicate transcript IDs with non-identical sequences were found. ",
                "Example(s): ",
                paste(
                    head(bad_dup, 10L),
                    collapse = ", "
                )
            ),
            call. = FALSE
        )
    }
}


transcripts <- unique(
    transcripts,
    by = "transcript_id"
)


genes_with_sequence <- unique(
    transcripts$gene_id
)


genes_without_sequence <- setdiff(
    shortlist_gene_ids,
    genes_with_sequence
)


log_step(
    "05",
    paste(
        "Recovered",
        format(nrow(transcripts), big.mark = ","),
        "unique transcripts for",
        format(length(genes_with_sequence), big.mark = ","),
        "of",
        format(length(shortlist_gene_ids), big.mark = ","),
        "shortlisted genes"
    )
)


if (length(genes_without_sequence) > 0L) {

    warning(
        paste0(
            length(genes_without_sequence),
            " shortlisted gene(s) had no cDNA/ncRNA sequence in release 63. ",
            "They will be retained with missing Step 05 editability metrics."
        ),
        call. = FALSE
    )
}


# -----------------------------------------------------------------------------
# 9. SpCas9 target scanner
# -----------------------------------------------------------------------------

gc_percent <- function(s) {

    s <- toupper(
        s
    )

    chars <- strsplit(
        s,
        "",
        fixed = TRUE
    )[[1]]

    if (length(chars) == 0L) {
        return(NA_real_)
    }

    100 *
        sum(
            chars %in%
                c(
                    "G",
                    "C"
                )
        ) /
        length(chars)
}


scan_one_orientation <- function(
        sequence,
        strand,
        transcript_id,
        gene_id,
        sequence_source
) {

    sequence <- toupper(
        sequence
    )


    n <- nchar(
        sequence
    )


    if (
        !is.finite(n) ||
        n <
        (GUIDE_LENGTH + PAM_LENGTH)
    ) {

        return(
            NULL
        )
    }


    # Scan the supplied orientation for:
    # 20 nt protospacer + NGG PAM.
    #
    # A look-ahead is used so overlapping candidate sites are retained.
    starts <- gregexpr(
        "(?=([ACGT]{20}[ACGT]GG))",
        sequence,
        perl = TRUE
    )[[1]]


    if (
        length(starts) == 1L &&
        starts[[1]] == -1L
    ) {

        return(
            NULL
        )
    }


    target23 <- substring(
        sequence,
        first = starts,
        last = starts + 22L
    )


    protospacer <- substr(
        target23,
        1L,
        GUIDE_LENGTH
    )


    pam <- substr(
        target23,
        GUIDE_LENGTH + 1L,
        GUIDE_LENGTH + PAM_LENGTH
    )


    if (strand == "+") {

        transcript_start <- starts
        transcript_end <- starts + 22L

    } else {

        # Coordinates are converted back to the original transcript orientation.
        transcript_start <-
            n -
            starts -
            21L

        transcript_end <-
            n -
            starts +
            1L
    }


    gc <- vapply(
        protospacer,
        gc_percent,
        numeric(1)
    )


    has_homopolymer4 <- grepl(
        paste0(
            "A{",
            HOMOPOLYMER_RUN,
            "}|C{",
            HOMOPOLYMER_RUN,
            "}|G{",
            HOMOPOLYMER_RUN,
            "}|T{",
            HOMOPOLYMER_RUN,
            "}"
        ),
        protospacer,
        perl = TRUE
    )


    has_tttt <- grepl(
        "TTTT",
        protospacer,
        fixed = TRUE
    )


    passes_gc_30_70 <-
        is.finite(gc) &
        gc >= MIN_GUIDE_GC &
        gc <= MAX_GUIDE_GC


    # Reference sequence-quality filter.
    #
    # IMPORTANT:
    # This is NOT a genome-wide off-target specificity filter and NOT a
    # population-conservation/DRA filter.
    passes_reference_sequence_quality <-
        passes_gc_30_70 &
        !has_homopolymer4 &
        !has_tttt


    # Simple sequence-only balance score:
    # 1.0 at 50% GC, declining linearly to 0 at <=30% or >=70%.
    gc_balance_score <- 1 -
        abs(
            gc -
            50
        ) /
        20


    gc_balance_score <- pmin(
        1,
        pmax(
            0,
            gc_balance_score
        )
    )


    data.table(
        gene_id = gene_id,
        transcript_id = transcript_id,
        sequence_source = sequence_source,
        strand = strand,
        transcript_start = as.integer(transcript_start),
        transcript_end = as.integer(transcript_end),
        protospacer_20nt = protospacer,
        pam = pam,
        target_23nt = target23,
        guide_gc_percent = gc,
        gc_balance_score = gc_balance_score,
        has_homopolymer_4plus = has_homopolymer4,
        has_TTTT = has_tttt,
        passes_gc_30_70 = passes_gc_30_70,
        passes_reference_sequence_quality =
            passes_reference_sequence_quality
    )
}


scan_transcript <- function(
        gene_id,
        transcript_id,
        sequence_source,
        sequence
) {

    sequence <- toupper(
        sequence
    )


    plus <- scan_one_orientation(
        sequence = sequence,
        strand = "+",
        transcript_id = transcript_id,
        gene_id = gene_id,
        sequence_source = sequence_source
    )


    rc <- as.character(
        Biostrings::reverseComplement(
            Biostrings::DNAString(
                sequence
            )
        )
    )


    minus <- scan_one_orientation(
        sequence = rc,
        strand = "-",
        transcript_id = transcript_id,
        gene_id = gene_id,
        sequence_source = sequence_source
    )


    rbindlist(
        list(
            plus,
            minus
        ),
        fill = TRUE,
        use.names = TRUE
    )
}


# -----------------------------------------------------------------------------
# 10. Scan all shortlisted transcripts
# -----------------------------------------------------------------------------

guide_list <- vector(
    "list",
    nrow(transcripts)
)


for (
    i in seq_len(
        nrow(transcripts)
    )
) {

    if (
        i == 1L ||
        i %% 100L == 0L ||
        i == nrow(transcripts)
    ) {

        log_step(
            "05",
            paste(
                "Scanning transcript",
                format(i, big.mark = ","),
                "of",
                format(nrow(transcripts), big.mark = ",")
            )
        )
    }


    guide_list[[i]] <- scan_transcript(
        gene_id = transcripts$gene_id[[i]],
        transcript_id = transcripts$transcript_id[[i]],
        sequence_source = transcripts$sequence_source[[i]],
        sequence = transcripts$sequence[[i]]
    )
}


guides_transcript <- rbindlist(
    guide_list,
    fill = TRUE,
    use.names = TRUE
)


if (nrow(guides_transcript) == 0L) {

    warning(
        "No canonical 20-nt + NGG SpCas9 targets were detected in the shortlisted transcripts.",
        call. = FALSE
    )

    guides_transcript <- data.table(
        gene_id = character(),
        transcript_id = character(),
        sequence_source = character(),
        strand = character(),
        transcript_start = integer(),
        transcript_end = integer(),
        protospacer_20nt = character(),
        pam = character(),
        target_23nt = character(),
        guide_gc_percent = numeric(),
        gc_balance_score = numeric(),
        has_homopolymer_4plus = logical(),
        has_TTTT = logical(),
        passes_gc_30_70 = logical(),
        passes_reference_sequence_quality = logical()
    )
}


# -----------------------------------------------------------------------------
# 11. Mark guide sharing across transcript isoforms
# -----------------------------------------------------------------------------

if (nrow(guides_transcript) > 0L) {

    guide_isoform_support <- guides_transcript[
        ,
        .(
            n_transcripts_supporting_guide =
                uniqueN(transcript_id)
        ),
        by = .(
            gene_id,
            protospacer_20nt,
            pam
        )
    ]


    guides_transcript <- merge(
        guides_transcript,
        guide_isoform_support,
        by = c(
            "gene_id",
            "protospacer_20nt",
            "pam"
        ),
        all.x = TRUE,
        sort = FALSE
    )


    # A gene-level unique-guide table avoids counting an identical guide once
    # for every transcript isoform carrying it.
    setorder(
        guides_transcript,
        gene_id,
        protospacer_20nt,
        pam,
        -passes_reference_sequence_quality,
        -n_transcripts_supporting_guide,
        transcript_id,
        strand,
        transcript_start
    )


    guides_unique <- guides_transcript[
        ,
        .SD[1L],
        by = .(
            gene_id,
            protospacer_20nt,
            pam
        )
    ]


    guides_unique[
        ,
        guide_sequence_key :=
            paste0(
                protospacer_20nt,
                pam
            )
    ]


    # Ranking among reference-design-ready guides only.
    guides_unique[
        ,
        reference_guide_quality_score :=
            data.table::fifelse(
                passes_reference_sequence_quality,
                gc_balance_score,
                NA_real_
            )
    ]


    data.table::setorder(
        guides_unique,
        gene_id,
        -passes_reference_sequence_quality,
        -reference_guide_quality_score,
        -n_transcripts_supporting_guide,
        protospacer_20nt,
        pam
    )


    guides_unique[
        ,
        reference_guide_rank_within_gene :=
            seq_len(.N),
        by = gene_id
    ]

} else {

    guides_unique <- copy(
        guides_transcript
    )

    guides_unique[
        ,
        `:=`(
            n_transcripts_supporting_guide = integer(),
            guide_sequence_key = character(),
            reference_guide_quality_score = numeric(),
            reference_guide_rank_within_gene = integer()
        )
    ]
}


# -----------------------------------------------------------------------------
# 12. Transcript-level summary
# -----------------------------------------------------------------------------

transcript_summary <- transcripts[
    ,
    .(
        n_transcripts_scanned =
            uniqueN(transcript_id),

        total_transcript_bp_scanned =
            sum(
                transcript_length_bp,
                na.rm = TRUE
            ),

        max_transcript_length_bp =
            max(
                transcript_length_bp,
                na.rm = TRUE
            ),

        n_cdna_transcripts =
            uniqueN(
                transcript_id[
                    sequence_source ==
                        "Ensembl_cDNA"
                ]
            ),

        n_ncrna_transcripts =
            uniqueN(
                transcript_id[
                    sequence_source ==
                        "Ensembl_ncRNA"
                ]
            )
    ),
    by = gene_id
]


if (nrow(guides_transcript) > 0L) {

    transcript_guide_summary <- guides_transcript[
        ,
        .(
            n_raw_ngg_guide_instances =
                .N,

            n_gc30_70_guide_instances =
                sum(
                    passes_gc_30_70,
                    na.rm = TRUE
                ),

            n_reference_quality_guide_instances =
                sum(
                    passes_reference_sequence_quality,
                    na.rm = TRUE
                ),

            n_transcripts_with_any_ngg =
                uniqueN(
                    transcript_id
                ),

            n_transcripts_with_reference_quality_guide =
                uniqueN(
                    transcript_id[
                        passes_reference_sequence_quality ==
                            TRUE
                    ]
                )
        ),
        by = gene_id
    ]

} else {

    transcript_guide_summary <- data.table(
        gene_id = character(),
        n_raw_ngg_guide_instances = integer(),
        n_gc30_70_guide_instances = integer(),
        n_reference_quality_guide_instances = integer(),
        n_transcripts_with_any_ngg = integer(),
        n_transcripts_with_reference_quality_guide = integer()
    )
}


# -----------------------------------------------------------------------------
# 13. Unique-guide gene-level summary
# -----------------------------------------------------------------------------

if (nrow(guides_unique) > 0L) {

    unique_guide_summary <- guides_unique[
        ,
        .(
            n_unique_ngg_guides =
                .N,

            n_unique_gc30_70_guides =
                sum(
                    passes_gc_30_70,
                    na.rm = TRUE
                ),

            n_unique_reference_quality_guides =
                sum(
                    passes_reference_sequence_quality,
                    na.rm = TRUE
                ),

            best_reference_guide_gc_balance_score =
                if (
                    any(
                        passes_reference_sequence_quality ==
                            TRUE
                    )
                ) {
                    max(
                        gc_balance_score[
                            passes_reference_sequence_quality ==
                                TRUE
                        ],
                        na.rm = TRUE
                    )
                } else {
                    NA_real_
                },

            best_reference_guide_gc_percent =
                if (
                    any(
                        passes_reference_sequence_quality ==
                            TRUE
                    )
                ) {

                    tmp <- guide_gc_percent[
                        passes_reference_sequence_quality ==
                            TRUE
                    ]

                    tmp[
                        which.min(
                            abs(
                                tmp -
                                50
                            )
                        )
                    ]

                } else {

                    NA_real_
                },

            median_reference_guide_gc_percent =
                if (
                    any(
                        passes_reference_sequence_quality ==
                            TRUE
                    )
                ) {
                    median(
                        guide_gc_percent[
                            passes_reference_sequence_quality ==
                                TRUE
                        ],
                        na.rm = TRUE
                    )
                } else {
                    NA_real_
                },

            max_transcript_isoform_support =
                if (.N > 0L) {
                    max(
                        n_transcripts_supporting_guide,
                        na.rm = TRUE
                    )
                } else {
                    NA_integer_
                }
        ),
        by = gene_id
    ]

} else {

    unique_guide_summary <- data.table(
        gene_id = character(),
        n_unique_ngg_guides = integer(),
        n_unique_gc30_70_guides = integer(),
        n_unique_reference_quality_guides = integer(),
        best_reference_guide_gc_balance_score = numeric(),
        best_reference_guide_gc_percent = numeric(),
        median_reference_guide_gc_percent = numeric(),
        max_transcript_isoform_support = integer()
    )
}


# -----------------------------------------------------------------------------
# 14. Assemble complete Step 05 gene-level table
# -----------------------------------------------------------------------------

editability <- merge(
    shortlist,
    transcript_summary,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


editability <- merge(
    editability,
    transcript_guide_summary,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


editability <- merge(
    editability,
    unique_guide_summary,
    by = "gene_id",
    all.x = TRUE,
    sort = FALSE
)


has_sequence <- editability$gene_id %in%
    genes_with_sequence


editability[
    ,
    reference_sequence_available :=
        has_sequence
]


count_columns <- c(
    "n_transcripts_scanned",
    "total_transcript_bp_scanned",
    "max_transcript_length_bp",
    "n_cdna_transcripts",
    "n_ncrna_transcripts",
    "n_raw_ngg_guide_instances",
    "n_gc30_70_guide_instances",
    "n_reference_quality_guide_instances",
    "n_transcripts_with_any_ngg",
    "n_transcripts_with_reference_quality_guide",
    "n_unique_ngg_guides",
    "n_unique_gc30_70_guides",
    "n_unique_reference_quality_guides"
)


for (
    col in count_columns
) {

    if (!col %in% names(editability)) {
        next
    }


    editability[
        reference_sequence_available == TRUE &
            is.na(
                get(col)
            ),
        (col) := 0
    ]
}


editability[
    ,
    fraction_transcripts_with_reference_quality_guide :=
        data.table::fifelse(
            reference_sequence_available &
                n_transcripts_scanned > 0L,
            n_transcripts_with_reference_quality_guide /
                n_transcripts_scanned,
            NA_real_
        )
]


editability[
    ,
    reference_quality_guide_instance_density_per_kb :=
        data.table::fifelse(
            reference_sequence_available &
                total_transcript_bp_scanned > 0,
            1000 *
                n_reference_quality_guide_instances /
                total_transcript_bp_scanned,
            NA_real_
        )
]


editability[
    ,
    has_reference_quality_guide :=
        data.table::fifelse(
            reference_sequence_available,
            n_unique_reference_quality_guides >= 1L,
            NA
        )
]


editability[
    ,
    has_at_least_2_reference_quality_guides :=
        data.table::fifelse(
            reference_sequence_available,
            n_unique_reference_quality_guides >=
                MULTIPLEX_TARGET_2,
            NA
        )
]


editability[
    ,
    has_at_least_4_reference_quality_guides :=
        data.table::fifelse(
            reference_sequence_available,
            n_unique_reference_quality_guides >=
                MULTIPLEX_TARGET_4,
            NA
        )
]


# -----------------------------------------------------------------------------
# 15. Explainable reference-editability score
# -----------------------------------------------------------------------------
#
# This score ranks only the Step 04 shortlist and is NOT a probability of
# successful gene editing.
#
# Components:
#   35% unique reference-quality guide abundance
#   25% reference-quality guide-instance density
#   20% multiplexability (saturates at 4 unique guides)
#   10% transcript coverage
#   10% best GC-balance score
#
# Off-target specificity and population conservation are deliberately absent.

editability[
    ,
    editability_abundance_component :=
        rank01(
            log1p(
                n_unique_reference_quality_guides
            )
        )
]


editability[
    ,
    editability_density_component :=
        rank01(
            log1p(
                reference_quality_guide_instance_density_per_kb
            )
        )
]


editability[
    ,
    editability_multiplex_component :=
        data.table::fifelse(
            reference_sequence_available,
            pmin(
                1,
                n_unique_reference_quality_guides /
                    MULTIPLEX_TARGET_4
            ),
            NA_real_
        )
]


editability[
    ,
    editability_transcript_coverage_component :=
        clip01(
            fraction_transcripts_with_reference_quality_guide
        )
]


editability[
    ,
    editability_gc_balance_component :=
        clip01(
            best_reference_guide_gc_balance_score
        )
]


editability_component_columns <- c(
    "editability_abundance_component",
    "editability_density_component",
    "editability_multiplex_component",
    "editability_transcript_coverage_component",
    "editability_gc_balance_component"
)


editability_component_weights <- c(
    0.35,
    0.25,
    0.20,
    0.10,
    0.10
)


component_matrix <- as.matrix(
    editability[
        ,
        ..editability_component_columns
    ]
)


storage.mode(
    component_matrix
) <- "double"


observed_components <- is.finite(
    component_matrix
)


weighted_numerator <- rowSums(
    sweep(
        ifelse(
            observed_components,
            component_matrix,
            0
        ),
        2,
        editability_component_weights,
        `*`
    )
)


observed_weight <- rowSums(
    sweep(
        observed_components * 1,
        2,
        editability_component_weights,
        `*`
    )
)


editability_score <- weighted_numerator /
    observed_weight


editability_score[
    observed_weight <= 0
] <- NA_real_


editability[
    ,
    crispr_reference_editability_score :=
        editability_score
]


editability[
    ,
    crispr_reference_editability_component_coverage :=
        observed_weight /
        sum(
            editability_component_weights
        )
]


editability[
    ,
    crispr_reference_editability_rank :=
        data.table::frank(
            -crispr_reference_editability_score,
            ties.method = "first",
            na.last = "keep"
        )
]


# -----------------------------------------------------------------------------
# 16. Status and interpretation fields
# -----------------------------------------------------------------------------

editability[
    ,
    step05_status :=
        data.table::fcase(
            !reference_sequence_available,
            "no_release63_transcript_sequence",

            n_unique_ngg_guides == 0L,
            "sequence_available_no_NGG",

            n_unique_reference_quality_guides == 0L,
            "NGG_present_no_reference_quality_guide",

            n_unique_reference_quality_guides >=
                MULTIPLEX_TARGET_4,
            "reference_editable_multiplex_4plus",

            n_unique_reference_quality_guides >=
                MULTIPLEX_TARGET_2,
            "reference_editable_multiplex_2plus",

            n_unique_reference_quality_guides >= 1L,
            "reference_editable_single_plus",

            default = "unclassified"
        )
]


editability[
    ,
    population_conservation_status :=
        "PENDING_TARGETED_AG1000G"
]


editability[
    ,
    genomewide_offtarget_status :=
        "NOT_EVALUATED_IN_STEP05"
]


# -----------------------------------------------------------------------------
# 17. Stable ordering
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


if (nrow(guides_transcript) > 0L) {

    data.table::setorder(
        guides_transcript,
        gene_id,
        -passes_reference_sequence_quality,
        -gc_balance_score,
        transcript_id,
        strand,
        transcript_start
    )
}


if (nrow(guides_unique) > 0L) {

    data.table::setorder(
        guides_unique,
        gene_id,
        -passes_reference_sequence_quality,
        reference_guide_rank_within_gene,
        protospacer_20nt
    )
}


# -----------------------------------------------------------------------------
# 18. Benchmark-specific table
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
# 19. QC
# -----------------------------------------------------------------------------

n_shortlist <- nrow(
    editability
)


n_with_sequence <- editability[
    reference_sequence_available == TRUE,
    .N
]


n_without_sequence <- editability[
    reference_sequence_available == FALSE,
    .N
]


n_with_ngg <- editability[
    reference_sequence_available == TRUE &
        n_unique_ngg_guides >= 1L,
    .N
]


n_with_gc30_70 <- editability[
    reference_sequence_available == TRUE &
        n_unique_gc30_70_guides >= 1L,
    .N
]


n_with_reference_quality <- editability[
    has_reference_quality_guide == TRUE,
    .N
]


n_multiplex2 <- editability[
    has_at_least_2_reference_quality_guides == TRUE,
    .N
]


n_multiplex4 <- editability[
    has_at_least_4_reference_quality_guides == TRUE,
    .N
]


qc <- data.table(
    metric = c(
        "Step 04 targeted shortlist genes",
        "Genes with release-63 transcript sequence",
        "Genes without release-63 transcript sequence",
        "Unique shortlisted transcripts scanned",
        "Transcript-level raw NGG guide instances",
        "Gene-level unique NGG guides",
        "Genes with >=1 canonical NGG guide",
        "Genes with >=1 guide passing 30-70% GC",
        "Genes with >=1 reference sequence-quality guide",
        "Genes with >=2 reference sequence-quality guides",
        "Genes with >=4 reference sequence-quality guides",
        "Population conservation evaluated",
        "Genome-wide off-target specificity evaluated"
    ),

    value = c(
        n_shortlist,
        n_with_sequence,
        n_without_sequence,
        nrow(transcripts),
        nrow(guides_transcript),
        nrow(guides_unique),
        n_with_ngg,
        n_with_gc30_70,
        n_with_reference_quality,
        n_multiplex2,
        n_multiplex4,
        0,
        0
    )
)


# -----------------------------------------------------------------------------
# 20. Provenance
# -----------------------------------------------------------------------------

provenance <- data.table(
    component = c(
        "Shortlist",
        "Reference transcript sequences",
        "SpCas9 target definition",
        "Schmidt-compatible GC filter",
        "Additional sequence-quality flags",
        "Guide deduplication",
        "Reference editability score",
        "Population conservation",
        "Off-target specificity"
    ),

    specification = c(
        paste0(
            "Step 04 targeted shortlist from ",
            SHORTLIST_FILE,
            "."
        ),

        paste0(
            "Ensembl Metazoa release ",
            ENSEMBL_GENOMES_RELEASE,
            ", Anopheles gambiae PEST assembly ",
            ENSEMBL_ASSEMBLY,
            "; release-pinned cDNA and ncRNA FASTA."
        ),

        paste0(
            "Canonical SpCas9 candidate = 20-nt protospacer followed by NGG PAM; ",
            "both transcript orientations scanned; overlapping sites retained."
        ),

        paste0(
            "Guide protospacer GC content between ",
            MIN_GUIDE_GC,
            "% and ",
            MAX_GUIDE_GC,
            "%, matching the reference-sequence GC criterion used by ",
            "Schmidt et al. 2020."
        ),

        paste0(
            "Reference sequence-quality guide requires 30-70% GC, no run of ",
            HOMOPOLYMER_RUN,
            "+ identical nucleotides, and no TTTT motif. These flags are ",
            "sequence-only and do not substitute for off-target analysis."
        ),

        paste0(
            "Identical protospacer+PAM sequences observed in multiple transcript ",
            "isoforms of the same gene are represented once in the gene-level ",
            "unique-guide table; isoform support count is retained."
        ),

        paste0(
            "Explainable shortlist-only score: 35% guide abundance, 25% guide ",
            "instance density, 20% multiplexability (saturates at 4 guides), ",
            "10% transcript coverage, 10% best GC balance. It is not a probability ",
            "of successful editing."
        ),

        paste0(
            "Not evaluated in Step 05. Population conservation and DRA frequency ",
            "remain PENDING_TARGETED_AG1000G."
        ),

        paste0(
            "Not evaluated in Step 05 because genome-wide off-target searching ",
            "requires the complete reference genome and a dedicated specificity ",
            "algorithm. No sequence is labelled off-target-safe here."
        )
    )
)


# -----------------------------------------------------------------------------
# 21. Write outputs
# -----------------------------------------------------------------------------

fwrite(
    guides_transcript,
    TRANSCRIPT_GUIDE_FILE,
    na = "NA"
)


fwrite(
    guides_unique,
    UNIQUE_GUIDE_FILE,
    na = "NA"
)


fwrite(
    editability,
    EDITABILITY_FILE,
    na = "NA"
)


fwrite(
    benchmark_editability,
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
    SHORTLIST_FILE,
    CDNA_FASTA,
    NCRNA_FASTA,
    TRANSCRIPT_GUIDE_FILE,
    UNIQUE_GUIDE_FILE,
    EDITABILITY_FILE,
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
# 22. Integrity checks
# -----------------------------------------------------------------------------

stopifnot(
    uniqueN(editability$gene_id) ==
        nrow(editability)
)


stopifnot(
    nrow(editability) ==
        nrow(shortlist)
)


stopifnot(
    all(
        editability[
            reference_sequence_available == TRUE,
            n_unique_reference_quality_guides
        ] >= 0,
        na.rm = TRUE
    )
)


if (nrow(guides_unique) > 0L) {

    stopifnot(
        uniqueN(
            guides_unique[
                ,
                paste(
                    gene_id,
                    protospacer_20nt,
                    pam,
                    sep = "|"
                )
            ]
        ) ==
        nrow(guides_unique)
    )
}


# -----------------------------------------------------------------------------
# 23. Console summary
# -----------------------------------------------------------------------------

cat(
    "\n",
    "============================================================\n",
    "MOSQEDIT-R STEP 05 COMPLETED SUCCESSFULLY\n",
    "Targeted SpCas9 reference-sequence editability\n",
    "============================================================\n",
    "Step 04 shortlist genes:                    ",
    format(n_shortlist, big.mark = ","),
    "\n",
    "Genes with release-63 transcript sequence:  ",
    format(n_with_sequence, big.mark = ","),
    "\n",
    "Genes without transcript sequence:          ",
    format(n_without_sequence, big.mark = ","),
    "\n",
    "Unique transcripts scanned:                 ",
    format(nrow(transcripts), big.mark = ","),
    "\n",
    "Raw transcript-level NGG guide instances:   ",
    format(nrow(guides_transcript), big.mark = ","),
    "\n",
    "Gene-level unique NGG guides:               ",
    format(nrow(guides_unique), big.mark = ","),
    "\n",
    "Genes with >=1 reference-quality guide:     ",
    format(n_with_reference_quality, big.mark = ","),
    "\n",
    "Genes with >=2 reference-quality guides:    ",
    format(n_multiplex2, big.mark = ","),
    "\n",
    "Genes with >=4 reference-quality guides:    ",
    format(n_multiplex4, big.mark = ","),
    "\n",
    "Population conservation evaluated:          NO\n",
    "Genome-wide off-target specificity:         NO\n",
    "Python used:                                NO\n",
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
    "\nTop 20 preliminary-ranked genes with Step 05 editability:\n"
)


display_columns <- c(
    "preliminary_rank",
    "gene_id",
    "benchmark_name",
    "n_transcripts_scanned",
    "n_unique_ngg_guides",
    "n_unique_reference_quality_guides",
    "fraction_transcripts_with_reference_quality_guide",
    "crispr_reference_editability_score",
    "crispr_reference_editability_rank",
    "step05_status"
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
        "\nBenchmark editability:\n"
    )


    benchmark_display_columns <- c(
        "gene_id",
        "benchmark_name",
        "preliminary_rank",
        "n_transcripts_scanned",
        "n_unique_ngg_guides",
        "n_unique_reference_quality_guides",
        "best_reference_guide_gc_percent",
        "crispr_reference_editability_score",
        "crispr_reference_editability_rank",
        "step05_status",
        "population_conservation_status"
    )


    benchmark_display_columns <- benchmark_display_columns[
        benchmark_display_columns %in%
            names(benchmark_editability)
    ]


    print(
        benchmark_editability[
            ,
            ..benchmark_display_columns
        ]
    )
}


log_step(
    "05",
    "Targeted SpCas9 reference-sequence editability completed successfully"
)

