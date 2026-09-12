# =============================================================================
# 03_acquire_orthology_flybase.R
#
# MosqEdit-R Manuscript 1
#
# Acquisition and harmonisation of:
#   1. Ensembl Metazoa release-63 Compara orthology between
#      Anopheles gambiae (AgamP4) and Drosophila melanogaster (BDGP6.54)
#   2. FlyBase FB2026_02 genotype-phenotype evidence
#
# IMPORTANT SOURCE DESIGN
# -----------------------
# Ensembl Genomes release 63 corresponds to Ensembl release 116.
# Orthology is downloaded from the official Compara TSV dumps rather than
# BioMart. This avoids archive-host DNS/BioMart HTTP failures and keeps the
# source genuinely release-pinned.
#
# Ensembl Compara species-specific homology files contain an arbitrary subset
# of orthologies involving that genome. Therefore BOTH the Anopheles gambiae
# and Drosophila melanogaster species-specific files are downloaded and merged.
# Protein and ncRNA collections are both included so that non-coding loci in
# the Step 01 universe are not silently excluded.
#
# Main outputs:
#   data_processed/03_ensembl_homology_file_manifest.csv
#   data_processed/03_anopheles_dmel_orthology.csv
#   data_processed/03_flybase_gene_phenotype_evidence.csv
#   data_processed/03_gene_orthology_phenotypes.csv
#   data_processed/03_orthology_phenotype_qc.csv
#   data_processed/03_orthology_provenance.csv
#   data_processed/03_flybase_phenotype_schema.csv
#   data_processed/03_flybase_allele_gene_schema.csv
#   logs/03_flybase_release.tsv
#   logs/03_ensembl_orthology_query.tsv
#   logs/03_sessionInfo.txt
#   logs/03_checksums.tsv
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Load helpers and packages
# -----------------------------------------------------------------------------

source("R/helpers.R")

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(curl)
})

# -----------------------------------------------------------------------------
# 2. Directories
# -----------------------------------------------------------------------------

dir.create("data_raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data_raw/EnsemblMetazoa63", recursive = TRUE, showWarnings = FALSE)
dir.create("data_raw/FlyBase", recursive = TRUE, showWarnings = FALSE)
dir.create("data_processed", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

log_step(
  "03",
  "Acquiring pinned Ensembl Compara orthology and FlyBase phenotype evidence"
)

# -----------------------------------------------------------------------------
# 3. Fixed source versions
# -----------------------------------------------------------------------------

ENSEMBL_GENOMES_RELEASE <- 63L
ENSEMBL_COMPARA_RELEASE <- 116L
ENSEMBL_DIVISION <- "metazoa"
ENSEMBL_AG_SPECIES <- "anopheles_gambiae"
ENSEMBL_DMEL_SPECIES <- "drosophila_melanogaster"
ENSEMBL_AG_ASSEMBLY <- "AgamP4"
ENSEMBL_DMEL_ASSEMBLY <- "BDGP6.54"

ENSEMBL_BASE_URLS <- c(
  "https://ftp.ebi.ac.uk/ensemblgenomes/pub",
  "https://ftp.ensemblgenomes.ebi.ac.uk/pub"
)

fb_release <- Sys.getenv(
  "MOSQEDIT_FLYBASE_RELEASE",
  unset = "FB2026_02"
)

if (!grepl("^FB[0-9]{4}_[0-9]{2}$", fb_release)) {
  stop(
    paste0(
      "Unexpected FlyBase release format: ", fb_release,
      ". Expected a value such as FB2026_02."
    ),
    call. = FALSE
  )
}

fb_tag <- tolower(sub("^FB", "", fb_release))

# -----------------------------------------------------------------------------
# 4. Load Step 01 AGAP universe
# -----------------------------------------------------------------------------

gene_annotation_file <- "data_processed/01_gene_annotation.csv"

if (!file.exists(gene_annotation_file)) {
  stop(
    paste0(
      "Step 01 output is missing: ", gene_annotation_file,
      "\nRun Step 01 before Step 03."
    ),
    call. = FALSE
  )
}

current_genes <- fread(gene_annotation_file)

if (!"gene_id" %in% names(current_genes)) {
  stop("Step 01 gene annotation does not contain gene_id.", call. = FALSE)
}

current_gene_ids <- unique(
  current_genes[
    !is.na(gene_id) & grepl("^AGAP", gene_id),
    gene_id
  ]
)

if (length(current_gene_ids) == 0L) {
  stop("Step 01 contains no AGAP gene identifiers.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 5. Utility functions
# -----------------------------------------------------------------------------

collapse_unique_nonempty <- function(x, sep = ";") {
  x <- unique(trimws(as.character(x)))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(NA_character_)
  paste(sort(x), collapse = sep)
}

safe_max_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(NA_real_)
  max(x, na.rm = TRUE)
}

parse_logical_or_numeric <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_real_, length(y))
  out[y %in% c("1", "true", "yes", "y")] <- 1
  out[y %in% c("0", "false", "no", "n")] <- 0
  numeric_hits <- suppressWarnings(as.numeric(y))
  replace <- is.na(out) & !is.na(numeric_hits)
  out[replace] <- numeric_hits[replace]
  out
}

# Download from multiple candidate URLs, keeping a valid cached file.
download_with_fallback <- function(
    urls,
    destfile,
    minimum_bytes = 1000L,
    label = basename(destfile)
) {
  if (
    file.exists(destfile) &&
      !is.na(file.info(destfile)$size) &&
      file.info(destfile)$size >= minimum_bytes
  ) {
    message("Using existing file: ", destfile)
    return("existing_local_file")
  }

  errors <- character()

  for (u in urls) {
    message("Trying: ", u)

    if (file.exists(destfile)) unlink(destfile)

    ok <- tryCatch(
      {
        curl::curl_download(
          u,
          destfile = destfile,
          quiet = FALSE
        )
        TRUE
      },
      error = function(e) {
        errors <<- c(errors, paste0(u, " -> ", conditionMessage(e)))
        FALSE
      }
    )

    if (
      ok &&
        file.exists(destfile) &&
        !is.na(file.info(destfile)$size) &&
        file.info(destfile)$size >= minimum_bytes
    ) {
      message(
        "Downloaded ", label, ": ",
        format(file.info(destfile)$size, big.mark = ","),
        " bytes"
      )
      return(u)
    }

    if (file.exists(destfile)) unlink(destfile)
  }

  stop(
    paste0(
      "Unable to download ", label, ".\n\nAttempts:\n",
      paste(errors, collapse = "\n")
    ),
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 6. Discover the actual release-63 Ensembl Compara homology directories/files
#
# IMPORTANT:
# Ensembl Metazoa release 63 contains multiple protein-tree collections
# (e.g. Metazoa, Protostomia, Insecta, Drosophilidae). Therefore the filename
# "protein_default" must NOT be assumed. We inspect the official FTP index,
# find the real species-specific directory for each genome, then identify all
# protein/ncRNA homology files shared by BOTH species.
# -----------------------------------------------------------------------------

ENSEMBL_HOMOLOGY_ROOTS <- paste0(
  ENSEMBL_BASE_URLS,
  "/",
  ENSEMBL_DIVISION,
  "/release-",
  ENSEMBL_GENOMES_RELEASE,
  "/tsv/ensembl-compara/homologies/"
)

fetch_index_links <- function(url) {
  h <- curl::new_handle()
  curl::handle_setopt(h, followlocation = TRUE, timeout = 120)

  res <- tryCatch(
    curl::curl_fetch_memory(url, handle = h),
    error = function(e) {
      stop(
        paste0(
          "Could not read Ensembl FTP directory index:\n",
          url,
          "\n",
          conditionMessage(e)
        ),
        call. = FALSE
      )
    }
  )

  if (res$status_code < 200L || res$status_code >= 300L) {
    stop(
      paste0(
        "Ensembl FTP directory returned HTTP ",
        res$status_code,
        ": ",
        url
      ),
      call. = FALSE
    )
  }

  txt <- rawToChar(res$content)
  hrefs <- stringr::str_match_all(
    txt,
    "href=[\\\"']([^\\\"']+)[\\\"']"
  )[[1]]

  if (nrow(hrefs) == 0L) return(character())

  links <- utils::URLdecode(hrefs[, 2])
  links <- links[!is.na(links) & nzchar(links)]
  links <- links[!links %in% c("../", "./")]
  unique(links)
}

get_working_homology_root <- function(urls) {
  errors <- character()

  for (u in urls) {
    message("Inspecting Ensembl homology index: ", u)

    out <- tryCatch(
      fetch_index_links(u),
      error = function(e) {
        errors <<- c(errors, paste0(u, " -> ", conditionMessage(e)))
        NULL
      }
    )

    if (!is.null(out) && length(out) > 0L) {
      return(list(root = u, links = out))
    }
  }

  stop(
    paste0(
      "Could not read any release-63 Ensembl homology directory index.\n\n",
      paste(errors, collapse = "\n")
    ),
    call. = FALSE
  )
}

join_url <- function(base, child) {
  if (grepl("^https?://", child, ignore.case = TRUE)) return(child)

  if (grepl("^/", child)) {
    origin <- sub("^(https?://[^/]+).*$", "\\1", base)
    return(paste0(origin, child))
  }

  paste0(sub("/+$", "", base), "/", sub("^/+", "", child))
}

find_species_homology_dir <- function(root, root_links, species) {
  # First try a direct species directory.
  direct <- root_links[
    grepl(paste0("(^|/)", species, "/?$"), root_links, ignore.case = TRUE)
  ]

  if (length(direct) > 0L) {
    return(join_url(root, direct[[1]]))
  }

  # Some genomes are stored inside core-database collection directories.
  collection_dirs <- root_links[
    grepl("_collection/?$", root_links, ignore.case = TRUE)
  ]

  for (d in collection_dirs) {
    d_url <- join_url(root, d)

    links <- tryCatch(
      fetch_index_links(d_url),
      error = function(e) character()
    )

    hit <- links[
      grepl(
        paste0("(^|/)", species, "([_/]|$)"),
        links,
        ignore.case = TRUE
      )
    ]

    hit <- hit[grepl("/$", hit)]

    if (length(hit) > 0L) {
      return(join_url(d_url, hit[[1]]))
    }
  }

  stop(
    paste0(
      "Could not discover the release-63 homology directory for species: ",
      species,
      "\nThe FTP structure may have changed; inspect: ",
      root
    ),
    call. = FALSE
  )
}

list_species_homology_files <- function(species_dir, species) {
  links <- fetch_index_links(species_dir)

  pattern <- paste0(
    "^Compara\\.",
    ENSEMBL_COMPARA_RELEASE,
    "\\.(protein|ncrna)_([^/]+)\\.homologies\\.tsv\\.gz$"
  )

  files <- basename(links)
  keep <- grepl(pattern, files)
  files <- files[keep]

  if (length(files) == 0L) {
    stop(
      paste0(
        "No Compara homology TSV files were found for ",
        species,
        " in:\n",
        species_dir
      ),
      call. = FALSE
    )
  }

  m <- stringr::str_match(files, pattern)

  data.table(
    species = species,
    member_type = m[, 2],
    tree_collection = m[, 3],
    remote_filename = files,
    remote_url = vapply(
      files,
      function(f) join_url(species_dir, f),
      character(1)
    )
  )
}

root_info <- get_working_homology_root(ENSEMBL_HOMOLOGY_ROOTS)
ENSEMBL_HOMOLOGY_ROOT <- root_info$root
root_links <- root_info$links

log_step(
  "03",
  paste("Using Ensembl homology root:", ENSEMBL_HOMOLOGY_ROOT)
)

ag_homology_dir <- find_species_homology_dir(
  root = ENSEMBL_HOMOLOGY_ROOT,
  root_links = root_links,
  species = ENSEMBL_AG_SPECIES
)

dmel_homology_dir <- find_species_homology_dir(
  root = ENSEMBL_HOMOLOGY_ROOT,
  root_links = root_links,
  species = ENSEMBL_DMEL_SPECIES
)

log_step("03", paste("A. gambiae homology directory:", ag_homology_dir))
log_step("03", paste("D. melanogaster homology directory:", dmel_homology_dir))

ag_available <- list_species_homology_files(
  ag_homology_dir,
  ENSEMBL_AG_SPECIES
)

dmel_available <- list_species_homology_files(
  dmel_homology_dir,
  ENSEMBL_DMEL_SPECIES
)

# Retain only gene-tree collections available for BOTH species.
shared_collections <- merge(
  unique(ag_available[, .(member_type, tree_collection)]),
  unique(dmel_available[, .(member_type, tree_collection)]),
  by = c("member_type", "tree_collection")
)

if (nrow(shared_collections) == 0L) {
  stop(
    paste0(
      "A. gambiae and D. melanogaster have no shared Compara homology ",
      "collections in the release-63 FTP index."
    ),
    call. = FALSE
  )
}

if (!any(shared_collections$member_type == "protein")) {
  stop(
    "No shared protein-tree homology collection was found.",
    call. = FALSE
  )
}

if (!any(shared_collections$member_type == "ncrna")) {
  warning(
    paste0(
      "No shared ncRNA-tree homology collection was found for A. gambiae ",
      "and D. melanogaster. Protein orthology will still be processed, but ",
      "ncRNA orthology coverage will be unavailable from this source."
    )
  )
}

cat("\nShared Ensembl Compara collections discovered:\n")
print(shared_collections[order(member_type, tree_collection)])

# Build the exact set of files that must be downloaded from BOTH species.
ensembl_files <- rbindlist(
  list(
    merge(
      ag_available,
      shared_collections,
      by = c("member_type", "tree_collection")
    ),
    merge(
      dmel_available,
      shared_collections,
      by = c("member_type", "tree_collection")
    )
  ),
  use.names = TRUE,
  fill = TRUE
)

ensembl_files[, local_file := file.path(
  "data_raw/EnsemblMetazoa63",
  paste0(
    species,
    "_",
    member_type,
    "_",
    tree_collection,
    "_homologies.tsv.gz"
  )
)]

ensembl_files[, source_url := remote_url]

# Save discovery manifest BEFORE downloads for reproducibility/debugging.
fwrite(
  ensembl_files[
    ,
    .(
      species,
      member_type,
      tree_collection,
      remote_filename,
      remote_url,
      local_file
    )
  ],
  "data_processed/03_ensembl_homology_file_manifest.csv",
  na = "NA"
)

# -----------------------------------------------------------------------------
# 7. Download all shared species-specific Compara files
# -----------------------------------------------------------------------------

for (i in seq_len(nrow(ensembl_files))) {
  sp <- ensembl_files$species[[i]]
  mt <- ensembl_files$member_type[[i]]
  tc <- ensembl_files$tree_collection[[i]]
  dest <- ensembl_files$local_file[[i]]
  u <- ensembl_files$remote_url[[i]]

  log_step(
    "03",
    paste(
      "Acquiring Ensembl Compara",
      mt,
      "collection",
      tc,
      "for",
      sp
    )
  )

  # The index has already confirmed the exact URL exists, so there is no
  # filename guessing here.
  ensembl_files$source_url[[i]] <- download_with_fallback(
    urls = u,
    destfile = dest,
    minimum_bytes = 1000L,
    label = paste(sp, mt, tc, "homology TSV")
  )
}

# -----------------------------------------------------------------------------
# 8. Read and validate Ensembl Compara homology TSVs
# -----------------------------------------------------------------------------

expected_compara_columns <- c(
  "gene_stable_id",
  "protein_stable_id",
  "species",
  "identity",
  "homology_type",
  "homology_gene_stable_id",
  "homology_protein_stable_id",
  "homology_species",
  "homology_identity",
  "dn",
  "ds",
  "goc_score",
  "wga_coverage",
  "is_high_confidence",
  "homology_id"
)

read_compara_tsv <- function(
    path,
    source_species,
    member_type,
    tree_collection
) {
  dt <- fread(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    fill = TRUE,
    data.table = TRUE,
    showProgress = FALSE
  )

  missing_cols <- setdiff(expected_compara_columns, names(dt))

  if (length(missing_cols) > 0L) {
    stop(
      paste0(
        "Unexpected Ensembl Compara schema in ", path, ".\n",
        "Missing columns: ", paste(missing_cols, collapse = ", "), "\n",
        "Observed columns: ", paste(names(dt), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  dt[, source_species_file := source_species]
  dt[, homology_collection := member_type]
  dt[, gene_tree_collection := tree_collection]
  dt
}

compara_list <- vector("list", nrow(ensembl_files))

for (i in seq_len(nrow(ensembl_files))) {
  compara_list[[i]] <- read_compara_tsv(
    path = ensembl_files$local_file[[i]],
    source_species = ensembl_files$species[[i]],
    member_type = ensembl_files$member_type[[i]],
    tree_collection = ensembl_files$tree_collection[[i]]
  )
}

compara_all <- rbindlist(
  compara_list,
  fill = TRUE,
  use.names = TRUE
)

if (nrow(compara_all) == 0L) {
  stop("Ensembl Compara files contained zero homology rows.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 9. Extract A. gambiae <-> D. melanogaster ORTHOLOGIES and orient consistently
# -----------------------------------------------------------------------------

pair_rows <- compara_all[
  (
    species == ENSEMBL_AG_SPECIES &
      homology_species == ENSEMBL_DMEL_SPECIES
  ) |
    (
      species == ENSEMBL_DMEL_SPECIES &
        homology_species == ENSEMBL_AG_SPECIES
    )
]

pair_rows <- pair_rows[
  grepl("ortholog", homology_type, ignore.case = TRUE)
]

if (nrow(pair_rows) == 0L) {
  stop(
    paste0(
      "No A. gambiae <-> D. melanogaster orthology rows were found in the ",
      "release-pinned Compara files."
    ),
    call. = FALSE
  )
}

forward <- pair_rows[species == ENSEMBL_AG_SPECIES]
reverse <- pair_rows[species == ENSEMBL_DMEL_SPECIES]

orth_forward <- data.table(
  gene_id = as.character(forward$gene_stable_id),
  dmel_ortholog = as.character(forward$homology_gene_stable_id),
  ag_percent_identity = suppressWarnings(as.numeric(forward$identity)),
  dmel_percent_identity = suppressWarnings(as.numeric(forward$homology_identity)),
  orthology_type = as.character(forward$homology_type),
  orthology_confidence = parse_logical_or_numeric(forward$is_high_confidence),
  goc_score = suppressWarnings(as.numeric(forward$goc_score)),
  wga_coverage = suppressWarnings(as.numeric(forward$wga_coverage)),
  homology_id = as.character(forward$homology_id),
  homology_collection = as.character(forward$homology_collection),
  gene_tree_collection = as.character(forward$gene_tree_collection),
  source_species_file = as.character(forward$source_species_file)
)

orth_reverse <- data.table(
  gene_id = as.character(reverse$homology_gene_stable_id),
  dmel_ortholog = as.character(reverse$gene_stable_id),
  ag_percent_identity = suppressWarnings(as.numeric(reverse$homology_identity)),
  dmel_percent_identity = suppressWarnings(as.numeric(reverse$identity)),
  orthology_type = as.character(reverse$homology_type),
  orthology_confidence = parse_logical_or_numeric(reverse$is_high_confidence),
  goc_score = suppressWarnings(as.numeric(reverse$goc_score)),
  wga_coverage = suppressWarnings(as.numeric(reverse$wga_coverage)),
  homology_id = as.character(reverse$homology_id),
  homology_collection = as.character(reverse$homology_collection),
  gene_tree_collection = as.character(reverse$gene_tree_collection),
  source_species_file = as.character(reverse$source_species_file)
)

orth_relationships <- rbindlist(
  list(orth_forward, orth_reverse),
  fill = TRUE,
  use.names = TRUE
)

orth_relationships <- orth_relationships[
  !is.na(gene_id) &
    grepl("^AGAP", gene_id) &
    !is.na(dmel_ortholog) &
    grepl("^FBgn", dmel_ortholog)
]

orth_relationships <- orth_relationships[
  gene_id %in% current_gene_ids
]

# The same homology can occur in both species-specific files. Deduplicate by
# biological relationship + collection, keeping the maximum available metrics.
orth_relationships <- orth_relationships[
  ,
  .(
    ag_percent_identity = safe_max_numeric(ag_percent_identity),
    dmel_percent_identity = safe_max_numeric(dmel_percent_identity),
    orthology_confidence = safe_max_numeric(orthology_confidence),
    goc_score = safe_max_numeric(goc_score),
    wga_coverage = safe_max_numeric(wga_coverage),
    source_species_file = collapse_unique_nonempty(source_species_file),
    homology_id = collapse_unique_nonempty(homology_id),
    gene_tree_collection = collapse_unique_nonempty(gene_tree_collection)
  ),
  by = .(
    gene_id,
    dmel_ortholog,
    orthology_type,
    homology_collection
  )
]

if (nrow(orth_relationships) == 0L) {
  stop(
    "No current Step 01 AGAP loci retained a D. melanogaster ortholog.",
    call. = FALSE
  )
}

setorder(orth_relationships, gene_id, dmel_ortholog, homology_collection)

# Dmel names are added after the FlyBase gene table is parsed below.
orth_relationships[, dmel_gene_name := NA_character_]

# -----------------------------------------------------------------------------
# 10. Record Ensembl source provenance immediately
# -----------------------------------------------------------------------------

ensembl_log <- rbindlist(
  list(
    data.table(
      key = c(
        "Ensembl_Genomes_release",
        "Ensembl_release_Compara",
        "division",
        "Anopheles_species_key",
        "Drosophila_species_key",
        "Anopheles_assembly",
        "Drosophila_assembly",
        "orthology_retrieval_method"
      ),
      value = c(
        as.character(ENSEMBL_GENOMES_RELEASE),
        as.character(ENSEMBL_COMPARA_RELEASE),
        ENSEMBL_DIVISION,
        ENSEMBL_AG_SPECIES,
        ENSEMBL_DMEL_SPECIES,
        ENSEMBL_AG_ASSEMBLY,
        ENSEMBL_DMEL_ASSEMBLY,
        "release-pinned Ensembl Compara species-specific TSV dumps"
      )
    ),
    ensembl_files[
      ,
      .(
        key = paste0("source_", species, "_", member_type, "_", tree_collection),
        value = source_url
      )
    ]
  )
)

fwrite(
  ensembl_log,
  "logs/03_ensembl_orthology_query.tsv",
  sep = "\t",
  na = "NA"
)

# -----------------------------------------------------------------------------
# 11. Pinned FlyBase bulk-data files
#
# FlyBase hosts its public archive in the AWS Open Data bucket
# s3ftp.flybase.org. Some networks cannot resolve the custom hostname
# s3ftp.flybase.org even though standard AWS S3 endpoints are reachable.
# Therefore each pinned file is attempted through:
#   1. canonical FlyBase hostname
#   2. AWS path-style global S3 endpoint
#   3. AWS path-style us-east-1 endpoint
#
# All URLs address the exact same public S3 bucket/object.
# -----------------------------------------------------------------------------

phen_file <- file.path(
  "data_raw/FlyBase",
  paste0("genotype_phenotype_data_fb_", fb_tag, ".tsv.gz")
)

map_file <- file.path(
  "data_raw/FlyBase",
  paste0("fbal_to_fbgn_fb_", fb_tag, ".tsv.gz")
)

gene_map_file <- file.path(
  "data_raw/FlyBase",
  paste0("gene_map_table_fb_", fb_tag, ".tsv.gz")
)

flybase_object_key <- function(subdir, filename) {
  paste0(
    "releases/",
    fb_release,
    "/precomputed_files/",
    subdir,
    "/",
    filename
  )
}

flybase_candidate_urls <- function(subdir, filename) {
  key <- flybase_object_key(subdir, filename)

  c(
    paste0("https://s3ftp.flybase.org/", key),
    paste0("https://s3.amazonaws.com/s3ftp.flybase.org/", key),
    paste0("https://s3.us-east-1.amazonaws.com/s3ftp.flybase.org/", key)
  )
}

phen_urls <- flybase_candidate_urls(
  "alleles",
  basename(phen_file)
)

map_urls <- flybase_candidate_urls(
  "alleles",
  basename(map_file)
)

gene_map_urls <- flybase_candidate_urls(
  "genes",
  basename(gene_map_file)
)

# Canonical URLs retained separately for reporting/citation.
phen_url <- phen_urls[[1]]
map_url <- map_urls[[1]]
gene_map_url <- gene_map_urls[[1]]

# -----------------------------------------------------------------------------
# 12. Download FlyBase files
# -----------------------------------------------------------------------------

log_step(
  "03",
  paste("Acquiring FlyBase phenotype file from pinned release", fb_release)
)

phen_source <- download_with_fallback(
  urls = phen_urls,
  destfile = phen_file,
  minimum_bytes = 10000L,
  label = "FlyBase genotype-phenotype table"
)

log_step(
  "03",
  paste("Acquiring FlyBase allele-to-gene file from pinned release", fb_release)
)

map_source <- download_with_fallback(
  urls = map_urls,
  destfile = map_file,
  minimum_bytes = 1000L,
  label = "FlyBase FBal-to-FBgn table"
)

log_step(
  "03",
  paste("Acquiring optional FlyBase gene-map file from pinned release", fb_release)
)

gene_map_source <- tryCatch(
  download_with_fallback(
    urls = gene_map_urls,
    destfile = gene_map_file,
    minimum_bytes = 1000L,
    label = "FlyBase gene map table"
  ),
  error = function(e) {
    warning(
      paste0(
        "FlyBase gene-map table could not be acquired. Dmel gene names will ",
        "remain unavailable; phenotype mapping can still proceed.\n",
        conditionMessage(e)
      )
    )
    NA_character_
  }
)

writeLines(
  c(
    paste0("FlyBase_release\t", fb_release),
    paste0("phenotype_canonical_url\t", phen_url),
    paste0("phenotype_resolved_source\t", phen_source),
    paste0("allele_gene_canonical_url\t", map_url),
    paste0("allele_gene_resolved_source\t", map_source),
    paste0("gene_map_canonical_url\t", gene_map_url),
    paste0("gene_map_resolved_source\t", ifelse(is.na(gene_map_source), "NA", gene_map_source)),
    paste0("phenotype_file\t", phen_file),
    paste0("allele_gene_file\t", map_file),
    paste0("gene_map_file\t", gene_map_file),
    paste0("retrieval_date\t", Sys.Date())
  ),
  "logs/03_flybase_release.tsv"
)

# -----------------------------------------------------------------------------
# 13. Robust reader for FlyBase commented TSV files
# -----------------------------------------------------------------------------

read_flybase_tsv <- function(path) {
  con <- gzfile(path, open = "rt")

  preview <- tryCatch(
    readLines(con, n = 500L, warn = FALSE),
    finally = close(con)
  )

  if (length(preview) == 0L) {
    stop(paste("FlyBase file is empty:", path), call. = FALSE)
  }

  # Select the first tab-delimited line that looks like a data header.
  tabbed <- which(grepl("\t", preview, fixed = TRUE) & nzchar(trimws(preview)))

  if (length(tabbed) == 0L) {
    stop(
      paste0("Could not identify a tab-delimited header in FlyBase file: ", path),
      call. = FALSE
    )
  }

  header_line <- tabbed[[1]]
  header_text <- sub("^#+[[:space:]]*", "", preview[[header_line]])
  header_names <- strsplit(header_text, "\t", fixed = TRUE)[[1]]
  header_names <- make.unique(trimws(header_names))

  out <- fread(
    path,
    sep = "\t",
    header = FALSE,
    skip = header_line,
    quote = "",
    fill = TRUE,
    col.names = header_names,
    data.table = TRUE,
    showProgress = FALSE
  )

  if (nrow(out) == 0L) {
    stop(paste("FlyBase file contains no data rows:", path), call. = FALSE)
  }

  out
}

ph <- read_flybase_tsv(phen_file)
am <- read_flybase_tsv(map_file)

# Optional gene-map table for Dmel display names/symbols.
gm <- NULL
if (file.exists(gene_map_file)) {
  gm <- tryCatch(
    read_flybase_tsv(gene_map_file),
    error = function(e) {
      warning(
        paste0(
          "Could not parse FlyBase gene-map table; continuing without Dmel names.\n",
          conditionMessage(e)
        )
      )
      NULL
    }
  )
}

# -----------------------------------------------------------------------------
# 14. Save actual FlyBase schemas used
# -----------------------------------------------------------------------------

fwrite(
  data.table(
    column_index = seq_along(names(ph)),
    column_name = names(ph),
    R_class = vapply(ph, function(x) class(x)[[1]], character(1))
  ),
  "data_processed/03_flybase_phenotype_schema.csv"
)

fwrite(
  data.table(
    column_index = seq_along(names(am)),
    column_name = names(am),
    R_class = vapply(am, function(x) class(x)[[1]], character(1))
  ),
  "data_processed/03_flybase_allele_gene_schema.csv"
)

# -----------------------------------------------------------------------------
# 15. Construct row text defensively
# -----------------------------------------------------------------------------

make_row_blob <- function(dt) {
  do.call(
    paste,
    c(
      lapply(
        dt,
        function(x) {
          y <- as.character(x)
          y[is.na(y)] <- ""
          y
        }
      ),
      sep = " | "
    )
  )
}

ph[, phenotype_row_id := .I]
ph[, phenotype_blob := make_row_blob(.SD), .SDcols = setdiff(names(ph), "phenotype_blob")]

am[, map_row_id := .I]
am[, map_blob := make_row_blob(.SD), .SDcols = setdiff(names(am), "map_blob")]

# -----------------------------------------------------------------------------
# 16. Expand all FBal identifiers from phenotype records
# -----------------------------------------------------------------------------

ph_long_list <- lapply(
  seq_len(nrow(ph)),
  function(i) {
    ids <- unique(
      stringr::str_extract_all(
        ph$phenotype_blob[[i]],
        "FBal[0-9]+"
      )[[1]]
    )

    if (length(ids) == 0L) return(NULL)

    data.table(
      phenotype_row_id = ph$phenotype_row_id[[i]],
      phenotype_blob = ph$phenotype_blob[[i]],
      fbal = ids
    )
  }
)

ph_long <- rbindlist(ph_long_list, fill = TRUE)

if (nrow(ph_long) == 0L) {
  stop(
    "No FBal identifiers were recovered from FlyBase phenotype records.",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 17. Parse FBal -> FBgn mapping with strict ambiguity QC
# -----------------------------------------------------------------------------

am[, fbal_ids := str_extract_all(map_blob, "FBal[0-9]+")]
am[, fbgn_ids := str_extract_all(map_blob, "FBgn[0-9]+")]
am[, n_fbal_ids := lengths(fbal_ids)]
am[, n_fbgn_ids := lengths(fbgn_ids)]

ambiguous_map_rows <- am[n_fbal_ids != 1L | n_fbgn_ids != 1L]

if (nrow(ambiguous_map_rows) > 0L) {
  fwrite(
    ambiguous_map_rows[
      ,
      .(map_row_id, n_fbal_ids, n_fbgn_ids, map_blob)
    ],
    "data_processed/03_flybase_ambiguous_allele_gene_rows.csv",
    na = "NA"
  )

  stop(
    paste0(
      "The FlyBase FBal-to-FBgn file contains ",
      nrow(ambiguous_map_rows),
      " row(s) that do not contain exactly one FBal and one FBgn. ",
      "Inspect data_processed/03_flybase_ambiguous_allele_gene_rows.csv."
    ),
    call. = FALSE
  )
}

am_map <- unique(
  data.table(
    fbal = vapply(am$fbal_ids, function(x) x[[1]], character(1)),
    fbgn = vapply(am$fbgn_ids, function(x) x[[1]], character(1))
  )
)

if (nrow(am_map) == 0L) {
  stop("Could not construct FlyBase FBal-to-FBgn mapping.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# 18. Optional FlyBase FBgn -> gene-name/symbol mapping
# -----------------------------------------------------------------------------

fbgn_name_map <- data.table(
  dmel_ortholog = character(),
  dmel_gene_name = character()
)

if (!is.null(gm) && nrow(gm) > 0L) {
  gm_blob <- make_row_blob(gm)
  gm_fbgn <- str_extract(gm_blob, "FBgn[0-9]+")

  normalized_names <- tolower(gsub("[^a-z0-9]+", "_", names(gm)))
  symbol_hits <- which(
    normalized_names %in% c(
      "symbol",
      "gene_symbol",
      "gene_name",
      "gene"
    )
  )

  if (length(symbol_hits) > 0L) {
    symbol_col <- names(gm)[symbol_hits[[1]]]
    symbol_values <- trimws(as.character(gm[[symbol_col]]))

    fbgn_name_map <- unique(
      data.table(
        dmel_ortholog = gm_fbgn,
        dmel_gene_name = symbol_values
      )[
        !is.na(dmel_ortholog) &
          nzchar(dmel_ortholog) &
          !is.na(dmel_gene_name) &
          nzchar(dmel_gene_name)
      ]
    )

    # Retain only unambiguous names for a given FBgn.
    fbgn_name_map <- fbgn_name_map[
      ,
      .(
        n_names = uniqueN(dmel_gene_name),
        dmel_gene_name = dmel_gene_name[[1]]
      ),
      by = dmel_ortholog
    ][n_names == 1L, .(dmel_ortholog, dmel_gene_name)]
  }
}

if (nrow(fbgn_name_map) > 0L) {
  orth_relationships[, dmel_gene_name := NULL]
  orth_relationships <- merge(
    orth_relationships,
    fbgn_name_map,
    by = "dmel_ortholog",
    all.x = TRUE
  )
}

setorder(orth_relationships, gene_id, dmel_ortholog, homology_collection)

fwrite(
  orth_relationships,
  "data_processed/03_ensembl_homology_file_manifest.csv",
  "data_processed/03_anopheles_dmel_orthology.csv",
  na = "NA"
)

# -----------------------------------------------------------------------------
# 19. Attach FlyBase phenotype records to genes
# -----------------------------------------------------------------------------

ph_gene <- merge(
  ph_long,
  am_map,
  by = "fbal",
  all = FALSE
)

if (nrow(ph_gene) == 0L) {
  stop(
    "Could not map FlyBase genotype-phenotype rows to FBgn genes.",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# 20. Phenotype vocabularies
#
# IMPORTANT:
# The strict female-sterility label contains explicit fertility/sterility terms.
# Broader developmental/reproductive terms such as maternal effect and
# oogenesis are NOT automatically treated as female sterility.
# -----------------------------------------------------------------------------

female_sterile_re <- paste(
  c(
    "female sterile",
    "female sterility",
    "female infertile",
    "female infertility",
    "reduced female fertility",
    "female fertility reduced",
    "reduced female fecundity",
    "female fecundity reduced"
  ),
  collapse = "|"
)

repro_re <- paste(
  c(
    female_sterile_re,
    "grandchildless",
    "maternal[- ]effect",
    "fertility",
    "fecundity",
    "ovary",
    "ovarian",
    "oogenesis",
    "oocyte",
    "egg[- ]laying",
    "egg laying",
    "germline",
    "germ cell",
    "gametogenesis"
  ),
  collapse = "|"
)

ph_gene[, female_sterile_hit := grepl(
  female_sterile_re,
  phenotype_blob,
  ignore.case = TRUE,
  perl = TRUE
)]

ph_gene[, reproductive_hit := grepl(
  repro_re,
  phenotype_blob,
  ignore.case = TRUE,
  perl = TRUE
)]

# -----------------------------------------------------------------------------
# 21. FlyBase evidence summarized per Dmel gene
# -----------------------------------------------------------------------------

fb_ev <- ph_gene[
  ,
  .(
    flybase_phenotype_record_count = uniqueN(phenotype_row_id),
    flybase_female_sterile = any(female_sterile_hit, na.rm = TRUE),
    flybase_female_sterile_record_count = uniqueN(
      phenotype_row_id[female_sterile_hit %in% TRUE]
    ),
    flybase_reproductive_phenotype_count = uniqueN(
      phenotype_row_id[reproductive_hit %in% TRUE]
    )
  ),
  by = fbgn
]

setorder(fb_ev, fbgn)

fwrite(
  fb_ev,
  "data_processed/03_flybase_gene_phenotype_evidence.csv",
  na = "NA"
)

# -----------------------------------------------------------------------------
# 22. Gene-level orthology summary for ALL current Step 01 loci
# -----------------------------------------------------------------------------

orth_summary <- orth_relationships[
  ,
  .(
    dmel_ortholog_count = uniqueN(dmel_ortholog),
    dmel_ortholog = collapse_unique_nonempty(dmel_ortholog),
    dmel_gene_name = collapse_unique_nonempty(dmel_gene_name),
    orthology_type = collapse_unique_nonempty(orthology_type),
    orthology_collection = collapse_unique_nonempty(homology_collection),
    orthology_confidence_max = safe_max_numeric(orthology_confidence),
    dmel_percent_identity_max = safe_max_numeric(dmel_percent_identity),
    ag_percent_identity_max = safe_max_numeric(ag_percent_identity),
    goc_score_max = safe_max_numeric(goc_score),
    wga_coverage_max = safe_max_numeric(wga_coverage)
  ),
  by = gene_id
]

# -----------------------------------------------------------------------------
# 23. Summarize FlyBase evidence across all Dmel orthologs of each AGAP locus
# -----------------------------------------------------------------------------

orth_to_ph <- merge(
  unique(
    orth_relationships[
      ,
      .(gene_id, dmel_ortholog)
    ]
  ),
  ph_gene[
    ,
    .(
      fbgn,
      phenotype_row_id,
      female_sterile_hit,
      reproductive_hit
    )
  ],
  by.x = "dmel_ortholog",
  by.y = "fbgn",
  all = FALSE
)

if (nrow(orth_to_ph) > 0L) {
  gene_ph_summary <- orth_to_ph[
    ,
    .(
      flybase_phenotype_record_count = uniqueN(phenotype_row_id),
      flybase_female_sterile = any(female_sterile_hit, na.rm = TRUE),
      flybase_female_sterile_record_count = uniqueN(
        phenotype_row_id[female_sterile_hit %in% TRUE]
      ),
      flybase_reproductive_phenotype_count = uniqueN(
        phenotype_row_id[reproductive_hit %in% TRUE]
      )
    ),
    by = gene_id
  ]
} else {
  gene_ph_summary <- data.table(
    gene_id = character(),
    flybase_phenotype_record_count = integer(),
    flybase_female_sterile = logical(),
    flybase_female_sterile_record_count = integer(),
    flybase_reproductive_phenotype_count = integer()
  )
}

# -----------------------------------------------------------------------------
# 24. Construct final Step 03 gene-level table
# -----------------------------------------------------------------------------

gene_ev <- data.table(gene_id = sort(current_gene_ids))

gene_ev <- merge(
  gene_ev,
  orth_summary,
  by = "gene_id",
  all.x = TRUE
)

gene_ev <- merge(
  gene_ev,
  gene_ph_summary,
  by = "gene_id",
  all.x = TRUE
)

gene_ev[, has_dmel_ortholog := !is.na(dmel_ortholog_count) & dmel_ortholog_count > 0L]
gene_ev[, has_flybase_phenotype_records := !is.na(flybase_phenotype_record_count) & flybase_phenotype_record_count > 0L]

# Absence of a mapped phenotype record is NOT treated as a biological negative.
gene_ev[
  !has_flybase_phenotype_records,
  `:=`(
    flybase_female_sterile = NA,
    flybase_female_sterile_record_count = NA_integer_,
    flybase_reproductive_phenotype_count = NA_integer_
  )
]

setorder(gene_ev, gene_id)

fwrite(
  gene_ev,
  "data_processed/03_gene_orthology_phenotypes.csv",
  na = "NA"
)

# -----------------------------------------------------------------------------
# 25. QC metrics
# -----------------------------------------------------------------------------

n_step01_genes <- length(current_gene_ids)
n_genes_with_ortholog <- gene_ev[has_dmel_ortholog == TRUE, .N]
n_genes_with_flybase_records <- gene_ev[has_flybase_phenotype_records == TRUE, .N]
n_genes_female_sterile <- gene_ev[flybase_female_sterile %in% TRUE, .N]
n_genes_reproductive <- gene_ev[
  !is.na(flybase_reproductive_phenotype_count) &
    flybase_reproductive_phenotype_count > 0L,
  .N
]

qc <- data.table(
  metric = c(
    "Step 01 current AGAP loci",
    "Raw pairwise orthology rows after orientation/deduplication",
    "Unique current AGAP loci with Dmel ortholog",
    "Orthology coverage percent",
    "Unique Dmel ortholog genes",
    "Protein orthology relationship rows",
    "ncRNA orthology relationship rows",
    "High-confidence orthology rows",
    "FlyBase phenotype source rows",
    "FlyBase allele-gene mapping rows",
    "FlyBase phenotype rows containing FBal",
    "FlyBase phenotype-to-gene relationship rows",
    "Current AGAP loci with mapped FlyBase phenotype records",
    "Current AGAP loci with explicit female-sterility evidence",
    "Current AGAP loci with broader reproductive phenotype evidence"
  ),
  value = c(
    n_step01_genes,
    nrow(orth_relationships),
    n_genes_with_ortholog,
    round(100 * n_genes_with_ortholog / n_step01_genes, 2),
    uniqueN(orth_relationships$dmel_ortholog),
    orth_relationships[homology_collection == "protein", .N],
    orth_relationships[homology_collection == "ncrna", .N],
    orth_relationships[orthology_confidence == 1, .N],
    nrow(ph),
    nrow(am_map),
    uniqueN(ph_long$phenotype_row_id),
    nrow(ph_gene),
    n_genes_with_flybase_records,
    n_genes_female_sterile,
    n_genes_reproductive
  )
)

fwrite(qc, "data_processed/03_orthology_phenotype_qc.csv")

# -----------------------------------------------------------------------------
# 26. Provenance
# -----------------------------------------------------------------------------

ensembl_resource <- paste0(
  "Ensembl Genomes release ", ENSEMBL_GENOMES_RELEASE,
  " / Ensembl Compara ", ENSEMBL_COMPARA_RELEASE,
  " species-specific homology TSV dumps"
)

provenance <- data.table(
  source_layer = c(
    "Ensembl Metazoa orthology",
    "FlyBase genotype phenotype",
    "FlyBase allele-to-gene map",
    "FlyBase Dmel gene map"
  ),
  release = c(
    paste0(
      "Ensembl Genomes ", ENSEMBL_GENOMES_RELEASE,
      " / Ensembl ", ENSEMBL_COMPARA_RELEASE
    ),
    fb_release,
    fb_release,
    fb_release
  ),
  resource = c(
    ensembl_resource,
    ifelse(identical(phen_source, "existing_local_file"), phen_url, phen_source),
    ifelse(identical(map_source, "existing_local_file"), map_url, map_source),
    ifelse(is.na(gene_map_source) || identical(gene_map_source, "existing_local_file"), gene_map_url, gene_map_source)
  ),
  local_file = c(
    "data_processed/03_ensembl_homology_file_manifest.csv",
    phen_file,
    map_file,
    gene_map_file
  ),
  notes = c(
    paste0(
      "Dynamically discovered and combined all shared species-specific protein/ncRNA Compara collections from the release-63 FTP index; ",
      "A. gambiae AgamP4 <-> D. melanogaster BDGP6.54; ortholog rows only"
    ),
    "Pinned FlyBase bulk genotype-phenotype table",
    "Pinned FlyBase FBal-to-FBgn mapping table",
    "Optional pinned FlyBase gene map used only for Dmel display names"
  ),
  retrieval_date = as.character(Sys.Date())
)

fwrite(
  provenance,
  "data_processed/03_orthology_provenance.csv",
  na = "NA"
)

# -----------------------------------------------------------------------------
# 27. Save session information
# -----------------------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = "logs/03_sessionInfo.txt"
)

# -----------------------------------------------------------------------------
# 28. Checksums
# -----------------------------------------------------------------------------

checksum_files <- c(
  "data_processed/03_ensembl_homology_file_manifest.csv",
  "data_processed/03_anopheles_dmel_orthology.csv",
  "data_processed/03_flybase_gene_phenotype_evidence.csv",
  "data_processed/03_gene_orthology_phenotypes.csv",
  "data_processed/03_orthology_phenotype_qc.csv",
  "data_processed/03_orthology_provenance.csv",
  "data_processed/03_flybase_phenotype_schema.csv",
  "data_processed/03_flybase_allele_gene_schema.csv",
  "logs/03_flybase_release.tsv",
  "logs/03_ensembl_orthology_query.tsv",
  "logs/03_sessionInfo.txt",
  ensembl_files$local_file,
  phen_file,
  map_file
)

if (file.exists(gene_map_file)) {
  checksum_files <- c(checksum_files, gene_map_file)
}

if (file.exists("data_processed/03_flybase_ambiguous_allele_gene_rows.csv")) {
  checksum_files <- c(
    checksum_files,
    "data_processed/03_flybase_ambiguous_allele_gene_rows.csv"
  )
}

checksum_files <- unique(checksum_files[file.exists(checksum_files)])

write_checksum(
  checksum_files,
  "logs/03_checksums.tsv"
)

# -----------------------------------------------------------------------------
# 29. Console validation
# -----------------------------------------------------------------------------

cat(
  "\n",
  "============================================================\n",
  "MOSQEDIT-R STEP 03 COMPLETED SUCCESSFULLY\n",
  "============================================================\n",
  "Ensembl Genomes release:        ", ENSEMBL_GENOMES_RELEASE, "\n",
  "Ensembl Compara release:        ", ENSEMBL_COMPARA_RELEASE, "\n",
  "Orthology source:               pinned Compara TSV dumps\n",
  "Anopheles assembly:             ", ENSEMBL_AG_ASSEMBLY, "\n",
  "Drosophila assembly:            ", ENSEMBL_DMEL_ASSEMBLY, "\n",
  "FlyBase release:                ", fb_release, "\n",
  "Shared Compara collections:      ", paste(paste0(shared_collections$member_type, ":", shared_collections$tree_collection), collapse = "; "), "\n",
  "Step 01 current AGAP loci:      ", format(n_step01_genes, big.mark = ","), "\n",
  "Orthology relationship rows:    ", format(nrow(orth_relationships), big.mark = ","), "\n",
  "AGAP loci with Dmel ortholog:   ", format(n_genes_with_ortholog, big.mark = ","), "\n",
  "Orthology coverage:             ", round(100 * n_genes_with_ortholog / n_step01_genes, 2), "%\n",
  "Unique Dmel ortholog genes:     ", format(uniqueN(orth_relationships$dmel_ortholog), big.mark = ","), "\n",
  "Protein orthology rows:         ", format(orth_relationships[homology_collection == "protein", .N], big.mark = ","), "\n",
  "ncRNA orthology rows:           ", format(orth_relationships[homology_collection == "ncrna", .N], big.mark = ","), "\n",
  "AGAP loci with FlyBase records: ", format(n_genes_with_flybase_records, big.mark = ","), "\n",
  "Female-sterility evidence loci: ", format(n_genes_female_sterile, big.mark = ","), "\n",
  "Reproductive evidence loci:     ", format(n_genes_reproductive, big.mark = ","), "\n",
  "============================================================\n",
  sep = ""
)

cat("\nQC summary:\n")
print(qc)

cat("\nOrthology types:\n")
print(
  orth_relationships[
    ,
    .N,
    by = .(orthology_type, homology_collection)
  ][order(-N)]
)

cat("\nFirst 10 gene-level Step 03 records:\n")
print(gene_ev[1:min(10L, nrow(gene_ev))])

cat("\nOutput files:\n")
cat(
  "  data_processed/03_anopheles_dmel_orthology.csv\n",
  "  data_processed/03_flybase_gene_phenotype_evidence.csv\n",
  "  data_processed/03_gene_orthology_phenotypes.csv\n",
  "  data_processed/03_orthology_phenotype_qc.csv\n",
  "  data_processed/03_orthology_provenance.csv\n",
  "  data_processed/03_flybase_phenotype_schema.csv\n",
  "  data_processed/03_flybase_allele_gene_schema.csv\n",
  "  logs/03_flybase_release.tsv\n",
  "  logs/03_ensembl_orthology_query.tsv\n",
  "  logs/03_sessionInfo.txt\n",
  "  logs/03_checksums.tsv\n",
  sep = ""
)

log_step(
  "03",
  "Orthology and FlyBase phenotype evidence acquisition completed successfully"
)

