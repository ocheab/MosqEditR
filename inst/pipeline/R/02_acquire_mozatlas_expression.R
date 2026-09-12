# =============================================================================
# 02_acquire_mozatlas_expression.R
#
# MosqEdit-R Manuscript 1
#
# Acquisition, harmonisation and feature engineering of the
# MozAtlas Anopheles gambiae adult tissue-expression dataset.
#
# Primary source:
#   NCBI Gene Expression Omnibus
#   Series:   GSE21689
#   Platform: GPL1321
#
# Study:
#   MozAtlas: A comprehensive gene expression atlas of sex- and
#   tissue-specificity in the malaria vector, Anopheles gambiae
#
# Experimental design:
#   60 samples
#   15 adult tissue/sex groups
#   4 biological replicates per group
#
# IMPORTANT:
#   ACPS_1â€“4 are male accessory-gland / Acps samples.
#   They are biological tissue samples and MUST NOT be removed as controls.
#
# Expression data:
#   VALUE              = GEO-supplied processed GC-RMA signal intensity
#   ABS_CALL           = MAS 5.0 detection call
#   DETECTION P-VALUE  = MAS 5.0 detection P-value
#
# MosqEdit-R transformation policy:
#   The GEO-supplied processed expression values are retained exactly as
#   retrieved. No additional log2 or inverse-log transformation is applied.
#
# Network strategy:
#   This script DOES NOT depend on ftp.ncbi.nlm.nih.gov.
#   GEO records are acquired through:
#
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi
#
# Main outputs:
#
#   data_raw/GEO/GSE21689_family.soft
#
#   data_processed/02_mozatlas_sample_manifest.csv
#   data_processed/02_mozatlas_probe_annotation.csv
#   data_processed/02_mozatlas_probe_gene_mapping.csv
#   data_processed/02_mozatlas_probe_expression.csv
#   data_processed/02_mozatlas_probe_detection.csv
#   data_processed/02_mozatlas_gene_expression.csv
#   data_processed/02_mozatlas_expression_qc.csv
#   data_processed/02_mozatlas_detection_qc.csv
#   data_processed/02_mozatlas_expression_provenance.csv
#
#   logs/02_expression_scale.tsv
#   logs/02_sessionInfo.txt
#   logs/02_checksums.tsv
#
# =============================================================================


# -----------------------------------------------------------------------------
# 1. Load helper functions and packages
# -----------------------------------------------------------------------------

source("R/helpers.R")


suppressPackageStartupMessages({
  
  library(GEOquery)
  
  library(Biobase)
  
  library(data.table)
  
  library(stringr)
})


# -----------------------------------------------------------------------------
# 2. Create directories
# -----------------------------------------------------------------------------

dir.create(
  "data_raw/GEO",
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
# 3. Fixed study identifiers
# -----------------------------------------------------------------------------

GEO_SERIES   <- "GSE21689"
GEO_PLATFORM <- "GPL1321"

SPECIES <- "Anopheles gambiae"


log_step(
  "02",
  paste(
    "Acquiring",
    GEO_SERIES,
    "MozAtlas expression data"
  )
)


# -----------------------------------------------------------------------------
# 4. Increase network timeout
#
# Do not use on.exit() at script top level.
# Individual helper functions restore their own timeout settings.
# -----------------------------------------------------------------------------

options(
  timeout = max(
    1800,
    getOption("timeout")
  )
)


# -----------------------------------------------------------------------------
# 5. Robust GEO download helper
#
# Uses www.ncbi.nlm.nih.gov rather than ftp.ncbi.nlm.nih.gov.
# -----------------------------------------------------------------------------

download_geo_www <- function(
    url,
    destfile,
    minimum_bytes = 1000L,
    timeout = 1800
) {
  
  old_timeout <- getOption("timeout")
  
  
  on.exit(
    options(
      timeout = old_timeout
    ),
    add = TRUE
  )
  
  
  options(
    timeout = max(
      timeout,
      old_timeout
    )
  )
  
  
  if (file.exists(destfile)) {
    
    unlink(
      destfile
    )
  }
  
  
  warning_messages <- character()
  
  
  status <- tryCatch(
    
    {
      
      withCallingHandlers(
        
        download.file(
          url      = url,
          destfile = destfile,
          mode     = "wb",
          method   = "libcurl",
          quiet    = FALSE
        ),
        
        warning = function(w) {
          
          warning_messages <<- c(
            warning_messages,
            conditionMessage(w)
          )
          
          
          invokeRestart(
            "muffleWarning"
          )
        }
      )
    },
    
    error = function(e) {
      
      stop(
        paste0(
          "\nUnable to download GEO record.\n\n",
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
    is.na(status) ||
    !identical(
      as.integer(status),
      0L
    )
  ) {
    
    stop(
      paste(
        "GEO download returned status:",
        status
      ),
      call. = FALSE
    )
  }
  
  
  if (!file.exists(destfile)) {
    
    stop(
      "GEO download completed without creating the expected file.",
      call. = FALSE
    )
  }
  
  
  downloaded_size <- file.info(
    destfile
  )$size
  
  
  if (
    is.na(downloaded_size) ||
    downloaded_size < minimum_bytes
  ) {
    
    stop(
      paste0(
        "Downloaded GEO file appears unexpectedly small: ",
        downloaded_size,
        " bytes."
      ),
      call. = FALSE
    )
  }
  
  
  message(
    "\nGEO download successful."
  )
  
  
  message(
    "Downloaded size: ",
    format(
      downloaded_size,
      big.mark = ","
    ),
    " bytes"
  )
  
  
  if (
    length(
      warning_messages
    ) > 0L
  ) {
    
    message(
      "\nWarnings captured during download:\n",
      paste(
        unique(
          warning_messages
        ),
        collapse = "\n"
      )
    )
  }
  
  
  invisible(
    url
  )
}


# -----------------------------------------------------------------------------
# 6. Test access to the ordinary NCBI GEO web endpoint
# -----------------------------------------------------------------------------

geo_test_url <- paste0(
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?",
  "acc=",
  GEO_SERIES,
  "&targ=self",
  "&view=brief",
  "&form=text"
)


geo_test_file <- file.path(
  "data_raw/GEO",
  paste0(
    GEO_SERIES,
    "_connectivity_test.txt"
  )
)


log_step(
  "02",
  "Testing NCBI GEO access through www.ncbi.nlm.nih.gov"
)


download_geo_www(
  url           = geo_test_url,
  destfile      = geo_test_file,
  minimum_bytes = 100L,
  timeout       = 900
)


test_lines <- readLines(
  geo_test_file,
  n    = 100L,
  warn = FALSE
)


if (
  !any(
    grepl(
      GEO_SERIES,
      test_lines,
      fixed = TRUE
    )
  )
) {
  
  stop(
    paste0(
      "NCBI GEO returned a file, but ",
      GEO_SERIES,
      " was not detected in the response."
    ),
    call. = FALSE
  )
}


log_step(
  "02",
  "NCBI GEO www endpoint is accessible"
)


# -----------------------------------------------------------------------------
# 7. Download complete GSE21689 family SOFT
#
# targ=all:
#   retrieve the Series, Samples and associated Platform.
#
# view=full:
#   include full sample/platform data tables.
#
# form=text:
#   return machine-readable SOFT text.
# -----------------------------------------------------------------------------

family_soft_url <- paste0(
  
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?",
  
  "acc=",
  GEO_SERIES,
  
  "&targ=all",
  
  "&view=full",
  
  "&form=text"
)


family_soft_file <- file.path(
  "data_raw/GEO",
  paste0(
    GEO_SERIES,
    "_family.soft"
  )
)


if (
  file.exists(
    family_soft_file
  ) &&
  !is.na(
    file.info(
      family_soft_file
    )$size
  ) &&
  file.info(
    family_soft_file
  )$size > 1000000
) {
  
  log_step(
    "02",
    paste(
      "Existing GEO family SOFT found:",
      family_soft_file
    )
  )
  
  
  family_soft_source <-
    "Existing local GSE21689 family SOFT"
  
} else {
  
  log_step(
    "02",
    "Downloading complete GSE21689 family SOFT from NCBI GEO"
  )
  
  
  download_geo_www(
    
    url =
      family_soft_url,
    
    destfile =
      family_soft_file,
    
    minimum_bytes =
      1000000L,
    
    timeout =
      3600
  )
  
  
  family_soft_source <-
    family_soft_url
}


# -----------------------------------------------------------------------------
# 8. Validate downloaded family SOFT
# -----------------------------------------------------------------------------

soft_preview <- readLines(
  family_soft_file,
  n    = 500L,
  warn = FALSE
)


if (
  !any(
    grepl(
      GEO_SERIES,
      soft_preview,
      fixed = TRUE
    )
  )
) {
  
  stop(
    paste0(
      "Downloaded family SOFT does not appear to contain ",
      GEO_SERIES,
      "."
    ),
    call. = FALSE
  )
}


log_step(
  "02",
  paste(
    "Validated family SOFT:",
    format(
      file.info(
        family_soft_file
      )$size,
      big.mark = ","
    ),
    "bytes"
  )
)


# -----------------------------------------------------------------------------
# 9. Parse local family SOFT
#
# GSEMatrix=FALSE ensures GEOquery parses the SOFT family rather than trying
# the Series Matrix route.
# -----------------------------------------------------------------------------

log_step(
  "02",
  "Parsing locally cached GSE21689 family SOFT"
)


geo_family <- tryCatch(
  
  GEOquery::getGEO(
    filename  = family_soft_file,
    GSEMatrix = FALSE
  ),
  
  error = function(e) {
    
    stop(
      paste0(
        "GSE21689 family SOFT downloaded successfully, ",
        "but GEOquery could not parse it.\n\n",
        "Original error:\n",
        conditionMessage(e)
      ),
      call. = FALSE
    )
  }
)


# -----------------------------------------------------------------------------
# 10. Resolve returned GSE object
# -----------------------------------------------------------------------------

if (inherits(geo_family, "GSE")) {
  
  gse_soft <- geo_family
  
} else if (
  is.list(geo_family) &&
  length(geo_family) >= 1L
) {
  
  gse_hits <- which(
    vapply(
      geo_family,
      function(x) inherits(x, "GSE"),
      logical(1)
    )
  )
  
  if (length(gse_hits) == 0L) {
    
    stop(
      paste0(
        "No GSE object was returned after parsing the family SOFT.\n",
        "Returned class: ",
        paste(
          class(geo_family),
          collapse = ", "
        )
      ),
      call. = FALSE
    )
  }
  
  # IMPORTANT: use [[...]] together
  gse_soft <- geo_family[[gse_hits[[1]]]]
  
} else {
  
  stop(
    paste0(
      "Unexpected object returned from GEOquery.\n",
      "Class: ",
      paste(
        class(geo_family),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


# Confirm what was recovered
cat(
  "\nRecovered GEO object class:\n"
)

print(
  class(gse_soft)
)


# -----------------------------------------------------------------------------
# 11. Extract GSM sample objects
# -----------------------------------------------------------------------------

gsm_list <- GEOquery::GSMList(
  gse_soft
)


if (
  length(
    gsm_list
  ) == 0L
) {
  
  stop(
    "No GSM samples were recovered from GSE21689.",
    call. = FALSE
  )
}


if (
  length(
    gsm_list
  ) != 60L
) {
  
  stop(
    paste0(
      "Unexpected number of GSE21689 samples: ",
      length(
        gsm_list
      ),
      ". Expected 60."
    ),
    call. = FALSE
  )
}


log_step(
  "02",
  paste(
    "Recovered",
    length(
      gsm_list
    ),
    "GSM samples"
  )
)


# -----------------------------------------------------------------------------
# 12. Metadata scalar helper
# -----------------------------------------------------------------------------


meta_scalar <- function(gsm, field) {
  
  metadata <- GEOquery::Meta(gsm)
  
  if (!field %in% names(metadata)) {
    return(NA_character_)
  }
  
  value <- metadata[[field]]
  
  if (length(value) == 0L) {
    return(NA_character_)
  }
  
  as.character(value[[1]])
}


# -----------------------------------------------------------------------------
# 13. Extract sample metadata
# -----------------------------------------------------------------------------

sample_accessions <- names(
  gsm_list
)


sample_titles <- vapply(
  
  gsm_list,
  
  meta_scalar,
  
  character(1),
  
  field = "title"
)


sample_platforms <- vapply(
  
  gsm_list,
  
  meta_scalar,
  
  character(1),
  
  field = "platform_id"
)


if (
  any(
    is.na(
      sample_titles
    )
  )
) {
  
  stop(
    "One or more GSE21689 samples lack a GEO title.",
    call. = FALSE
  )
}


if (
  any(
    !is.na(
      sample_platforms
    ) &
    sample_platforms !=
    GEO_PLATFORM
  )
) {
  
  stop(
    paste0(
      "One or more GSE21689 samples use a platform other than ",
      GEO_PLATFORM,
      "."
    ),
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 14. Parse MozAtlas sample titles
#
# Typical titles:
#
#   ACPS_1
#   Female_Body_1
#   Female_Carcass_1
#   Female_Head_1
#   Female_Malpighian_1
#   Female_Midgut_1
#   Female_Ovary_1
#   Female_Salivary_1
#   Male_Body_1
#   Male_Carcass_1
#   Male_Head_1
#   Male_Malpighian_1
#   Male_Midgut_1
#   Male_Salivary_1
#   Male_Testis_1
#
# ACPS = male accessory-gland samples.
# -----------------------------------------------------------------------------

parse_mozatlas_title <- function(x) {
  
  x <- as.character(
    x
  )
  
  
  # ACPS
  if (
    grepl(
      "^ACPS_[0-9]+$",
      x,
      ignore.case = TRUE
    )
  ) {
    
    replicate_number <- as.integer(
      str_extract(
        x,
        "[0-9]+$"
      )
    )
    
    
    return(
      data.table(
        
        sample_title =
          x,
        
        sex =
          "male",
        
        tissue =
          "accessory_gland",
        
        group =
          "male_accessory_gland",
        
        replicate =
          replicate_number
      )
    )
  }
  
  
  hit <- str_match(
    
    x,
    
    regex(
      "^(Female|Male)_([A-Za-z]+)_([0-9]+)$",
      ignore_case = TRUE
    )
  )
  
  
  if (
    is.na(
      hit[1, 1]
    )
  ) {
    
    stop(
      paste(
        "Unable to parse MozAtlas sample title:",
        x
      ),
      call. = FALSE
    )
  }
  
  
  sex <- tolower(
    hit[1, 2]
  )
  
  
  tissue <- tolower(
    hit[1, 3]
  )
  
  
  replicate_number <- as.integer(
    hit[1, 4]
  )
  
  
  data.table(
    
    sample_title =
      x,
    
    sex =
      sex,
    
    tissue =
      tissue,
    
    group =
      paste(
        sex,
        tissue,
        sep = "_"
      ),
    
    replicate =
      replicate_number
  )
}


sample_manifest <- rbindlist(
  
  lapply(
    sample_titles,
    parse_mozatlas_title
  )
)


sample_manifest[
  ,
  sample_accession :=
    sample_accessions
]


setcolorder(
  
  sample_manifest,
  
  c(
    "sample_accession",
    "sample_title",
    "sex",
    "tissue",
    "group",
    "replicate"
  )
)


# -----------------------------------------------------------------------------
# 15. Validate expected tissue/sex groups
# -----------------------------------------------------------------------------

expected_groups <- c(
  
  "male_accessory_gland",
  
  "female_body",
  "male_body",
  
  "female_carcass",
  "male_carcass",
  
  "female_head",
  "male_head",
  
  "female_malpighian",
  "male_malpighian",
  
  "female_midgut",
  "male_midgut",
  
  "female_ovary",
  
  "female_salivary",
  "male_salivary",
  
  "male_testis"
)


observed_groups <- unique(
  sample_manifest$group
)


missing_groups <- setdiff(
  expected_groups,
  observed_groups
)


unexpected_groups <- setdiff(
  observed_groups,
  expected_groups
)


if (
  length(
    missing_groups
  ) > 0L
) {
  
  stop(
    paste0(
      "Missing expected MozAtlas groups: ",
      paste(
        missing_groups,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


if (
  length(
    unexpected_groups
  ) > 0L
) {
  
  stop(
    paste0(
      "Unexpected MozAtlas groups: ",
      paste(
        unexpected_groups,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


group_counts <- sample_manifest[
  ,
  .N,
  by = group
]


group_counts[
  ,
  sort_order :=
    match(
      group,
      expected_groups
    )
]


setorder(
  group_counts,
  sort_order
)


group_counts[
  ,
  sort_order := NULL
]


if (
  any(
    group_counts$N !=
    4L
  )
) {
  
  stop(
    paste0(
      "Expected exactly 4 biological replicates per MozAtlas group.\n",
      paste(
        paste0(
          group_counts$group,
          "=",
          group_counts$N
        ),
        collapse = "; "
      )
    ),
    call. = FALSE
  )
}


log_step(
  "02",
  "Validated 15 tissue/sex groups with 4 biological replicates each"
)


fwrite(
  sample_manifest,
  "data_processed/02_mozatlas_sample_manifest.csv",
  na = "NA"
)


# -----------------------------------------------------------------------------
# 16. Extract each GSM data table
#
# Expected columns:
#
#   ID_REF
#   VALUE
#   ABS_CALL
#   DETECTION P-VALUE
# -----------------------------------------------------------------------------

gsm_tables <- vector(
  "list",
  length(
    gsm_list
  )
)


names(
  gsm_tables
) <- sample_accessions


for (
  i in seq_along(
    gsm_list
  )
) {
  
  dt <- as.data.table(
    GEOquery::Table(
      gsm_list[[i]]
    )
  )
  
  
  if (
    !all(
      c(
        "ID_REF",
        "VALUE"
      ) %in%
      names(
        dt
      )
    )
  ) {
    
    stop(
      paste0(
        "Sample ",
        sample_accessions[[i]],
        " does not contain required ID_REF and VALUE columns."
      ),
      call. = FALSE
    )
  }
  
  
  dt[
    ,
    ID_REF :=
      as.character(
        ID_REF
      )
  ]
  
  
  dt[
    ,
    VALUE :=
      suppressWarnings(
        as.numeric(
          VALUE
        )
      )
  ]
  
  
  gsm_tables[[i]] <- dt
}


# -----------------------------------------------------------------------------
# 17. Establish common probe universe
# -----------------------------------------------------------------------------

reference_probe_ids <- gsm_tables[[1]]$ID_REF


if (
  length(
    reference_probe_ids
  ) == 0L
) {
  
  stop(
    "First GSE21689 sample contains zero probe rows.",
    call. = FALSE
  )
}


for (
  i in seq_along(
    gsm_tables
  )
) {
  
  current_ids <- gsm_tables[[i]]$ID_REF
  
  
  if (
    length(
      current_ids
    ) !=
    length(
      reference_probe_ids
    )
  ) {
    
    stop(
      paste0(
        "Probe count differs for sample ",
        sample_accessions[[i]],
        "."
      ),
      call. = FALSE
    )
  }
  
  
  if (
    !setequal(
      current_ids,
      reference_probe_ids
    )
  ) {
    
    stop(
      paste0(
        "Probe universe differs for sample ",
        sample_accessions[[i]],
        "."
      ),
      call. = FALSE
    )
  }
  
  
  reorder_index <- match(
    reference_probe_ids,
    current_ids
  )
  
  
  if (
    anyNA(
      reorder_index
    )
  ) {
    
    stop(
      paste0(
        "Probe ordering failed for sample ",
        sample_accessions[[i]],
        "."
      ),
      call. = FALSE
    )
  }
  
  
  gsm_tables[[i]] <- gsm_tables[[i]][
    reorder_index
  ]
}


n_probes <- length(
  reference_probe_ids
)


n_samples <- length(
  gsm_tables
)


log_step(
  "02",
  paste(
    "Validated common probe universe:",
    format(
      n_probes,
      big.mark = ","
    ),
    "probe sets"
  )
)


if (
  n_probes !=
  22769L
) {
  
  warning(
    paste0(
      "GSE21689 contains ",
      n_probes,
      " probe sets; historical GPL1321 size is 22,769."
    )
  )
}


# -----------------------------------------------------------------------------
# 18. Build expression matrix
# -----------------------------------------------------------------------------

expr <- matrix(
  
  NA_real_,
  
  nrow =
    n_probes,
  
  ncol =
    n_samples,
  
  dimnames = list(
    
    reference_probe_ids,
    
    sample_accessions
  )
)


for (
  i in seq_along(
    gsm_tables
  )
) {
  
  expr[
    ,
    i
  ] <- gsm_tables[[i]]$VALUE
}


if (
  all(
    is.na(
      expr
    )
  )
) {
  
  stop(
    "Expression matrix contains only missing values.",
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 19. Build ABS_CALL detection matrix
# -----------------------------------------------------------------------------

detection_call <- matrix(
  
  NA_character_,
  
  nrow =
    n_probes,
  
  ncol =
    n_samples,
  
  dimnames = list(
    
    reference_probe_ids,
    
    sample_accessions
  )
)


has_abs_call <- vapply(
  
  gsm_tables,
  
  function(x) {
    
    "ABS_CALL" %in%
      names(
        x
      )
  },
  
  logical(1)
)


if (
  all(
    has_abs_call
  )
) {
  
  for (
    i in seq_along(
      gsm_tables
    )
  ) {
    
    detection_call[
      ,
      i
    ] <- toupper(
      as.character(
        gsm_tables[[i]]$ABS_CALL
      )
    )
  }
  
} else {
  
  warning(
    "ABS_CALL is not available in every GSE21689 sample table."
  )
}


# -----------------------------------------------------------------------------
# 20. Locate and build detection-P matrix
# -----------------------------------------------------------------------------

find_detection_p_column <- function(
    dt
) {
  
  normalized_names <- toupper(
    gsub(
      "[^A-Z0-9]",
      "",
      names(
        dt
      )
    )
  )
  
  
  hit <- which(
    normalized_names %in%
      c(
        "DETECTIONPVALUE",
        "DETECTIONPVAL"
      )
  )
  
  
  if (
    length(
      hit
    ) == 0L
  ) {
    
    return(
      NA_character_
    )
  }
  
  
  names(
    dt
  )[
    hit[[1]]
  ]
}


detection_p_columns <- vapply(
  
  gsm_tables,
  
  find_detection_p_column,
  
  character(1)
)


detection_p <- matrix(
  
  NA_real_,
  
  nrow =
    n_probes,
  
  ncol =
    n_samples,
  
  dimnames = list(
    
    reference_probe_ids,
    
    sample_accessions
  )
)


if (
  all(
    !is.na(
      detection_p_columns
    )
  )
) {
  
  for (
    i in seq_along(
      gsm_tables
    )
  ) {
    
    detection_p[
      ,
      i
    ] <- suppressWarnings(
      as.numeric(
        gsm_tables[[i]][[detection_p_columns[[i]]]]
      )
    )
  }
  
} else {
  
  warning(
    "Detection P-value column is not available in every sample."
  )
}


# -----------------------------------------------------------------------------
# 21. Expression-scale audit
#
# No transformation is performed.
# -----------------------------------------------------------------------------

expr_quantiles <- quantile(
  
  as.numeric(
    expr
  ),
  
  probs = c(
    0,
    0.01,
    0.05,
    0.25,
    0.50,
    0.75,
    0.95,
    0.99,
    1
  ),
  
  na.rm = TRUE
)


expression_scale_decision <- data.table(
  
  property = c(
    
    "GEO_series",
    
    "GEO_platform",
    
    "GEO_VALUE_description",
    
    "MosqEditR_transformation",
    
    paste0(
      "expression_quantile_",
      names(
        expr_quantiles
      )
    )
  ),
  
  value = c(
    
    GEO_SERIES,
    
    GEO_PLATFORM,
    
    paste0(
      "GEO-supplied normalized GC-RMA signal intensity; ",
      "retained exactly as submitted"
    ),
    
    "NONE",
    
    as.character(
      as.numeric(
        expr_quantiles
      )
    )
  )
)


fwrite(
  expression_scale_decision,
  "logs/02_expression_scale.tsv",
  sep = "\t"
)


# -----------------------------------------------------------------------------
# 22. Recover GPL1321 from the downloaded family SOFT
# -----------------------------------------------------------------------------

gpl_list <- tryCatch(
  
  GEOquery::GPLList(
    gse_soft
  ),
  
  error = function(e) {
    
    list()
  }
)


gpl <- NULL


if (
  length(
    gpl_list
  ) > 0L
) {
  
  if (
    GEO_PLATFORM %in%
    names(
      gpl_list
    )
  ) {
    
    gpl <- gpl_list[[GEO_PLATFORM]]
    
  } else {
    
    gpl_hits <- which(
      vapply(
        gpl_list,
        inherits,
        logical(1),
        what = "GPL"
      )
    )
    
    
    if (
      length(
        gpl_hits
      ) > 0L
    ) {
      
      gpl <- gpl_list[[gpl_hits[[1]]]]
    }
  }
}


# -----------------------------------------------------------------------------
# 23. Fallback: download GPL1321 separately if family SOFT did not contain it
# -----------------------------------------------------------------------------

if (
  is.null(
    gpl
  )
) {
  
  log_step(
    "02",
    "GPL1321 not recovered from family SOFT; downloading platform separately"
  )
  
  
  gpl_url <- paste0(
    
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?",
    
    "acc=",
    GEO_PLATFORM,
    
    "&targ=self",
    
    "&view=full",
    
    "&form=text"
  )
  
  
  gpl_soft_file <- file.path(
    
    "data_raw/GEO",
    
    paste0(
      GEO_PLATFORM,
      "_full.soft"
    )
  )
  
  
  if (
    !file.exists(
      gpl_soft_file
    ) ||
    is.na(
      file.info(
        gpl_soft_file
      )$size
    ) ||
    file.info(
      gpl_soft_file
    )$size < 10000
  ) {
    
    download_geo_www(
      
      url =
        gpl_url,
      
      destfile =
        gpl_soft_file,
      
      minimum_bytes =
        10000L,
      
      timeout =
        1800
    )
  }
  
  
  gpl_object <- GEOquery::getGEO(
    filename  = gpl_soft_file,
    GSEMatrix = FALSE
  )
  
  
  if (
    inherits(
      gpl_object,
      "GPL"
    )
  ) {
    
    gpl <- gpl_object
    
  } else if (
    is.list(
      gpl_object
    )
  ) {
    
    gpl_hits <- which(
      vapply(
        gpl_object,
        inherits,
        logical(1),
        what = "GPL"
      )
    )
    
    
    if (
      length(
        gpl_hits
      ) > 0L
    ) {
      
      gpl <- gpl_object[[gpl_hits[[1]]]]
    }
  }
}


if (
  is.null(
    gpl
  ) ||
  !inherits(
    gpl,
    "GPL"
  )
) {
  
  stop(
    "Unable to recover GPL1321 platform annotation.",
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 24. Extract platform table
# -----------------------------------------------------------------------------

platform_dt <- as.data.table(
  GEOquery::Table(
    gpl
  )
)


if (
  nrow(
    platform_dt
  ) == 0L
) {
  
  stop(
    "GPL1321 platform annotation contains zero rows.",
    call. = FALSE
  )
}


# Identify platform probe-ID column
if (
  "ID" %in%
  names(
    platform_dt
  )
) {
  
  setnames(
    platform_dt,
    "ID",
    "probe_id"
  )
  
} else if (
  "ID_REF" %in%
  names(
    platform_dt
  )
) {
  
  setnames(
    platform_dt,
    "ID_REF",
    "probe_id"
  )
  
} else {
  
  stop(
    paste0(
      "GPL1321 annotation does not contain an ID or ID_REF column.\n",
      "Available columns:\n",
      paste(
        names(
          platform_dt
        ),
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}


platform_dt[
  ,
  probe_id :=
    as.character(
      probe_id
    )
]


# Restrict platform annotation to probes represented in expression data
platform_dt <- platform_dt[
  probe_id %in%
    reference_probe_ids
]


if (
  nrow(
    platform_dt
  ) == 0L
) {
  
  stop(
    "No GPL1321 platform probes overlap the GSE21689 expression matrix.",
    call. = FALSE
  )
}


fwrite(
  platform_dt,
  "data_processed/02_mozatlas_probe_annotation.csv",
  na = "NA"
)


log_step(
  "02",
  paste(
    "Parsed GPL1321 platform annotation:",
    format(
      nrow(
        platform_dt
      ),
      big.mark = ","
    ),
    "probe sets"
  )
)


# -----------------------------------------------------------------------------
# 25. Create searchable annotation text for every platform probe
# -----------------------------------------------------------------------------

annotation_columns <- setdiff(
  names(
    platform_dt
  ),
  "probe_id"
)


if (
  length(
    annotation_columns
  ) == 0L
) {
  
  stop(
    "GPL1321 contains no annotation columns beyond probe ID.",
    call. = FALSE
  )
}


annotation_blob <- apply(
  
  platform_dt[
    ,
    ..annotation_columns
  ],
  
  1,
  
  function(z) {
    
    paste(
      z,
      collapse = " | "
    )
  }
)


# -----------------------------------------------------------------------------
# 26. Direct probe-to-AGAP mapping
#
# Use ALL AGAP identifiers present in platform annotation rather than extracting
# only the first one.
# -----------------------------------------------------------------------------

agap_list <- stringr::str_extract_all(
  
  annotation_blob,
  
  "AGAP[0-9]{6}"
)


direct_map_list <- lapply(
  
  seq_along(
    agap_list
  ),
  
  function(i) {
    
    ids <- unique(
      agap_list[[i]]
    )
    
    
    if (
      length(
        ids
      ) == 0L
    ) {
      
      return(
        NULL
      )
    }
    
    
    data.table(
      
      probe_id =
        platform_dt$probe_id[
          i
        ],
      
      gene_id =
        ids,
      
      mapping_method =
        "direct_AGAP_from_GPL1321"
    )
  }
)


direct_map <- rbindlist(
  direct_map_list,
  fill = TRUE
)


if (
  nrow(
    direct_map
  ) == 0L
) {
  
  direct_map <- data.table(
    
    probe_id =
      character(),
    
    gene_id =
      character(),
    
    mapping_method =
      character()
  )
}


# -----------------------------------------------------------------------------
# 27. Load current release-63 AgamP4 universe from Step 01
# -----------------------------------------------------------------------------

gene_annotation_file <-
  "data_processed/01_gene_annotation.csv"


if (
  !file.exists(
    gene_annotation_file
  )
) {
  
  stop(
    paste0(
      "Step 01 output is missing:\n",
      gene_annotation_file,
      "\nRun Step 01 before Step 02."
    ),
    call. = FALSE
  )
}


current_genes <- fread(
  gene_annotation_file
)


if (
  !"gene_id" %in%
  names(
    current_genes
  )
) {
  
  stop(
    "Step 01 annotation does not contain gene_id.",
    call. = FALSE
  )
}


current_gene_ids <- unique(
  current_genes$gene_id
)


# -----------------------------------------------------------------------------
# 28. Direct mapping QC
# -----------------------------------------------------------------------------

if (
  nrow(
    direct_map
  ) > 0L
) {
  
  direct_map[
    ,
    current_AgamP4 :=
      gene_id %in%
      current_gene_ids
  ]
  
  
  direct_map[
    ,
    n_gene_ids_per_probe :=
      uniqueN(
        gene_id
      ),
    by = probe_id
  ]
  
} else {
  
  direct_map[
    ,
    current_AgamP4 :=
      logical()
  ]
  
  
  direct_map[
    ,
    n_gene_ids_per_probe :=
      integer()
  ]
}


# -----------------------------------------------------------------------------
# 29. Conservative unique gene-symbol fallback
#
# The fallback is permitted only when:
#
#   - no direct AGAP mapping exists;
#   - platform annotation explicitly mentions Anopheles;
#   - exactly one gene symbol is supplied;
#   - that symbol uniquely identifies exactly one current AgamP4 locus.
#
# This prevents accidental mapping of Plasmodium probes.
# -----------------------------------------------------------------------------

normalized_platform_names <- tolower(
  gsub(
    "[^a-z0-9]+",
    "_",
    names(
      platform_dt
    )
  )
)


symbol_candidates <- which(
  normalized_platform_names %in%
    c(
      "gene_symbol",
      "gene_symbols",
      "symbol"
    )
)


symbol_map <- data.table(
  
  probe_id =
    character(),
  
  gene_id =
    character(),
  
  mapping_method =
    character(),
  
  current_AgamP4 =
    logical(),
  
  n_gene_ids_per_probe =
    integer()
)


if (
  length(
    symbol_candidates
  ) > 0L &&
  "gene_symbol" %in%
  names(
    current_genes
  )
) {
  
  symbol_col <- names(
    platform_dt
  )[
    symbol_candidates[[1]]
  ]
  
  
  current_symbol_map <- current_genes[
    !is.na(
      gene_symbol
    ) &
      nzchar(
        trimws(
          gene_symbol
        )
      ),
    .(
      n_genes =
        uniqueN(
          gene_id
        ),
      
      gene_id =
        first(
          gene_id
        )
    ),
    by = .(
      gene_symbol =
        trimws(
          gene_symbol
        )
    )
  ][
    n_genes == 1L
  ]
  
  
  platform_symbols <- data.table(
    
    probe_id =
      platform_dt$probe_id,
    
    gene_symbol =
      trimws(
        as.character(
          platform_dt[[symbol_col]]
        )
      ),
    
    annotation_text =
      annotation_blob
  )
  
  
  platform_symbols <- platform_symbols[
    
    !is.na(
      gene_symbol
    ) &
      
      nzchar(
        gene_symbol
      ) &
      
      grepl(
        "Anopheles",
        annotation_text,
        ignore.case = TRUE
      ) &
      
      !grepl(
        "///|//|;|\\|",
        gene_symbol
      )
  ]
  
  
  symbol_map_temp <- merge(
    
    platform_symbols,
    
    current_symbol_map[
      ,
      .(
        gene_symbol,
        gene_id
      )
    ],
    
    by =
      "gene_symbol",
    
    all =
      FALSE
  )
  
  
  if (
    nrow(
      direct_map
    ) > 0L
  ) {
    
    symbol_map_temp <- symbol_map_temp[
      !probe_id %in%
        direct_map$probe_id
    ]
  }
  
  
  if (
    nrow(
      symbol_map_temp
    ) > 0L
  ) {
    
    symbol_map <- symbol_map_temp[
      ,
      .(
        probe_id,
        gene_id,
        
        mapping_method =
          "unique_Anopheles_gene_symbol",
        
        current_AgamP4 =
          TRUE,
        
        n_gene_ids_per_probe =
          1L
      )
    ]
  }
}


# -----------------------------------------------------------------------------
# 30. Combine mapping strategies
# -----------------------------------------------------------------------------

mapping_all <- rbindlist(
  
  list(
    
    direct_map[
      ,
      .(
        probe_id,
        gene_id,
        mapping_method,
        current_AgamP4,
        n_gene_ids_per_probe
      )
    ],
    
    symbol_map
  ),
  
  fill = TRUE
)


mapping_all <- unique(
  mapping_all
)


if (
  nrow(
    mapping_all
  ) == 0L
) {
  
  stop(
    paste0(
      "No probe-to-AgamP4 mappings could be constructed.\n",
      "Inspect:\n",
      "data_processed/02_mozatlas_probe_annotation.csv"
    ),
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 31. Restrict primary analysis to unambiguous current AgamP4 mappings
# -----------------------------------------------------------------------------

mapping_primary <- mapping_all[
  current_AgamP4 == TRUE &
    n_gene_ids_per_probe == 1L
]


mapping_primary <- unique(
  
  mapping_primary[
    ,
    .(
      probe_id,
      gene_id,
      mapping_method
    )
  ]
)


# Combined mapping may theoretically yield different genes for the same probe.
mapping_primary[
  ,
  current_gene_count :=
    uniqueN(
      gene_id
    ),
  by = probe_id
]


mapping_primary <- mapping_primary[
  current_gene_count == 1L
]


mapping_primary[
  ,
  current_gene_count := NULL
]


n_mapped_probes <- uniqueN(
  mapping_primary$probe_id
)


n_mapped_genes <- uniqueN(
  mapping_primary$gene_id
)


log_step(
  "02",
  paste(
    "Mapped",
    format(
      n_mapped_probes,
      big.mark = ","
    ),
    "probe sets to",
    format(
      n_mapped_genes,
      big.mark = ","
    ),
    "current AgamP4 loci"
  )
)


fwrite(
  mapping_all,
  "data_processed/02_mozatlas_probe_gene_mapping.csv",
  na = "NA"
)


if (
  n_mapped_genes <
  1000L
) {
  
  stop(
    paste0(
      "Only ",
      n_mapped_genes,
      " AgamP4 genes were mapped. ",
      "This is too low for genome-wide MozAtlas feature engineering."
    ),
    call. = FALSE
  )
}


if (
  n_mapped_genes <
  8000L
) {
  
  warning(
    paste0(
      "Only ",
      n_mapped_genes,
      " current AgamP4 genes were mapped from GPL1321. ",
      "Review mapping coverage before model fitting."
    )
  )
}


# -----------------------------------------------------------------------------
# 32. Align sample columns to sample manifest
# -----------------------------------------------------------------------------

sample_manifest[
  ,
  matrix_column :=
    match(
      sample_accession,
      colnames(
        expr
      )
    )
]


if (
  anyNA(
    sample_manifest$matrix_column
  )
) {
  
  stop(
    "One or more sample metadata rows cannot be matched to expression columns.",
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 33. Aggregate expression replicates at probe level
# -----------------------------------------------------------------------------

probe_expression_matrix <- sapply(
  
  expected_groups,
  
  function(grp) {
    
    sample_indices <- sample_manifest[
      group == grp,
      matrix_column
    ]
    
    
    rowMeans(
      expr[
        ,
        sample_indices,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  }
)


if (
  is.null(
    dim(
      probe_expression_matrix
    )
  )
) {
  
  probe_expression_matrix <- matrix(
    
    probe_expression_matrix,
    
    ncol =
      length(
        expected_groups
      )
  )
}


rownames(
  probe_expression_matrix
) <- reference_probe_ids


colnames(
  probe_expression_matrix
) <- paste0(
  "expr_",
  expected_groups
)


probe_expr_dt <- as.data.table(
  
  probe_expression_matrix,
  
  keep.rownames =
    "probe_id"
)


# -----------------------------------------------------------------------------
# 34. Aggregate detection calls at probe level
#
# present_frac_* = fraction of four biological replicates with ABS_CALL == "P".
# -----------------------------------------------------------------------------

present_fraction_matrix <- matrix(
  
  NA_real_,
  
  nrow =
    n_probes,
  
  ncol =
    length(
      expected_groups
    ),
  
  dimnames = list(
    
    reference_probe_ids,
    
    paste0(
      "present_frac_",
      expected_groups
    )
  )
)


if (
  !all(
    is.na(
      detection_call
    )
  )
) {
  
  for (
    j in seq_along(
      expected_groups
    )
  ) {
    
    grp <- expected_groups[[j]]
    
    
    sample_indices <- sample_manifest[
      group == grp,
      matrix_column
    ]
    
    
    present_matrix <- detection_call[
      ,
      sample_indices,
      drop = FALSE
    ] == "P"
    
    
    present_fraction_matrix[
      ,
      j
    ] <- rowMeans(
      present_matrix,
      na.rm = TRUE
    )
    
    
    present_fraction_matrix[
      !is.finite(
        present_fraction_matrix[
          ,
          j
        ]
      ),
      j
    ] <- NA_real_
  }
}


# -----------------------------------------------------------------------------
# 35. Aggregate detection P-values at probe level
# -----------------------------------------------------------------------------

detection_p_matrix_group <- matrix(
  
  NA_real_,
  
  nrow =
    n_probes,
  
  ncol =
    length(
      expected_groups
    ),
  
  dimnames = list(
    
    reference_probe_ids,
    
    paste0(
      "detp_mean_",
      expected_groups
    )
  )
)


if (
  !all(
    is.na(
      detection_p
    )
  )
) {
  
  for (
    j in seq_along(
      expected_groups
    )
  ) {
    
    grp <- expected_groups[[j]]
    
    
    sample_indices <- sample_manifest[
      group == grp,
      matrix_column
    ]
    
    
    detection_p_matrix_group[
      ,
      j
    ] <- rowMeans(
      detection_p[
        ,
        sample_indices,
        drop = FALSE
      ],
      na.rm = TRUE
    )
    
    
    detection_p_matrix_group[
      !is.finite(
        detection_p_matrix_group[
          ,
          j
        ]
      ),
      j
    ] <- NA_real_
  }
}


probe_detection_dt <- data.table(
  
  probe_id =
    reference_probe_ids
)


probe_detection_dt <- cbind(
  
  probe_detection_dt,
  
  as.data.table(
    present_fraction_matrix
  ),
  
  as.data.table(
    detection_p_matrix_group
  )
)


fwrite(
  probe_detection_dt,
  "data_processed/02_mozatlas_probe_detection.csv",
  na = "NA"
)


# -----------------------------------------------------------------------------
# 36. Harmonise reproductive-tissue column names
#
# Female ovary, male testis and male ACPS have no opposite-sex equivalents.
# -----------------------------------------------------------------------------

if (
  "expr_female_ovary" %in%
  names(
    probe_expr_dt
  )
) {
  
  setnames(
    probe_expr_dt,
    "expr_female_ovary",
    "expr_ovary"
  )
}


if (
  "expr_male_testis" %in%
  names(
    probe_expr_dt
  )
) {
  
  setnames(
    probe_expr_dt,
    "expr_male_testis",
    "expr_testis"
  )
}


if (
  "expr_male_accessory_gland" %in%
  names(
    probe_expr_dt
  )
) {
  
  setnames(
    probe_expr_dt,
    "expr_male_accessory_gland",
    "expr_accessory_gland"
  )
}


# Detection columns
if (
  "present_frac_female_ovary" %in%
  names(
    probe_detection_dt
  )
) {
  
  setnames(
    probe_detection_dt,
    "present_frac_female_ovary",
    "present_frac_ovary"
  )
}


if (
  "present_frac_male_testis" %in%
  names(
    probe_detection_dt
  )
) {
  
  setnames(
    probe_detection_dt,
    "present_frac_male_testis",
    "present_frac_testis"
  )
}


if (
  "present_frac_male_accessory_gland" %in%
  names(
    probe_detection_dt
  )
) {
  
  setnames(
    probe_detection_dt,
    "present_frac_male_accessory_gland",
    "present_frac_accessory_gland"
  )
}


if (
  "detp_mean_female_ovary" %in%
  names(
    probe_detection_dt
  )
) {
  
  setnames(
    probe_detection_dt,
    "detp_mean_female_ovary",
    "detp_mean_ovary"
  )
}


if (
  "detp_mean_male_testis" %in%
  names(
    probe_detection_dt
  )
) {
  
  setnames(
    probe_detection_dt,
    "detp_mean_male_testis",
    "detp_mean_testis"
  )
}


if (
  "detp_mean_male_accessory_gland" %in%
  names(
    probe_detection_dt
  )
) {
  
  setnames(
    probe_detection_dt,
    "detp_mean_male_accessory_gland",
    "detp_mean_accessory_gland"
  )
}


# -----------------------------------------------------------------------------
# 37. Merge expression, detection and gene mapping at probe level
# -----------------------------------------------------------------------------

probe_dt <- merge(
  
  probe_expr_dt,
  
  probe_detection_dt,
  
  by =
    "probe_id",
  
  all.x =
    TRUE
)


probe_dt <- merge(
  
  probe_dt,
  
  mapping_primary,
  
  by =
    "probe_id",
  
  all =
    FALSE
)


setorder(
  probe_dt,
  gene_id,
  probe_id
)


fwrite(
  probe_dt,
  "data_processed/02_mozatlas_probe_expression.csv",
  na = "NA"
)


# -----------------------------------------------------------------------------
# 38. Collapse probes to gene level
#
# Multiple unambiguously mapped probes for the same gene are averaged.
# -----------------------------------------------------------------------------

expression_columns <- grep(
  
  "^expr_",
  
  names(
    probe_dt
  ),
  
  value = TRUE
)


present_columns <- grep(
  
  "^present_frac_",
  
  names(
    probe_dt
  ),
  
  value = TRUE
)


detp_columns <- grep(
  
  "^detp_mean_",
  
  names(
    probe_dt
  ),
  
  value = TRUE
)


aggregate_columns <- c(
  
  expression_columns,
  
  present_columns,
  
  detp_columns
)


gene_expr <- probe_dt[
  ,
  c(
    
    list(
      
      n_probes_mapped =
        uniqueN(
          probe_id
        )
    ),
    
    lapply(
      
      .SD,
      
      function(x) {
        
        result <- mean(
          x,
          na.rm = TRUE
        )
        
        
        if (
          !is.finite(
            result
          )
        ) {
          
          return(
            NA_real_
          )
        }
        
        
        result
      }
    )
  ),
  
  by =
    gene_id,
  
  .SDcols =
    aggregate_columns
]


setorder(
  gene_expr,
  gene_id
)


# -----------------------------------------------------------------------------
# 39. Safe Tau tissue-specificity function
#
# Tau approximately ranges:
#
#   0 = broadly expressed
#   1 = strongly tissue-specific
#
# Calculated using the GEO-supplied processed values AS RETRIEVED.
# -----------------------------------------------------------------------------

tau_index_safe <- function(x) {
  
  x <- as.numeric(
    x
  )
  
  
  x <- x[
    is.finite(
      x
    )
  ]
  
  
  if (
    length(
      x
    ) < 2L
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  xmax <- max(
    x
  )
  
  
  if (
    !is.finite(
      xmax
    ) ||
    xmax <= 0
  ) {
    
    return(
      NA_real_
    )
  }
  
  
  sum(
    1 -
      x /
      xmax
  ) /
    (
      length(
        x
      ) -
        1L
    )
}


# -----------------------------------------------------------------------------
# 40. Calculate tissue-specificity features
# -----------------------------------------------------------------------------

expression_columns <- grep(
  
  "^expr_",
  
  names(
    gene_expr
  ),
  
  value = TRUE
)


gene_expr[
  ,
  tau_tissue_specificity :=
    apply(
      .SD,
      1,
      tau_index_safe
    ),
  .SDcols =
    expression_columns
]


# -----------------------------------------------------------------------------
# 41. Expression range
# -----------------------------------------------------------------------------

gene_expr[
  ,
  expression_range :=
    apply(
      
      .SD,
      
      1,
      
      function(x) {
        
        x <- as.numeric(
          x
        )
        
        
        x <- x[
          is.finite(
            x
          )
        ]
        
        
        if (
          length(
            x
          ) == 0L
        ) {
          
          return(
            NA_real_
          )
        }
        
        
        max(
          x
        ) -
          min(
            x
          )
      }
    ),
  .SDcols =
    expression_columns
]


# -----------------------------------------------------------------------------
# 42. Identify tissue/group with maximum expression
# -----------------------------------------------------------------------------

gene_expr[
  ,
  max_expression_group :=
    apply(
      
      .SD,
      
      1,
      
      function(x) {
        
        x <- as.numeric(
          x
        )
        
        
        if (
          all(
            !is.finite(
              x
            )
          )
        ) {
          
          return(
            NA_character_
          )
        }
        
        
        idx <- which.max(
          ifelse(
            is.finite(
              x
            ),
            x,
            -Inf
          )
        )
        
        
        sub(
          "^expr_",
          "",
          expression_columns[
            idx
          ]
        )
      }
    ),
  .SDcols =
    expression_columns
]


gene_expr[
  ,
  max_expression_value :=
    apply(
      
      .SD,
      
      1,
      
      function(x) {
        
        x <- as.numeric(
          x
        )
        
        
        x <- x[
          is.finite(
            x
          )
        ]
        
        
        if (
          length(
            x
          ) == 0L
        ) {
          
          return(
            NA_real_
          )
        }
        
        
        max(
          x
        )
      }
    ),
  .SDcols =
    expression_columns
]


# -----------------------------------------------------------------------------
# 43. Matched female-versus-male tissue features
#
# Only tissues represented in BOTH sexes are compared.
# -----------------------------------------------------------------------------

matched_tissues <- c(
  "body",
  "carcass",
  "head",
  "malpighian",
  "midgut",
  "salivary"
)


matched_tissues <- matched_tissues[
  
  paste0(
    "expr_female_",
    matched_tissues
  ) %in%
    names(
      gene_expr
    ) &
    
    paste0(
      "expr_male_",
      matched_tissues
    ) %in%
    names(
      gene_expr
    )
]


female_cols <- paste0(
  "expr_female_",
  matched_tissues
)


male_cols <- paste0(
  "expr_male_",
  matched_tissues
)


if (
  length(
    female_cols
  ) > 0L &&
  length(
    male_cols
  ) > 0L
) {
  
  gene_expr[
    ,
    female_mean_matched :=
      rowMeans(
        .SD,
        na.rm = TRUE
      ),
    .SDcols =
      female_cols
  ]
  
  
  gene_expr[
    ,
    male_mean_matched :=
      rowMeans(
        .SD,
        na.rm = TRUE
      ),
    .SDcols =
      male_cols
  ]
  
  
  gene_expr[
    ,
    female_male_delta :=
      female_mean_matched -
      male_mean_matched
  ]
}


# -----------------------------------------------------------------------------
# 44. Ovary-specificity feature
# -----------------------------------------------------------------------------

if (
  "expr_ovary" %in%
  names(
    gene_expr
  )
) {
  
  female_nonovary <- paste0(
    "expr_female_",
    matched_tissues
  )
  
  
  female_nonovary <- female_nonovary[
    female_nonovary %in%
      names(
        gene_expr
      )
  ]
  
  
  if (
    length(
      female_nonovary
    ) > 0L
  ) {
    
    gene_expr[
      ,
      female_nonovary_mean :=
        rowMeans(
          .SD,
          na.rm = TRUE
        ),
      .SDcols =
        female_nonovary
    ]
    
    
    gene_expr[
      ,
      ovary_specificity_delta :=
        expr_ovary -
        female_nonovary_mean
    ]
  }
}


# -----------------------------------------------------------------------------
# 45. Testis-specificity feature
# -----------------------------------------------------------------------------

if (
  "expr_testis" %in%
  names(
    gene_expr
  )
) {
  
  male_nontestis <- c(
    
    paste0(
      "expr_male_",
      matched_tissues
    ),
    
    "expr_accessory_gland"
  )
  
  
  male_nontestis <- male_nontestis[
    male_nontestis %in%
      names(
        gene_expr
      )
  ]
  
  
  if (
    length(
      male_nontestis
    ) > 0L
  ) {
    
    gene_expr[
      ,
      male_nontestis_mean :=
        rowMeans(
          .SD,
          na.rm = TRUE
        ),
      .SDcols =
        male_nontestis
    ]
    
    
    gene_expr[
      ,
      testis_specificity_delta :=
        expr_testis -
        male_nontestis_mean
    ]
  }
}


# -----------------------------------------------------------------------------
# 46. Male accessory-gland specificity
# -----------------------------------------------------------------------------

if (
  "expr_accessory_gland" %in%
  names(
    gene_expr
  )
) {
  
  male_nonaccessory <- c(
    
    paste0(
      "expr_male_",
      matched_tissues
    ),
    
    "expr_testis"
  )
  
  
  male_nonaccessory <- male_nonaccessory[
    male_nonaccessory %in%
      names(
        gene_expr
      )
  ]
  
  
  if (
    length(
      male_nonaccessory
    ) > 0L
  ) {
    
    gene_expr[
      ,
      male_nonaccessory_mean :=
        rowMeans(
          .SD,
          na.rm = TRUE
        ),
      .SDcols =
        male_nonaccessory
    ]
    
    
    gene_expr[
      ,
      accessory_gland_specificity_delta :=
        expr_accessory_gland -
        male_nonaccessory_mean
    ]
  }
}


# -----------------------------------------------------------------------------
# 47. Maximum reproductive-tissue expression
# -----------------------------------------------------------------------------

reproductive_cols <- intersect(
  
  c(
    "expr_ovary",
    "expr_testis",
    "expr_accessory_gland"
  ),
  
  names(
    gene_expr
  )
)


if (
  length(
    reproductive_cols
  ) > 0L
) {
  
  gene_expr[
    ,
    reproductive_expression_max :=
      apply(
        
        .SD,
        
        1,
        
        function(x) {
          
          x <- as.numeric(
            x
          )
          
          
          x <- x[
            is.finite(
              x
            )
          ]
          
          
          if (
            length(
              x
            ) == 0L
          ) {
            
            return(
              NA_real_
            )
          }
          
          
          max(
            x
          )
        }
      ),
    .SDcols =
      reproductive_cols
  ]
}


# -----------------------------------------------------------------------------
# 48. Detection-breadth features
#
# present fraction >= 0.75 means at least 3 of 4 biological replicates were
# called Present by MAS 5.0.
# -----------------------------------------------------------------------------

present_columns <- grep(
  
  "^present_frac_",
  
  names(
    gene_expr
  ),
  
  value = TRUE
)


if (
  length(
    present_columns
  ) > 0L
) {
  
  gene_expr[
    ,
    detection_breadth_3of4 :=
      apply(
        
        .SD,
        
        1,
        
        function(x) {
          
          x <- as.numeric(
            x
          )
          
          
          if (
            all(
              is.na(
                x
              )
            )
          ) {
            
            return(
              NA_integer_
            )
          }
          
          
          sum(
            x >=
              0.75,
            na.rm = TRUE
          )
        }
      ),
    .SDcols =
      present_columns
  ]
}


# -----------------------------------------------------------------------------
# 49. Current AgamP4 membership sanity check
# -----------------------------------------------------------------------------

gene_expr[
  ,
  current_AgamP4 :=
    gene_id %in%
    current_gene_ids
]


n_noncurrent <- gene_expr[
  current_AgamP4 == FALSE,
  .N
]


if (
  n_noncurrent >
  0L
) {
  
  warning(
    paste(
      n_noncurrent,
      "MozAtlas gene IDs are not present in the current Step 01 universe."
    )
  )
}


gene_expr <- gene_expr[
  current_AgamP4 == TRUE
]


gene_expr[
  ,
  current_AgamP4 := NULL
]


setorder(
  gene_expr,
  gene_id
)


# -----------------------------------------------------------------------------
# 50. Save final gene-level expression feature table
# -----------------------------------------------------------------------------

fwrite(
  gene_expr,
  "data_processed/02_mozatlas_gene_expression.csv",
  na = "NA"
)


# -----------------------------------------------------------------------------
# 51. Detection-call QC
# -----------------------------------------------------------------------------

detection_values <- as.vector(
  detection_call
)


detection_values <- detection_values[
  !is.na(
    detection_values
  ) &
    nzchar(
      detection_values
    )
]


if (
  length(
    detection_values
  ) > 0L
) {
  
  detection_qc <- as.data.table(
    table(
      detection_values
    )
  )
  
  
  setnames(
    detection_qc,
    c(
      "detection_values",
      "N"
    ),
    c(
      "ABS_CALL",
      "count"
    )
  )
  
  
  detection_qc[
    ,
    proportion :=
      count /
      sum(
        count
      )
  ]
  
} else {
  
  detection_qc <- data.table(
    
    ABS_CALL =
      character(),
    
    count =
      integer(),
    
    proportion =
      numeric()
  )
}


fwrite(
  detection_qc,
  "data_processed/02_mozatlas_detection_qc.csv",
  na = "NA"
)


# -----------------------------------------------------------------------------
# 52. Main QC summary
# -----------------------------------------------------------------------------

step01_gene_count <- uniqueN(
  current_genes$gene_id
)


expression_coverage_percent <- round(
  
  100 *
    nrow(
      gene_expr
    ) /
    step01_gene_count,
  
  2
)


qc <- data.table(
  
  metric = c(
    
    "GEO samples",
    
    "Expected tissue/sex groups",
    
    "Observed tissue/sex groups",
    
    "Replicates per group",
    
    "GPL1321 probe sets in expression tables",
    
    "GPL1321 platform rows retained",
    
    "Primary mapped probe sets",
    
    "Current AgamP4 loci with MozAtlas mapping",
    
    "Step 01 total AgamP4 loci",
    
    "AgamP4 expression coverage percent",
    
    "Expression missing values",
    
    "ABS_CALL values available",
    
    "Detection P-values available"
  ),
  
  value = c(
    
    n_samples,
    
    length(
      expected_groups
    ),
    
    uniqueN(
      sample_manifest$group
    ),
    
    4L,
    
    n_probes,
    
    nrow(
      platform_dt
    ),
    
    n_mapped_probes,
    
    nrow(
      gene_expr
    ),
    
    step01_gene_count,
    
    expression_coverage_percent,
    
    sum(
      is.na(
        expr
      )
    ),
    
    as.integer(
      !all(
        is.na(
          detection_call
        )
      )
    ),
    
    as.integer(
      !all(
        is.na(
          detection_p
        )
      )
    )
  )
)


fwrite(
  qc,
  "data_processed/02_mozatlas_expression_qc.csv"
)


# -----------------------------------------------------------------------------
# 53. Save provenance
# -----------------------------------------------------------------------------

provenance <- data.table(
  
  source_database =
    "NCBI Gene Expression Omnibus",
  
  GEO_series =
    GEO_SERIES,
  
  GEO_platform =
    GEO_PLATFORM,
  
  species =
    SPECIES,
  
  experiment_type =
    "Expression profiling by array",
  
  biological_design =
    "15 adult tissue/sex groups; 4 biological replicates per group",
  
  total_samples =
    n_samples,
  
  total_probe_sets =
    n_probes,
  
  acquisition_endpoint =
    "www.ncbi.nlm.nih.gov/geo/query/acc.cgi",
  
  family_SOFT_source =
    family_soft_source,
  
  expression_measure =
    "GEO-supplied normalized GC-RMA signal intensity",
  
  MosqEditR_expression_transformation =
    "None",
  
  replicate_aggregation =
    "Arithmetic mean across four biological replicates within each group",
  
  ABS_CALL_use =
    paste0(
      "MAS 5.0 detection calls retained; present fraction calculated ",
      "within each tissue/sex group"
    ),
  
  detection_P_value_use =
    "Mean detection P-value calculated within each tissue/sex group",
  
  probe_mapping_strategy =
    paste0(
      "Direct AGAP stable identifiers from GPL1321; conservative unique ",
      "Anopheles gene-symbol fallback; ambiguous mappings excluded"
    ),
  
  multiple_probe_gene_aggregation =
    "Arithmetic mean across uniquely mapped probe sets",
  
  ACPS_interpretation =
    "Male accessory-gland/Acps biological samples; retained",
  
  tissue_specificity_metric =
    paste0(
      "Tau calculated on GEO processed expression values as retrieved; ",
      "no additional scale transformation"
    ),
  
  mapped_current_AgamP4_loci =
    nrow(
      gene_expr
    ),
  
  current_AgamP4_expression_coverage_percent =
    expression_coverage_percent,
  
  retrieval_date =
    as.character(
      Sys.Date()
    )
)


fwrite(
  provenance,
  "data_processed/02_mozatlas_expression_provenance.csv",
  na = "NA"
)


# -----------------------------------------------------------------------------
# 54. Save session information
# -----------------------------------------------------------------------------

capture.output(
  
  sessionInfo(),
  
  file =
    "logs/02_sessionInfo.txt"
)


# -----------------------------------------------------------------------------
# 55. Checksums
# -----------------------------------------------------------------------------

checksum_files <- c(
  
  geo_test_file,
  
  family_soft_file,
  
  "data_processed/02_mozatlas_sample_manifest.csv",
  
  "data_processed/02_mozatlas_probe_annotation.csv",
  
  "data_processed/02_mozatlas_probe_gene_mapping.csv",
  
  "data_processed/02_mozatlas_probe_expression.csv",
  
  "data_processed/02_mozatlas_probe_detection.csv",
  
  "data_processed/02_mozatlas_gene_expression.csv",
  
  "data_processed/02_mozatlas_expression_qc.csv",
  
  "data_processed/02_mozatlas_detection_qc.csv",
  
  "data_processed/02_mozatlas_expression_provenance.csv",
  
  "logs/02_expression_scale.tsv",
  
  "logs/02_sessionInfo.txt"
)


# Add separately downloaded GPL file if fallback was required.
possible_gpl_file <- file.path(
  "data_raw/GEO",
  paste0(
    GEO_PLATFORM,
    "_full.soft"
  )
)


if (
  file.exists(
    possible_gpl_file
  )
) {
  
  checksum_files <- c(
    
    checksum_files,
    
    possible_gpl_file
  )
}


checksum_files <- unique(
  checksum_files[
    file.exists(
      checksum_files
    )
  ]
)


write_checksum(
  checksum_files,
  "logs/02_checksums.tsv"
)


# -----------------------------------------------------------------------------
# 56. Console validation
# -----------------------------------------------------------------------------

cat(
  "\n",
  "============================================================\n",
  "MOSQEDIT-R STEP 02 COMPLETED SUCCESSFULLY\n",
  "============================================================\n",
  "GEO Series:                  ", GEO_SERIES, "\n",
  "Platform:                    ", GEO_PLATFORM, "\n",
  "Species:                     ", SPECIES, "\n",
  "Samples:                     ",
  format(
    n_samples,
    big.mark = ","
  ),
  "\n",
  "Tissue/sex groups:           ",
  uniqueN(
    sample_manifest$group
  ),
  "\n",
  "Replicates per group:        4\n",
  "Array probe sets:            ",
  format(
    n_probes,
    big.mark = ","
  ),
  "\n",
  "Mapped probe sets:           ",
  format(
    n_mapped_probes,
    big.mark = ","
  ),
  "\n",
  "Mapped current AgamP4 loci:  ",
  format(
    nrow(
      gene_expr
    ),
    big.mark = ","
  ),
  "\n",
  "AgamP4 expression coverage:  ",
  expression_coverage_percent,
  "%\n",
  "Expression transformation:   NONE\n",
  "ACPS samples:                RETAINED\n",
  "ABS_CALL data:               ",
  ifelse(
    all(
      is.na(
        detection_call
      )
    ),
    "NOT AVAILABLE",
    "AVAILABLE"
  ),
  "\n",
  "Detection P-values:          ",
  ifelse(
    all(
      is.na(
        detection_p
      )
    ),
    "NOT AVAILABLE",
    "AVAILABLE"
  ),
  "\n",
  "============================================================\n",
  sep = ""
)


cat(
  "\nMozAtlas sample groups:\n"
)


print(
  group_counts
)


cat(
  "\nExpression quantiles as retrieved from GEO:\n"
)


print(
  expr_quantiles
)


cat(
  "\nProbe-to-gene mapping summary:\n"
)


print(
  mapping_all[
    ,
    .(
      probe_gene_links =
        .N,
      
      unique_probes =
        uniqueN(
          probe_id
        ),
      
      unique_genes =
        uniqueN(
          gene_id
        )
    ),
    by =
      mapping_method
  ]
)


cat(
  "\nDetection-call summary:\n"
)


print(
  detection_qc
)


cat(
  "\nFirst 10 gene-level MozAtlas records:\n"
)


print(
  gene_expr[
    1:min(
      10L,
      nrow(
        gene_expr
      )
    )
  ]
)


cat(
  "\nGenes by tissue/group of maximum expression:\n"
)


print(
  gene_expr[
    ,
    .N,
    by =
      max_expression_group
  ][
    order(
      -N
    )
  ]
)


cat(
  "\nQC summary:\n"
)


print(
  qc
)


cat(
  "\nOutput files:\n",
  "  ",
  family_soft_file,
  "\n",
  "  data_processed/02_mozatlas_sample_manifest.csv\n",
  "  data_processed/02_mozatlas_probe_annotation.csv\n",
  "  data_processed/02_mozatlas_probe_gene_mapping.csv\n",
  "  data_processed/02_mozatlas_probe_expression.csv\n",
  "  data_processed/02_mozatlas_probe_detection.csv\n",
  "  data_processed/02_mozatlas_gene_expression.csv\n",
  "  data_processed/02_mozatlas_expression_qc.csv\n",
  "  data_processed/02_mozatlas_detection_qc.csv\n",
  "  data_processed/02_mozatlas_expression_provenance.csv\n",
  "  logs/02_expression_scale.tsv\n",
  "  logs/02_sessionInfo.txt\n",
  "  logs/02_checksums.tsv\n",
  sep = ""
)


log_step(
  "02",
  "MozAtlas expression acquisition and feature engineering completed successfully"
)
