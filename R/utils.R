#' Tissue-specificity Tau index
#'
#' Compute the Tau tissue-specificity index from a numeric expression vector.
#' Higher values indicate greater concentration of expression in a subset of
#' tissues. Missing and non-finite values are removed.
#'
#' @param x Numeric vector of expression values.
#' @return Numeric scalar in approximately [0, 1], or NA when not estimable.
#' @export
mosqedit_tau <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) < 2L || max(x, na.rm = TRUE) <= 0) return(NA_real_)
  z <- x / max(x, na.rm = TRUE)
  sum(1 - z) / (length(z) - 1)
}

#' Rank-normalize values to percentiles
#'
#' @param x Numeric vector.
#' @param higher_is_better Logical; if FALSE, smaller values receive larger
#'   percentiles.
#' @return Numeric vector in [0, 1] with NA preserved for non-finite values.
#' @export
mosqedit_percentile <- function(x, higher_is_better = TRUE) {
  x <- suppressWarnings(as.numeric(x))
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n_ok <- sum(ok)
  if (n_ok == 0L) return(out)
  if (n_ok == 1L) {
    out[ok] <- 0.5
    return(out)
  }
  z <- x[ok]
  if (!isTRUE(higher_is_better)) z <- -z
  out[ok] <- (rank(z, ties.method = "average") - 1) / (n_ok - 1)
  out
}

#' Deterministic MosqEdit-R ranking
#'
#' Ranks finite scores from highest to lowest, using gene identifiers as a
#' deterministic tie-breaker.
#'
#' @param score Numeric vector of scores.
#' @param gene_id Character vector of unique identifiers of the same length.
#' @return Integer rank vector; non-finite scores receive NA.
#' @export
mosqedit_rank <- function(score, gene_id) {
  score <- suppressWarnings(as.numeric(score))
  gene_id <- as.character(gene_id)
  if (length(score) != length(gene_id)) {
    stop("score and gene_id must have the same length.", call. = FALSE)
  }
  if (anyDuplicated(gene_id)) {
    stop("gene_id must be unique for deterministic ranking.", call. = FALSE)
  }
  idx <- seq_along(score)
  sortable <- ifelse(is.finite(score), score, -Inf)
  ord <- order(-sortable, gene_id, na.last = TRUE)
  out <- rep(NA_integer_, length(score))
  out[ord] <- seq_along(ord)
  out[!is.finite(score)] <- NA_integer_
  out
}

#' Integrate biology, PU support, and CRISPR tractability
#'
#' Implements the primary MosqEdit-R integration rule. Missing domains are not
#' interpreted as negative evidence: weights are renormalized over observed
#' domains and a modest completeness adjustment is then applied.
#'
#' @param biology Numeric biology-evidence percentile.
#' @param pu Numeric positive-unlabeled support percentile.
#' @param crispr Numeric CRISPR tractability percentile.
#' @param weights Named numeric vector with elements biology, pu, crispr.
#' @param completeness_floor Lower multiplicative floor for incomplete evidence.
#' @return A data.frame with raw_score, adjusted_score, domain_coverage, and
#'   observed_weight.
#' @export
mosqedit_integrated_score <- function(
    biology,
    pu,
    crispr,
    weights = c(biology = 0.50, pu = 0.30, crispr = 0.20),
    completeness_floor = 0.80) {

  if (!all(c("biology", "pu", "crispr") %in% names(weights))) {
    stop("weights must be named biology, pu, and crispr.", call. = FALSE)
  }
  weights <- as.numeric(weights[c("biology", "pu", "crispr")])
  if (any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("weights must be finite, non-negative, and have positive sum.", call. = FALSE)
  }
  if (!is.finite(completeness_floor) || completeness_floor < 0 || completeness_floor > 1) {
    stop("completeness_floor must lie in [0, 1].", call. = FALSE)
  }

  n <- max(length(biology), length(pu), length(crispr))
  biology <- rep_len(suppressWarnings(as.numeric(biology)), n)
  pu <- rep_len(suppressWarnings(as.numeric(pu)), n)
  crispr <- rep_len(suppressWarnings(as.numeric(crispr)), n)
  components <- cbind(biology = biology, pu = pu, crispr = crispr)
  observed <- is.finite(components)
  wmat <- matrix(weights, nrow = n, ncol = 3L, byrow = TRUE)
  numerator <- rowSums(ifelse(observed, components, 0) * wmat)
  observed_weight <- rowSums(observed * wmat)
  raw <- numerator / observed_weight
  raw[observed_weight <= 0] <- NA_real_
  coverage <- rowMeans(observed)
  adjusted <- raw * (completeness_floor + (1 - completeness_floor) * coverage)

  data.frame(
    raw_score = raw,
    adjusted_score = adjusted,
    domain_coverage = coverage,
    observed_weight = observed_weight
  )
}

# Internal helpers used by copied manuscript scripts.
.mosqedit_first_nonempty <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x)) x[[1]] else NA_character_
}

.mosqedit_log_step <- function(step, msg, log_dir = "logs") {
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  line <- sprintf(
    "%s\t%s\t%s",
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    step,
    msg
  )
  cat(line, "\n")
  cat(line, "\n", file = file.path(log_dir, "pipeline.log"), append = TRUE)
  invisible(line)
}

.mosqedit_write_checksum <- function(files, outfile) {
  files <- files[file.exists(files)]
  if (!length(files)) return(invisible(NULL))
  hashes <- tools::md5sum(files)
  # Prefer SHA-256 when digest is installed; retain a working base-R fallback.
  if (requireNamespace("digest", quietly = TRUE)) {
    hashes <- vapply(
      files,
      digest::digest,
      character(1),
      algo = "sha256",
      file = TRUE
    )
  }
  tab <- data.frame(file = files, checksum = unname(hashes), stringsAsFactors = FALSE)
  utils::write.table(tab, outfile, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(tab)
}

