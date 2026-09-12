log_step <- function(step, msg) {
  dir.create("logs", showWarnings=FALSE, recursive=TRUE)
  line <- sprintf("%s\t%s\t%s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), step, msg)
  cat(line, "\n")
  cat(line, "\n", file="logs/pipeline.log", append=TRUE)
}
first_nonempty <- function(x) { x <- x[!is.na(x) & nzchar(x)]; if (length(x)) x[[1]] else NA_character_ }
tau_index <- function(x) {
  x <- as.numeric(x); x <- x[is.finite(x)]
  if (length(x) < 2 || max(x, na.rm=TRUE) <= 0) return(NA_real_)
  z <- x/max(x, na.rm=TRUE)
  sum(1-z)/(length(z)-1)
}
write_checksum <- function(files, outfile) {
  files <- files[file.exists(files)]
  if (!length(files)) return(invisible(NULL))
  h <- data.frame(file=files, sha256=vapply(files, digest::digest, character(1), algo="sha256", file=TRUE))
  write.table(h, outfile, sep="\t", row.names=FALSE, quote=FALSE)
}

