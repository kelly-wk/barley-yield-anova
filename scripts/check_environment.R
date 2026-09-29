#!/usr/bin/env Rscript

minimum_r <- numeric_version("4.2.0")
if (getRversion() < minimum_r) {
  stop("R >= 4.2.0 is required; found ", getRversion(), call. = FALSE)
}
if (!requireNamespace("lattice", quietly = TRUE)) {
  stop("The recommended R package 'lattice' is required.", call. = FALSE)
}

cat("R:", R.version.string, "\n")
cat("lattice:", as.character(utils::packageVersion("lattice")), "\n")
cat("Environment check passed.\n")
