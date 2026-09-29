#!/usr/bin/env Rscript

parse_output <- function(args) {
  value <- sub("^--output=", "", args[grepl("^--output=", args)])
  if (length(value) == 0L) "data/barley.csv" else value[[1L]]
}

output_path <- parse_output(commandArgs(trailingOnly = TRUE))
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("lattice", quietly = TRUE)) {
  stop("The recommended R package 'lattice' is required.", call. = FALSE)
}

source(file.path("R", "analysis_helpers.R"), local = TRUE)
environment <- new.env(parent = emptyenv())
utils::data("barley", package = "lattice", envir = environment)
data <- get("barley", envir = environment)
data <- coerce_barley(data)
validate_barley(data)

data <- data[order(data$year, data$variety, data$site), ]
data_out <- transform(
  data,
  variety = as.character(variety),
  year = as.character(year),
  site = as.character(site)
)
utils::write.csv(data_out, output_path, row.names = FALSE, na = "")

checksum <- unname(tools::md5sum(output_path))
writeLines(paste(checksum, basename(output_path)), paste0(output_path, ".md5"))
cat("Wrote", nrow(data_out), "rows to", output_path, "\n")
cat("MD5:", checksum, "\n")
