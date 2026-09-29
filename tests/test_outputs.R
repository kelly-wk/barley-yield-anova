#!/usr/bin/env Rscript

required_files <- c(
  "REPORT.md",
  "data/barley.csv",
  "data/barley.csv.md5",
  "results/anova.csv",
  "results/selected_model_anova.csv",
  "results/marginal_means.csv",
  "results/pairwise_holm.csv",
  "results/selected_model_marginal_means.csv",
  "results/selected_model_pairwise_holm.csv",
  "results/year_by_site_contrasts.csv",
  "results/effect_sizes.csv",
  "results/model_comparison.csv",
  "results/diagnostics.csv",
  "results/selected_model_diagnostics.csv",
  "results/influence.csv",
  "results/coefficients_ols_hc3.csv",
  "results/year_effect_robustness.csv",
  "results/cell_predictions.csv",
  "results/sensitivity_morris_swap.csv",
  "results/model_selection_sensitivity.csv",
  "results/residuals.csv",
  "results/summary.json",
  "results/session_info.txt",
  "figures/model_diagnostics.png",
  "figures/marginal_means.png",
  "figures/variety_site_interaction.png",
  "figures/year_site_interaction.png"
)

missing <- required_files[!file.exists(required_files)]
if (length(missing) > 0L) {
  stop("Missing generated files: ", paste(missing, collapse = ", "), call. = FALSE)
}
stopifnot(all(file.info(required_files)$size > 0L))

anova_table <- utils::read.csv("results/anova.csv", check.names = FALSE)
pairwise <- utils::read.csv("results/pairwise_holm.csv", check.names = FALSE)
comparison <- utils::read.csv("results/model_comparison.csv", check.names = FALSE)
sensitivity <- utils::read.csv("results/sensitivity_morris_swap.csv", check.names = FALSE)

stopifnot(nrow(anova_table) == 5L)
stopifnot(nrow(pairwise) == 61L) # 1 year + 45 variety + 15 site comparisons
stopifnot(nrow(comparison) == 12L)
stopifnot(nrow(sensitivity) == 2L)
stopifnot(all(c("published_lattice", "morris_years_swapped") %in% sensitivity$scenario))
stopifnot(comparison$AICc_rank[1L] == 1L)
stopifnot(all(diff(comparison$AICc) >= 0))

report <- paste(readLines("REPORT.md", warn = FALSE), collapse = "\n")
stopifnot(grepl("Holm", report, fixed = TRUE))
stopifnot(grepl("Morris", report, fixed = TRUE))
stopifnot(!grepl("\\b[0-9]{8}\\b", report, perl = TRUE))

cat("Output tests passed.\n")
