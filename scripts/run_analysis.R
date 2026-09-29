#!/usr/bin/env Rscript

script_arguments <- commandArgs(trailingOnly = TRUE)

argument_value <- function(name, default) {
  prefix <- paste0("--", name, "=")
  matches <- script_arguments[startsWith(script_arguments, prefix)]
  if (length(matches) == 0L) default else sub(prefix, "", matches[[1L]], fixed = TRUE)
}

data_path <- argument_value("data", "data/barley.csv")
results_dir <- argument_value("results", "results")
figures_dir <- argument_value("figures", "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("R", "analysis_helpers.R"), local = TRUE)

write_csv <- function(data, filename) {
  utils::write.csv(data, file.path(results_dir, filename), row.names = FALSE, na = "")
}

markdown_table <- function(data, digits = 3L) {
  formatted <- data
  for (name in names(formatted)) {
    if (is.numeric(formatted[[name]])) {
      formatted[[name]] <- ifelse(
        is.na(formatted[[name]]), "",
        formatC(formatted[[name]], digits = digits, format = "f")
      )
    }
    formatted[[name]] <- gsub("\\|", "\\\\|", as.character(formatted[[name]]))
  }
  header <- paste0("| ", paste(names(formatted), collapse = " | "), " |")
  rule <- paste0("| ", paste(rep("---", ncol(formatted)), collapse = " | "), " |")
  rows <- apply(formatted, 1L, function(row) {
    paste0("| ", paste(row, collapse = " | "), " |")
  })
  paste(c(header, rule, rows), collapse = "\n")
}

data <- read_barley(data_path)
models <- candidate_models(data)
model <- primary_model(data)

comparison <- model_comparison(models)
anova_table <- anova_effects(model)
diagnostics <- diagnostic_summary(model)
influence <- influence_table(model, data)

# The year-by-site model is the strongest small-sample information-criterion
# candidate. It is reported as a model-selection sensitivity analysis, while
# the source-motivated variety-by-site model remains the planned test above.
selected_model <- models$year_by_site
selected_anova <- anova_effects(selected_model)
selected_diagnostics <- diagnostic_summary(selected_model)

factor_names <- c("year", "variety", "site")
emm_objects <- lapply(factor_names, function(name) {
  estimated_marginal_means(model, name)
})
names(emm_objects) <- factor_names
marginal_means <- do.call(rbind, lapply(emm_objects, `[[`, "table"))
pairwise <- do.call(rbind, lapply(factor_names, function(name) {
  pairwise_from_emm(model, emm_objects[[name]], name)
}))
rownames(marginal_means) <- NULL
rownames(pairwise) <- NULL

selected_emm_objects <- lapply(factor_names, function(name) {
  estimated_marginal_means(selected_model, name)
})
names(selected_emm_objects) <- factor_names
selected_marginal_means <- do.call(rbind, lapply(selected_emm_objects, `[[`, "table"))
selected_pairwise <- do.call(rbind, lapply(factor_names, function(name) {
  pairwise_from_emm(selected_model, selected_emm_objects[[name]], name)
}))
selected_year_site <- year_contrasts_by_site(selected_model)
rownames(selected_marginal_means) <- NULL
rownames(selected_pairwise) <- NULL

year_ols <- pairwise[pairwise$factor == "year", ]
year_hc3_emm <- estimated_marginal_means(model, "year", covariance = hc3_vcov(model))
year_hc3 <- pairwise_from_emm(model, year_hc3_emm, "year")
year_robustness <- rbind(
  data.frame(method = "OLS covariance", year_ols[, c(
    "contrast", "estimate", "standard_error", "df", "t_value", "p_value",
    "conf_low", "conf_high", "standardized_effect"
  )]),
  data.frame(method = "HC3 covariance", year_hc3[, c(
    "contrast", "estimate", "standard_error", "df", "t_value", "p_value",
    "conf_low", "conf_high", "standardized_effect"
  )])
)

coefficient_names <- names(stats::coef(model))
ols_se <- sqrt(diag(stats::vcov(model)))
hc3_se <- sqrt(diag(hc3_vcov(model)))
coefficients <- data.frame(
  term = coefficient_names,
  estimate = stats::coef(model),
  OLS_standard_error = ols_se,
  OLS_t = stats::coef(model) / ols_se,
  OLS_p = 2 * stats::pt(abs(stats::coef(model) / ols_se),
                       df = stats::df.residual(model), lower.tail = FALSE),
  HC3_standard_error = hc3_se,
  HC3_t = stats::coef(model) / hc3_se,
  HC3_p = 2 * stats::pt(abs(stats::coef(model) / hc3_se),
                       df = stats::df.residual(model), lower.tail = FALSE),
  stringsAsFactors = FALSE,
  row.names = NULL
)

prediction_grid <- expand.grid(
  year = levels(data$year),
  variety = levels(data$variety),
  site = levels(data$site),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
for (name in factor_names) {
  prediction_grid[[name]] <- factor(prediction_grid[[name]], levels = levels(data[[name]]))
}
prediction_interval <- stats::predict(model, newdata = prediction_grid,
                                      interval = "confidence", level = 0.95)
cell_predictions <- cbind(
  transform(
    prediction_grid,
    year = as.character(year), variety = as.character(variety), site = as.character(site)
  ),
  data.frame(
    estimate = prediction_interval[, "fit"],
    conf_low = prediction_interval[, "lwr"],
    conf_high = prediction_interval[, "upr"]
  )
)

sensitivity_row <- function(scenario, scenario_data) {
  scenario_model <- primary_model(scenario_data)
  scenario_anova <- anova_effects(scenario_model)
  scenario_emm <- estimated_marginal_means(scenario_model, "year")
  scenario_year <- pairwise_from_emm(scenario_model, scenario_emm, "year")
  scenario_diag <- diagnostic_summary(scenario_model)
  interaction_row <- scenario_anova[scenario_anova$term == "variety:site", ]
  data.frame(
    scenario = scenario,
    year_contrast = scenario_year$contrast,
    year_estimate = scenario_year$estimate,
    year_conf_low = scenario_year$conf_low,
    year_conf_high = scenario_year$conf_high,
    year_p_value = scenario_year$p_value,
    interaction_F = interaction_row$F_value,
    interaction_p_value = interaction_row$p_value,
    AIC = stats::AIC(scenario_model),
    residual_sigma = summary(scenario_model)$sigma,
    BP_p_value = extract_metric(scenario_diag, "bp_p"),
    max_cooks_distance = extract_metric(scenario_diag, "max_cooks_distance"),
    stringsAsFactors = FALSE
  )
}

swapped_data <- swap_morris_years(data)
sensitivity <- rbind(
  sensitivity_row("published_lattice", data),
  sensitivity_row("morris_years_swapped", swapped_data)
)

swapped_year_site <- year_contrasts_by_site(
  candidate_models(swapped_data)$year_by_site
)
year_site_sensitivity <- rbind(
  transform(selected_year_site, scenario = "published_lattice"),
  transform(swapped_year_site, scenario = "morris_years_swapped")
)
year_site_sensitivity <- year_site_sensitivity[, c(
  "scenario", setdiff(names(year_site_sensitivity), "scenario")
)]

selection_sensitivity <- rbind(
  transform(model_comparison(candidate_models(data)), scenario = "published_lattice"),
  transform(model_comparison(candidate_models(swapped_data)), scenario = "morris_years_swapped")
)
selection_sensitivity <- selection_sensitivity[, c("scenario", setdiff(
  names(selection_sensitivity), "scenario"
))]

residual_data <- data.frame(
  observation = seq_len(nrow(data)),
  year = as.character(data$year),
  variety = as.character(data$variety),
  site = as.character(data$site),
  yield = data$yield,
  fitted = stats::fitted(model),
  residual = stats::residuals(model),
  standardized_residual = stats::rstandard(model),
  studentized_residual = stats::rstudent(model),
  stringsAsFactors = FALSE
)

write_csv(anova_table, "anova.csv")
write_csv(selected_anova, "selected_model_anova.csv")
write_csv(marginal_means, "marginal_means.csv")
write_csv(pairwise, "pairwise_holm.csv")
write_csv(selected_marginal_means, "selected_model_marginal_means.csv")
write_csv(selected_pairwise, "selected_model_pairwise_holm.csv")
write_csv(year_site_sensitivity, "year_by_site_contrasts.csv")
write_csv(pairwise[, c(
  "factor", "contrast", "estimate", "conf_low", "conf_high",
  "standardized_effect", "standardized_conf_low", "standardized_conf_high"
)], "effect_sizes.csv")
write_csv(comparison, "model_comparison.csv")
write_csv(diagnostics, "diagnostics.csv")
write_csv(selected_diagnostics, "selected_model_diagnostics.csv")
write_csv(influence, "influence.csv")
write_csv(coefficients, "coefficients_ols_hc3.csv")
write_csv(year_robustness, "year_effect_robustness.csv")
write_csv(cell_predictions, "cell_predictions.csv")
write_csv(sensitivity, "sensitivity_morris_swap.csv")
write_csv(selection_sensitivity, "model_selection_sensitivity.csv")
write_csv(residual_data, "residuals.csv")
writeLines(capture.output(utils::sessionInfo()), file.path(results_dir, "session_info.txt"))

# Diagnostic panels.
grDevices::png(file.path(figures_dir, "model_diagnostics.png"),
               width = 1800, height = 1800, res = 180)
old_par <- graphics::par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1.5),
                         oma = c(0, 0, 2, 0))
graphics::plot(stats::fitted(model), stats::residuals(model),
               pch = 19, col = grDevices::adjustcolor("#1b6ca8", 0.7),
               xlab = "Fitted yield", ylab = "Residual", main = "Residuals vs fitted")
graphics::abline(h = 0, lty = 2, col = "#555555")
stats::qqnorm(stats::rstudent(model), pch = 19,
              col = grDevices::adjustcolor("#1b6ca8", 0.7),
              main = "Studentized-residual Q-Q")
stats::qqline(stats::rstudent(model), col = "#c44e52", lwd = 2)
graphics::plot(stats::fitted(model), sqrt(abs(stats::rstandard(model))),
               pch = 19, col = grDevices::adjustcolor("#1b6ca8", 0.7),
               xlab = "Fitted yield", ylab = expression(sqrt("|standardized residual|")),
               main = "Scale-location")
graphics::plot(seq_len(nrow(data)), stats::cooks.distance(model),
               pch = 19, col = grDevices::adjustcolor("#1b6ca8", 0.7),
               xlab = "Observation", ylab = "Cook's distance", main = "Influence")
graphics::abline(h = 4 / nrow(data), lty = 2, col = "#c44e52")
graphics::mtext("Primary model diagnostics", outer = TRUE, cex = 1.15, font = 2)
graphics::par(old_par)
grDevices::dev.off()

# Site-specific year contrasts from the AICc-leading year-by-site model.
grDevices::png(file.path(figures_dir, "year_site_interaction.png"),
               width = 1800, height = 1200, res = 180)
site_year_means <- stats::aggregate(yield ~ site + year, data = data, FUN = mean)
site_year_matrix <- xtabs(yield ~ site + year, data = site_year_means)
old_par <- graphics::par(mar = c(8, 5, 3, 1))
graphics::matplot(
  seq_len(nrow(site_year_matrix)), site_year_matrix,
  type = "b", lty = 1, lwd = 2, pch = c(16, 17),
  col = c("#1b6ca8", "#c44e52"), xaxt = "n", xlab = "Site",
  ylab = "Mean yield (bushels/acre)", main = "Year-by-site pattern in published data"
)
graphics::axis(1, at = seq_len(nrow(site_year_matrix)),
               labels = rownames(site_year_matrix), las = 2, cex.axis = 0.8)
graphics::legend("topright", legend = colnames(site_year_matrix),
                 col = c("#1b6ca8", "#c44e52"), pch = c(16, 17),
                 lty = 1, lwd = 2, bty = "n")
graphics::par(old_par)
grDevices::dev.off()

# Model-derived marginal means and 95% confidence intervals.
grDevices::png(file.path(figures_dir, "marginal_means.png"),
               width = 1800, height = 1800, res = 180)
old_par <- graphics::par(mfrow = c(1, 3), mar = c(5, 8, 3, 1))
for (name in factor_names) {
  table <- emm_objects[[name]]$table
  y <- seq_len(nrow(table))
  xlim <- range(c(table$conf_low, table$conf_high))
  graphics::plot(table$estimate, y, xlim = xlim,
                 ylim = c(0.5, nrow(table) + 0.5), yaxt = "n",
                 pch = 19, col = "#1b6ca8", xlab = "Estimated marginal mean",
                 ylab = "", main = name)
  graphics::segments(table$conf_low, y, table$conf_high, y,
                     col = "#1b6ca8", lwd = 2)
  graphics::axis(2, at = y, labels = table$level, las = 2, cex.axis = 0.72)
  graphics::grid(nx = NA, ny = NULL, col = "#eeeeee")
  graphics::points(table$estimate, y, pch = 19, col = "#1b6ca8")
}
graphics::par(old_par)
grDevices::dev.off()

# Variety-by-site patterns averaged across the two years.
interaction_means <- stats::aggregate(yield ~ variety + site, data = data, FUN = mean)
interaction_matrix <- xtabs(yield ~ variety + site, data = interaction_means)
grDevices::png(file.path(figures_dir, "variety_site_interaction.png"),
               width = 1800, height = 1200, res = 180)
old_par <- graphics::par(mar = c(8, 5, 3, 9), xpd = NA)
palette <- c("#1b6ca8", "#c44e52", "#55a868", "#8172b2", "#ccb974", "#64b5cd")
graphics::matplot(
  seq_len(nrow(interaction_matrix)), interaction_matrix,
  type = "b", lty = 1, lwd = 1.8, pch = seq_len(ncol(interaction_matrix)) + 14,
  col = palette, xaxt = "n", xlab = "Variety", ylab = "Mean yield (bushels/acre)",
  main = "Variety-by-site yield patterns"
)
graphics::axis(1, at = seq_len(nrow(interaction_matrix)),
               labels = rownames(interaction_matrix), las = 2, cex.axis = 0.75)
graphics::legend("topright", inset = c(-0.32, 0), legend = colnames(interaction_matrix),
                 col = palette, lty = 1, pch = seq_len(ncol(interaction_matrix)) + 14,
                 bty = "n", cex = 0.8)
graphics::par(old_par)
grDevices::dev.off()

year_result <- year_ols[1L, ]
interaction_result <- anova_table[anova_table$term == "variety:site", ]
year_site_result <- selected_anova[selected_anova$term == "year:site", ]
top_model <- comparison[1L, ]
published <- sensitivity[sensitivity$scenario == "published_lattice", ]
swapped <- sensitivity[sensitivity$scenario == "morris_years_swapped", ]
significant_variety <- sum(pairwise$factor == "variety" & pairwise$significant_holm_0_05)
significant_site <- sum(pairwise$factor == "site" & pairwise$significant_holm_0_05)

summary_json <- c(
  "{",
  sprintf("  \"n\": %d,", nrow(data)),
  sprintf("  \"primary_model\": \"%s\",", "yield ~ year + variety + site + variety:site"),
  sprintf("  \"year_contrast\": \"%s\",", year_result$contrast),
  sprintf("  \"year_estimate\": %.8f,", year_result$estimate),
  sprintf("  \"year_conf_low\": %.8f,", year_result$conf_low),
  sprintf("  \"year_conf_high\": %.8f,", year_result$conf_high),
  sprintf("  \"year_p_value\": %.12g,", year_result$p_value),
  sprintf("  \"variety_site_interaction_p\": %.12g,", interaction_result$p_value),
  sprintf("  \"best_AICc_model\": \"%s\",", top_model$model),
  sprintf("  \"best_AICc\": %.8f,", top_model$AICc),
  sprintf("  \"morris_swap_year_estimate\": %.8f", swapped$year_estimate),
  "}"
)
writeLines(summary_json, file.path(results_dir, "summary.json"))

anova_report <- anova_table[anova_table$term != "Residuals", c(
  "term", "df", "F_value", "p_value", "partial_eta_squared", "omega_squared"
)]
names(anova_report) <- c("Term", "df", "F", "p", "partial eta2", "omega2")
model_report <- head(comparison[, c(
  "model", "parameters", "AICc", "delta_AICc", "AICc_weight", "BIC", "LOOCV_RMSE"
)], 6L)
names(model_report) <- c("Model", "k", "AICc", "delta AICc", "weight", "BIC", "LOOCV RMSE")
robust_report <- year_robustness[, c(
  "method", "estimate", "standard_error", "conf_low", "conf_high", "p_value"
)]
names(robust_report) <- c("Covariance", "Estimate", "SE", "95% low", "95% high", "p")
year_site_report <- selected_year_site[, c(
  "site", "estimate", "conf_low", "conf_high", "p_holm", "significant_holm_0_05"
)]
names(year_site_report) <- c("Site", "1932 - 1931", "95% low", "95% high", "Holm p", "Holm significant")

report <- c(
  "# 大麦収量の要因効果分析 / Factorial Analysis of Barley Yield",
  "",
  "## Executive summary",
  "",
  sprintf(
    "The analysis uses all %d published observations and a pre-specified fixed-effects model, `yield ~ year + variety + site + variety:site`. The model-estimated %s contrast was %.2f bushels/acre (95%% CI %.2f to %.2f; p %s).",
    nrow(data), year_result$contrast, year_result$estimate,
    year_result$conf_low, year_result$conf_high, format_p(year_result$p_value)
  ),
  "",
  sprintf(
    "The variety-by-site interaction was %s (F(%d, %d) = %.2f, p %s). Holm-adjusted model contrasts identified %d of 45 variety pairs and %d of 15 site pairs at alpha = 0.05. These comparisons use the fitted factorial model; they are not pooled raw pairwise t tests.",
    if (interaction_result$p_value < 0.05) "detectable" else "not detectable",
    interaction_result$df, stats::df.residual(model), interaction_result$F_value,
    format_p(interaction_result$p_value), significant_variety, significant_site
  ),
  "",
  sprintf(
    "AICc ranked `%s` first (AICc %.2f), so the interaction model is retained as the pre-specified inferential model rather than claimed as the uniquely best predictive model. Full rankings and leave-one-out errors are machine-readable.",
    top_model$model, top_model$AICc
  ),
  "",
  "## Data and design",
  "",
  "The public `lattice::barley` data contain mean yield in bushels per acre for 10 varieties at six Minnesota sites in 1931 and 1932. Every year-by-variety-by-site combination appears once. Each published value is already averaged across three randomized blocks. The original plot-level replicates are unavailable, so block variance and the original randomized-block model cannot be recovered.",
  "",
  "The official R documentation is the source of record: <https://search.r-project.org/CRAN/refmans/lattice/html/barley.html>. It cites Immer, Hayes, and Powers (1934) and Wright (2013).",
  "",
  "## Methods",
  "",
  "The primary model includes fixed main effects for year, variety, and site plus variety-by-site interaction. Because the data are a complete balanced factorial, the reported sequential ANOVA terms are orthogonal for this specification. The three-way interaction is not estimable separately from error after block averaging and is not fitted.",
  "",
  "Estimated marginal means are linear combinations of the model coefficients, averaged equally over the other design factors. Pairwise contrasts use the model residual variance and degrees of freedom. P-values are Holm-adjusted separately within the year, variety, and site families. Reported 95% contrast intervals are unadjusted model intervals; the standardized intervals divide those limits by the residual standard deviation and are approximate because uncertainty in that standard deviation is not propagated.",
  "",
  "Model-selection sensitivity compares nested scientific candidates with AIC, small-sample AICc, BIC, adjusted R-squared, and exact linear-model leave-one-out RMSE. Diagnostics include residual plots, Shapiro-Wilk, a fitted-value Koenker-Breusch-Pagan auxiliary regression, leverage, externally studentized residuals, Cook's distance, and HC3 covariance sensitivity.",
  "",
  "## Factorial ANOVA and effect magnitudes",
  "",
  markdown_table(anova_report),
  "",
  "Partial eta-squared describes each term relative to that term plus residual sum of squares. Omega-squared is a less biased variance-share estimate; both should be interpreted descriptively because sites and varieties are fixed and the data are aggregated.",
  "",
  "## Year effect and heteroskedasticity sensitivity",
  "",
  markdown_table(robust_report),
  "",
  sprintf(
    "The fitted-value heteroskedasticity check gave p %s. The Shapiro-Wilk residual check gave p %s. The largest Cook's distance was %.3f; %d observations exceeded the screening threshold 4/n = %.3f. These thresholds are diagnostics, not automatic deletion rules; no observations were removed.",
    format_p(extract_metric(diagnostics, "bp_p")),
    format_p(extract_metric(diagnostics, "shapiro_p")),
    extract_metric(diagnostics, "max_cooks_distance"),
    as.integer(extract_metric(diagnostics, "count_high_cooks")),
    extract_metric(diagnostics, "cooks_flag_threshold")
  ),
  "",
  "## Model-selection sensitivity",
  "",
  markdown_table(model_report),
  "",
  sprintf(
    "The leading model contains a year-by-site interaction. Its joint interaction test was F(%d, %d) = %.2f (p %s). Site-specific model contrasts are below; Holm adjustment covers the six site contrasts.",
    year_site_result$df, stats::df.residual(selected_model), year_site_result$F_value,
    format_p(year_site_result$p_value)
  ),
  "",
  markdown_table(year_site_report),
  "",
  "AIC-based selection addresses relative expected predictive fit, whereas ANOVA tests address a pre-specified effect. They need not return the same wording of a conclusion. Accordingly, the report keeps the variety-by-site model for the planned interaction test and treats the year-by-site model as strong sensitivity evidence. Results from both models are exported.",
  "",
  "## Morris year-label sensitivity",
  "",
  sprintf(
    "The package documentation notes Cleveland's suggestion that Morris may have its two years switched. The published data give a %s estimate of %.2f (95%% CI %.2f to %.2f). Swapping only the Morris year labels gives %.2f (95%% CI %.2f to %.2f). The interaction p-value changes from %s to %s. This scenario is a sensitivity analysis, not a correction to the canonical data.",
    published$year_contrast, published$year_estimate, published$year_conf_low,
    published$year_conf_high, swapped$year_estimate, swapped$year_conf_low,
    swapped$year_conf_high, format_p(published$interaction_p_value),
    format_p(swapped$interaction_p_value)
  ),
  "",
  "The year-by-site model has exactly the same likelihood after a within-Morris year-label swap: it can absorb that relabeling by reversing the Morris simple effect. By contrast, the additive model improves sharply after the swap and comes within about 2 AICc units of the year-by-site model. This shows that the published Morris pattern drives much, but not necessarily all, of the apparent year-by-site heterogeneity.",
  "",
  "## Conclusions",
  "",
  "Year and site differences are large in this aggregated historical experiment, while variety differences are smaller. The model explicitly asks whether varieties behave differently by site; the joint interaction result and all model-based contrasts are reported without replacing the design with independent raw t tests. Conclusions are limited to these observed years, sites, and varieties and should not be presented as a modern agronomic forecast.",
  "",
  "## Reproducibility record",
  "",
  "Run `make all` from this directory. It rebuilds the data snapshot from `lattice`, runs the analysis, regenerates this report and every table/figure, and executes the tests. The verified runtime is recorded in `results/session_info.txt`; exact R package pins are in `renv.lock`.",
  ""
)
writeLines(report, "REPORT.md")

cat("Analysis complete.\n")
cat("Primary year contrast:", year_result$contrast, sprintf("%.3f", year_result$estimate), "\n")
cat("Top AICc model:", top_model$model, "\n")
cat("Report: REPORT.md\n")
