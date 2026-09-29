#!/usr/bin/env Rscript

source(file.path("R", "analysis_helpers.R"), local = TRUE)

data <- read_barley(file.path("data", "barley.csv"))
stopifnot(nrow(data) == 120L)
stopifnot(nlevels(data$year) == 2L)
stopifnot(nlevels(data$variety) == 10L)
stopifnot(nlevels(data$site) == 6L)
stopifnot(all(xtabs(~ year + variety + site, data = data) == 1L))

model <- primary_model(data)
stopifnot(stats::df.residual(model) == 59L)
stopifnot(length(stats::coef(model)) == 61L)
stopifnot(!anyNA(stats::coef(model)))

anova_table <- anova_effects(model)
stopifnot(setequal(
  anova_table$term,
  c("year", "variety", "site", "variety:site", "Residuals")
))
stopifnot(anova_table$p_value[anova_table$term == "year"] < 0.001)
stopifnot(anova_table$p_value[anova_table$term == "site"] < 0.001)

for (name in c("year", "variety", "site")) {
  emm <- estimated_marginal_means(model, name)
  expected_levels <- nlevels(data[[name]])
  stopifnot(nrow(emm$table) == expected_levels)
  pairwise <- pairwise_from_emm(model, emm, name)
  stopifnot(nrow(pairwise) == choose(expected_levels, 2L))
  stopifnot(all(pairwise$p_holm >= pairwise$p_value))
  stopifnot(all(pairwise$p_holm >= 0 & pairwise$p_holm <= 1))
}

swapped <- swap_morris_years(data)
stopifnot(sum(as.character(swapped$year) != as.character(data$year)) == 20L)
stopifnot(all(xtabs(~ year + variety + site, data = swapped) == 1L))

hc3 <- hc3_vcov(model)
stopifnot(all(dim(hc3) == c(61L, 61L)))
stopifnot(max(abs(hc3 - t(hc3))) < 1e-8)
stopifnot(all(diag(hc3) > 0))

selected_model <- candidate_models(data)$year_by_site
site_year <- year_contrasts_by_site(selected_model)
stopifnot(nrow(site_year) == 6L)
stopifnot(all(site_year$contrast == "1932 - 1931"))
stopifnot(all(site_year$p_holm >= site_year$p_value))

cat("Helper tests passed.\n")
