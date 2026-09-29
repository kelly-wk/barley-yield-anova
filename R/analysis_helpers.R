abort <- function(message) {
  stop(message, call. = FALSE)
}

assert_that <- function(condition, message) {
  if (!isTRUE(condition)) abort(message)
}

barley_levels <- function() {
  list(
    year = c("1931", "1932"),
    variety = c(
      "Svansota", "No. 462", "Manchuria", "No. 475", "Velvet",
      "Peatland", "Glabron", "No. 457", "Wisconsin No. 38", "Trebi"
    ),
    site = c(
      "Grand Rapids", "Duluth", "University Farm", "Morris",
      "Crookston", "Waseca"
    )
  )
}

coerce_barley <- function(data) {
  required <- c("yield", "variety", "year", "site")
  assert_that(setequal(names(data), required),
              "Data must contain exactly yield, variety, year, and site.")
  data <- data[, required]
  lev <- barley_levels()
  data$yield <- as.numeric(data$yield)
  data$year <- factor(as.character(data$year), levels = lev$year)
  data$variety <- factor(as.character(data$variety), levels = lev$variety)
  data$site <- factor(as.character(data$site), levels = lev$site)
  data
}

validate_barley <- function(data) {
  assert_that(nrow(data) == 120L, "Expected exactly 120 observations.")
  assert_that(!anyNA(data), "The barley data must not contain missing values.")
  assert_that(all(is.finite(data$yield)), "Yield must be finite and numeric.")
  assert_that(all(data$yield > 0), "All yields must be positive.")
  cell_counts <- xtabs(~ year + variety + site, data = data)
  assert_that(all(cell_counts == 1L),
              "Expected one observation in every year-by-variety-by-site cell.")
  invisible(TRUE)
}

read_barley <- function(path) {
  assert_that(file.exists(path), paste("Data file not found:", path))
  data <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  data <- coerce_barley(data)
  validate_barley(data)
  data
}

candidate_models <- function(data) {
  list(
    intercept = stats::lm(yield ~ 1, data = data),
    year = stats::lm(yield ~ year, data = data),
    variety = stats::lm(yield ~ variety, data = data),
    site = stats::lm(yield ~ site, data = data),
    year_variety = stats::lm(yield ~ year + variety, data = data),
    year_site = stats::lm(yield ~ year + site, data = data),
    variety_site = stats::lm(yield ~ variety + site, data = data),
    additive = stats::lm(yield ~ year + variety + site, data = data),
    site_by_variety = stats::lm(
      yield ~ year + variety + site + variety:site, data = data
    ),
    year_by_site = stats::lm(
      yield ~ year + variety + site + year:site, data = data
    ),
    year_by_variety = stats::lm(
      yield ~ year + variety + site + year:variety, data = data
    ),
    all_two_way = stats::lm(
      yield ~ (year + variety + site)^2, data = data
    )
  )
}

primary_model <- function(data) {
  stats::lm(yield ~ year + variety + site + variety:site, data = data)
}

model_comparison <- function(models) {
  rows <- lapply(names(models), function(name) {
    model <- models[[name]]
    n <- stats::nobs(model)
    k <- attr(stats::logLik(model), "df")
    aic <- stats::AIC(model)
    aicc <- if (n > k + 1) aic + (2 * k * (k + 1)) / (n - k - 1) else Inf
    leverage <- stats::hatvalues(model)
    press_residual <- stats::residuals(model) / pmax(1 - leverage, .Machine$double.eps)
    data.frame(
      model = name,
      formula = paste(deparse(stats::formula(model)), collapse = " "),
      n = n,
      parameters = k,
      residual_df = stats::df.residual(model),
      log_likelihood = as.numeric(stats::logLik(model)),
      AIC = aic,
      AICc = aicc,
      BIC = stats::BIC(model),
      adjusted_r_squared = summary(model)$adj.r.squared,
      LOOCV_RMSE = sqrt(mean(press_residual^2)),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$delta_AICc <- result$AICc - min(result$AICc)
  relative <- exp(-0.5 * result$delta_AICc)
  result$AICc_weight <- relative / sum(relative)
  result <- result[order(result$AICc, result$LOOCV_RMSE), ]
  result$AICc_rank <- seq_len(nrow(result))
  rownames(result) <- NULL
  result
}

anova_effects <- function(model) {
  table <- stats::anova(model)
  result <- data.frame(
    term = rownames(table),
    df = table[["Df"]],
    sum_sq = table[["Sum Sq"]],
    mean_sq = table[["Mean Sq"]],
    F_value = table[["F value"]],
    p_value = table[["Pr(>F)"]],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  residual_ss <- tail(result$sum_sq, 1L)
  residual_ms <- tail(result$mean_sq, 1L)
  total_ss <- sum(result$sum_sq)
  result$partial_eta_squared <- ifelse(
    result$term == "Residuals", NA_real_,
    result$sum_sq / (result$sum_sq + residual_ss)
  )
  result$omega_squared <- ifelse(
    result$term == "Residuals", NA_real_,
    pmax(0, (result$sum_sq - result$df * residual_ms) / (total_ss + residual_ms))
  )
  result
}

hc3_vcov <- function(model) {
  x <- stats::model.matrix(model)
  e <- stats::residuals(model)
  h <- stats::hatvalues(model)
  bread <- chol2inv(qr.R(qr(x)))
  scaled <- e^2 / pmax((1 - h)^2, .Machine$double.eps)
  meat <- crossprod(x, x * scaled)
  bread %*% meat %*% bread
}

breusch_pagan_fitted <- function(model) {
  e2 <- stats::residuals(model)^2
  fitted <- stats::fitted(model)
  auxiliary <- stats::lm(e2 ~ fitted)
  statistic <- length(e2) * summary(auxiliary)$r.squared
  data.frame(
    test = "Koenker-Breusch-Pagan (fitted-value auxiliary regression)",
    statistic = unname(statistic),
    df = 1,
    p_value = stats::pchisq(statistic, df = 1, lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}

reference_grid <- function(model, target) {
  factor_names <- c("year", "variety", "site")
  assert_that(target %in% factor_names, paste("Unknown target factor:", target))
  other <- setdiff(factor_names, target)
  base <- expand.grid(
    lapply(model$xlevels[other], identity),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  names(base) <- other
  target_levels <- model$xlevels[[target]]
  grids <- lapply(target_levels, function(level) {
    grid <- base
    grid[[target]] <- level
    for (name in factor_names) {
      grid[[name]] <- factor(grid[[name]], levels = model$xlevels[[name]])
    }
    grid[, factor_names, drop = FALSE]
  })
  names(grids) <- target_levels
  grids
}

estimated_marginal_means <- function(model, target, covariance = stats::vcov(model),
                                     conf_level = 0.95) {
  grids <- reference_grid(model, target)
  beta <- stats::coef(model)
  assert_that(!anyNA(beta), "Model has aliased coefficients; EMMs are undefined.")
  term_object <- stats::delete.response(stats::terms(model))
  l_matrix <- do.call(rbind, lapply(grids, function(grid) {
    mm <- stats::model.matrix(
      term_object, data = grid, contrasts.arg = model$contrasts,
      xlev = model$xlevels
    )
    mm <- mm[, names(beta), drop = FALSE]
    colMeans(mm)
  }))
  estimates <- as.vector(l_matrix %*% beta)
  standard_errors <- sqrt(diag(l_matrix %*% covariance %*% t(l_matrix)))
  critical <- stats::qt(1 - (1 - conf_level) / 2, df = stats::df.residual(model))
  table <- data.frame(
    factor = target,
    level = rownames(l_matrix),
    estimate = estimates,
    standard_error = standard_errors,
    df = stats::df.residual(model),
    conf_low = estimates - critical * standard_errors,
    conf_high = estimates + critical * standard_errors,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  list(table = table, L = l_matrix, covariance = covariance)
}

pairwise_from_emm <- function(model, emm, target, conf_level = 0.95) {
  levels <- rownames(emm$L)
  pairs <- utils::combn(seq_along(levels), 2L)
  sigma <- summary(model)$sigma
  critical <- stats::qt(1 - (1 - conf_level) / 2, df = stats::df.residual(model))
  rows <- lapply(seq_len(ncol(pairs)), function(index) {
    first <- pairs[1L, index]
    second <- pairs[2L, index]
    contrast_vector <- emm$L[second, ] - emm$L[first, ]
    estimate <- as.numeric(contrast_vector %*% stats::coef(model))
    standard_error <- sqrt(as.numeric(
      contrast_vector %*% emm$covariance %*% contrast_vector
    ))
    statistic <- estimate / standard_error
    conf_low <- estimate - critical * standard_error
    conf_high <- estimate + critical * standard_error
    data.frame(
      factor = target,
      contrast = paste(levels[second], "-", levels[first]),
      level_1 = levels[first],
      level_2 = levels[second],
      estimate = estimate,
      standard_error = standard_error,
      df = stats::df.residual(model),
      t_value = statistic,
      p_value = 2 * stats::pt(abs(statistic), df = stats::df.residual(model),
                              lower.tail = FALSE),
      conf_low = conf_low,
      conf_high = conf_high,
      standardized_effect = estimate / sigma,
      standardized_conf_low = conf_low / sigma,
      standardized_conf_high = conf_high / sigma,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$p_holm <- stats::p.adjust(result$p_value, method = "holm")
  result$significant_holm_0_05 <- result$p_holm < 0.05
  rownames(result) <- NULL
  result
}

year_contrasts_by_site <- function(model, covariance = stats::vcov(model),
                                   conf_level = 0.95) {
  beta <- stats::coef(model)
  assert_that(!anyNA(beta), "Model has aliased coefficients; contrasts are undefined.")
  assert_that(all(c("year", "variety", "site") %in% names(model$xlevels)),
              "The model must include year, variety, and site.")
  assert_that(setequal(model$xlevels$year, c("1931", "1932")),
              "Expected year levels 1931 and 1932.")
  term_object <- stats::delete.response(stats::terms(model))
  df <- stats::df.residual(model)
  critical <- stats::qt(1 - (1 - conf_level) / 2, df = df)
  sigma <- summary(model)$sigma

  rows <- lapply(model$xlevels$site, function(site_level) {
    grids <- lapply(model$xlevels$year, function(year_level) {
      grid <- data.frame(
        year = factor(year_level, levels = model$xlevels$year),
        variety = factor(model$xlevels$variety, levels = model$xlevels$variety),
        site = factor(site_level, levels = model$xlevels$site)
      )
      mm <- stats::model.matrix(
        term_object, data = grid, contrasts.arg = model$contrasts,
        xlev = model$xlevels
      )
      colMeans(mm[, names(beta), drop = FALSE])
    })
    contrast_vector <- grids[[2L]] - grids[[1L]]
    estimate <- as.numeric(contrast_vector %*% beta)
    standard_error <- sqrt(as.numeric(
      contrast_vector %*% covariance %*% contrast_vector
    ))
    statistic <- estimate / standard_error
    conf_low <- estimate - critical * standard_error
    conf_high <- estimate + critical * standard_error
    data.frame(
      site = site_level,
      contrast = "1932 - 1931",
      estimate = estimate,
      standard_error = standard_error,
      df = df,
      t_value = statistic,
      p_value = 2 * stats::pt(abs(statistic), df = df, lower.tail = FALSE),
      conf_low = conf_low,
      conf_high = conf_high,
      standardized_effect = estimate / sigma,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$p_holm <- stats::p.adjust(result$p_value, method = "holm")
  result$significant_holm_0_05 <- result$p_holm < 0.05
  rownames(result) <- NULL
  result
}

diagnostic_summary <- function(model) {
  studentized <- stats::rstudent(model)
  leverage <- stats::hatvalues(model)
  cooks <- stats::cooks.distance(model)
  n <- stats::nobs(model)
  p <- length(stats::coef(model))
  shapiro <- stats::shapiro.test(stats::residuals(model))
  bp <- breusch_pagan_fitted(model)
  data.frame(
    metric = c(
      "n", "parameters", "residual_df", "residual_sigma",
      "shapiro_W", "shapiro_p", "bp_statistic", "bp_df", "bp_p",
      "max_abs_studentized_residual", "count_abs_studentized_gt_3",
      "max_leverage", "leverage_flag_threshold", "count_high_leverage",
      "max_cooks_distance", "cooks_flag_threshold", "count_high_cooks"
    ),
    value = c(
      n, p, stats::df.residual(model), summary(model)$sigma,
      unname(shapiro$statistic), shapiro$p.value,
      bp$statistic, bp$df, bp$p_value,
      max(abs(studentized)), sum(abs(studentized) > 3),
      max(leverage), 2 * p / n, sum(leverage > 2 * p / n),
      max(cooks), 4 / n, sum(cooks > 4 / n)
    ),
    stringsAsFactors = FALSE
  )
}

influence_table <- function(model, data) {
  result <- data.frame(
    observation = seq_len(nrow(data)),
    year = as.character(data$year),
    variety = as.character(data$variety),
    site = as.character(data$site),
    yield = data$yield,
    fitted = stats::fitted(model),
    residual = stats::residuals(model),
    studentized_residual = stats::rstudent(model),
    leverage = stats::hatvalues(model),
    cooks_distance = stats::cooks.distance(model),
    stringsAsFactors = FALSE
  )
  result$flag_abs_studentized_gt_3 <- abs(result$studentized_residual) > 3
  result$flag_cooks_gt_4_over_n <- result$cooks_distance > 4 / nrow(data)
  result[order(result$cooks_distance, decreasing = TRUE), ]
}

swap_morris_years <- function(data) {
  result <- data
  index <- result$site == "Morris"
  swapped <- ifelse(as.character(result$year[index]) == "1931", "1932", "1931")
  result$year[index] <- factor(swapped, levels = levels(result$year))
  validate_barley(result)
  result
}

extract_metric <- function(diagnostics, name) {
  diagnostics$value[match(name, diagnostics$metric)]
}

format_p <- function(value) {
  if (is.na(value)) return("")
  if (value < 0.001) return("<0.001")
  sprintf("%.3f", value)
}
