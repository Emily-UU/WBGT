#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
} else {
  getwd()
}
eval(parse(file = file.path(script_dir, "00_config.R")))

analysis <- readRDS(survival_rds)

fit_linear_iqr <- function(df, outcome, exposure) {
  d <- make_survival_df(df, outcome) %>% filter(!is.na(.data[[exposure]]))
  iqr <- exposure_iqr(d[[exposure]])
  if (nrow(d) == 0 || sum(d$analysis_event) < 20 || is.na(iqr)) return(tibble())

  prep <- prep_covariates(d, clinical_covars)
  d <- prep$data
  d$exposure_iqr <- d[[exposure]] / iqr
  fit <- coxph(
    as.formula(paste0(
      "Surv(analysis_time_day, analysis_event) ~ exposure_iqr + ",
      paste(prep$terms, collapse = " + ")
    )),
    data = d,
    ties = "efron"
  )

  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "exposure_iqr") %>%
    transmute(
      outcome = outcome, exposure = exposure, analysis_type = "baseline_iqr_linear",
      n = nrow(d), events = sum(d$analysis_event), exposure_iqr = iqr,
      hazard_ratio = estimate, conf.low = conf.low, conf.high = conf.high,
      p.value = p.value, hr_ci = fmt_ci(hazard_ratio, conf.low, conf.high)
    )
}

fit_quartile <- function(df, outcome, exposure) {
  d <- make_survival_df(df, outcome) %>% filter(!is.na(.data[[exposure]]))
  if (nrow(d) == 0 || sum(d$analysis_event) < 20) return(tibble())

  qs <- unique(quantile(d[[exposure]], probs = seq(0, 1, 0.25), na.rm = TRUE, names = FALSE))
  if (length(qs) < 5) return(tibble())
  d$exposure_q <- cut(d[[exposure]], breaks = qs, include.lowest = TRUE, labels = paste0("Q", 1:4))
  d$exposure_q <- relevel(factor(d$exposure_q), ref = "Q1")
  d$exposure_q_score <- as.numeric(d$exposure_q)

  prep <- prep_covariates(d, clinical_covars)
  d <- prep$data
  fit <- coxph(
    as.formula(paste0(
      "Surv(analysis_time_day, analysis_event) ~ exposure_q + ",
      paste(prep$terms, collapse = " + ")
    )),
    data = d,
    ties = "efron"
  )
  trend_fit <- coxph(
    as.formula(paste0(
      "Surv(analysis_time_day, analysis_event) ~ exposure_q_score + ",
      paste(prep$terms, collapse = " + ")
    )),
    data = d,
    ties = "efron"
  )
  p_trend <- broom::tidy(trend_fit)$p.value[broom::tidy(trend_fit)$term == "exposure_q_score"][1]
  counts <- d %>% group_by(exposure_q) %>% summarise(category_n = n(), category_events = sum(analysis_event), .groups = "drop")

  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(startsWith(term, "exposure_q")) %>%
    mutate(quartile = str_remove(term, "exposure_q")) %>%
    left_join(counts, by = c("quartile" = "exposure_q")) %>%
    transmute(
      outcome, exposure, quartile, n = nrow(d), events = sum(d$analysis_event),
      category_n, category_events, hazard_ratio = estimate, conf.low, conf.high,
      p.value, p_trend, hr_ci = fmt_ci(hazard_ratio, conf.low, conf.high)
    )
}

fit_high_threshold <- function(df, outcome, exposure, probs = c(0.90, 0.95)) {
  d0 <- make_survival_df(df, outcome) %>% filter(!is.na(.data[[exposure]]))
  if (nrow(d0) == 0 || sum(d0$analysis_event) < 20) return(tibble())

  safe_bind(lapply(probs, function(prob) {
    cutoff <- as.numeric(quantile(d0[[exposure]], probs = prob, na.rm = TRUE, names = FALSE))
    d <- d0
    d$high_wbgt <- as.integer(d[[exposure]] >= cutoff)
    if (length(unique(d$high_wbgt)) < 2 || sum(d$analysis_event[d$high_wbgt == 1]) < 5) return(tibble())

    prep <- prep_covariates(d, clinical_covars)
    d <- prep$data
    fit <- coxph(
      as.formula(paste0(
        "Surv(analysis_time_day, analysis_event) ~ high_wbgt + ",
        paste(prep$terms, collapse = " + ")
      )),
      data = d,
      ties = "efron"
    )

    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "high_wbgt") %>%
      transmute(
        outcome, exposure, threshold = paste0("p", prob * 100),
        cutoff, n = nrow(d), events = sum(d$analysis_event),
        high_n = sum(d$high_wbgt == 1), high_events = sum(d$analysis_event[d$high_wbgt == 1]),
        hazard_ratio = estimate, conf.low, conf.high, p.value,
        hr_ci = fmt_ci(hazard_ratio, conf.low, conf.high)
      )
  }))
}

fit_spline_test <- function(df, outcome, exposure) {
  d <- make_survival_df(df, outcome) %>% filter(!is.na(.data[[exposure]]))
  iqr <- exposure_iqr(d[[exposure]])
  if (nrow(d) == 0 || sum(d$analysis_event) < 20 || is.na(iqr) || length(unique(d[[exposure]])) < 5) {
    return(tibble())
  }

  prep <- prep_covariates(d, clinical_covars)
  d <- prep$data
  d$exposure_linear <- d[[exposure]] / iqr
  covar_rhs <- paste(prep$terms, collapse = " + ")
  fit_linear <- coxph(as.formula(paste0("Surv(analysis_time_day, analysis_event) ~ exposure_linear + ", covar_rhs)), data = d, ties = "efron")
  fit_spline <- coxph(as.formula(paste0("Surv(analysis_time_day, analysis_event) ~ ns(", exposure, ", df = 3) + ", covar_rhs)), data = d, ties = "efron")
  lrt <- anova(fit_linear, fit_spline, test = "LRT")

  tibble(
    outcome, exposure, n = nrow(d), events = sum(d$analysis_event),
    linear_loglik = as.numeric(logLik(fit_linear)),
    spline_loglik = as.numeric(logLik(fit_spline)),
    p_nonlinearity = lrt$`Pr(>|Chi|)`[2]
  )
}

fit_spline_curve <- function(df, outcome, exposure) {
  d <- make_survival_df(df, outcome) %>% filter(!is.na(.data[[exposure]]))
  iqr <- exposure_iqr(d[[exposure]])
  if (nrow(d) == 0 || sum(d$analysis_event) < 20 || is.na(iqr) || length(unique(d[[exposure]])) < 5) {
    return(tibble())
  }

  prep <- prep_covariates(d, clinical_covars)
  d <- prep$data
  fit <- coxph(
    as.formula(paste0(
      "Surv(analysis_time_day, analysis_event) ~ ns(", exposure, ", df = 3) + ",
      paste(prep$terms, collapse = " + ")
    )),
    data = d,
    ties = "efron"
  )

  xgrid <- seq(
    quantile(d[[exposure]], 0.05, na.rm = TRUE, names = FALSE),
    quantile(d[[exposure]], 0.95, na.rm = TRUE, names = FALSE),
    length.out = 80
  )
  xref <- median(d[[exposure]], na.rm = TRUE)
  pred <- predict(fit, newdata = reference_newdata(d, prep$terms, exposure, xgrid), type = "lp", se.fit = TRUE)
  ref_lp <- as.numeric(predict(fit, newdata = reference_newdata(d, prep$terms, exposure, xref), type = "lp"))

  tibble(
    outcome, exposure, exposure_value = xgrid, reference_value = xref,
    n = nrow(d), events = sum(d$analysis_event),
    hazard_ratio = exp(pred$fit - ref_lp),
    conf.low = exp(pred$fit - ref_lp - 1.96 * pred$se.fit),
    conf.high = exp(pred$fit - ref_lp + 1.96 * pred$se.fit)
  )
}

fit_subgroup <- function(df, outcome, exposure, subgroup_var) {
  if (!subgroup_var %in% names(df)) return(tibble())
  levels <- sort(unique(na.omit(as.character(df[[subgroup_var]]))))
  if (length(levels) < 2) return(tibble())

  safe_bind(lapply(levels, function(level) {
    fit_linear_iqr(
      df %>% filter(as.character(.data[[subgroup_var]]) == level),
      outcome,
      exposure
    ) %>%
      mutate(subgroup = subgroup_var, level = level)
  }))
}

pairs <- tidyr::expand_grid(outcome = outcomes, exposure = baseline_exposures)

linear_results <- safe_bind(lapply(seq_len(nrow(pairs)), function(i) {
  fit_linear_iqr(analysis, pairs$outcome[i], pairs$exposure[i])
}))
quartile_results <- safe_bind(lapply(seq_len(nrow(pairs)), function(i) {
  fit_quartile(analysis, pairs$outcome[i], pairs$exposure[i])
}))
threshold_results <- safe_bind(lapply(seq_len(nrow(pairs)), function(i) {
  fit_high_threshold(analysis, pairs$outcome[i], pairs$exposure[i])
}))
spline_tests <- safe_bind(lapply(seq_len(nrow(pairs)), function(i) {
  fit_spline_test(analysis, pairs$outcome[i], pairs$exposure[i])
}))
spline_curves <- safe_bind(lapply(seq_len(nrow(pairs)), function(i) {
  fit_spline_curve(analysis, pairs$outcome[i], pairs$exposure[i])
}))
subgroup_results <- safe_bind(lapply(seq_len(nrow(pairs)), function(i) {
  safe_bind(lapply(subgroup_vars, function(sg) {
    fit_subgroup(analysis, pairs$outcome[i], pairs$exposure[i], sg)
  }))
}))

if (nrow(linear_results) > 0) linear_results <- linear_results %>% mutate(fdr_p = p.adjust(p.value, method = "BH"))
if (nrow(threshold_results) > 0) threshold_results <- threshold_results %>% mutate(fdr_p = p.adjust(p.value, method = "BH"))

write.csv(linear_results, file.path(results_dir, "04_linear_iqr_results.csv"), row.names = FALSE)
write.csv(quartile_results, file.path(results_dir, "04_quartile_results.csv"), row.names = FALSE)
write.csv(threshold_results, file.path(results_dir, "04_high_threshold_results.csv"), row.names = FALSE)
write.csv(spline_tests, file.path(results_dir, "04_spline_tests.csv"), row.names = FALSE)
write.csv(spline_curves, file.path(results_dir, "04_spline_curves.csv"), row.names = FALSE)
write.csv(subgroup_results, file.path(results_dir, "04_subgroup_results.csv"), row.names = FALSE)
message("Wrote focused survival analysis results.")
