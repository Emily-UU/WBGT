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

fit_baseline_cox <- function(df, outcome, exposure, covars, model_name) {
  d <- make_survival_df(df, outcome) %>%
    filter(!is.na(.data[[exposure]]))

  if (nrow(d) == 0 || sum(d$analysis_event) < 10) {
    return(tibble(
      outcome = outcome, exposure = exposure, model = model_name,
      n = nrow(d), events = sum(d$analysis_event), term = "per 1 deg C higher WBGT",
      hazard_ratio = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
      p.value = NA_real_, note = "Model skipped: too few records or events."
    ))
  }

  prep <- prep_covariates(d, covars)
  d <- prep$data
  xvar <- paste0(exposure, "_per_1c")
  d[[xvar]] <- d[[exposure]]
  fit <- coxph(
    as.formula(paste0(
      "Surv(analysis_time_day, analysis_event) ~ ",
      paste(c(xvar, prep$terms), collapse = " + ")
    )),
    data = d,
    ties = "efron"
  )

  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == xvar) %>%
    transmute(
      outcome = outcome,
      exposure = exposure,
      model = model_name,
      n = nrow(d),
      events = sum(d$analysis_event),
      term = "per 1 deg C higher WBGT",
      hazard_ratio = estimate,
      conf.low = conf.low,
      conf.high = conf.high,
      p.value = p.value,
      hr_ci = fmt_ci(hazard_ratio, conf.low, conf.high)
    )
}

results <- safe_bind(lapply(outcomes, function(outcome) {
  safe_bind(lapply(baseline_exposures, function(exposure) {
    safe_bind(list(
      fit_baseline_cox(analysis, outcome, exposure, base_covars, "Model 1"),
      fit_baseline_cox(analysis, outcome, exposure, clinical_covars, "Model 2")
    ))
  }))
}))

results <- results %>%
  group_by(outcome, model) %>%
  mutate(fdr_p = p.adjust(p.value, method = "BH")) %>%
  ungroup()

write.csv(results, file.path(results_dir, "02_baseline_cox_results.csv"), row.names = FALSE)
message("Wrote baseline Cox model results.")
