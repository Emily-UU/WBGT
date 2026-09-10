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
annual <- read.csv(annual_exposure_file, stringsAsFactors = FALSE, check.names = FALSE)

for (v in intersect(c("index_year", "exposure_year", annual_exposures), names(annual))) {
  annual[[v]] <- to_num(annual[[v]])
}

annual <- annual %>%
  mutate(index_date = as.Date(index_date), exposure_year = as.integer(exposure_year)) %>%
  select(participant_id, index_date, index_year, exposure_year, all_of(intersect(annual_exposures, names(.)))) %>%
  distinct()

make_time_updated_data <- function(df, outcome) {
  d <- make_survival_df(df, outcome) %>%
    mutate(
      end_day = analysis_time_day,
      event_day = if_else(analysis_event == 1, analysis_time_day, NA_real_),
      analysis_event_person = analysis_event
    )

  d %>%
    select(
      participant_id, index_date, index_year, end_day, event_day,
      analysis_event_person, all_of(intersect(c(clinical_covars, subgroup_vars), names(d)))
    ) %>%
    left_join(annual, by = c("participant_id", "index_date", "index_year")) %>%
    filter(!is.na(exposure_year), exposure_year >= as.integer(format(index_date, "%Y"))) %>%
    mutate(
      year_start_date = as.Date(paste0(exposure_year, "-01-01")),
      year_end_date = as.Date(paste0(exposure_year + 1, "-01-01")),
      interval_start_date = pmax(index_date, year_start_date),
      interval_start = pmax(as.numeric(interval_start_date - index_date), 0),
      interval_stop = pmin(as.numeric(year_end_date - index_date), end_day),
      analysis_event = if_else(
        analysis_event_person == 1 & !is.na(event_day) &
          event_day > interval_start & event_day <= interval_stop,
        1L, 0L
      )
    ) %>%
    filter(!is.na(interval_start), !is.na(interval_stop), interval_stop > interval_start) %>%
    arrange(participant_id, interval_start)
}

fit_time_updated_cox <- function(tv, outcome, exposure, covars, model_name) {
  tv <- tv %>% filter(!is.na(.data[[exposure]]))

  if (nrow(tv) == 0 || sum(tv$analysis_event) < 20) {
    return(tibble(
      outcome = outcome, exposure = exposure, model = model_name,
      n_persons = length(unique(tv$participant_id)), n_intervals = nrow(tv),
      events = sum(tv$analysis_event), term = "per 1 deg C higher annual WBGT",
      hazard_ratio = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
      p.value = NA_real_, note = "Model skipped: too few intervals or events."
    ))
  }

  prep <- prep_covariates(tv, covars)
  tv <- prep$data
  xvar <- paste0(exposure, "_per_1c")
  tv[[xvar]] <- tv[[exposure]]
  fit <- coxph(
    as.formula(paste0(
      "Surv(interval_start, interval_stop, analysis_event) ~ ",
      paste(c(xvar, prep$terms, "cluster(participant_id)"), collapse = " + ")
    )),
    data = tv,
    ties = "efron"
  )

  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == xvar) %>%
    transmute(
      outcome = outcome,
      exposure = exposure,
      model = model_name,
      n_persons = length(unique(tv$participant_id)),
      n_intervals = nrow(tv),
      events = sum(tv$analysis_event),
      term = "per 1 deg C higher annual WBGT",
      hazard_ratio = estimate,
      conf.low = conf.low,
      conf.high = conf.high,
      p.value = p.value,
      hr_ci = fmt_ci(hazard_ratio, conf.low, conf.high)
    )
}

tv_by_outcome <- setNames(lapply(outcomes, function(x) make_time_updated_data(analysis, x)), outcomes)

results <- safe_bind(lapply(outcomes, function(outcome) {
  tv <- tv_by_outcome[[outcome]]
  safe_bind(lapply(annual_exposures, function(exposure) {
    safe_bind(list(
      fit_time_updated_cox(tv, outcome, exposure, base_covars, "Model 1"),
      fit_time_updated_cox(tv, outcome, exposure, clinical_covars, "Model 2")
    ))
  }))
}))

results <- results %>%
  group_by(outcome, model) %>%
  mutate(fdr_p = p.adjust(p.value, method = "BH")) %>%
  ungroup()

qc <- safe_bind(lapply(names(tv_by_outcome), function(outcome) {
  tv <- tv_by_outcome[[outcome]]
  tibble(
    outcome = outcome,
    n_persons = dplyr::n_distinct(tv$participant_id),
    n_intervals = nrow(tv),
    events = sum(tv$analysis_event),
    first_exposure_year = min(tv$exposure_year, na.rm = TRUE),
    last_exposure_year = max(tv$exposure_year, na.rm = TRUE)
  )
}))

write.csv(results, file.path(results_dir, "03_time_updated_cox_results.csv"), row.names = FALSE)
write.csv(qc, file.path(results_dir, "03_time_updated_qc.csv"), row.names = FALSE)
message("Wrote annual time-updated Cox model results.")
