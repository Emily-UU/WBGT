#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE))
} else {
  getwd()
}
eval(parse(file = file.path(script_dir, "00_config.R")))

suppressPackageStartupMessages({
  library(sandwich)
  library(lmtest)
})

monthly <- read.csv(monthly_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  mutate(
    participant_id = as.character(participant_id),
    measurement_date = as.Date(measurement_date),
    sex = factor(sex),
    cohort_group = factor(cohort_group),
    gap_days_scaled = gap_days / 30.4375,
    egfr_decline_fraction = if_else(
      !is.na(prior_egfr) & prior_egfr > 0,
      (prior_egfr - egfr) / prior_egfr,
      NA_real_
    ),
    creatinine_ratio = creatinine / prior_creatinine,
    creatinine_increase_1p5x = as.integer(
      !is.na(prior_creatinine) & prior_creatinine > 0 & creatinine_ratio >= 1.5
    ),
    log_crp = log1p(crp)
  ) %>%
  filter(!is.na(monthly_wbgt), !is.na(gap_days), gap_days > 0)

impute_with_indicator <- function(dat, vars) {
  out <- dat
  for (v in vars) {
    missing_var <- paste0(v, "_missing")
    imputed_var <- paste0(v, "_imp")
    out[[missing_var]] <- as.integer(is.na(out[[v]]))
    med <- median(out[[v]], na.rm = TRUE)
    if (is.na(med) || is.infinite(med)) med <- 0
    out[[imputed_var]] <- out[[v]]
    out[[imputed_var]][is.na(out[[imputed_var]])] <- med
  }
  out
}

monthly <- impute_with_indicator(monthly, c("bmi", "hba1c", "ldl", "log_crp"))

fit_cluster_robust_logit <- function(dat, outcome, model_name, threshold = NA_real_, prior_term = "prior_egfr") {
  dat <- dat %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(monthly_wbgt),
      !is.na(.data[[prior_term]]),
      !is.na(age_at_measurement),
      !is.na(sex),
      !is.na(gap_days),
      gap_days > 0
    )

  base <- tibble(
    outcome = outcome,
    decline_threshold_percent = threshold,
    model = model_name,
    n_rows = nrow(dat),
    n_ids = dplyr::n_distinct(dat$participant_id),
    events = sum(dat[[outcome]], na.rm = TRUE),
    event_percent = 100 * mean(dat[[outcome]], na.rm = TRUE),
    term = "monthly_wbgt",
    monthly_wbgt_iqr = IQR(dat$monthly_wbgt, na.rm = TRUE)
  )

  if (nrow(dat) < 100 || length(unique(dat[[outcome]])) < 2) {
    return(base %>%
      mutate(
        beta = NA_real_, robust_se = NA_real_,
        OR_per_1C = NA_real_, CI_low_per_1C = NA_real_, CI_high_per_1C = NA_real_,
        OR_per_IQR = NA_real_, CI_low_per_IQR = NA_real_, CI_high_per_IQR = NA_real_,
        p_value = NA_real_
      ))
  }

  covars <- c(
    "monthly_wbgt", prior_term, "age_at_measurement", "sex", "gap_days_scaled",
    "bmi_imp", "bmi_missing", "nsaid_use", "acei_arb_use", "glp1_use", "sglt2i_use",
    "hba1c_imp", "hba1c_missing", "ldl_imp", "ldl_missing",
    "log_crp_imp", "log_crp_missing"
  )
  if (length(unique(dat$cohort_group)) > 1) covars <- c(covars, "cohort_group")

  fit <- glm(as.formula(paste(outcome, "~", paste(covars, collapse = " + "))), family = binomial(), data = dat)
  robust_vcov <- sandwich::vcovCL(fit, cluster = dat$participant_id, type = "HC0")
  ct <- lmtest::coeftest(fit, vcov. = robust_vcov)
  beta <- unname(ct["monthly_wbgt", "Estimate"])
  se <- unname(ct["monthly_wbgt", "Std. Error"])
  p <- unname(ct["monthly_wbgt", "Pr(>|z|)"])
  iqr <- base$monthly_wbgt_iqr[1]

  base %>%
    mutate(
      beta = beta,
      robust_se = se,
      OR_per_1C = exp(beta),
      CI_low_per_1C = exp(beta - 1.96 * se),
      CI_high_per_1C = exp(beta + 1.96 * se),
      OR_per_IQR = exp(beta * iqr),
      CI_low_per_IQR = exp((beta - 1.96 * se) * iqr),
      CI_high_per_IQR = exp((beta + 1.96 * se) * iqr),
      p_value = p
    )
}

thresholds <- c(15, 20, 30, 40, 50)
for (threshold in thresholds) {
  monthly[[paste0("egfr_decline_", threshold)]] <-
    as.integer(!is.na(monthly$egfr_decline_fraction) & monthly$egfr_decline_fraction > threshold / 100)
}

egfr_results <- safe_bind(lapply(thresholds, function(threshold) {
  fit_cluster_robust_logit(
    monthly,
    paste0("egfr_decline_", threshold),
    "Repeated monthly kidney-function analysis",
    threshold,
    "prior_egfr"
  )
}))

creatinine_results <- fit_cluster_robust_logit(
  monthly,
  "creatinine_increase_1p5x",
  "Repeated monthly kidney-function analysis",
  NA_real_,
  "prior_creatinine"
)

qc <- tibble(
  metric = c(
    "analysis_rows",
    "analysis_unique_ids",
    "bmi_missing_percent",
    "hba1c_available_percent",
    "ldl_available_percent",
    "crp_available_percent"
  ),
  value = c(
    nrow(monthly),
    dplyr::n_distinct(monthly$participant_id),
    100 * mean(is.na(monthly$bmi)),
    100 * mean(!is.na(monthly$hba1c)),
    100 * mean(!is.na(monthly$ldl)),
    100 * mean(!is.na(monthly$crp))
  )
)

write.csv(egfr_results, file.path(results_dir, "05_monthly_egfr_decline_results.csv"), row.names = FALSE)
write.csv(creatinine_results, file.path(results_dir, "05_monthly_creatinine_increase_results.csv"), row.names = FALSE)
write.csv(qc, file.path(results_dir, "05_monthly_model_qc.csv"), row.names = FALSE)
message("Wrote monthly kidney-function model results.")
