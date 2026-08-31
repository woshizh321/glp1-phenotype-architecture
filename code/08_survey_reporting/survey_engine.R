# Design-based survey functions used by the publication workflow.
# It does not read NHANES data or execute weighted estimates.

phase1_load_config <- function(path) {
  yaml::read_yaml(path)
}

phase1_build_design <- function(data, weight_variable) {
  stopifnot(weight_variable %in% c("WTMECPRP", "WTINTPRP"))
  required <- c("RIDAGEYR", "SDMVPSU", "SDMVSTRA", weight_variable)
  stopifnot(all(required %in% names(data)))
  options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
  eligible <- is.finite(data[[weight_variable]]) & data[[weight_variable]] > 0 &
    !is.na(data$SDMVPSU) & !is.na(data$SDMVSTRA)
  full_design <- survey::svydesign(
    ids = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = stats::as.formula(paste0("~", weight_variable)),
    nest = TRUE,
    data = data[eligible, , drop = FALSE]
  )
  subset(full_design, RIDAGEYR >= 20)
}

phase1_subset_target <- function(design, domain_id) {
  target_column <- paste0("domain_target_flag__", domain_id)
  stopifnot(target_column %in% names(design$variables))
  keep <- design$variables[[target_column]] == "TRUE"
  design[!is.na(keep) & keep, ]
}

phase1_status_proportion <- function(target_design, status_label, status_column) {
  stopifnot(status_column %in% names(target_design$variables))
  indicator <- target_design$variables[[status_column]] == status_label
  target_design$variables$.phase1_indicator <- indicator
  estimate <- survey::svyciprop(
    ~.phase1_indicator,
    target_design,
    method = "beta",
    level = 0.95,
    df = survey::degf(target_design)
  )
  mean_for_deff <- survey::svymean(~.phase1_indicator, target_design, deff = TRUE, na.rm = FALSE)
  list(estimate = estimate, mean_for_deff = mean_for_deff)
}

phase1_reliability <- function(p, lower, upper, denominator_n, event_n, design_effect, df) {
  effective_n <- if (is.finite(design_effect) && design_effect > 0) {
    min(denominator_n, denominator_n / design_effect)
  } else if (event_n %in% c(0, denominator_n)) {
    denominator_n
  } else {
    NA_real_
  }
  width <- upper - lower
  relative_width <- if (p > 0) width / p else Inf
  if (denominator_n < 30) return("SUPPRESS_DENOMINATOR_N_LT30")
  if (!is.finite(effective_n) || effective_n < 30) return("SUPPRESS_EFFECTIVE_N_LT30")
  if (event_n %in% c(0, denominator_n)) return("SUPPRESS_ZERO_OR_FULL_EVENT_AUTOMATED_REVIEW_RULE")
  if (!is.finite(width) || width >= 0.30) return("SUPPRESS_ABSOLUTE_CI_WIDTH_GE0_30")
  if (width > 0.05 && width < 0.30 && relative_width > 1.30) return("SUPPRESS_RELATIVE_CI_WIDTH_GT130_PERCENT")
  if (df < 8) return("SUPPRESS_DF_LT8_AUTOMATED_REVIEW_RULE")
  "PRESENT"
}
