#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survey)
  library(readr)
  library(dplyr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
input_csv <- args[[1]]
out_root <- args[[2]]
repo_root <- Sys.getenv("PUBLIC_REPO_ROOT")
stopifnot(nzchar(repo_root))
source(file.path(repo_root, "code/08_survey_reporting/survey_engine.R"))
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE, warn = 1)

d <- readr::read_csv(input_csv, show_col_types = FALSE, progress = FALSE)
des <- phase1_build_design(d, "WTMECPRP")
internal <- file.path(out_root, "intermediate")
dir.create(internal, recursive = TRUE, showWarnings = FALSE)

trials <- c("ELIXA", "LEADER", "SUSTAIN-6", "EXSCEL", "HARMONY Outcomes", "PIONEER 6", "REWIND", "AMPLITUDE-O", "SOUL", "SURPASS-CVOT")
states <- c("DIRECT", "POTENTIAL", "NOT_ALIGNED", "INDETERMINATE")
cv_states <- c("DIRECT", "CONDITIONAL", "NO_DIRECT", "UNCLASSIFIED")
k_states <- c("DIRECT_COMPATIBLE", "CONDITIONAL", "NO_DIRECT", "UNCLASSIFIED")

safe_num <- function(x) if (length(x) == 0 || !is.finite(x[[1]])) NA_real_ else as.numeric(x[[1]])

estimate_binary <- function(indicator, denominator, module, estimand, category = "", denominator_label = "", extra = list()) {
  indicator <- as.logical(indicator); denominator <- as.logical(denominator)
  keep <- !is.na(denominator) & denominator
  sd <- des[keep, ]
  z <- indicator[keep]
  stopifnot(length(z) == nrow(sd$variables), !anyNA(z))
  if (length(z) == 0) {
    base <- tibble(
      module = module, estimand = estimand, category = category, denominator = denominator_label,
      unweighted_denominator_n = 0L, internal_unweighted_numerator = 0L,
      weighted_proportion = NA_real_, ci95_lower = NA_real_, ci95_upper = NA_real_, design_se = NA_real_,
      design_df = NA_real_, design_effect = NA_real_, effective_n = NA_real_,
      reliability_status = "SUPPRESS_DENOMINATOR_N_LT30", publication_display_status = "SUPPRESSED"
    )
    if (length(extra)) for (nm in names(extra)) base[[nm]] <- extra[[nm]]
    return(base)
  }
  sd$variables$.pa05_indicator <- z
  n <- length(z); events <- sum(z); df <- survey::degf(sd)
  mo <- survey::svymean(~.pa05_indicator, sd, deff = TRUE, na.rm = FALSE)
  p <- safe_num(coef(mo)); se <- safe_num(SE(mo)); deff <- suppressWarnings(safe_num(survey::deff(mo)))
  ci <- c(NA_real_, NA_real_)
  cp <- tryCatch(survey::svyciprop(~.pa05_indicator, sd, method = "beta", level = 0.95, df = df), error = function(e) NULL)
  if (!is.null(cp)) {
    p <- safe_num(coef(cp)); se <- safe_num(SE(cp)); ci <- as.numeric(confint(cp)[1, ])
  }
  rel <- phase1_reliability(p, ci[1], ci[2], n, events, deff, df)
  eff <- if (is.finite(deff) && deff > 0) min(n, n / deff) else if (events %in% c(0, n)) n else NA_real_
  base <- tibble(
    module = module, estimand = estimand, category = category, denominator = denominator_label,
    unweighted_denominator_n = n, internal_unweighted_numerator = events,
    weighted_proportion = p, ci95_lower = ci[1], ci95_upper = ci[2], design_se = se,
    design_df = df, design_effect = deff, effective_n = eff,
    reliability_status = rel,
    publication_display_status = ifelse(rel == "PRESENT", "PRESENT", "SUPPRESSED")
  )
  if (length(extra)) for (nm in names(extra)) base[[nm]] <- extra[[nm]]
  base
}

t2d <- des$variables$pa05_t2d
bmi25 <- des$variables$pa05_bmi25

# Module A: trial, portfolio, and original-versus-augmented matrix.
a_trial <- bind_rows(lapply(trials, function(trial) bind_rows(lapply(states, function(st) {
  estimate_binary(des$variables[[paste0("trial__", trial)]] == st, t2d, "A", "TRIAL_LEVEL_ALIGNMENT", st, "DIAGNOSED_T2D", list(trial = trial))
})))) %>% select(trial, everything())

a_port <- bind_rows(c(lapply(states, function(st) {
  estimate_binary(des$variables$portfolio_status == st, t2d, "A", "PORTFOLIO_ALIGNMENT", st, "DIAGNOSED_T2D")
}), list(estimate_binary(des$variables$portfolio_envelope, t2d, "A", "PORTFOLIO_ALIGNMENT", "DIRECT_OR_POTENTIAL_ENVELOPE", "DIAGNOSED_T2D"))))

a_matrix <- bind_rows(lapply(cv_states, function(cv) bind_rows(lapply(states, function(st) {
  estimate_binary(des$variables$original_cv == cv & des$variables$portfolio_status == st, t2d, "A", "ORIGINAL_VS_AUGMENTED_JOINT", paste(cv, st, sep = "__"), "DIAGNOSED_T2D", list(original_cv = cv, augmented_status = st))
})))) %>% select(original_cv, augmented_status, everything())
for (cv in cv_states) {
  for (st in states) {
    a_matrix <- bind_rows(a_matrix, estimate_binary(des$variables$portfolio_status == st, t2d & des$variables$original_cv == cv, "A", "ORIGINAL_VS_AUGMENTED_ROW_CONDITIONAL", paste(cv, st, sep = "__"), paste0("ORIGINAL_CV_", cv), list(original_cv = cv, augmented_status = st)) %>% select(original_cv, augmented_status, everything()))
  }
}

# Module B: depth, exact depth, accrual, breadth/depth, LOTO.
depth_bin <- function(x) ifelse(x == 0, "0", ifelse(x == 1, "1", ifelse(x <= 4, "2-4", ">=5")))
b_depth <- bind_rows(lapply(c("DIRECT", "DIRECT_OR_POTENTIAL"), function(measure) {
  x <- if (measure == "DIRECT") des$variables$n_direct else des$variables$n_direct_or_potential
  bind_rows(lapply(c("0", "1", "2-4", ">=5"), function(bin) estimate_binary(depth_bin(x) == bin, t2d, "B", "DEPTH_GROUP", bin, "DIAGNOSED_T2D", list(depth_measure = measure))))
})) %>% select(depth_measure, everything())

b_exact <- bind_rows(lapply(0:10, function(k) estimate_binary(des$variables$n_direct == k, t2d, "B", "EXACT_DIRECT_DEPTH", as.character(k), "DIAGNOSED_T2D", list(exact_depth = k)))) %>% select(exact_depth, everything())

years <- c(2015, 2016, 2017, 2018, 2019, 2021, 2025)
b_accrual <- bind_rows(lapply(c("DIRECT", "DIRECT_OR_POTENTIAL"), function(measure) bind_rows(lapply(years, function(y) {
  col <- if (measure == "DIRECT") paste0("cum_direct_", y) else paste0("cum_dp_", y)
  cur <- des$variables[[col]]
  prev <- if (y == years[1]) rep(FALSE, length(cur)) else des$variables[[if (measure == "DIRECT") paste0("cum_direct_", years[match(y, years)-1]) else paste0("cum_dp_", years[match(y, years)-1])]]
  bind_rows(
    estimate_binary(cur, t2d, "B", "CUMULATIVE_ACCRUAL", as.character(y), "DIAGNOSED_T2D", list(accrual_measure = measure, milestone_year = y, result_type = "CUMULATIVE")),
    estimate_binary(cur & !prev, t2d, "B", "INCREMENTAL_BREADTH", as.character(y), "DIAGNOSED_T2D", list(accrual_measure = measure, milestone_year = y, result_type = "INCREMENT_PERCENTAGE_POINTS"))
  )
})))) %>% select(accrual_measure, milestone_year, result_type, everything())

b_bd <- bind_rows(lapply(years, function(y) {
  cur_trials <- strsplit(readr::read_csv(file.path(repo_root, "metadata/evidence_milestone_map.csv"), show_col_types = FALSE)$trials_entering_bundle[match(y, years)], ";", fixed = TRUE)[[1]]
  cur_cols <- paste0("trial__", cur_trials)
  cur_any <- des$variables[[paste0("cum_direct_", y)]]
  prev_any <- if (y == years[1]) rep(FALSE, length(cur_any)) else des$variables[[paste0("cum_direct_", years[match(y, years)-1])]]
  new_trial_direct <- apply(des$variables[, cur_cols, drop = FALSE] == "DIRECT", 1, any)
  bind_rows(
    estimate_binary(cur_any & !prev_any, t2d, "B", "BREADTH_DEPTH_DECOMPOSITION", "NEW_BREADTH", "DIAGNOSED_T2D", list(milestone_year = y)),
    estimate_binary(prev_any & new_trial_direct, t2d, "B", "BREADTH_DEPTH_DECOMPOSITION", "ADDED_DEPTH", "DIAGNOSED_T2D", list(milestone_year = y))
  )
})) %>% select(milestone_year, everything())

b_loto <- bind_rows(lapply(trials, function(omit) bind_rows(lapply(c("DIRECT", "DIRECT_OR_POTENTIAL"), function(measure) {
  cols <- paste0("trial__", setdiff(trials, omit))
  all_col <- if (measure == "DIRECT") des$variables$n_direct >= 1 else des$variables$n_direct_or_potential >= 1
  minus <- if (measure == "DIRECT") {
    apply(des$variables[, cols, drop = FALSE] == "DIRECT", 1, any)
  } else {
    apply(sapply(des$variables[, cols, drop = FALSE], function(x) x %in% c("DIRECT", "POTENTIAL")), 1, any)
  }
  estimate_binary(all_col & !minus, t2d, "B", "LEAVE_ONE_TRIAL_OUT_INFLUENCE", omit, "DIAGNOSED_T2D", list(trial_omitted = omit, breadth_measure = measure))
})))) %>% select(trial_omitted, breadth_measure, everything())

# Module C.
c_cv <- bind_rows(lapply(c("CV-C1", "CV-C2", "CV-C3", "CV-C4"), function(cl) bind_rows(
  estimate_binary(des$variables$cv_conditional_class == cl, t2d & des$variables$original_cv == "CONDITIONAL", "C", "CV_CONDITIONAL_DECOMPOSITION", cl, "ORIGINAL_CV_CONDITIONAL", list(denominator_scope = "CONDITIONAL")),
  estimate_binary(des$variables$cv_conditional_class == cl, t2d, "C", "CV_CONDITIONAL_DECOMPOSITION", cl, "DIAGNOSED_T2D", list(denominator_scope = "FULL_T2D"))
))) %>% select(denominator_scope, everything())

c_context <- bind_rows(
  estimate_binary(des$variables$ascvd_group == "POSITIVE", t2d, "C", "CV_PORTFOLIO_CONTEXT", "CV1_ASCVD_ANCHOR", "CV0_DIAGNOSED_T2D"),
  estimate_binary(des$variables$portfolio_status == "DIRECT", t2d, "C", "CV_PORTFOLIO_CONTEXT", "CV2_DIRECT_PORTFOLIO", "CV0_DIAGNOSED_T2D"),
  estimate_binary(des$variables$portfolio_envelope, t2d, "C", "CV_PORTFOLIO_CONTEXT", "CV3_DIRECT_OR_POTENTIAL", "CV0_DIAGNOSED_T2D")
)

c_select_branch <- bind_rows(lapply(paste0("SELECT-C", 1:4), function(br) {
  estimate_binary(des$variables$select_branch == br, des$variables$status__OBESITY_CARDIOVASCULAR_OUTCOME == "CONDITIONAL_DOMAIN_REPRESENTATION", "C", "SELECT_CONDITIONAL_BRANCH", br, "ORIGINAL_OBESITY_CV_CONDITIONAL")
}))

c_select_funnel <- bind_rows(
  estimate_binary(des$variables$select_D1, des$variables$select_D0, "C", "SELECT_FUNNEL", "D1/D0", "D0_BMI_GE25"),
  estimate_binary(des$variables$select_D2, des$variables$select_D0, "C", "SELECT_FUNNEL", "D2/D0", "D0_BMI_GE25"),
  estimate_binary(des$variables$select_D2, des$variables$select_D1, "C", "SELECT_FUNNEL", "D2/D1", "D1_AGE45_BMI27"),
  estimate_binary(des$variables$select_D3, des$variables$select_D0, "C", "SELECT_FUNNEL", "D3/D0", "D0_BMI_GE25"),
  estimate_binary(des$variables$select_D3, des$variables$select_D2, "C", "SELECT_FUNNEL", "D3/D2", "D2_NO_DIAGNOSED_DIABETES"),
  estimate_binary(des$variables$select_D3_PAD, des$variables$select_D0, "C", "SELECT_FUNNEL", "D3_PAD/D0", "D0_BMI_GE25"),
  estimate_binary(des$variables$select_D3_ENVELOPE, des$variables$select_D0, "C", "SELECT_FUNNEL", "D3_OR_D3_PAD/D0", "D0_BMI_GE25")
)

k_reasons <- c("DIRECT_COMPATIBLE", "NO_DIRECT", "CONDITIONAL_REASON_UACR_AT_OR_ABOVE_5000", "CONDITIONAL_REASON_EGFR_15_TO_LT25", "CONDITIONAL_REASON_EGFR_GT75_WITH_CKD_MARKER", "CONDITIONAL_REASON_ALBUMINURIA_BELOW_BRANCH_THRESHOLD")
c_kidney <- bind_rows(lapply(k_reasons, function(reason) bind_rows(
  estimate_binary(des$variables$kidney_reason_class %in% reason, t2d & des$variables$original_kidney == "CONDITIONAL", "C", "KIDNEY_CONDITIONAL_REASON", reason, "ORIGINAL_KIDNEY_CONDITIONAL", list(denominator_scope = "CONDITIONAL")),
  estimate_binary(des$variables$kidney_reason_class %in% reason, t2d, "C", "KIDNEY_CONDITIONAL_REASON", reason, "DIAGNOSED_T2D", list(denominator_scope = "FULL_T2D"))
))) %>% select(denominator_scope, everything())

c_flow <- bind_rows(
  estimate_binary(des$variables$flow_K1, des$variables$flow_K0, "C", "FLOW_FUNNEL", "K1/K0", "K0_DIAGNOSED_T2D"),
  estimate_binary(des$variables$flow_K2, des$variables$flow_K0, "C", "FLOW_FUNNEL", "K2/K0", "K0_DIAGNOSED_T2D"),
  estimate_binary(des$variables$flow_K2, des$variables$flow_K1, "C", "FLOW_FUNNEL", "K2/K1", "K1_OBSERVABLE_CKD"),
  estimate_binary(des$variables$flow_K3, des$variables$flow_K0, "C", "FLOW_FUNNEL", "K3/K0", "K0_DIAGNOSED_T2D"),
  estimate_binary(des$variables$flow_K3, des$variables$flow_K1, "C", "FLOW_FUNNEL", "K3/K1", "K1_OBSERVABLE_CKD"),
  estimate_binary(des$variables$flow_K3, des$variables$flow_K2, "C", "FLOW_FUNNEL", "K3/K2", "K2_BROAD_FLOW_RANGE")
)

# Module D.
d_orig <- bind_rows(lapply(cv_states, function(cv) bind_rows(lapply(k_states, function(ks) {
  estimate_binary(des$variables$original_cv == cv & des$variables$original_kidney == ks, t2d, "D", "ORIGINAL_CV_KIDNEY_MATRIX", paste(cv, ks, sep = "__"), "DIAGNOSED_T2D", list(cv_state = cv, kidney_state = ks))
})))) %>% select(cv_state, kidney_state, everything())

d_selected <- bind_rows(
  estimate_binary(des$variables$original_kidney == "NO_DIRECT", t2d & des$variables$original_cv == "DIRECT", "D", "SELECTED_ESTIMAND", "P_KIDNEY_NO_DIRECT_GIVEN_CV_DIRECT", "CV_DIRECT"),
  estimate_binary(des$variables$original_kidney == "DIRECT_COMPATIBLE", t2d & des$variables$original_cv == "DIRECT", "D", "SELECTED_ESTIMAND", "P_KIDNEY_DIRECT_GIVEN_CV_DIRECT", "CV_DIRECT"),
  estimate_binary(des$variables$original_kidney == "NO_DIRECT", t2d & des$variables$original_cv %in% c("DIRECT", "CONDITIONAL"), "D", "SELECTED_ESTIMAND", "P_KIDNEY_NO_DIRECT_GIVEN_CV_DIRECT_OR_CONDITIONAL", "CV_DIRECT_OR_CONDITIONAL"),
  estimate_binary(des$variables$original_cv %in% c("DIRECT", "CONDITIONAL") & des$variables$original_kidney == "NO_DIRECT", t2d, "D", "SELECTED_ESTIMAND", "P_CV_DIRECT_OR_CONDITIONAL_AND_KIDNEY_NO_DIRECT", "DIAGNOSED_T2D")
)

d_aug <- bind_rows(lapply(states, function(st) bind_rows(lapply(k_states, function(ks) {
  estimate_binary(des$variables$portfolio_status == st & des$variables$original_kidney == ks, t2d, "D", "AUGMENTED_CV_KIDNEY_MATRIX", paste(st, ks, sep = "__"), "DIAGNOSED_T2D", list(portfolio_state = st, kidney_state = ks))
})))) %>% select(portfolio_state, kidney_state, everything())

outputs <- list(
  MODULE_A_TRIAL_LEVEL_ALIGNMENT_RESULTS = a_trial,
  MODULE_A_PORTFOLIO_ALIGNMENT_RESULTS = a_port,
  MODULE_A_ORIGINAL_VS_AUGMENTED_MATRIX = a_matrix,
  MODULE_B_DEPTH_RESULTS = b_depth,
  MODULE_B_EXACT_DEPTH_INTERNAL = b_exact,
  MODULE_B_CUMULATIVE_ACCRUAL_RESULTS = b_accrual,
  MODULE_B_BREADTH_DEPTH_DECOMPOSITION = b_bd,
  MODULE_B_LEAVE_ONE_TRIAL_OUT_RESULTS = b_loto,
  MODULE_C_CV_CONDITIONAL_DECOMPOSITION = c_cv,
  MODULE_C_CV_PORTFOLIO_CONTEXT = c_context,
  MODULE_C_SELECT_BRANCH_RESULTS = c_select_branch,
  MODULE_C_SELECT_FUNNEL_RESULTS = c_select_funnel,
  MODULE_C_KIDNEY_CONDITIONAL_REASON_RESULTS = c_kidney,
  MODULE_C_FLOW_FUNNEL_RESULTS = c_flow,
  MODULE_D_ORIGINAL_CV_KIDNEY_MATRIX = d_orig,
  MODULE_D_SELECTED_ESTIMANDS = d_selected,
  MODULE_D_AUGMENTED_CV_KIDNEY_MATRIX = d_aug
)
for (nm in names(outputs)) readr::write_csv(outputs[[nm]], file.path(internal, paste0(nm, ".csv")), na = "")

diag <- tibble(
  metric = c("adult_mec_frame_n", "diagnosed_t2d_n", "bmi_ge25_n", "design_df", "strata_n", "psu_n"),
  value = c(nrow(des$variables), sum(t2d), sum(bmi25), survey::degf(des), length(unique(des$variables$SDMVSTRA)), nrow(unique(des$variables[, c("SDMVSTRA", "SDMVPSU")])) )
)
readr::write_csv(diag, file.path(internal, "PA05_SURVEY_DESIGN_DIAGNOSTICS_INTERNAL.csv"), na = "")
