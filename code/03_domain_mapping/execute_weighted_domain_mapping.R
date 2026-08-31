#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(survey)
  library(readr)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
input_csv <- args[[1]]
out_root <- args[[2]]
repo_root <- Sys.getenv("PUBLIC_REPO_ROOT")
stopifnot(nzchar(repo_root))
source(file.path(repo_root, "code/08_survey_reporting/survey_engine.R"))

options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
options(warn = 1)
d <- readr::read_csv(input_csv, show_col_types = FALSE, progress = FALSE)

domains <- c("T2D_CARDIOVASCULAR_OUTCOME", "OBESITY_CARDIOVASCULAR_OUTCOME", "T2D_KIDNEY_OUTCOME")
statuses <- c("DIRECT_DOMAIN_REPRESENTATION", "CONDITIONAL_DOMAIN_REPRESENTATION", "NO_DIRECT_COMPLETED_DOMAIN_REPRESENTATION", "UNCLASSIFIED_DATA_LIMITATION")
target_states <- c("TRUE", "FALSE", "UNCLASSIFIED")

dir.create(file.path(out_root, "estimates"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "sensitivities"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "subgroups"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "qc"), recursive = TRUE, showWarnings = FALSE)

design_mec <- phase1_build_design(d, "WTMECPRP")
design_int <- phase1_build_design(d, "WTINTPRP")

safe_num <- function(x) if (length(x) == 0 || !is.finite(x[[1]])) NA_real_ else as.numeric(x[[1]])

estimate_prop <- function(des, indicator, domain, estimand, category, denominator_label) {
  indicator <- as.logical(indicator)
  stopifnot(length(indicator) == nrow(des$variables), !anyNA(indicator))
  des$variables$.phase1_indicator <- indicator
  den_n <- length(indicator)
  event_n <- sum(indicator)
  df <- survey::degf(des)
  mean_obj <- survey::svymean(~.phase1_indicator, des, deff = TRUE, na.rm = FALSE)
  p <- safe_num(coef(mean_obj))
  se <- safe_num(SE(mean_obj))
  deff <- suppressWarnings(safe_num(survey::deff(mean_obj)))
  ci <- c(NA_real_, NA_real_)
  cp <- tryCatch(
    survey::svyciprop(~.phase1_indicator, des, method = "beta", level = 0.95, df = df),
    error = function(e) NULL
  )
  if (!is.null(cp)) {
    p <- safe_num(coef(cp)); se <- safe_num(SE(cp)); ci <- as.numeric(confint(cp)[1, ])
  }
  reliability <- phase1_reliability(p, ci[1], ci[2], den_n, event_n, deff, df)
  effective_n <- if (is.finite(deff) && deff > 0) min(den_n, den_n / deff) else if (event_n %in% c(0, den_n)) den_n else NA_real_
  tibble(
    domain = domain, estimand = estimand, denominator = denominator_label,
    category = category, unweighted_denominator_n = den_n,
    unweighted_event_n = event_n, weighted_proportion = p, design_se = se,
    ci95_lower = ci[1], ci95_upper = ci[2], design_df = df,
    design_effect = deff, effective_n = effective_n,
    reliability_flag = reliability,
    suppression_flag = ifelse(reliability == "PRESENT", "PRESENT", "NOT_FOR_MANUSCRIPT_REPORTING"),
    suppression_reason = ifelse(reliability == "PRESENT", "", reliability)
  )
}

estimate_mean <- function(des, values, domain, display_row, denominator_label) {
  valid <- is.finite(values)
  full_n <- length(values)
  valid_n <- sum(valid)
  missing_n <- full_n - valid_n
  subdes <- des[valid, ]
  subdes$variables$.phase1_value <- values[valid]
  df <- survey::degf(subdes)
  obj <- survey::svymean(~.phase1_value, subdes, deff = TRUE, na.rm = FALSE)
  est <- safe_num(coef(obj)); se <- safe_num(SE(obj)); deff <- suppressWarnings(safe_num(survey::deff(obj)))
  ci <- est + c(-1, 1) * stats::qt(0.975, df = df) * se
  eff <- if (is.finite(deff) && deff > 0) min(valid_n, valid_n / deff) else NA_real_
  rel <- if (valid_n < 30) "SUPPRESS_DENOMINATOR_N_LT30" else if (!is.finite(eff) || eff < 30) "SUPPRESS_EFFECTIVE_N_LT30" else if (df < 8) "SUPPRESS_DF_LT8_AUTOMATED_REVIEW_RULE" else if (!all(is.finite(c(est, se, ci)))) "SUPPRESS_NONFINITE_CONTINUOUS_ESTIMATE" else "PRESENT"
  miss <- estimate_prop(des, !valid, domain, paste0("TABLE2_MISSINGNESS_", display_row), "MISSING_OR_INVALID", denominator_label)
  tibble(
    domain = domain, display_row = display_row, row_type = "CONTINUOUS_MEAN", category = "MEAN",
    denominator = denominator_label, unweighted_denominator_n = full_n,
    unweighted_valid_n = valid_n, unweighted_missing_n = missing_n,
    unweighted_event_n = NA_integer_, weighted_estimate = est, design_se = se,
    ci95_lower = ci[1], ci95_upper = ci[2], design_df = df,
    design_effect = deff, effective_n = eff, weighted_missing_proportion = miss$weighted_proportion,
    reliability_flag = rel, suppression_flag = ifelse(rel == "PRESENT", "PRESENT", "NOT_FOR_MANUSCRIPT_REPORTING"),
    suppression_reason = ifelse(rel == "PRESENT", "", rel), p_value = NA_real_
  )
}

estimate_categorical <- function(des, values, categories, domain, display_row, denominator_label) {
  bind_rows(lapply(categories, function(cat) {
    z <- estimate_prop(des, values == cat, domain, paste0("TABLE2_", display_row), cat, denominator_label)
    transmute(z,
      domain = domain, display_row = display_row, row_type = "CATEGORICAL_PROPORTION", category = cat,
      denominator = denominator, unweighted_denominator_n = unweighted_denominator_n,
      unweighted_valid_n = NA_integer_, unweighted_missing_n = NA_integer_, unweighted_event_n = unweighted_event_n,
      weighted_estimate = weighted_proportion, design_se = design_se, ci95_lower = ci95_lower,
      ci95_upper = ci95_upper, design_df = design_df, design_effect = design_effect,
      effective_n = effective_n, weighted_missing_proportion = NA_real_, reliability_flag = reliability_flag,
      suppression_flag = suppression_flag, suppression_reason = suppression_reason, p_value = NA_real_)
  }))
}

target_rows <- list(); status_rows <- list(); reason_rows <- list()
for (domain in domains) {
  tcol <- paste0("domain_target_flag__", domain)
  scol <- paste0("status__", domain)
  rcol <- paste0("reason__", domain)
  for (state in target_states) {
    target_rows[[length(target_rows)+1]] <- estimate_prop(design_mec, design_mec$variables[[tcol]] == state, domain, "TARGET_ASCERTAINMENT", state, "AGE20_MEC_ANALYTIC_FRAME")
  }
  td <- phase1_subset_target(design_mec, domain)
  for (status in statuses) {
    status_rows[[length(status_rows)+1]] <- estimate_prop(td, td$variables[[scol]] == status, domain, "PRIMARY_EVIDENCE_STATUS", status, "DOMAIN_TARGET_TRUE")
  }
  cell_column <- if (domain == "T2D_KIDNEY_OUTCOME") "kidney_cell_id__T2D_KIDNEY_OUTCOME" else rcol
  reasons <- if (domain == "T2D_KIDNEY_OUTCOME") {
    readr::read_csv(file.path(repo_root, "metadata/kidney_cell_registry.csv"), show_col_types=FALSE)$kidney_cell_id
  } else sort(unique(td$variables[[cell_column]]))
  reasons <- reasons[!is.na(reasons) & reasons != ""]
  for (reason in reasons) {
    cell_indicator <- !is.na(td$variables[[cell_column]]) & td$variables[[cell_column]] == reason
    reason_rows[[length(reason_rows)+1]] <- estimate_prop(td, cell_indicator, domain, ifelse(domain == "T2D_KIDNEY_OUTCOME", "KIDNEY_30_CELL", "PHENOTYPE_REASON_CELL"), reason, "DOMAIN_TARGET_TRUE")
  }
}
target_out <- bind_rows(target_rows)
status_out <- bind_rows(status_rows)
cell_out <- bind_rows(reason_rows)

readr::write_csv(target_out, file.path(out_root, "estimates/PHASE1_TARGET_ASCERTAINMENT_ESTIMATES.csv"), na = "")
readr::write_csv(status_out, file.path(out_root, "estimates/PHASE1_EVIDENCE_STATUS_ESTIMATES.csv"), na = "")
readr::write_csv(cell_out, file.path(out_root, "estimates/PHASE1_PHENOTYPE_CELL_ESTIMATES.csv"), na = "")

# Prespecified Table 2 characteristics.
char_rows <- list()
age_cats <- c("20-44", "45-64", "65-74", "75+", "UNCLASSIFIED")
sex_cats <- c("Male", "Female", "UNCLASSIFIED")
race_cats <- c("Mexican American", "Other Hispanic", "Non-Hispanic White", "Non-Hispanic Black", "Non-Hispanic Asian", "Other race including multiracial", "UNCLASSIFIED")
diab_cats <- c("DIAGNOSED_DIABETES", "NO_DIAGNOSED_DIABETES", "UNCLASSIFIED")
ascvd_cats <- c("POSITIVE", "NEGATIVE", "UNCLASSIFIED")
for (domain in domains) {
  td <- phase1_subset_target(design_mec, domain)
  label <- paste0(domain, " target TRUE")
  char_rows[[length(char_rows)+1]] <- estimate_mean(td, td$variables$age_valid, domain, "age_mean", label)
  char_rows[[length(char_rows)+1]] <- estimate_categorical(td, td$variables$age_group, age_cats, domain, "age_group", label)
  char_rows[[length(char_rows)+1]] <- estimate_categorical(td, td$variables$sex_group, sex_cats, domain, "sex", label)
  char_rows[[length(char_rows)+1]] <- estimate_categorical(td, td$variables$race_group, race_cats, domain, "race_ethnicity", label)
  char_rows[[length(char_rows)+1]] <- estimate_mean(td, td$variables$bmi_valid, domain, "bmi_mean", label)
  char_rows[[length(char_rows)+1]] <- estimate_categorical(td, td$variables$diabetes_group, diab_cats, domain, "diabetes_phenotype", label)
  if (domain != "T2D_KIDNEY_OUTCOME") {
    char_rows[[length(char_rows)+1]] <- estimate_categorical(td, td$variables$ascvd_group, ascvd_cats, domain, "ascvd_proxy", label)
  } else {
    char_rows[[length(char_rows)+1]] <- estimate_mean(td, td$variables$egfr_2021, domain, "egfr_2021_mean", label)
    char_rows[[length(char_rows)+1]] <- estimate_mean(td, td$variables$uacr_valid, domain, "uacr_mean", label)
  }
}
char_out <- bind_rows(char_rows)
readr::write_csv(char_out, file.path(out_root, "estimates/PHASE1_DOMAIN_CHARACTERISTICS.csv"), na = "")

# S1 and S2: prespecified reclassifications, alongside primary estimates.
run_sensitivity <- function(prefix, domains_use, target_prefix, status_prefix, out_file) {
  rows <- list()
  for (domain in domains_use) {
    p_tcol <- paste0("domain_target_flag__", domain); p_scol <- paste0("status__", domain)
    s_tcol <- paste0(target_prefix, domain); s_scol <- paste0(status_prefix, domain)
    for (analysis in c("PRIMARY", prefix)) {
      tcol <- if (analysis == "PRIMARY") p_tcol else s_tcol
      scol <- if (analysis == "PRIMARY") p_scol else s_scol
      for (state in target_states) {
        z <- estimate_prop(design_mec, design_mec$variables[[tcol]] == state, domain, "TARGET_ASCERTAINMENT", state, "AGE20_MEC_ANALYTIC_FRAME")
        z$analysis <- analysis; rows[[length(rows)+1]] <- z
      }
      keep <- design_mec$variables[[tcol]] == "TRUE"; td <- design_mec[!is.na(keep) & keep, ]
      for (status in statuses) {
        z <- estimate_prop(td, td$variables[[scol]] == status, domain, "EVIDENCE_STATUS", status, "DOMAIN_TARGET_TRUE")
        z$analysis <- analysis; rows[[length(rows)+1]] <- z
      }
    }
  }
  out <- bind_rows(rows) %>% select(analysis, everything())
  primary <- out %>% filter(analysis == "PRIMARY") %>% select(domain, estimand, category, primary_proportion = weighted_proportion)
  out <- out %>% left_join(primary, by = c("domain", "estimand", "category")) %>% mutate(difference_percentage_points = 100 * (weighted_proportion - primary_proportion))
  readr::write_csv(out, out_file, na = "")
}
run_sensitivity("S1", domains, "s1_target__", "s1_status__", file.path(out_root, "sensitivities/PHASE1_S1_DIABETES_SENSITIVITY.csv"))
run_sensitivity("S2", c("T2D_CARDIOVASCULAR_OUTCOME", "OBESITY_CARDIOVASCULAR_OUTCOME"), "domain_target_flag__", "s2_status__", file.path(out_root, "sensitivities/PHASE1_S2_ASCVD_SENSITIVITY.csv"))

# S4: same phenotype/classifier, separate MEC and interview designs.
s4_rows <- list()
for (analysis in c("PRIMARY_MEC_WTMECPRP", "S4_INTERVIEW_WTINTPRP")) {
  des <- if (analysis == "PRIMARY_MEC_WTMECPRP") design_mec else design_int
  domain <- "T2D_CARDIOVASCULAR_OUTCOME"; tcol <- paste0("domain_target_flag__", domain); scol <- paste0("status__", domain)
  for (state in target_states) {
    z <- estimate_prop(des, des$variables[[tcol]] == state, domain, "TARGET_ASCERTAINMENT", state, ifelse(analysis == "PRIMARY_MEC_WTMECPRP", "AGE20_MEC_FRAME", "AGE20_INTERVIEW_FRAME")); z$analysis <- analysis; s4_rows[[length(s4_rows)+1]] <- z
  }
  keep <- des$variables[[tcol]] == "TRUE"; td <- des[!is.na(keep) & keep, ]
  for (status in statuses) {
    z <- estimate_prop(td, td$variables[[scol]] == status, domain, "EVIDENCE_STATUS", status, ifelse(analysis == "PRIMARY_MEC_WTMECPRP", "MEC_TARGET_TRUE", "INTERVIEW_TARGET_TRUE")); z$analysis <- analysis; s4_rows[[length(s4_rows)+1]] <- z
  }
}
s4 <- bind_rows(s4_rows) %>% select(analysis, everything())
s4p <- s4 %>% filter(analysis == "PRIMARY_MEC_WTMECPRP") %>% select(estimand, category, primary_proportion = weighted_proportion)
s4 <- s4 %>% left_join(s4p, by = c("estimand", "category")) %>% mutate(difference_percentage_points = 100 * (weighted_proportion - primary_proportion))
readr::write_csv(s4, file.path(out_root, "sensitivities/PHASE1_S4_SURVEY_FRAME_SENSITIVITY.csv"), na = "")

# S3: paired complete T2D records, weighted eGFR-category/gateway concordance and threshold-near estimates.
t2d_keep <- design_mec$variables$domain_target_flag__T2D_KIDNEY_OUTCOME == "TRUE"
t2d_des <- design_mec[!is.na(t2d_keep) & t2d_keep, ]
s3_rows <- list(); reclass_rows <- list()
for (variant in c("2009_NO_RACE", "2009_HISTORICAL_RACE")) {
  c09 <- if (variant == "2009_NO_RACE") "egfr_2009_no_race" else "egfr_2009_historical_race"
  cat09 <- if (variant == "2009_NO_RACE") "egfr_cat_2009_no_race" else "egfr_cat_2009_historical_race"
  gate09 <- if (variant == "2009_NO_RACE") "flow_gateway_2009_no_race" else "flow_gateway_2009_historical_race"
  status09 <- if (variant == "2009_NO_RACE") "kidney_status_2009_no_race" else "kidney_status_2009_historical_race"
  complete <- is.finite(t2d_des$variables$egfr_2021) & is.finite(t2d_des$variables[[c09]]) & is.finite(t2d_des$variables$uacr_valid)
  pd <- t2d_des[complete, ]
  # Mean paired difference.
  pd$variables$.diff <- pd$variables$egfr_2021 - pd$variables[[c09]]
  mo <- survey::svymean(~.diff, pd, deff = TRUE); df <- survey::degf(pd); est <- safe_num(coef(mo)); se <- safe_num(SE(mo)); ci <- est + c(-1,1)*qt(.975,df)*se
  s3_rows[[length(s3_rows)+1]] <- tibble(variant=variant, record_type="CONTINUOUS_DIFFERENCE_2021_MINUS_2009", category_2021="", category_2009="", threshold=NA_real_, equation="2021_MINUS_2009", unweighted_denominator_n=nrow(pd$variables), unweighted_event_n=NA_integer_, weighted_estimate=est, design_se=se, ci95_lower=ci[1], ci95_upper=ci[2], design_df=df, reliability_flag="PRESENT", suppression_flag="PRESENT")
  for (a in sort(unique(pd$variables$egfr_cat_2021))) for (b in sort(unique(pd$variables[[cat09]]))) {
    z <- estimate_prop(pd, pd$variables$egfr_cat_2021 == a & pd$variables[[cat09]] == b, "T2D_KIDNEY_OUTCOME", "S3_EGFR_CATEGORY_CROSSTAB", paste(a,b,sep="__"), "PAIRED_COMPLETE_T2D")
    s3_rows[[length(s3_rows)+1]] <- transmute(z, variant=variant, record_type="EGFR_CATEGORY_CROSSTAB", category_2021=a, category_2009=b, threshold=NA_real_, equation="", unweighted_denominator_n, unweighted_event_n, weighted_estimate=weighted_proportion, design_se, ci95_lower, ci95_upper, design_df, reliability_flag, suppression_flag)
  }
  for (a in c(FALSE, TRUE)) for (b in c(FALSE, TRUE)) {
    z <- estimate_prop(pd, pd$variables$flow_gateway_2021 == a & pd$variables[[gate09]] == b, "T2D_KIDNEY_OUTCOME", "S3_FLOW_GATEWAY_CONCORDANCE", paste(a,b,sep="__"), "PAIRED_COMPLETE_T2D")
    reclass_rows[[length(reclass_rows)+1]] <- transmute(z, variant=variant, record_type="FLOW_GATEWAY", primary_2021=as.character(a), sensitivity_2009=as.character(b), direction=case_when(a==b ~ "CONCORDANT", !a & b ~ "2009_DIRECT_ONLY", a & !b ~ "2021_DIRECT_ONLY"), unweighted_denominator_n, unweighted_event_n, weighted_proportion, design_se, ci95_lower, ci95_upper, design_df, reliability_flag, suppression_flag)
  }
  for (a in statuses) for (b in statuses) {
    z <- estimate_prop(pd, pd$variables$status__T2D_KIDNEY_OUTCOME == a & pd$variables[[status09]] == b, "T2D_KIDNEY_OUTCOME", "S3_STATUS_RECLASSIFICATION", paste(a,b,sep="__"), "PAIRED_COMPLETE_T2D")
    reclass_rows[[length(reclass_rows)+1]] <- transmute(z, variant=variant, record_type="STATUS", primary_2021=a, sensitivity_2009=b, direction=ifelse(a==b,"CONCORDANT","RECLASSIFIED"), unweighted_denominator_n, unweighted_event_n, weighted_proportion, design_se, ci95_lower, ci95_upper, design_df, reliability_flag, suppression_flag)
  }
}
for (eq in c("2021", "2009_NO_RACE", "2009_HISTORICAL_RACE")) {
  for (threshold in c(15,25,50,60,75)) {
    col <- paste0("near_", tolower(eq), "_", threshold)
    complete <- !is.na(t2d_des$variables[[col]])
    pd <- t2d_des[complete, ]
    z <- estimate_prop(pd, pd$variables[[col]], "T2D_KIDNEY_OUTCOME", "S3_NEAR_THRESHOLD", as.character(threshold), "T2D_COMPLETE_EGFR")
    s3_rows[[length(s3_rows)+1]] <- transmute(z, variant="ALL", record_type="NEAR_THRESHOLD", category_2021="", category_2009="", threshold=threshold, equation=eq, unweighted_denominator_n, unweighted_event_n, weighted_estimate=weighted_proportion, design_se, ci95_lower, ci95_upper, design_df, reliability_flag, suppression_flag)
  }
}
readr::write_csv(bind_rows(s3_rows), file.path(out_root, "sensitivities/PHASE1_S3_EGFR_CONCORDANCE.csv"), na = "")
readr::write_csv(bind_rows(reclass_rows), file.path(out_root, "sensitivities/PHASE1_S3_FLOW_RECLASSIFICATION.csv"), na = "")

# Prespecified descriptive subgroup status distributions.
run_subgroups <- function(group_col, group_values, filename) {
  rows <- list()
  for (domain in domains) {
    tcol <- paste0("domain_target_flag__", domain); scol <- paste0("status__", domain)
    for (g in group_values) {
      keep <- design_mec$variables[[tcol]] == "TRUE" & design_mec$variables[[group_col]] == g
      sd <- design_mec[!is.na(keep) & keep, ]
      for (status in statuses) {
        z <- estimate_prop(sd, sd$variables[[scol]] == status, domain, paste0("SUBGROUP_", group_col), status, paste0("TARGET_TRUE_AND_",group_col,"=",g)); z$subgroup <- g; rows[[length(rows)+1]] <- z
      }
    }
  }
  readr::write_csv(bind_rows(rows) %>% select(subgroup, everything()), filename, na = "")
}
run_subgroups("age_group", c("20-44","45-64","65-74","75+"), file.path(out_root,"subgroups/PHASE1_SUBGROUP_AGE.csv"))
run_subgroups("sex_group", c("Male","Female"), file.path(out_root,"subgroups/PHASE1_SUBGROUP_SEX.csv"))
run_subgroups("race_group", race_cats[1:6], file.path(out_root,"subgroups/PHASE1_SUBGROUP_RACE_ETHNICITY.csv"))

# Publication-facing tables retain suppression flags; rendering blanks is deferred to figure/table formatting.
table2_public <- char_out %>% mutate(
  weighted_estimate = ifelse(suppression_flag == "PRESENT", weighted_estimate, NA_real_),
  design_se = ifelse(suppression_flag == "PRESENT", design_se, NA_real_),
  ci95_lower = ifelse(suppression_flag == "PRESENT", ci95_lower, NA_real_),
  ci95_upper = ifelse(suppression_flag == "PRESENT", ci95_upper, NA_real_),
  weighted_missing_proportion = ifelse(suppression_flag == "PRESENT", weighted_missing_proportion, NA_real_)
)
table3_public <- status_out %>% mutate(
  weighted_proportion = ifelse(suppression_flag == "PRESENT", weighted_proportion, NA_real_),
  design_se = ifelse(suppression_flag == "PRESENT", design_se, NA_real_),
  ci95_lower = ifelse(suppression_flag == "PRESENT", ci95_lower, NA_real_),
  ci95_upper = ifelse(suppression_flag == "PRESENT", ci95_upper, NA_real_),
  design_effect = ifelse(suppression_flag == "PRESENT", design_effect, NA_real_),
  effective_n = ifelse(suppression_flag == "PRESENT", effective_n, NA_real_)
)
readr::write_csv(table2_public, file.path(out_root, "tables/TABLE2_DOMAIN_CHARACTERISTICS.csv"), na = "")
readr::write_csv(table3_public, file.path(out_root, "tables/TABLE3_EVIDENCE_STATUS.csv"), na = "")

all_reliability <- bind_rows(
  target_out %>% mutate(source_file="PHASE1_TARGET_ASCERTAINMENT_ESTIMATES.csv"),
  status_out %>% mutate(source_file="PHASE1_EVIDENCE_STATUS_ESTIMATES.csv"),
  cell_out %>% mutate(source_file="PHASE1_PHENOTYPE_CELL_ESTIMATES.csv")
) %>% select(source_file, everything())
readr::write_csv(all_reliability, file.path(out_root, "estimates/PHASE1_RELIABILITY_FLAGS.csv"), na = "")

diag <- tibble(
  metric=c("release_n","age20_n","mec_frame_n","interview_frame_n","mec_strata_n","mec_psu_stratum_n","mec_design_df","mec_lonely_strata_n"),
  value=c(nrow(d),sum(d$RIDAGEYR>=20,na.rm=TRUE),nrow(design_mec$variables),nrow(design_int$variables),length(unique(design_mec$variables$SDMVSTRA)),nrow(unique(design_mec$variables[c("SDMVSTRA","SDMVPSU")])),survey::degf(design_mec),sum(table(design_mec$variables$SDMVSTRA)==1))
)
readr::write_csv(diag, file.path(out_root,"qc/PHASE1_SURVEY_DESIGN_DIAGNOSTICS.csv"), na="")

cat("PHASE1_WEIGHTED_R_EXECUTION_COMPLETE\n")
