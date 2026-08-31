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

safe_num <- function(x) if (length(x) == 0 || !is.finite(x[[1]])) NA_real_ else as.numeric(x[[1]])

estimate_binary <- function(indicator, denominator, analysis, estimand, category, denominator_label) {
  indicator <- as.logical(indicator); denominator <- as.logical(denominator)
  keep <- !is.na(denominator) & denominator
  sd <- des[keep, ]
  z <- indicator[keep]
  stopifnot(length(z) == nrow(sd$variables), !anyNA(z), length(z) > 0)
  sd$variables$.cfwse_indicator <- z
  n <- length(z); events <- sum(z); df <- survey::degf(sd)
  mo <- survey::svymean(~.cfwse_indicator, sd, deff = TRUE, na.rm = FALSE)
  p <- safe_num(coef(mo)); se <- safe_num(SE(mo)); deff <- suppressWarnings(safe_num(survey::deff(mo)))
  cp <- tryCatch(survey::svyciprop(~.cfwse_indicator, sd, method = "beta", level = 0.95, df = df), error = function(e) NULL)
  ci <- c(NA_real_, NA_real_)
  if (!is.null(cp)) {
    p <- safe_num(coef(cp)); se <- safe_num(SE(cp)); ci <- as.numeric(confint(cp)[1, ])
  }
  rel <- phase1_reliability(p, ci[1], ci[2], n, events, deff, df)
  eff <- if (is.finite(deff) && deff > 0) min(n, n / deff) else if (events %in% c(0, n)) n else NA_real_
  tibble(
    analysis = analysis, estimand = estimand, category = category, denominator = denominator_label,
    unweighted_denominator_n = n, internal_unweighted_numerator = events,
    weighted_proportion = p, ci95_lower = ci[1], ci95_upper = ci[2], design_se = se,
    design_df = df, design_effect = deff, effective_n = eff,
    reliability_status = rel,
    isolated_reporting_status = ifelse(rel == "PRESENT", "PRESENT", "SUPPRESSED"),
    result_layer = "INTERNAL_AUDIT_ONLY"
  )
}

t2d <- des$variables$pa05_t2d
existing_ge5 <- t2d & des$variables$cfwse_existing_n_direct_ge5

bin_a <- function(x) ifelse(x == 0, "0", ifelse(x == 1, "1", ifelse(x <= 4, "2-4", ">=5")))
a <- bind_rows(lapply(c("0", "1", "2-4", ">=5"), function(k) {
  estimate_binary(
    bin_a(des$variables$cfwse_n_direct_minus_exscel) == k, t2d,
    "A_EXSCEL_REMOVED", "DIRECT_DEPTH_MINUS_EXSCEL", k, "FROZEN_DIAGNOSED_T2D"
  )
}))

bin_primary <- function(x) ifelse(x == 1, "1", ifelse(x == 2, "2", ">=3"))
b_primary <- bind_rows(lapply(c("1", "2", ">=3"), function(k) {
  estimate_binary(
    bin_primary(des$variables$cfwse_n_unique_direct_signatures) == k, existing_ge5,
    "B_PRIMARY_UNIQUE_SIGNATURE", "UNIQUE_SIGNATURE_DEPTH_WITHIN_EXISTING_GE5_TRIALS", k,
    "FROZEN_EXISTING_N_DIRECT_GE5"
  )
}))

b_exact <- bind_rows(lapply(0:6, function(k) {
  estimate_binary(
    des$variables$cfwse_n_unique_direct_signatures == k, t2d,
    "B_SECONDARY_EXACT_SIGNATURE", "EXACT_UNIQUE_SIGNATURE_DEPTH", as.character(k),
    "FROZEN_DIAGNOSED_T2D"
  )
}))

readr::write_csv(a, file.path(internal, "EXSCEL_REMOVED_DIRECT_DEPTH_RESULTS_INTERNAL.csv"), na = "")
readr::write_csv(b_primary, file.path(internal, "UNIQUE_SIGNATURE_PRIMARY_RESULTS_INTERNAL.csv"), na = "")
readr::write_csv(b_exact, file.path(internal, "UNIQUE_SIGNATURE_EXACT_DEPTH_RESULTS_INTERNAL.csv"), na = "")

checks <- bind_rows(
  tibble(analysis="A_EXSCEL_REMOVED", partition_sum=sum(a$weighted_proportion), tolerance=1e-12),
  tibble(analysis="B_PRIMARY_UNIQUE_SIGNATURE", partition_sum=sum(b_primary$weighted_proportion), tolerance=1e-12),
  tibble(analysis="B_SECONDARY_EXACT_SIGNATURE", partition_sum=sum(b_exact$weighted_proportion), tolerance=1e-12)
) %>% mutate(absolute_discrepancy=abs(partition_sum-1), status=ifelse(absolute_discrepancy<=tolerance,"PASS","FAIL"))
readr::write_csv(checks, file.path(internal, "WEIGHTED_PARTITION_ARITHMETIC_CHECKS_INTERNAL.csv"), na = "")

diag <- tibble(
  metric=c("adult_mec_frame_n","diagnosed_t2d_n","existing_n_direct_ge5_n","design_df","strata_n","psu_n"),
  value=c(nrow(des$variables),sum(t2d),sum(existing_ge5),survey::degf(des),length(unique(des$variables$SDMVSTRA)),nrow(unique(des$variables[,c("SDMVSTRA","SDMVPSU")])) )
)
readr::write_csv(diag, file.path(internal, "WEIGHTED_SENSITIVITY_DESIGN_DIAGNOSTICS_INTERNAL.csv"), na = "")

session <- capture.output(sessionInfo())
writeLines(session, file.path(internal, "R_SESSION_INFO.txt"))
