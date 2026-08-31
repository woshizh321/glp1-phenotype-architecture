#!/usr/bin/env python3
"""Build publication participant-state inputs without weighted estimation."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ["PUBLIC_REPO_ROOT"])
SOURCE = Path(os.environ["PHASE1_INPUT_CSV"])
R1 = ROOT / "metadata"
PRE = ROOT / "metadata"

STATUSES = ["DIRECT", "POTENTIAL", "NOT_ALIGNED", "INDETERMINATE"]


def tri_age(s: pd.Series, threshold: float) -> pd.Series:
    return pd.Series(np.where(s.notna(), np.where(s >= threshold, "T", "F"), "U"), index=s.index)


def tri_binary_yes_no(s: pd.Series) -> pd.Series:
    return pd.Series(np.where(s.eq(1), "T", np.where(s.eq(2), "F", "U")), index=s.index)


def tri_ascvd(d: pd.DataFrame) -> pd.Series:
    x = d[["MCQ160C", "MCQ160E", "MCQ160F"]]
    return pd.Series(np.where(x.eq(1).any(axis=1), "T", np.where(x.eq(2).all(axis=1), "F", "U")), index=d.index)


def evaluate_criterion(d: pd.DataFrame, row: pd.Series) -> pd.Series:
    if row.observability_class == "UNAVAILABLE":
        return pd.Series("X", index=d.index)
    cid = str(row.criterion_id)
    desc = str(row.criterion_description)
    if "Diagnosed type 2 diabetes" in desc:
        return tri_binary_yes_no(d.DIQ010)
    if "Established coronary or cerebrovascular" in desc:
        return tri_ascvd(d)
    if desc.startswith("Age ≥"):
        threshold = float(desc.split("≥", 1)[1].split()[0])
        return tri_age(d.RIDAGEYR, threshold)
    if "BMI ≥25" in desc:
        valid = d.BMXBMI.notna() & d.BMXBMI.between(10, 100)
        return pd.Series(np.where(valid, np.where(d.BMXBMI >= 25, "T", "F"), "U"), index=d.index)
    if "eGFR ≥25 and <60" in desc:
        valid = d.egfr_2021.notna()
        return pd.Series(np.where(valid, np.where((d.egfr_2021 >= 25) & (d.egfr_2021 < 60), "T", "F"), "U"), index=d.index)
    if "eGFR <60" in desc:
        valid = d.egfr_2021.notna()
        return pd.Series(np.where(valid, np.where(d.egfr_2021 < 60, "T", "F"), "U"), index=d.index)
    if "Male age ≥50 or female age ≥55" in desc:
        valid = d.RIAGENDR.isin([1, 2]) & d.RIDAGEYR.notna()
        passed = (d.RIAGENDR.eq(1) & d.RIDAGEYR.ge(50)) | (d.RIAGENDR.eq(2) & d.RIDAGEYR.ge(55))
        return pd.Series(np.where(valid, np.where(passed, "T", "F"), "U"), index=d.index)
    raise ValueError(f"No evaluator for {cid}: {desc}")


def collapse_branch(values: pd.DataFrame) -> pd.Series:
    # R1 order: contradiction; unresolved observed/proxy; structural unavailable; direct.
    return pd.Series(
        np.where(values.eq("F").any(axis=1), "NOT_ALIGNED",
        np.where(values.eq("U").any(axis=1), "INDETERMINATE",
        np.where(values.eq("X").any(axis=1), "POTENTIAL", "DIRECT"))),
        index=values.index,
    )


def collapse_trial(values: pd.DataFrame) -> pd.Series:
    return pd.Series(
        np.where(values.eq("DIRECT").any(axis=1), "DIRECT",
        np.where(values.eq("POTENTIAL").any(axis=1), "POTENTIAL",
        np.where(values.eq("INDETERMINATE").any(axis=1), "INDETERMINATE", "NOT_ALIGNED"))),
        index=values.index,
    )


def norm_original_cv(s: pd.Series) -> pd.Series:
    m = {
        "DIRECT_DOMAIN_REPRESENTATION": "DIRECT",
        "CONDITIONAL_DOMAIN_REPRESENTATION": "CONDITIONAL",
        "NO_DIRECT_COMPLETED_DOMAIN_REPRESENTATION": "NO_DIRECT",
        "UNCLASSIFIED_DATA_LIMITATION": "UNCLASSIFIED",
    }
    return s.map(m)


def norm_original_kidney(s: pd.Series) -> pd.Series:
    m = {
        "DIRECT_DOMAIN_REPRESENTATION": "DIRECT_COMPATIBLE",
        "CONDITIONAL_DOMAIN_REPRESENTATION": "CONDITIONAL",
        "NO_DIRECT_COMPLETED_DOMAIN_REPRESENTATION": "NO_DIRECT",
        "UNCLASSIFIED_DATA_LIMITATION": "UNCLASSIFIED",
    }
    return s.map(m)


def build(out_root: Path) -> None:
    d = pd.read_csv(SOURCE)
    crosswalk = pd.read_csv(R1 / "trial_phenotype_crosswalk_public.csv")
    branch_logic = pd.read_csv(R1 / "trial_branch_logic_public.csv")
    milestones = pd.read_csv(PRE / "evidence_milestone_map.csv")
    kidney_map = pd.read_csv(PRE / "kidney_conditional_reason_map.csv")

    required = {"SEQN", "RIDAGEYR", "WTMECPRP", "SDMVSTRA", "SDMVPSU", "DIQ010", "BMXBMI", "egfr_2021", "uacr_valid"}
    if not required.issubset(d.columns):
        raise RuntimeError(f"Missing required columns: {sorted(required - set(d.columns))}")

    adult_design = d.RIDAGEYR.ge(20) & d.WTMECPRP.gt(0) & d.SDMVSTRA.notna() & d.SDMVPSU.notna()
    t2d = adult_design & d.DIQ010.eq(1)
    if int(t2d.sum()) != 1324:
        raise RuntimeError(f"Diagnosed-T2D identity failed: {int(t2d.sum())}")
    d["pa05_t2d"] = t2d
    d["pa05_bmi25"] = adult_design & d.BMXBMI.notna() & d.BMXBMI.ge(25)

    branch_status = {}
    criterion_cache = {}
    for branch_id in branch_logic.branch_id:
        rows = crosswalk[crosswalk.phenotype_branch_id.eq(branch_id)].drop_duplicates("criterion_id")
        if rows.empty:
            raise RuntimeError(f"No crosswalk rows for {branch_id}")
        vals = {}
        for _, row in rows.iterrows():
            criterion_cache.setdefault(row.criterion_id, evaluate_criterion(d, row))
            vals[row.criterion_id] = criterion_cache[row.criterion_id]
        branch_status[branch_id] = collapse_branch(pd.DataFrame(vals))

    branch_df = pd.DataFrame(branch_status)
    trial_status = {}
    trials = list(dict.fromkeys(branch_logic.trial.tolist()))
    if len(trials) != 10:
        raise RuntimeError(f"Expected ten trials, found {len(trials)}")
    for trial in trials:
        branches = branch_logic.loc[branch_logic.trial.eq(trial), "branch_id"].tolist()
        trial_status[trial] = collapse_trial(branch_df[branches])
        d[f"trial__{trial}"] = trial_status[trial]

    trial_df = pd.DataFrame(trial_status)
    d["portfolio_status"] = collapse_trial(trial_df)
    d["n_direct"] = trial_df.eq("DIRECT").sum(axis=1)
    d["n_direct_or_potential"] = trial_df.isin(["DIRECT", "POTENTIAL"]).sum(axis=1)
    d["portfolio_envelope"] = d.n_direct_or_potential.ge(1)

    d["original_cv"] = norm_original_cv(d["status__T2D_CARDIOVASCULAR_OUTCOME"])
    d["original_kidney"] = norm_original_kidney(d["status__T2D_KIDNEY_OUTCOME"])
    d["cv_conditional_class"] = np.where(d.original_cv.ne("CONDITIONAL"), "OUTSIDE",
        np.where(d.portfolio_status.eq("DIRECT"), "CV-C1",
        np.where(d.portfolio_status.eq("POTENTIAL"), "CV-C2",
        np.where(d.portfolio_status.eq("INDETERMINATE"), "CV-C4", "CV-C3"))))

    obesity_status = d["status__OBESITY_CARDIOVASCULAR_OUTCOME"]
    obesity_reason = d["reason__OBESITY_CARDIOVASCULAR_OUTCOME"]
    reason_to_branch = {
        "OBESCV_ASCVD_AGE20_LT45": "SELECT-C1",
        "OBESCV_ASCVD_BMI25_LT27": "SELECT-C2",
        "OBESCV_ASCVD_WITH_DIAGNOSED_DIABETES": "SELECT-C3",
        "PAD_STATUS_UNOBSERVABLE": "SELECT-C4",
    }
    d["select_branch"] = np.where(obesity_status.eq("CONDITIONAL_DOMAIN_REPRESENTATION"), obesity_reason.map(reason_to_branch), "OUTSIDE")
    d0 = adult_design & d.BMXBMI.notna() & d.BMXBMI.ge(25)
    d1 = d0 & d.RIDAGEYR.ge(45) & d.BMXBMI.ge(27)
    d2 = d1 & d.DIQ010.eq(2)
    ascvd = tri_ascvd(d)
    d["select_D0"] = d0
    d["select_D1"] = d1
    d["select_D2"] = d2
    d["select_D3"] = d2 & ascvd.eq("T")
    d["select_D3_PAD"] = d2 & ascvd.eq("F")
    d["select_D3_ENVELOPE"] = d.select_D3 | d.select_D3_PAD

    cell_to_reason = kidney_map.set_index("kidney_cell_id")["scientific_reason_class"]
    d["kidney_reason_class"] = d["kidney_cell_id__T2D_KIDNEY_OUTCOME"].map(cell_to_reason)
    valid_kidney = d.egfr_2021.notna() & d.uacr_valid.notna()
    d["flow_K0"] = t2d
    d["flow_K1"] = t2d & valid_kidney & (d.egfr_2021.lt(60) | d.uacr_valid.ge(30))
    d["flow_K2"] = d.flow_K1 & d.egfr_2021.between(25, 75, inclusive="both") & d.uacr_valid.lt(5000)
    d["flow_K3"] = t2d & d.flow_gateway_2021.fillna(False).astype(bool)
    if (d.flow_K3 & ~d.flow_K2).any() or (d.flow_K2 & ~d.flow_K1).any() or (d.flow_K1 & ~d.flow_K0).any():
        raise RuntimeError("FLOW_FUNNEL_EMPIRICAL_NESTING failed")

    year_trials = {}
    entered = []
    for row in milestones.itertuples(index=False):
        entered += str(row.trials_entering_bundle).split(";")
        year_trials[int(row.milestone_year)] = list(entered)
        cols = [f"trial__{t}" for t in entered]
        d[f"cum_direct_{int(row.milestone_year)}"] = d[cols].eq("DIRECT").any(axis=1)
        d[f"cum_dp_{int(row.milestone_year)}"] = d[cols].isin(["DIRECT", "POTENTIAL"]).any(axis=1)

    out_root.mkdir(parents=True, exist_ok=True)
    internal = out_root / "intermediate"
    internal.mkdir(exist_ok=True)
    branch_long = branch_df.loc[t2d].assign(SEQN=d.loc[t2d, "SEQN"].values).melt(id_vars="SEQN", var_name="branch_id", value_name="branch_status")
    branch_long.to_csv(internal / "PA05_PARTICIPANT_BRANCH_STATUS_INTERNAL.csv", index=False)
    d.to_csv(internal / "PA05_SURVEY_INPUT_INTERNAL.csv", index=False)
    pd.DataFrame({"milestone_year": list(year_trials), "trials_available": [";".join(x) for x in year_trials.values()]}).to_csv(internal / "PA05_MILESTONE_TRIALS_INTERNAL.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out_root", type=Path)
    build(ap.parse_args().out_root)
