#!/usr/bin/env python3
"""Derive the study analytical frame from a harmonized official NHANES CSV."""

from __future__ import annotations

import argparse
import importlib.util
import math
from pathlib import Path

import pandas as pd

SOURCE_VARIABLES = [
    "SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH3", "WTMECPRP", "WTINTPRP",
    "SDMVPSU", "SDMVSTRA", "DIQ010", "BMXBMI", "MCQ160C", "MCQ160D",
    "MCQ160E", "MCQ160F", "LBXSCR", "URDACT", "LBXGH",
]
DOMAINS = [
    "T2D_CARDIOVASCULAR_OUTCOME",
    "OBESITY_CARDIOVASCULAR_OUTCOME",
    "T2D_KIDNEY_OUTCOME",
]


def load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def egfr_category(value):
    if value is None or not math.isfinite(float(value)) or value <= 0:
        return None
    if value < 15: return "G0_LT15"
    if value < 25: return "G1_15_LT25"
    if value < 50: return "G2_25_LT50"
    if value < 60: return "G3_50_LT60"
    if value <= 75: return "G4_60_LE75"
    return "G5_GT75"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    args = parser.parse_args()
    mapping_dir = Path(__file__).resolve().parents[1] / "03_domain_mapping"
    classifier = load(mapping_dir / "classifier.py", "publication_classifier")
    contract = load(mapping_dir / "domain_contract.py", "publication_domain_contract")
    source = pd.read_csv(args.source_csv, low_memory=False)
    missing = sorted(set(SOURCE_VARIABLES) - set(source.columns))
    if missing:
        raise ValueError(f"Missing required NHANES variables: {missing}")

    rows = []
    for raw in source[SOURCE_VARIABLES].where(pd.notnull(source), None).to_dict("records"):
        row = dict(raw)
        row["age_valid"] = raw["RIDAGEYR"] if raw["RIDAGEYR"] is not None and 20 <= raw["RIDAGEYR"] <= 80 else None
        row["bmi_valid"] = raw["BMXBMI"] if raw["BMXBMI"] is not None and 11.9 <= raw["BMXBMI"] <= 92.3 else None
        row["uacr_valid"] = raw["URDACT"] if raw["URDACT"] is not None and 0.27 <= raw["URDACT"] <= 11676.92 else None
        row["age_group"] = contract.age_group(raw["RIDAGEYR"])
        row["sex_group"] = contract.sex_group(raw["RIAGENDR"])
        row["race_group"] = contract.race_ethnicity_group(raw["RIDRETH3"])
        row["diabetes_group"] = "DIAGNOSED_DIABETES" if raw["DIQ010"] == 1 else ("NO_DIAGNOSED_DIABETES" if raw["DIQ010"] == 2 else "UNCLASSIFIED")
        ascvd, _ = classifier.ascvd_proxy(raw["MCQ160C"], raw["MCQ160E"], raw["MCQ160F"])
        row["ascvd_group"] = "POSITIVE" if ascvd is True else ("NEGATIVE" if ascvd is False else "UNCLASSIFIED")
        e21 = classifier.ckd_epi_2021(raw["RIDAGEYR"], raw["RIAGENDR"], raw["LBXSCR"])
        e09n = classifier.ckd_epi_2009(raw["RIDAGEYR"], raw["RIAGENDR"], raw["LBXSCR"], False)
        e09r = classifier.ckd_epi_2009(raw["RIDAGEYR"], raw["RIAGENDR"], raw["LBXSCR"], raw["RIDRETH3"] == 4)
        row.update(egfr_2021=e21, egfr_2009_no_race=e09n, egfr_2009_historical_race=e09r,
                   egfr_cat_2021=egfr_category(e21), egfr_cat_2009_no_race=egfr_category(e09n),
                   egfr_cat_2009_historical_race=egfr_category(e09r))
        for label, value in [("2021", e21), ("2009_no_race", e09n), ("2009_historical_race", e09r)]:
            result = classifier.classify_kidney_values(1, value, raw["URDACT"]) if value is not None else None
            row[f"flow_gateway_{label}"] = None if result is None or raw["URDACT"] is None else result["evidence_status"] == classifier.DIRECT
            row[f"kidney_status_{label}"] = None if raw["DIQ010"] != 1 or result is None else result["evidence_status"]
            for threshold in contract.FLOW_EGFR_THRESHOLDS:
                row[f"near_{label}_{int(threshold)}"] = contract.flow_egfr_near_threshold(value, threshold)
        sensitivity_diabetes = contract.sensitivity_diabetes_row(raw)
        sensitivity_ascvd = contract.ascvd_angina_sensitivity(raw["MCQ160C"], raw["MCQ160D"], raw["MCQ160E"], raw["MCQ160F"])
        sensitivity_ascvd_row = dict(raw)
        sensitivity_ascvd_row.update(MCQ160C=1 if sensitivity_ascvd is True else (2 if sensitivity_ascvd is False else None),
                                     MCQ160E=2 if sensitivity_ascvd is not None else None,
                                     MCQ160F=2 if sensitivity_ascvd is not None else None)
        direct = classifier.classify_participant(raw)
        for index, domain in enumerate(DOMAINS):
            primary = contract.classify_with_domain_contract(raw, domain)
            s1 = contract.classify_with_domain_contract(sensitivity_diabetes, domain)
            s2 = contract.classify_with_domain_contract(sensitivity_ascvd_row, domain)
            row[f"domain_target_flag__{domain}"] = primary["domain_target_flag"]
            row[f"status__{domain}"] = primary["evidence_status"]
            row[f"reason__{domain}"] = primary["triggered_branch"]
            row[f"kidney_cell_id__{domain}"] = direct[index].get("kidney_cell_id")
            row[f"s1_target__{domain}"] = s1["domain_target_flag"]
            row[f"s1_status__{domain}"] = s1["evidence_status"]
            row[f"s2_status__{domain}"] = s2["evidence_status"]
        rows.append(row)
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    main()
