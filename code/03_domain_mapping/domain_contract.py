#!/usr/bin/env python3
"""Publication denominator wrapper and sensitivity definitions.

This module performs classification logic only. It never computes survey-weighted
estimates. The publication classifier remains the evidence-status authority; this wrapper
separates analytic-frame membership and domain-target membership from status.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
from typing import Any


R2R_CLASSIFIER_PATH = Path(__file__).with_name("classifier.py")
_spec = importlib.util.spec_from_file_location("public_domain_classifier", R2R_CLASSIFIER_PATH)
if _spec is None or _spec.loader is None:
    raise ImportError(f"Cannot load frozen publication classifier: {R2R_CLASSIFIER_PATH}")
_classifier = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_classifier)

CLASSIFIER_VERSION = _classifier.CLASSIFIER_VERSION
DIRECT = _classifier.DIRECT
CONDITIONAL = _classifier.CONDITIONAL
NO_DIRECT = _classifier.NO_DIRECT
UNCLASSIFIED = _classifier.UNCLASSIFIED
T2D_CV = _classifier.T2D_CV
OBESITY_CV = _classifier.OBESITY_CV
T2D_KIDNEY = _classifier.T2D_KIDNEY

TARGET_TRUE = "TRUE"
TARGET_FALSE = "FALSE"
TARGET_UNCLASSIFIED = "UNCLASSIFIED"

FLOW_EGFR_THRESHOLDS = (15.0, 25.0, 50.0, 60.0, 75.0)
FLOW_EGFR_NEAR_THRESHOLD_WINDOW = 1.0


def _number(value: Any) -> float | None:
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _code(value: Any, allowed: set[int]) -> int | None:
    number = _number(value)
    if number is None or not number.is_integer() or int(number) not in allowed:
        return None
    return int(number)


def analytic_frame_flag(row: dict[str, Any], weight: str = "WTMECPRP") -> bool:
    """Primary/interview frame flag: age >=20 plus positive weight and design."""
    age = _number(row.get("RIDAGEYR"))
    survey_weight = _number(row.get(weight))
    psu = _number(row.get("SDMVPSU"))
    stratum = _number(row.get("SDMVSTRA"))
    return bool(age is not None and age >= 20 and survey_weight is not None and survey_weight > 0 and psu is not None and stratum is not None)


def domain_target_flag(row: dict[str, Any], domain: str) -> str:
    """Return the target-membership state, independent of evidence status."""
    if domain in {T2D_CV, T2D_KIDNEY}:
        diabetes = _code(row.get("DIQ010"), {1, 2})
        if diabetes is None:
            return TARGET_UNCLASSIFIED
        return TARGET_TRUE if diabetes == 1 else TARGET_FALSE
    if domain == OBESITY_CV:
        bmi = _number(row.get("BMXBMI"))
        if bmi is None or not 11.9 <= bmi <= 92.3:
            return TARGET_UNCLASSIFIED
        return TARGET_TRUE if bmi >= 25 else TARGET_FALSE
    raise ValueError(f"Unknown quantitative domain: {domain}")


def classify_with_domain_contract(row: dict[str, Any], domain: str) -> dict[str, Any]:
    """Assign status only inside target TRUE; retain explicit target uncertainty."""
    target = domain_target_flag(row, domain)
    output = {
        "domain": domain,
        "domain_target_flag": target,
        "evidence_status": None,
        "triggered_branch": None,
        "alignment_type": None,
        "alignment_reason": None,
        "classifier_version": CLASSIFIER_VERSION,
    }
    if target != TARGET_TRUE:
        return output
    if domain == T2D_CV:
        result = _classifier.classify_t2d_cv(row)
    elif domain == OBESITY_CV:
        result = _classifier.classify_obesity_cv(row)
    else:
        result = _classifier.classify_kidney(row)
    output.update(
        evidence_status=result["evidence_status"],
        triggered_branch=result["triggered_branch"],
        alignment_type=result.get("alignment_type") or None,
        alignment_reason=result.get("alignment_reason") or None,
    )
    return output


def augmented_diabetes_code(row: dict[str, Any]) -> int | None:
    """S1: diagnosed diabetes OR HbA1c >=6.5%; false requires both negatives."""
    diagnosed = _code(row.get("DIQ010"), {1, 2})
    hba1c = _number(row.get("LBXGH"))
    hba1c_valid = hba1c is not None and 2.0 <= hba1c <= 20.0
    if diagnosed == 1 or (hba1c_valid and hba1c >= 6.5):
        return 1
    if diagnosed == 2 and hba1c_valid and hba1c < 6.5:
        return 2
    return None


def sensitivity_diabetes_row(row: dict[str, Any]) -> dict[str, Any]:
    result = dict(row)
    result["DIQ010"] = augmented_diabetes_code(row)
    return result


def ascvd_angina_sensitivity(c: Any, d: Any, e: Any, f: Any) -> bool | None:
    """S2: four-component three-valued OR; PAD remains unobserved."""
    values = [_code(value, {1, 2}) for value in (c, d, e, f)]
    if 1 in values:
        return True
    if values == [2, 2, 2, 2]:
        return False
    return None


def flow_egfr_near_threshold(egfr: Any, threshold: Any) -> bool | None:
    """S3: exact prespecified near-threshold flag; None means invalid eGFR."""
    egfr_value = _number(egfr)
    threshold_value = _number(threshold)
    if threshold_value not in FLOW_EGFR_THRESHOLDS:
        raise ValueError(f"S3 threshold must be one of {FLOW_EGFR_THRESHOLDS}")
    if egfr_value is None or egfr_value <= 0:
        return None
    return abs(egfr_value - threshold_value) <= FLOW_EGFR_NEAR_THRESHOLD_WINDOW


def flow_egfr_near_threshold_flags(egfr: Any) -> dict[str, bool | None]:
    """Return the five prespecified S3 threshold-window flags for one eGFR value."""
    return {
        str(int(threshold)): flow_egfr_near_threshold(egfr, threshold)
        for threshold in FLOW_EGFR_THRESHOLDS
    }


def age_group(age: Any) -> str:
    value = _number(age)
    if value is None or value < 20:
        return "UNCLASSIFIED"
    if value < 45:
        return "20-44"
    if value < 65:
        return "45-64"
    if value < 75:
        return "65-74"
    return "75+"


def sex_group(code: Any) -> str:
    return {1: "Male", 2: "Female"}.get(_code(code, {1, 2}), "UNCLASSIFIED")


def race_ethnicity_group(code: Any) -> str:
    labels = {
        1: "Mexican American",
        2: "Other Hispanic",
        3: "Non-Hispanic White",
        4: "Non-Hispanic Black",
        6: "Non-Hispanic Asian",
        7: "Other race including multiracial",
    }
    return labels.get(_code(code, set(labels)), "UNCLASSIFIED")


def all_domain_assignments(row: dict[str, Any]) -> list[dict[str, Any]]:
    return [classify_with_domain_contract(row, d) for d in (T2D_CV, OBESITY_CV, T2D_KIDNEY)]
