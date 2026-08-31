#!/usr/bin/env python3
"""Minimal publication boundary repair over the base classifier.

Only the authorized SELECT/PAD status boundary and FLOW alignment qualifier are
changed. All other R2 classifier behavior is inherited from the base
source and remains diagnostic-only and unweighted.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
from typing import Any


R2_CLASSIFIER = Path(__file__).with_name("classifier_base.py")
_spec = importlib.util.spec_from_file_location("public_base_classifier", R2_CLASSIFIER)
if _spec is None or _spec.loader is None:
    raise ImportError(f"Cannot load base classifier: {R2_CLASSIFIER}")
_base = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_base)


CLASSIFIER_VERSION = "PUBLIC_DOMAIN_CLASSIFIER_1.0.0"

DIRECT = _base.DIRECT
CONDITIONAL = _base.CONDITIONAL
NO_DIRECT = _base.NO_DIRECT
UNCLASSIFIED = _base.UNCLASSIFIED
T2D_CV = _base.T2D_CV
OBESITY_CV = _base.OBESITY_CV
T2D_KIDNEY = _base.T2D_KIDNEY

ALIGNMENT_COMPATIBLE_PROXY = "COMPATIBLE_PROXY"
ALIGNMENT_EGFR_MISMATCH = "EGFR_EQUATION_MISMATCH_2009_VS_2021"

_number = _base._number
_integer_code = _base._integer_code
ascvd_proxy = _base.ascvd_proxy
ckd_epi_2021 = _base.ckd_epi_2021
kidney_cell = _base.kidney_cell
classify_t2d_cv = _base.classify_t2d_cv


def _versioned(result: dict[str, Any]) -> dict[str, Any]:
    result = dict(result)
    result["classifier_version"] = CLASSIFIER_VERSION
    result.setdefault("alignment_type", "")
    result.setdefault("alignment_reason", "")
    return result


def _result(
    domain: str,
    status: str,
    branch: str,
    missing: bool,
    kidney_cell_id: str | None = None,
) -> dict[str, Any]:
    return {
        "domain": domain,
        "evidence_status": status,
        "triggered_branch": branch,
        "missingness_flag": "YES" if missing else "NO",
        "kidney_cell_id": kidney_cell_id,
        "classifier_version": CLASSIFIER_VERSION,
        "alignment_type": "",
        "alignment_reason": "",
    }


def classify_obesity_cv(row: dict[str, Any]) -> dict[str, Any]:
    """Apply the authorized conservative PAD-unobservable boundary policy."""

    age = _number(row.get("RIDAGEYR"))
    bmi = _number(row.get("BMXBMI"))
    if age is None or not 20 <= age <= 80:
        return _result(OBESITY_CV, UNCLASSIFIED, "OBESCV_AGE_INVALID", True)
    if bmi is None or not 11.9 <= bmi <= 92.3:
        return _result(OBESITY_CV, UNCLASSIFIED, "OBESCV_BMI_MISSING_OR_OUTSIDE_RELEASED_RANGE", True)
    if bmi < 25:
        return _result(OBESITY_CV, NO_DIRECT, "OBESCV_BELOW_DOMAIN_BMI25", False)

    ascvd, branch = ascvd_proxy(row.get("MCQ160C"), row.get("MCQ160E"), row.get("MCQ160F"))
    if ascvd is None:
        return _result(OBESITY_CV, UNCLASSIFIED, f"OBESCV_{branch}", True)

    diabetes = _integer_code(row.get("DIQ010"), {1, 2})

    # PAD is not observed. If every other SELECT-defining observable condition
    # is met, an all-negative CHD/MI/stroke proxy cannot establish no CVD.
    if age >= 45 and bmi >= 27 and diabetes == 2 and ascvd is False:
        result = _result(OBESITY_CV, CONDITIONAL, "PAD_STATUS_UNOBSERVABLE", False)
        result["alignment_type"] = "SELECT_OBSERVABLE_CVD_PROXY_INCOMPLETE"
        result["alignment_reason"] = "SYMPTOMATIC_PAD_NOT_OBSERVED"
        return result

    # With the observable CVD proxy negative, unresolved diabetes cannot be
    # short-circuited: a valid completion of 2 invokes PAD uncertainty whereas
    # a completion of 1 clearly fails SELECT's no-diabetes gateway.
    if ascvd is False:
        if age >= 45 and bmi >= 27 and diabetes is None:
            return _result(
                OBESITY_CV,
                UNCLASSIFIED,
                "OBESCV_DIABETES_UNRESOLVED_WITH_PAD_UNOBSERVABLE",
                True,
            )
        return _result(OBESITY_CV, NO_DIRECT, "OBESCV_NO_ASCVD_PROXY", False)
    if diabetes is None:
        return _result(OBESITY_CV, UNCLASSIFIED, "OBESCV_DIABETES_UNRESOLVED", True)
    if age >= 45 and bmi >= 27 and diabetes == 2:
        return _result(OBESITY_CV, DIRECT, "SELECT_PROXY_AGE45_BMI27_ASCVD_NO_DIAGNOSED_DIABETES", False)
    if diabetes == 1:
        return _result(OBESITY_CV, CONDITIONAL, "OBESCV_ASCVD_WITH_DIAGNOSED_DIABETES", False)
    if bmi < 27:
        return _result(OBESITY_CV, CONDITIONAL, "OBESCV_ASCVD_BMI25_LT27", False)
    return _result(OBESITY_CV, CONDITIONAL, "OBESCV_ASCVD_AGE20_LT45", False)


def ckd_epi_2009(
    age: Any,
    sex: Any,
    serum_creatinine: Any,
    apply_black_coefficient: bool = False,
) -> float | None:
    """Return 2009 CKD-EPI creatinine eGFR (mL/min/1.73 m2).

    ``apply_black_coefficient`` applies the historical multiplicative 1.159
    term. It is used only for the prespecified equation-concordance sensitivity,
    never for the frozen 2021-equation primary classifier.
    """

    age_n = _number(age)
    scr = _number(serum_creatinine)
    sex_code = _integer_code(sex, {1, 2})
    if age_n is None or not 0 <= age_n <= 80 or scr is None or not 0.25 <= scr <= 14.97 or sex_code is None:
        return None
    female = sex_code == 2
    kappa = 0.7 if female else 0.9
    alpha = -0.329 if female else -0.411
    ratio = scr / kappa
    value = (
        141.0
        * min(ratio, 1.0) ** alpha
        * max(ratio, 1.0) ** -1.209
        * 0.993 ** age_n
        * (1.018 if female else 1.0)
        * (1.159 if apply_black_coefficient else 1.0)
    )
    return value if math.isfinite(value) and value > 0 else None


def classify_kidney_values(diabetes_code: Any, egfr: Any, uacr: Any) -> dict[str, Any]:
    result = _versioned(_base.classify_kidney_values(diabetes_code, egfr, uacr))
    if result["evidence_status"] == DIRECT:
        result["alignment_type"] = ALIGNMENT_COMPATIBLE_PROXY
        result["alignment_reason"] = ALIGNMENT_EGFR_MISMATCH
    return result


def classify_kidney(row: dict[str, Any]) -> dict[str, Any]:
    diabetes = _integer_code(row.get("DIQ010"), {1, 2})
    if diabetes is None:
        return _result(T2D_KIDNEY, UNCLASSIFIED, "KIDNEY_DIABETES_UNRESOLVED", True)
    if diabetes == 2:
        return _result(T2D_KIDNEY, NO_DIRECT, "KIDNEY_NO_DIAGNOSED_DIABETES", False)
    egfr = ckd_epi_2021(row.get("RIDAGEYR"), row.get("RIAGENDR"), row.get("LBXSCR"))
    result = classify_kidney_values(diabetes, egfr, row.get("URDACT"))
    result["derived_egfr"] = egfr
    return result


def classify_participant(row: dict[str, Any]) -> list[dict[str, Any]]:
    t2d = _versioned(classify_t2d_cv(row))
    return [t2d, classify_obesity_cv(row), classify_kidney(row)]
