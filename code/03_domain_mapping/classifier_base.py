#!/usr/bin/env python3
"""Deterministic participant evidence-status classifier.

This module implements diagnostic phenotype classification only. It does not
estimate trial eligibility, treatment effects, treatment recommendations, or
survey-weighted population quantities.
"""

from __future__ import annotations

import math
from typing import Any


CLASSIFIER_VERSION = "PUBLIC_DOMAIN_CLASSIFIER_1.0.0"

DIRECT = "DIRECT_DOMAIN_REPRESENTATION"
CONDITIONAL = "CONDITIONAL_DOMAIN_REPRESENTATION"
NO_DIRECT = "NO_DIRECT_COMPLETED_DOMAIN_REPRESENTATION"
UNCLASSIFIED = "UNCLASSIFIED_DATA_LIMITATION"

T2D_CV = "T2D_CARDIOVASCULAR_OUTCOME"
OBESITY_CV = "OBESITY_CARDIOVASCULAR_OUTCOME"
T2D_KIDNEY = "T2D_KIDNEY_OUTCOME"


def _number(value: Any) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _integer_code(value: Any, valid: set[int]) -> int | None:
    number = _number(value)
    if number is None or not number.is_integer():
        return None
    code = int(number)
    return code if code in valid else None


def ascvd_proxy(c: Any, e: Any, f: Any) -> tuple[bool | None, str]:
    """Evaluate CHD/MI/stroke as a three-valued OR.

    Any observed yes is sufficient for a positive proxy. A negative proxy
    requires all three observed no responses. Every other combination is
    unresolved and cannot be silently treated as disease absence.
    """

    codes = [_integer_code(v, {1, 2}) for v in (c, e, f)]
    if 1 in codes:
        return True, "ASCVD_OR_POSITIVE"
    if codes == [2, 2, 2]:
        return False, "ASCVD_OR_ALL_NEGATIVE"
    return None, "ASCVD_OR_UNRESOLVED"


def ckd_epi_2021(age: Any, sex: Any, serum_creatinine: Any) -> float | None:
    """Return 2021 CKD-EPI race-free eGFRcr (mL/min/1.73 m2)."""

    age_n = _number(age)
    scr = _number(serum_creatinine)
    sex_code = _integer_code(sex, {1, 2})
    if age_n is None or not 0 <= age_n <= 80 or scr is None or not 0.25 <= scr <= 14.97 or sex_code is None:
        return None
    female = sex_code == 2
    kappa = 0.7 if female else 0.9
    alpha = -0.241 if female else -0.302
    ratio = scr / kappa
    value = (
        142.0
        * min(ratio, 1.0) ** alpha
        * max(ratio, 1.0) ** -1.200
        * 0.9938 ** age_n
        * (1.012 if female else 1.0)
    )
    return value if math.isfinite(value) and value > 0 else None


def kidney_cell(egfr: Any, uacr: Any) -> tuple[str | None, str | None, str | None]:
    """Assign one mutually exclusive and exhaustive eGFR x UACR cell."""

    egfr_n = _number(egfr)
    uacr_n = _number(uacr)
    if egfr_n is None or egfr_n <= 0:
        return None, None, "INVALID_OR_MISSING_EGFR"
    if uacr_n is None or not 0.27 <= uacr_n <= 11676.92:
        return None, None, "INVALID_OR_MISSING_UACR"

    if egfr_n < 15:
        g = "G0_LT15"
    elif egfr_n < 25:
        g = "G1_15_LT25"
    elif egfr_n < 50:
        g = "G2_25_LT50"
    elif egfr_n < 60:
        g = "G3_50_LT60"
    elif egfr_n <= 75:
        g = "G4_60_LE75"
    else:
        g = "G5_GT75"

    if uacr_n < 30:
        a = "A0_LT30"
    elif uacr_n <= 100:
        a = "A1_30_LE100"
    elif uacr_n <= 300:
        a = "A2_GT100_LE300"
    elif uacr_n < 5000:
        a = "A3_GT300_LT5000"
    else:
        a = "A4_GE5000"
    return f"{g}__{a}", f"{g}|{a}", None


def classify_kidney_values(diabetes_code: Any, egfr: Any, uacr: Any) -> dict[str, Any]:
    """Classify the FLOW-like gateway from already derived eGFR and UACR."""

    diabetes = _integer_code(diabetes_code, {1, 2})
    if diabetes is None:
        return _result(T2D_KIDNEY, UNCLASSIFIED, "KIDNEY_DIABETES_UNRESOLVED", True)
    if diabetes == 2:
        return _result(T2D_KIDNEY, NO_DIRECT, "KIDNEY_NO_DIAGNOSED_DIABETES", False)

    cell_id, branch_pair, error = kidney_cell(egfr, uacr)
    if error:
        return _result(T2D_KIDNEY, UNCLASSIFIED, error, True)
    assert cell_id is not None and branch_pair is not None
    g, a = branch_pair.split("|")

    # Exact FLOW branches: eGFR >=25 to <50 with UACR >100 to <5000,
    # OR eGFR >=50 to <=75 with UACR >300 to <5000.
    if g == "G2_25_LT50" and a in {"A2_GT100_LE300", "A3_GT300_LT5000"}:
        return _result(T2D_KIDNEY, DIRECT, "FLOW_BRANCH_EGFR25_LT50_UACR_GT100_LT5000", False, cell_id)
    if g in {"G3_50_LT60", "G4_60_LE75"} and a == "A3_GT300_LT5000":
        return _result(T2D_KIDNEY, DIRECT, "FLOW_BRANCH_EGFR50_LE75_UACR_GT300_LT5000", False, cell_id)

    egfr_n = float(egfr)
    uacr_n = float(uacr)
    if egfr_n < 15:
        return _result(T2D_KIDNEY, NO_DIRECT, "KIDNEY_FAILURE_OUTSIDE_FLOW_GATEWAY", False, cell_id)
    if egfr_n < 60 or uacr_n >= 30:
        return _result(T2D_KIDNEY, CONDITIONAL, "BROADER_CKD_OUTSIDE_FLOW_GATEWAY", False, cell_id)
    return _result(T2D_KIDNEY, NO_DIRECT, "NO_OBSERVABLE_CKD_MARKER", False, cell_id)


def classify_t2d_cv(row: dict[str, Any]) -> dict[str, Any]:
    diabetes = _integer_code(row.get("DIQ010"), {1, 2})
    if diabetes is None:
        return _result(T2D_CV, UNCLASSIFIED, "T2DCV_DIABETES_UNRESOLVED", True)
    if diabetes == 2:
        return _result(T2D_CV, NO_DIRECT, "T2DCV_NO_DIAGNOSED_DIABETES", False)
    ascvd, branch = ascvd_proxy(row.get("MCQ160C"), row.get("MCQ160E"), row.get("MCQ160F"))
    if ascvd is True:
        return _result(T2D_CV, DIRECT, "T2DCV_DIAGNOSED_DIABETES_AND_ASCVD_PROXY", False)
    if ascvd is False:
        return _result(T2D_CV, CONDITIONAL, "T2DCV_DIAGNOSED_DIABETES_WITHOUT_ASCVD_PROXY", False)
    return _result(T2D_CV, UNCLASSIFIED, f"T2DCV_{branch}", True)


def classify_obesity_cv(row: dict[str, Any]) -> dict[str, Any]:
    age = _number(row.get("RIDAGEYR"))
    bmi = _number(row.get("BMXBMI"))
    if age is None or not 20 <= age <= 80:
        return _result(OBESITY_CV, UNCLASSIFIED, "OBESCV_AGE_INVALID", True)
    if bmi is None or not 11.9 <= bmi <= 92.3:
        return _result(OBESITY_CV, UNCLASSIFIED, "OBESCV_BMI_MISSING_OR_OUTSIDE_RELEASED_RANGE", True)
    if bmi < 25:
        return _result(OBESITY_CV, NO_DIRECT, "OBESCV_BELOW_DOMAIN_BMI25", False)

    ascvd, branch = ascvd_proxy(row.get("MCQ160C"), row.get("MCQ160E"), row.get("MCQ160F"))
    if ascvd is False:
        return _result(OBESITY_CV, NO_DIRECT, "OBESCV_NO_ASCVD_PROXY", False)
    if ascvd is None:
        return _result(OBESITY_CV, UNCLASSIFIED, f"OBESCV_{branch}", True)

    diabetes = _integer_code(row.get("DIQ010"), {1, 2})
    if diabetes is None:
        return _result(OBESITY_CV, UNCLASSIFIED, "OBESCV_DIABETES_UNRESOLVED", True)
    if age >= 45 and bmi >= 27 and diabetes == 2:
        return _result(OBESITY_CV, DIRECT, "SELECT_PROXY_AGE45_BMI27_ASCVD_NO_DIAGNOSED_DIABETES", False)
    if diabetes == 1:
        return _result(OBESITY_CV, CONDITIONAL, "OBESCV_ASCVD_WITH_DIAGNOSED_DIABETES", False)
    if bmi < 27:
        return _result(OBESITY_CV, CONDITIONAL, "OBESCV_ASCVD_BMI25_LT27", False)
    return _result(OBESITY_CV, CONDITIONAL, "OBESCV_ASCVD_AGE20_LT45", False)


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
    """Return exactly one status for each of the three quantitative domains."""

    return [classify_t2d_cv(row), classify_obesity_cv(row), classify_kidney(row)]


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
    }
