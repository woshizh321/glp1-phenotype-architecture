# GLP-1 phenotype architecture

Reproducibility code and public metadata for mapping clinical phenotype representation across GLP-1–based cardiometabolic outcome trials using NHANES and ClinicalTrials.gov/AACT.

The study asks: (1) Which clinical phenotypes are represented across GLP-1 outcome programs? (2) Which T2D-CV phenotypes are uniquely or repeatedly anchored? (3) How do cardiovascular and kidney representations differ within the same adults?

Quantitative analysis uses the nationally representative US NHANES 2017–March 2020 pre-pandemic release. Trial metadata use the AACT snapshot dated 2026-05-01. Obtain these official public sources and configure local paths; no source or row-level data are redistributed.

## Scientific boundary

Representation is observable phenotype mapping—not exact trial eligibility, comparative trial quality, treatment recommendation, or treatment-effect estimation. STEP-HFpEF, STEP-HFpEF DM, and SUMMIT form a qualitative HF expansion layer and do not contribute to quantitative NHANES estimates.

Observable-signature equivalence refers to the NHANES projection used in this study and does not imply identical source trial eligibility criteria. Ten T2D-CV trials map to seven observable Direct signatures; LEADER, SUSTAIN-6, PIONEER 6, and REWIND share one projected signature.

## Run

```bash
cp config/paths.example.yaml config/paths.yaml
python run_pipeline.py --config config/paths.yaml --validate-inputs
python run_pipeline.py --config config/paths.yaml
```

Python/R versions and consequential packages are in `environment/`. Figure/table mappings are in `docs/FIGURE_TABLE_REPRODUCTION_MAP.md`.

## Citation and license

Please cite the associated Cardiovascular Diabetology manuscript when available; no DOI has been assigned. No code license has yet been assigned, so reuse permission is not implied until the repository owner selects one.
