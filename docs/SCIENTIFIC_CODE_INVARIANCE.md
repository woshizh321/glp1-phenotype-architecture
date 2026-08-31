# Scientific code invariance

| Authoritative source | Public path | Publication edit | Scientific logic changed |
|---|---|---|---|
| `src/preflight00r2/classifier.py` | `code/03_domain_mapping/classifier_base.py` | publication naming | NO |
| `src/preflight00r2r/classifier.py` | `code/03_domain_mapping/classifier.py` | relative module path; publication naming | NO |
| `src/phase1_protocol_lock/protocol_lock.py` | `code/03_domain_mapping/domain_contract.py` | relative module path; publication naming | NO |
| `src/phase1_protocol_lock/phase1_survey_engine.R` | `code/08_survey_reporting/survey_engine.R` | comment cleanup | NO |
| `src/phase1_execution/run_phase1_restart.py` derivation | `code/02_nhanes_cohort/build_analysis_input.py` | extracted data-derivation entry point; private provenance checks omitted | NO |
| `src/phase1_execution/execute_weighted_phase1.R` | `code/03_domain_mapping/execute_weighted_domain_mapping.R` | repository-relative dependencies | NO |
| `src/portfolio_augmentation_pa05/build_pa05_inputs.py` | `code/04_trial_specific_mapping/build_trial_mapping.py` | configured paths and local-output naming | NO |
| `src/portfolio_augmentation_pa05/execute_pa05_weighted.R` | `code/05_portfolio_architecture/execute_portfolio_analysis.R` | repository-relative dependencies | NO |
| `src/clinician_first_weighted_sensitivity/execute_weighted.R` | `code/06_signature_analysis/execute_signature_sensitivity.R` | repository-relative dependencies | NO |
| final V5R3 asset-generation logic | `code/09_figures_tables/generate_publication_assets.py` | aggregate-input public renderer | NO; rendering only |

The non-weighted NHANES derivation was re-executed locally against the authoritative 15,560-row analytical input: all 72 comparable derived/source columns matched exactly. Weighted analyses were not rerun during curation.

`SCIENTIFIC_LOGIC_CHANGED = NO`
