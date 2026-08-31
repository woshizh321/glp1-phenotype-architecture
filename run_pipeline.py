#!/usr/bin/env python3
import argparse, os, subprocess, sys
from pathlib import Path
import pandas as pd
REQ={"SEQN","RIDAGEYR","RIAGENDR","RIDRETH3","WTMECPRP","WTINTPRP","SDMVPSU","SDMVSTRA","DIQ010","BMXBMI","MCQ160C","MCQ160D","MCQ160E","MCQ160F","LBXSCR","URDACT","LBXGH"}
p=argparse.ArgumentParser(); p.add_argument("--config",type=Path,required=True); p.add_argument("--validate-inputs",action="store_true"); p.add_argument("--dry-run",action="store_true"); a=p.parse_args(); root=Path(__file__).resolve().parent; cfg={line.split(":",1)[0].strip():line.split(":",1)[1].strip() for line in a.config.read_text().splitlines() if line.strip() and not line.lstrip().startswith("#")}; src=Path(cfg["nhanes_source_csv"]); out=Path(cfg["output_root"]); analysis=out/"phase1/nhanes_analysis_input.csv"
if not src.exists(): raise SystemExit(f"Missing input: {src}")
miss=sorted(REQ-set(pd.read_csv(src,nrows=0).columns))
if miss: raise SystemExit(f"Missing variables: {miss}")
out.mkdir(parents=True,exist_ok=True); print("Input validation passed")
if not a.validate_inputs:
 env=os.environ.copy(); env.update(PUBLIC_REPO_ROOT=str(root),PHASE1_INPUT_CSV=str(analysis)); cmds=[[sys.executable,root/"code/02_nhanes_cohort/build_analysis_input.py",src,analysis],["Rscript",root/"code/03_domain_mapping/execute_weighted_domain_mapping.R",analysis,out/"phase1"],[sys.executable,root/"code/04_trial_specific_mapping/build_trial_mapping.py",out/"portfolio"],["Rscript",root/"code/05_portfolio_architecture/execute_portfolio_analysis.R",out/"portfolio/intermediate/PA05_SURVEY_INPUT_INTERNAL.csv",out/"portfolio"],[sys.executable,root/"code/06_signature_analysis/prepare_signature_inputs.py",out/"portfolio/intermediate/PA05_SURVEY_INPUT_INTERNAL.csv",root/"metadata/signature_map_public.csv",out/"signature/signature_input.csv"],["Rscript",root/"code/06_signature_analysis/execute_signature_sensitivity.R",out/"signature/signature_input.csv",out/"signature"],[sys.executable,root/"code/09_figures_tables/generate_publication_assets.py","--results",out/"portfolio","--metadata",root/"metadata","--output",out/"figures"]]
 for cmd in cmds:
  print("+"," ".join(map(str,cmd)))
  if not a.dry_run: subprocess.run(list(map(str,cmd)),check=True,env=env)
