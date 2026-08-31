#!/usr/bin/env python3
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
def save(fig,out,stem):
 out.mkdir(parents=True,exist_ok=True)
 for ext in ["png","pdf","svg"]: fig.savefig(out/f"{stem}.{ext}",dpi=300 if ext=="png" else None,bbox_inches="tight")
 plt.close(fig)
def bar(path,out,stem):
 d=pd.read_csv(path); d=d[d.weighted_proportion.notna()]; fig,ax=plt.subplots(figsize=(8,max(4,.34*len(d)))); ax.barh(d.category.astype(str),100*d.weighted_proportion,color="#2f6f8f"); ax.set_xlabel("Weighted proportion (%)"); save(fig,out,stem)
p=argparse.ArgumentParser(); p.add_argument("--results",type=Path,required=True); p.add_argument("--metadata",type=Path,required=True); p.add_argument("--output",type=Path,required=True); a=p.parse_args()
d=pd.read_csv(a.metadata/"trial_portfolio_public.csv").dropna(subset=["publication_year"]); fig,ax=plt.subplots(figsize=(10,6))
for i,r in d.sort_values(["portfolio_status","publication_year","trial"]).reset_index().iterrows():
 ax.scatter(r.publication_year,i,marker="D" if r.portfolio_status=="QUALITATIVE_HF_LAYER" else "o"); ax.text(r.publication_year+.08,i,r.trial,va="center",fontsize=8)
ax.set_yticks([]); ax.set_xlabel("Primary results year"); ax.set_title("Outcome-evidence portfolio\nDiamonds: qualitative HF layer; excluded from quantitative NHANES mapping"); save(fig,a.output,"figure1_trial_portfolio")
for rel,stem in [("estimates/PHASE1_EVIDENCE_STATUS_ESTIMATES.csv","figure2_domain_representation"),("MODULE_B_DEPTH_RESULTS.csv","figure4_portfolio_depth"),("MODULE_D_SELECTED_ESTIMANDS.csv","figure5_cv_kidney")]:
 q=a.results/rel
 if q.exists(): bar(q,a.output,stem)
