#!/usr/bin/env python3
import argparse, numpy as np, pandas as pd
from pathlib import Path
NUMERIC=["internal_unweighted_numerator","weighted_proportion","ci95_lower","ci95_upper","design_se","design_effect","effective_n"]
p=argparse.ArgumentParser(); p.add_argument("input_csv",type=Path); p.add_argument("output_csv",type=Path); a=p.parse_args(); d=pd.read_csv(a.input_csv)
col="publication_display_status" if "publication_display_status" in d else "reliability_status"; hidden=~d[col].eq("PRESENT")
for c in NUMERIC:
    if c in d: d.loc[hidden,c]=np.nan
d["publication_reporting_status"]=np.where(hidden,"NR","PRESENT"); a.output_csv.parent.mkdir(parents=True,exist_ok=True); d.to_csv(a.output_csv,index=False)
