#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
p=argparse.ArgumentParser(); p.add_argument("input_csv",type=Path); p.add_argument("signature_map",type=Path); p.add_argument("output_csv",type=Path); a=p.parse_args()
d=pd.read_csv(a.input_csv,low_memory=False); m=pd.read_csv(a.signature_map); groups=m.groupby("signature_id").trial.apply(list).to_dict()
for sid,trials in groups.items(): d[f"signature__{sid}"]=d[[f"trial__{t}" for t in trials]].eq("DIRECT").any(axis=1)
d["cfwse_n_unique_direct_signatures"]=d[[f"signature__{s}" for s in groups]].sum(axis=1); d["cfwse_n_direct_minus_exscel"]=d.n_direct-d["trial__EXSCEL"].eq("DIRECT").astype(int)
a.output_csv.parent.mkdir(parents=True,exist_ok=True); d.to_csv(a.output_csv,index=False)
