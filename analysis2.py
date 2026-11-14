import pandas as pd, numpy as np, json
from lifelines import KaplanMeierFitter, CoxPHFitter

d = pd.read_csv("out/merged_raw.csv")
R2 = {}

# ---------- multiple-diagnosis frequency ----------
dx = pd.read_csv("data/data_timeline_diagnosis.txt", sep="\t", dtype=str)
dx["START_DATE"] = pd.to_numeric(dx["START_DATE"], errors="coerce")
ndx = dx.groupby("PATIENT_ID").size()
# gap between earliest and latest diagnosis (months) for those with >1
span = dx.groupby("PATIENT_ID")["START_DATE"].agg(lambda s: (s.max()-s.min()))
multi_ids = set(ndx[ndx>1].index)
R2["pct_multi_dx"] = round(100*len(multi_ids)/ndx.shape[0],1)
R2["n_multi_dx"] = int(len(multi_ids))
R2["n_with_dx"] = int(ndx.shape[0])

# ---------- analysis frame ----------
d = d[d["delay_months"].notna() & d["OS_MONTHS"].notna()].copy()
d = d[d["delay_months"]>=0].copy()
d["T_dx"]=d["delay_months"]+d["OS_MONTHS"]; d["entry"]=d["delay_months"]
d["single_dx"] = ~d["PATIENT_ID"].isin(multi_ids)

def med(dur,ev,entry=None):
    k=KaplanMeierFitter(); k.fit(dur,ev,entry=entry); return float(k.median_survival_time_)

# ---------- sensitivity: single-diagnosis only ----------
s = d[d["single_dx"]]
R2["sens_single_dx"] = {
  "n": int(len(s)),
  "median_naive": round(med(s["T_dx"],s["deceased"]),1),
  "median_corrected": round(med(s["T_dx"],s["deceased"],entry=s["entry"]),1),
}
R2["full_ref"] = {"median_naive": round(med(d["T_dx"],d["deceased"]),1),
                  "median_corrected": round(med(d["T_dx"],d["deceased"],entry=d["entry"]),1)}

json.dump(R2, open("out/results2.json","w"), indent=2)
print(json.dumps(R2, indent=2))
