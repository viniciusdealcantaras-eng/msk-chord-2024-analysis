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

# ---------- Cox adjusted for cancer type (drop age) ----------
cx = d[d["STAGE_HIGHEST_RECORDED"].isin(["Stage 1-3","Stage 4"])].copy()
cx["stage4"] = (cx["STAGE_HIGHEST_RECORDED"]=="Stage 4").astype(int)
dummies = pd.get_dummies(cx["cancer_type"], prefix="ct", drop_first=True).astype(int)
base = pd.concat([cx[["T_dx","OS_MONTHS","deceased","entry","stage4"]].reset_index(drop=True),
                  dummies.reset_index(drop=True)], axis=1)

def cox_stage(df, dur, entry=None):
    keep = [dur,"deceased","stage4"]+[c for c in df.columns if c.startswith("ct_")]
    if entry: keep = keep+[entry]
    use = df[keep].copy()
    cph=CoxPHFitter()
    cph.fit(use, duration_col=dur, event_col="deceased", entry_col=(entry if entry else None))
    return {"HR":round(float(np.exp(cph.params_["stage4"])),2),
            "lo":round(float(np.exp(cph.confidence_intervals_.loc["stage4"].iloc[0])),2),
            "hi":round(float(np.exp(cph.confidence_intervals_.loc["stage4"].iloc[1])),2),
            "p":float(cph.summary.loc["stage4","p"])}

# need entry col present in df for left-truncated fit
base_lt = base.copy()
R2["cox_stage_ct"] = {
  "naive":   cox_stage(base, "T_dx"),
  "corrected": cox_stage(base_lt.assign(), "T_dx", entry="entry"),
  "from_seq": cox_stage(base.rename(columns={}), "OS_MONTHS"),
}

json.dump(R2, open("out/results2.json","w"), indent=2)
print(json.dumps(R2, indent=2))
