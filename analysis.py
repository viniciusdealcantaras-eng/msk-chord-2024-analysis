import pandas as pd, numpy as np, json
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.utils import median_survival_times

df = pd.read_csv("out/merged_raw.csv")
R = {}

# ---- clean ----
d = df[df["delay_months"].notna() & df["OS_MONTHS"].notna()].copy()
n_anom = (d["delay_months"] < 0).sum()
d = d[d["delay_months"] >= 0].copy()
d["T_dx"] = d["delay_months"] + d["OS_MONTHS"]     # OS from diagnosis
d["entry"] = d["delay_months"]                      # left-truncation entry time
d["T_seq"] = d["OS_MONTHS"]                         # OS from sequencing (as released)
R["n_total"] = int(len(df))
R["n_analysis"] = int(len(d))
R["n_anomalous_neg_delay_dropped"] = int(n_anom)
R["n_deaths"] = int(d["deceased"].sum())

# ---- 1. delay distribution ----
dl = d["delay_months"]
R["delay_median"] = float(dl.median())
R["delay_iqr"] = [float(dl.quantile(.25)), float(dl.quantile(.75))]
R["delay_p90"] = float(dl.quantile(.90))
R["delay_mean"] = float(dl.mean())
R["delay_gt12_pct"] = float((dl>12).mean()*100)
R["delay_gt24_pct"] = float((dl>24).mean()*100)
R["delay_gt60_pct"] = float((dl>60).mean()*100)
# by stage
R["delay_by_stage"] = {}
for s in ["Stage 1-3","Stage 4"]:
    sub = d.loc[d["STAGE_HIGHEST_RECORDED"]==s,"delay_months"]
    R["delay_by_stage"][s] = {"n":int(len(sub)),"median":float(sub.median()),
                              "q25":float(sub.quantile(.25)),"q75":float(sub.quantile(.75))}
# by top cancer types
R["delay_by_type"] = {}
for c in d["cancer_type"].value_counts().head(5).index:
    sub = d.loc[d["cancer_type"]==c,"delay_months"]
    R["delay_by_type"][c] = {"n":int(len(sub)),"median":float(sub.median()),
                             "q25":float(sub.quantile(.25)),"q75":float(sub.quantile(.75))}

# ---- 2. immortal person-time ----
immortal = d["delay_months"].sum()
total_dx = d["T_dx"].sum()
R["immortal_persontime_pct"] = float(immortal/total_dx*100)
R["immortal_person_years"] = float(immortal/12)

# ---- 3. core KM: from-diagnosis naive vs left-truncated, + from-seq ----
def km_summary(dur, ev, entry=None):
    kmf = KaplanMeierFitter()
    kmf.fit(dur, ev, entry=entry)
    med = kmf.median_survival_time_
    ci = median_survival_times(kmf.confidence_interval_)
    lo = ci.iloc[0,0]; hi = ci.iloc[0,1]
    def S(t):
        try: return float(kmf.survival_function_at_times(t).iloc[0])
        except Exception: return float("nan")
    return {"median":float(med),"median_ci":[float(lo),float(hi)],
            "S24":S(24),"S60":S(60),"S120":S(120)}, kmf

naive, kmf_naive = km_summary(d["T_dx"], d["deceased"])                 # ignores truncation
corr,  kmf_corr  = km_summary(d["T_dx"], d["deceased"], entry=d["entry"])# left-truncated
seq,   kmf_seq   = km_summary(d["T_seq"], d["deceased"])                # as released (from seq)
R["km_from_dx_naive"] = naive
R["km_from_dx_corrected"] = corr
R["km_from_seq"] = seq
R["median_overest_months"] = float(naive["median"] - corr["median"])
R["S60_overest_abs"] = float(naive["S60"] - corr["S60"])

# ---- 4. spurious advantage of 'late-sequenced' (immortal time by design) ----
med_delay = d["delay_months"].median()
d["seq_timing"] = np.where(d["delay_months"]>=med_delay,"late (>=median delay)","early (<median delay)")
R["spurious_late_vs_early"] = {}
for grp in ["early (<median delay)","late (>=median delay)"]:
    g = d[d["seq_timing"]==grp]
    s_naive,_ = km_summary(g["T_dx"], g["deceased"])
    s_corr,_  = km_summary(g["T_dx"], g["deceased"], entry=g["entry"])
    s_seq,_   = km_summary(g["T_seq"], g["deceased"])
    R["spurious_late_vs_early"][grp] = {"n":int(len(g)),
        "median_from_dx_naive":s_naive["median"],
        "median_from_dx_corrected":s_corr["median"],
        "median_from_seq":s_seq["median"]}

# ---- 5. Cox: stage effect, naive vs left-truncated ----
cx = d[d["STAGE_HIGHEST_RECORDED"].isin(["Stage 1-3","Stage 4"])].copy()
cx["stage4"] = (cx["STAGE_HIGHEST_RECORDED"]=="Stage 4").astype(int)
cx["age_z"] = (cx["CURRENT_AGE_DEID"]-cx["CURRENT_AGE_DEID"].mean())/cx["CURRENT_AGE_DEID"].std()
cx = cx.dropna(subset=["age_z"])

def cox_hr(dframe, entry_col=None):
    cph = CoxPHFitter()
    cols = ["T_dx","deceased","stage4","age_z"]
    if entry_col: cols = cols+[entry_col]
    cph.fit(dframe[cols], duration_col="T_dx", event_col="deceased", entry_col=entry_col)
    out={}
    for v in ["stage4","age_z"]:
        out[v]={"HR":float(np.exp(cph.params_[v])),
                "ci":[float(np.exp(cph.confidence_intervals_.loc[v].iloc[0])),
                      float(np.exp(cph.confidence_intervals_.loc[v].iloc[1]))],
                "p":float(cph.summary.loc[v,"p"])}
    return out
R["cox_from_dx_naive"] = cox_hr(cx)
R["cox_from_dx_corrected"] = cox_hr(cx, entry_col="entry")
# reference: cox from sequencing (as released)
def cox_hr_seq(dframe):
    cph=CoxPHFitter(); cph.fit(dframe[["T_seq","deceased","stage4","age_z"]],
                               duration_col="T_seq",event_col="deceased")
    return {v:{"HR":float(np.exp(cph.params_[v])),"p":float(cph.summary.loc[v,"p"])} for v in ["stage4","age_z"]}
R["cox_from_seq"] = cox_hr_seq(cx)

json.dump(R, open("out/results.json","w"), indent=2)
print(json.dumps(R, indent=2))
