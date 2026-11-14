# -*- coding: utf-8 -*-
# Landmark analysis at 6 and 12 months, conditional Kendall's tau for
# quasi-independence (subsample test), Schoenfeld test and stage counts.
import pandas as pd, numpy as np, json
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import proportional_hazard_test

d = pd.read_csv("out/merged_raw.csv")
d = d[d["delay_months"].notna() & d["OS_MONTHS"].notna()].copy()
d = d[d["delay_months"] >= 0].copy()
d["T_dx"] = d["delay_months"] + d["OS_MONTHS"]     # OS from diagnosis
d["entry"] = d["delay_months"]                      # left-truncation entry
R = {}

def km_median(dur, ev, entry=None):
    k = KaplanMeierFitter(); k.fit(dur, ev, entry=entry)
    return float(k.median_survival_time_)

# reference marginal medians
R["ref_naive_median"] = round(km_median(d["T_dx"], d["deceased"]), 1)
R["ref_ltrunc_median"] = round(km_median(d["T_dx"], d["deceased"], entry=d["entry"]), 1)
R["ref_fromseq_median"] = round(km_median(d["OS_MONTHS"], d["deceased"]), 1)

# ---------- landmark analysis at 6 and 12 months from diagnosis ----------
# include patients sequenced by t0 (entry<=t0) AND alive at t0 (T_dx>t0);
# estimate survival measured from diagnosis, conditional on reaching t0 in-cohort.
R["landmark"] = {}
for t0 in [6, 12]:
    g = d[(d["entry"] <= t0) & (d["T_dx"] > t0)].copy()
    # KM measured from diagnosis among this landmark cohort (no immortal time: all truly at risk from t0)
    kmf = KaplanMeierFitter(); kmf.fit(g["T_dx"], g["deceased"], entry=np.full(len(g), t0))
    med_from_dx = float(kmf.median_survival_time_)
    # stage HR in landmark cohort (Cox, delayed entry at t0), stratified by cancer type
    cx = g[g["STAGE_HIGHEST_RECORDED"].isin(["Stage 1-3", "Stage 4"])].copy()
    cx["stage4"] = (cx["STAGE_HIGHEST_RECORDED"] == "Stage 4").astype(int)
    cx["t0"] = float(t0)
    cph = CoxPHFitter()
    cph.fit(cx[["T_dx", "deceased", "stage4", "cancer_type", "t0"]],
            duration_col="T_dx", event_col="deceased", entry_col="t0", strata=["cancer_type"])
    hr = float(np.exp(cph.params_["stage4"]))
    lo = float(np.exp(cph.confidence_intervals_.loc["stage4"].iloc[0]))
    hi = float(np.exp(cph.confidence_intervals_.loc["stage4"].iloc[1]))
    R["landmark"][f"t0_{t0}m"] = {"n": int(len(g)), "deaths": int(g["deceased"].sum()),
                                  "median_OS_from_dx": round(med_from_dx, 1),
                                  "stageIV_HR": round(hr, 2), "ci": [round(lo, 2), round(hi, 2)]}

# ---------- conditional Kendall's tau for quasi-independence (Martin-Betensky) ----------
# subsample for tractability (O(n^2) over ~25k is 6e8 pairs); fixed seed.
rng = np.random.RandomState(20260728)
sub = d.sample(n=5000, random_state=rng) if len(d) > 5000 else d
L = sub["entry"].values.astype(float)
X = sub["T_dx"].values.astype(float)
del_ = sub["deceased"].values.astype(int)
n = len(sub)
conc = disc = comp = 0
for i in range(n):
    Li, Xi, di = L[i], X[i], del_[i]
    Lj, Xj, dj = L[i+1:], X[i+1:], del_[i+1:]
    # comparable (orderable) pair: max entry < min exit  -> both at risk together
    orderable = np.maximum(Li, Lj) < np.minimum(Xi, Xj)
    # need the smaller X to be an observed event (uncensored) to order the pair
    xi_smaller = Xi < Xj
    smaller_is_event = np.where(xi_smaller, di == 1, dj == 1)
    use = orderable & smaller_is_event
    if use.any():
        s = np.sign((Li - Lj[use]) * (Xi - Xj[use]))
        conc += int((s > 0).sum()); disc += int((s < 0).sum()); comp += int(use.sum())
tau_c = (conc - disc) / comp if comp else float("nan")
# approximate SE via normal approximation on the concordant fraction (subsample)
p_conc = conc / comp
se = np.sqrt(4 * p_conc * (1 - p_conc) / comp)  # rough
z = tau_c / se if se > 0 else float("nan")
from math import erfc, sqrt
pval = erfc(abs(z) / sqrt(2))
R["cond_kendall_tau"] = {"n_subsample": int(n), "comparable_pairs": int(comp),
                         "tau_c": round(tau_c, 3), "z": round(float(z), 1),
                         "p": float(pval), "seed": 20260728}

# ---------- Schoenfeld PH test for stage IV (stratified by cancer type) ----------
cx = d[d["STAGE_HIGHEST_RECORDED"].isin(["Stage 1-3", "Stage 4"])].copy()
cx["stage4"] = (cx["STAGE_HIGHEST_RECORDED"] == "Stage 4").astype(int)
# PH residuals not implemented for delayed-entry models in lifelines; test on the
# from-sequencing clock (no entry), stratified by cancer type, same stage effect.
cph = CoxPHFitter()
cph.fit(cx[["OS_MONTHS", "deceased", "stage4", "cancer_type"]],
        duration_col="OS_MONTHS", event_col="deceased", strata=["cancer_type"])
ph = proportional_hazard_test(cph, cx[["OS_MONTHS", "deceased", "stage4", "cancer_type"]], time_transform="rank")
R["ph_test_stage4"] = {"test_stat": round(float(ph.summary["test_statistic"].iloc[0]), 2),
                       "p": float(ph.summary["p"].iloc[0])}

# ---------- stage counts in the analysis set ----------
vc = d["STAGE_HIGHEST_RECORDED"].value_counts(dropna=False)
R["stage_counts"] = {str(k): int(v) for k, v in vc.items()}
n13 = int((d["STAGE_HIGHEST_RECORDED"] == "Stage 1-3").sum())
n4 = int((d["STAGE_HIGHEST_RECORDED"] == "Stage 4").sum())
R["stage_pct_of_total"] = {"Stage 1-3": round(100*n13/len(d), 1), "Stage 4": round(100*n4/len(d), 1),
                           "denominator": int(len(d))}

json.dump(R, open("out/results4_landmark_tau.json", "w"), indent=2)
print(json.dumps(R, indent=2))
