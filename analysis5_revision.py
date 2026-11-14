"""Analyses added in revision R1: stage at index diagnosis, long intervals,
early risk sets, per-type medians and full-cohort conditional Kendall's tau.
Input: out/analytic_dataset.csv (make_analytic_dataset.py). Output: out/results5_revision.json"""
import sys, json
import numpy as np, pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter

f = pd.read_csv("out/analytic_dataset.csv")
d = f[f.included == 1].copy()
d["T"], d["L"], d["S"], d["E"] = d.os_months_from_diagnosis, d.dx_to_seq_months, d.os_months_from_sequencing, d.os_status
R = {}

# ---- stage at index diagnosis (earliest diagnosis event) ----
d["stage_dx"] = d.stage_at_index_diagnosis
ct = pd.crosstab(d.stage_dx, d.stage_highest_recorded, margins=True)
R["stage_crosstab"] = ct.to_dict()
up = ((d.stage_dx == "I-III") & (d.stage_highest_recorded == "IV")).sum()
R["upstaged_n"] = int(up); R["upstaged_pct_of_IV_highest"] = round(100 * up / (d.stage_highest_recorded == "IV").sum(), 1)
R["stage_dx_IV_n"] = int((d.stage_dx == "IV").sum())
R["delay_by_stage_dx"] = {s: [round(x, 1) for x in d.loc[d.stage_dx == s, "L"].quantile([.5, .25, .75])] for s in ["I-III", "IV"]}

def cox(df, stagecol, dur, entry=None):
    c = pd.DataFrame({"t": df[dur], "e": df.E, "x": (df[stagecol] == "IV").astype(int), "ct": df.cancer_type})
    cols = ["t", "e", "x", "ct"]
    if entry: c["en"] = df[entry]; cols.append("en")
    m = CoxPHFitter().fit(c[cols], "t", "e", strata=["ct"], entry_col=("en" if entry else None))
    ci = np.exp(m.confidence_intervals_.loc["x"].values)
    return [round(float(np.exp(m.params_["x"])), 2), round(float(ci[0]), 2), round(float(ci[1]), 2)]
for col in ["stage_highest_recorded", "stage_dx"]:
    R[f"HR_{col}"] = {"naive": cox(d, col, "T"), "lt": cox(d, col, "T", "L"), "seq": cox(d, col, "S")}

def km(x, lt):
    return KaplanMeierFitter().fit(x["T"], x["E"], entry=(x["L"] if lt else None))
def summ(x):
    a, b = km(x, False), km(x, True)
    return {"n": len(x), "med_naive": round(a.median_survival_time_, 1), "med_lt": round(b.median_survival_time_, 1),
            "S60_naive": round(100 * float(a.predict(60)), 1), "S60_lt": round(100 * float(b.predict(60)), 1)}

# ---- per cancer type medians ----
R["per_type"] = {k: summ(g) for k, g in d.groupby("cancer_type")}

# ---- sensitivity: exclude very long intervals ----
for cut in [120, 240]:
    s = d[d.L <= cut]
    R[f"sens_L_le_{cut}"] = {**summ(s), "excluded": int((d.L > cut).sum()),
                             "HR": {"naive": cox(s, "stage_highest_recorded", "T"), "lt": cox(s, "stage_highest_recorded", "T", "L")}}
R["L_gt240_by_type"] = d[d.L > 240].cancer_type.value_counts().to_dict()
R["L_zero_n"] = int((d.L == 0).sum())

# ---- early risk sets of the left-truncated estimator ----
R["risk_t0_by_type"] = {k: int((g.L <= 0).sum()) for k, g in d.groupby("cancer_type")}
# conditional on surviving 1 month after diagnosis (start where risk set is non-trivial)
def cond(x, t0):
    k = km(x, True); s0 = float(k.predict(t0))
    sf = k.survival_function_.iloc[:, 0]
    sc = sf[sf.index >= t0] / s0
    med = sc[sc <= 0.5].index.min()
    return round(float(med), 1)
R["lt_median_conditional_1mo"] = cond(d, 1.0)
R["lt_median_conditional_3mo"] = cond(d, 3.0)

# ---- landmark check (as described in manuscript) ----
for lm in [6, 12]:
    s = d[d["T"] >= lm].copy()
    s_seq = s[s.L <= lm]
    R[f"landmark_{lm}"] = {"alive_at_lm_any": int(len(s)), "alive_and_sequenced_by_lm": int(len(s_seq)),
                           "deaths_seqbylm": int(s_seq.E.sum())}

# ---- conditional Kendall's tau (Martin & Betensky 2005), full cohort ----
L, X, D = d.L.values, d["T"].values, d.E.values.astype(bool)
n = len(L); num = 0.0; den = 0
for i0 in range(0, n, 400):
    i = np.arange(i0, min(i0 + 400, n))[:, None]
    j = np.arange(n)[None, :]
    Li, Lj, Xi, Xj, Di, Dj = L[i], L[j], X[i], X[j], D[i], D[j]
    comp = (np.maximum(Li, Lj) <= np.minimum(Xi, Xj)) & (j > i)
    # orderable: the smaller observed time must be an event
    comp &= np.where(Xi < Xj, Di, np.where(Xj < Xi, Dj, Di & Dj))
    sgn = np.sign((Li - Lj) * (Xi - Xj))
    num += float((sgn * comp).sum()); den += int(comp.sum())
R["kendall_tau_full"] = round(num / den, 3); R["kendall_pairs"] = den

json.dump(R, open("out/results5_revision.json", "w"), indent=1, default=str)
print(json.dumps(R, indent=1, default=str))
