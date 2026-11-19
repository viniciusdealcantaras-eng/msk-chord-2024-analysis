"""Landmark with delayed entry, conditional Kendall's tau by stratum, and results by cancer type.
1. Landmark at 6 and 12 months with delayed entry for patients sequenced after the
   landmark (all patients alive at the landmark are kept), next to the original
   restricted version (only patients sequenced by the landmark).
2. Conditional Kendall's tau (Martin & Betensky 2005) in the full cohort, by cancer
   type and by cancer type x stage, with a U-statistic (Hoeffding projection)
   standard error computed on all patients, replacing the subsample test.
3. Descriptive table by cancer type: n, deaths, stage IV, interval, median follow-up
   from sequencing (reverse Kaplan-Meier), left-truncated stage IV HR per type.
4. Left-truncated median conditional on 1-month survival, by cancer type.
Input: out/analytic_dataset.csv. Output: out/results6_landmark_tau_by_type.json"""
import json
from math import erfc, sqrt
import numpy as np, pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter

f = pd.read_csv("out/analytic_dataset.csv")
d = f[f.included == 1].copy()
d["T"], d["L"], d["S"], d["E"] = d.os_months_from_diagnosis, d.dx_to_seq_months, d.os_months_from_sequencing, d.os_status
d["x"] = (d.stage_highest_recorded == "IV").astype(int)
R = {"n": int(len(d)), "stage_missing": int(d.stage_highest_recorded.isna().sum())}

def km(t, e, entry=None):
    return KaplanMeierFitter().fit(t, e, entry=entry)

def hr(df, dur, entry=None, strata=True):
    cols = {"t": df[dur].values, "e": df.E.values, "x": df.x.values}
    if strata: cols["ct"] = df.cancer_type.values
    if entry is not None: cols["en"] = np.asarray(entry, float)
    c = pd.DataFrame(cols)
    m = CoxPHFitter().fit(c, "t", "e", strata=(["ct"] if strata else None),
                          entry_col=("en" if entry is not None else None))
    ci = np.exp(m.confidence_intervals_.loc["x"].values)
    return [round(float(np.exp(m.params_["x"])), 2), round(float(ci[0]), 2), round(float(ci[1]), 2)]

# ---------------- 1. landmark ----------------
lt = km(d["T"], d.E, d.L)
R["landmark"] = {}
for t0 in [6, 12]:
    a = d[d["T"] > t0].copy()                      # alive at the landmark
    a["en"] = np.maximum(a.L, t0)                  # delayed entry if sequenced later
    k = km(a["T"], a.E, a.en)
    r = d[(d.L <= t0) & (d["T"] > t0)].copy()      # original: sequenced by the landmark
    kr = km(r["T"], r.E, np.full(len(r), float(t0)))
    # check: conditional median read off the marginal left-truncated curve
    s0 = float(lt.predict(t0)); tl = lt.survival_function_.index.values; sv = lt.survival_function_.iloc[:, 0].values
    cond_med = float(tl[np.argmax(sv <= 0.5 * s0)])
    R["landmark"][t0] = {
        "delayed_entry": {"n": int(len(a)), "n_seq_after": int((a.L > t0).sum()), "deaths": int(a.E.sum()),
                          "median": round(float(k.median_survival_time_), 1),
                          "median_ci": [round(float(v), 1) for v in
                                        __import__("lifelines").utils.median_survival_times(k.confidence_interval_).values[0]],
                          "HR": hr(a, "T", a.en)},
        "restricted_original": {"n": int(len(r)), "deaths": int(r.E.sum()),
                                "median": round(float(kr.median_survival_time_), 1), "HR": hr(r, "T", np.full(len(r), float(t0)))},
        "marginal_LT_conditional_median": round(cond_med, 1)}
    # stage composition of those excluded in the original version
    ex = a[a.L > t0]
    R["landmark"][t0]["excluded_pct_stageI_III"] = round(100 * (ex.x == 0).mean(), 1)
    R["landmark"][t0]["kept_pct_stageI_III"] = round(100 * (r.x == 0).mean(), 1)

# ---------------- 2. conditional Kendall's tau with U-statistic SE ----------------
def ctau(df):
    L, X, D = df.L.values, df["T"].values, df.E.values.astype(bool)
    n = len(L); a = np.zeros(n); b = np.zeros(n)
    for i0 in range(0, n, 400):
        i = np.arange(i0, min(i0 + 400, n))[:, None]; j = np.arange(n)[None, :]
        Li, Lj, Xi, Xj, Di, Dj = L[i], L[j], X[i], X[j], D[i], D[j]
        c = (np.maximum(Li, Lj) <= np.minimum(Xi, Xj)) & (j != i)
        c &= np.where(Xi < Xj, Di, np.where(Xj < Xi, Dj, Di & Dj))
        s = np.sign((Li - Lj) * (Xi - Xj)) * c
        a[i[:, 0]] = s.sum(1); b[i[:, 0]] = c.sum(1)
    tau = a.sum() / b.sum()                         # each pair counted twice, ratio unchanged
    # Hoeffding projection of the ratio of U-statistics (delta method)
    h = (a - tau * b) / (n - 1); ub = b.sum() / (n * (n - 1))
    se = 2 * np.std(h, ddof=1) / sqrt(n) / ub
    z = tau / se
    return {"n": int(n), "pairs": int(b.sum() // 2), "tau": round(float(tau), 3),
            "ci": [round(float(tau - 1.96 * se), 3), round(float(tau + 1.96 * se), 3)], "p": float(erfc(abs(z) / sqrt(2)))}

R["tau_full"] = ctau(d)
R["tau_by_type"] = {ct: ctau(g) for ct, g in d.groupby("cancer_type")}
R["tau_by_type_stage"] = {f"{ct} | {st}": ctau(g) for (ct, st), g in d.groupby(["cancer_type", "stage_highest_recorded"])}

# ---------------- 3. descriptive table by cancer type + per-type HR ----------------
def fu(df):  # reverse Kaplan-Meier median follow-up from sequencing
    k = km(df.S, 1 - df.E); return round(float(k.median_survival_time_), 1)
def q(v): return [round(float(x), 1) for x in np.quantile(v, [.5, .25, .75])]
tab = {}
for ct, g in list(d.groupby("cancer_type")) + [("All", d)]:
    tab[ct] = {"n": int(len(g)), "deaths": int(g.E.sum()), "stageIV_n": int(g.x.sum()),
               "stageIV_pct": round(100 * g.x.mean(), 1), "interval_med_iqr": q(g.L),
               "fu_seq_reverseKM": fu(g),
               "HR_lt": hr(g, "T", g.L, strata=(ct == "All")),
               "HR_naive": hr(g, "T", None, strata=(ct == "All")),
               "lt_median_cond_1mo": round(float(km(g[g["T"] > 1]["T"], g[g["T"] > 1].E,
                                                     np.maximum(g[g["T"] > 1].L, 1.0)).median_survival_time_), 1)}
R["by_type"] = tab

json.dump(R, open("out/results6_landmark_tau_by_type.json", "w"), indent=1, default=str)
print(json.dumps(R, indent=1, default=str))
