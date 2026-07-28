import pandas as pd, numpy as np
from lifelines import KaplanMeierFitter

d = pd.read_csv("out/merged_raw.csv")
d = d[d["delay_months"].notna() & d["OS_MONTHS"].notna()].copy()
d = d[d["delay_months"]>=0].copy()
d["T_dx"]=d["delay_months"]+d["OS_MONTHS"]; d["entry"]=d["delay_months"]

# ---- independent left-truncated KM (risk set with entry/exit) ----
entry=d["entry"].values; exitt=d["T_dx"].values; ev=d["deceased"].values
# death times
dt=np.sort(np.unique(exitt[ev==1]))
S=1.0; med=None; surv={}
for t in dt:
    at_risk=np.sum((entry< t) & (exitt>=t))          # entered before t, still under obs at t
    deaths=np.sum((exitt==t) & (ev==1))
    if at_risk>0:
        S*= (1-deaths/at_risk)
    surv[t]=S
    if med is None and S<=0.5: med=t
# S at 24/60
def Sat(tt):
    ks=[t for t in dt if t<=tt]
    return surv[ks[-1]] if ks else 1.0
print("[MANUAL left-trunc KM]  median=%.2f  S24=%.4f  S60=%.4f"%(med,Sat(24),Sat(60)))

kmf=KaplanMeierFitter(); kmf.fit(d["T_dx"],d["deceased"],entry=d["entry"])
print("[lifelines left-trunc ] median=%.2f  S24=%.4f  S60=%.4f"%(
    kmf.median_survival_time_,
    kmf.survival_function_at_times(24).iloc[0],
    kmf.survival_function_at_times(60).iloc[0]))

# naive manual (entry=0 for all)
S=1.0; medN=None
for t in dt:
    at_risk=np.sum(exitt>=t); deaths=np.sum((exitt==t)&(ev==1))
    S*=(1-deaths/at_risk)
    if medN is None and S<=0.5: medN=t
print("[MANUAL naive KM      ] median=%.2f"%medN)

# immortal person-time cross-check
print("immortal%% = %.2f"%(d["delay_months"].sum()/d["T_dx"].sum()*100))
# sanity: everyone's T_dx >= entry
print("all T_dx>=entry:", bool((d["T_dx"]>=d["entry"]).all()))
print("deaths:",int(d['deceased'].sum()),"N:",len(d))
