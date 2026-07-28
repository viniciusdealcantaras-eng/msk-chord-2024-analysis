import pandas as pd, numpy as np, json
from lifelines import CoxPHFitter

d = pd.read_csv("out/merged_raw.csv")
d = d[d["delay_months"].notna() & d["OS_MONTHS"].notna()].copy()
d = d[d["delay_months"]>=0].copy()
d["T_dx"]=d["delay_months"]+d["OS_MONTHS"]; d["entry"]=d["delay_months"]
cx = d[d["STAGE_HIGHEST_RECORDED"].isin(["Stage 1-3","Stage 4"])].copy()
cx["stage4"] = (cx["STAGE_HIGHEST_RECORDED"]=="Stage 4").astype(int)

def strat_cox(dur, entry=None):
    cols=[dur,"deceased","stage4","cancer_type"]+([entry] if entry else [])
    cph=CoxPHFitter()
    cph.fit(cx[cols], duration_col=dur, event_col="deceased",
            strata=["cancer_type"], entry_col=(entry if entry else None))
    return {"HR":round(float(np.exp(cph.params_["stage4"])),2),
            "lo":round(float(np.exp(cph.confidence_intervals_.loc["stage4"].iloc[0])),2),
            "hi":round(float(np.exp(cph.confidence_intervals_.loc["stage4"].iloc[1])),2),
            "p":float(cph.summary.loc["stage4","p"])}

out={"naive":strat_cox("T_dx"),
     "corrected":strat_cox("T_dx",entry="entry"),
     "from_seq":strat_cox("OS_MONTHS")}
json.dump(out, open("out/results3_strat.json","w"), indent=2)
print(json.dumps(out, indent=2))
