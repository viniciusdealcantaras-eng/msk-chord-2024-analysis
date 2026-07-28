import pandas as pd, numpy as np, json
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts

BLU="#0072B2"; ORA="#E69F00"; GRN="#009E73"; VER="#D55E00"; PUR="#CC79A7"; SKY="#56B4E9"; GREY="#7A7A7A"
plt.rcParams.update({"font.size":11,"font.family":"DejaVu Sans","axes.spines.top":False,
    "axes.spines.right":False,"axes.grid":True,"grid.color":"#E6E6E6","grid.linewidth":0.8,
    "axes.axisbelow":True,"figure.dpi":150})

d = pd.read_csv("out/merged_raw.csv")
d = d[d["delay_months"].notna() & d["OS_MONTHS"].notna()]; d=d[d["delay_months"]>=0].copy()
d["T_dx"]=d["delay_months"]+d["OS_MONTHS"]; d["entry"]=d["delay_months"]; d["T_seq"]=d["OS_MONTHS"]
R=json.load(open("out/results.json"))

def km_curve(dur,ev,entry=None):
    k=KaplanMeierFitter(); k.fit(dur,ev,entry=entry)
    return k.survival_function_.index.values, k.survival_function_.iloc[:,0].values, k

# ================= FIG 1 =================
fig,ax=plt.subplots(1,2,figsize=(11,4.3))
dl=d["delay_months"].clip(lower=0.03)
bins=np.logspace(np.log10(0.03),np.log10(dl.max()),45)
ax[0].hist(dl,bins=bins,color=BLU,alpha=0.85,edgecolor="white",linewidth=0.4)
ax[0].set_xscale("log")
ax[0].axvline(R["delay_median"],color=VER,lw=2,ls="--")
ax[0].text(R["delay_median"]*1.15,ax[0].get_ylim()[1]*0.9,f"median\n{R['delay_median']:.1f} mo",
           color=VER,fontsize=9.5,va="top")
ax[0].axvline(R["delay_p90"],color=GREY,lw=1.5,ls=":")
ax[0].text(R["delay_p90"]*1.05,ax[0].get_ylim()[1]*0.55,f"P90\n{R['delay_p90']:.0f} mo",color=GREY,fontsize=9)
ax[0].set_xlabel("Diagnosis-to-sequencing interval (months, log scale)")
ax[0].set_ylabel("Number of patients")
ax[0].set_title("(a) Immortal time per patient",loc="left",fontweight="bold",fontsize=11)

types=list(R["delay_by_type"].keys())
labels=[t.replace(" Cancer","").replace("Non-Small Cell Lung","NSCLC") for t in types]
meds=[R["delay_by_type"][t]["median"] for t in types]
q1=[R["delay_by_type"][t]["q25"] for t in types]; q3=[R["delay_by_type"][t]["q75"] for t in types]
order=np.argsort(meds)
labels=[labels[i] for i in order]; meds=[meds[i] for i in order]
q1=[q1[i] for i in order]; q3=[q3[i] for i in order]
y=np.arange(len(labels))
ax[1].hlines(y,q1,q3,color=SKY,lw=6,alpha=0.6)
ax[1].plot(meds,y,"o",color=BLU,ms=8,zorder=5)
for yi,m in zip(y,meds): ax[1].text(m,yi+0.18,f"{m:.1f}",fontsize=8.5,ha="center",color=BLU)
ax[1].set_yticks(y); ax[1].set_yticklabels(labels)
ax[1].set_xlabel("Median delay (months); bar = IQR")
ax[1].set_title("(b) Delay by cancer type",loc="left",fontweight="bold",fontsize=11)
ax[1].grid(axis="y",visible=False)
plt.tight_layout(); plt.savefig("out/fig1_immortal_time.png",bbox_inches="tight"); plt.close()

# ================= FIG 2 =================
fig,ax=plt.subplots(figsize=(8.0,5.6))
_,_,kn=km_curve(d["T_dx"],d["deceased"])
_,_,kc=km_curve(d["T_dx"],d["deceased"],entry=d["entry"])
_,_,ks=km_curve(d["T_seq"],d["deceased"])
for k,c,lw,ls,labl in [(kn,VER,2.4,"-",f"From diagnosis, naive (median {R['km_from_dx_naive']['median']:.0f} mo)"),
                       (kc,BLU,2.4,"-",f"From diagnosis, left-truncated (median {R['km_from_dx_corrected']['median']:.0f} mo)"),
                       (ks,GREY,1.8,"--",f"From sequencing, published (median {R['km_from_seq']['median']:.0f} mo)")]:
    ax.step(k.survival_function_.index.values,k.survival_function_.iloc[:,0].values,where="post",color=c,lw=lw,ls=ls,label=labl)
ax.axhline(0.5,color="#BBBBBB",lw=1,ls=":")
ax.set_xlim(0,150); ax.set_ylim(0,1)
ax.set_xlabel("Months"); ax.set_ylabel("Estimated overall survival")
ax.set_title("Immortal time bias in survival from diagnosis",loc="left",fontweight="bold",fontsize=11.5)
ax.annotate("",xy=(R['km_from_dx_corrected']['median'],0.5),xytext=(R['km_from_dx_naive']['median'],0.5),
            arrowprops=dict(arrowstyle="<->",color="#333333",lw=1.3))
ax.text((R['km_from_dx_corrected']['median']+R['km_from_dx_naive']['median'])/2,0.53,
        f"+{R['median_overest_months']:.0f} months",ha="center",fontsize=10,fontweight="bold")
ax.legend(loc="upper right",fontsize=9.2,frameon=False)
add_at_risk_counts(kn,kc,ks,ax=ax,labels=["Naive","Left-truncated","From sequencing"],rows_to_show=["At risk"])
plt.tight_layout(); plt.savefig("out/fig2_km_core.png",bbox_inches="tight"); plt.close()

# ================= FIG 3 =================
med_delay=d["delay_months"].median()
d["late"]=d["delay_months"]>=med_delay
fig,ax=plt.subplots(1,2,figsize=(12,5.6),sharey=True)
fitsA=[]; fitsB=[]; labs=[]
for grp,col,lab in [(False,GRN,"Sequenced early (< median)"),(True,PUR,"Sequenced late (>= median)")]:
    g=d[d["late"]==grp]
    _,_,ka=km_curve(g["T_dx"],g["deceased"])
    _,_,kb=km_curve(g["T_dx"],g["deceased"],entry=g["entry"])
    ax[0].step(ka.survival_function_.index.values,ka.survival_function_.iloc[:,0].values,where="post",color=col,lw=2.3,label=lab)
    ax[1].step(kb.survival_function_.index.values,kb.survival_function_.iloc[:,0].values,where="post",color=col,lw=2.3,label=lab)
    fitsA.append(ka); fitsB.append(kb); labs.append("Early" if not grp else "Late")
for a in ax:
    a.axhline(0.5,color="#BBBBBB",lw=1,ls=":"); a.set_xlim(0,150); a.set_ylim(0,1); a.set_xlabel("Months since diagnosis")
ax[0].set_ylabel("Estimated overall survival")
sp=R["spurious_late_vs_early"]
ax[0].set_title(f"(a) Naive: 'late' appears to live longer\n(medians {sp['late (>=median delay)']['median_from_dx_naive']:.0f} vs {sp['early (<median delay)']['median_from_dx_naive']:.0f} mo)",loc="left",fontsize=10.5,fontweight="bold")
ax[1].set_title(f"(b) Corrected: the advantage reverses\n(medians {sp['late (>=median delay)']['median_from_dx_corrected']:.0f} vs {sp['early (<median delay)']['median_from_dx_corrected']:.0f} mo)",loc="left",fontsize=10.5,fontweight="bold")
ax[0].legend(loc="lower left",fontsize=9,frameon=False)
add_at_risk_counts(*fitsA,ax=ax[0],labels=labs,rows_to_show=["At risk"])
add_at_risk_counts(*fitsB,ax=ax[1],labels=labs,rows_to_show=["At risk"])
plt.tight_layout(); plt.savefig("out/fig3_spurious.png",bbox_inches="tight"); plt.close()

# ================= FIG 4 (Stage IV HR, adjusted for cancer type) =================
cs=json.load(open("out/results3_strat.json"))
fig,ax=plt.subplots(figsize=(7.6,2.8))
items=[("From diagnosis, naive",cs["naive"],VER),
       ("From diagnosis, left-truncated",cs["corrected"],BLU),
       ("From sequencing, published",cs["from_seq"],GREY)]
for i,(name,r,col) in enumerate(items):
    yy=len(items)-1-i
    ax.plot([r["lo"],r["hi"]],[yy,yy],color=col,lw=2.4)
    ax.plot(r["HR"],yy,"o",color=col,ms=8)
    ax.text(r["hi"]+0.04,yy,f'{r["HR"]:.2f} ({r["lo"]:.2f}–{r["hi"]:.2f})',va="center",fontsize=9,color=col)
    ax.text(-0.02,yy,name,va="center",ha="right",fontsize=9.2,transform=ax.get_yaxis_transform())
ax.axvline(1,color="#999999",lw=1,ls="--")
ax.set_yticks([]); ax.set_ylim(-0.5,len(items)-0.5)
ax.set_xlabel("Hazard ratio for stage IV vs I–III (95% CI)"); ax.set_xlim(0.9,3.4)
ax.set_title("Immortal time inflates the stage-IV hazard ratio",loc="left",fontweight="bold",fontsize=11)
ax.spines["left"].set_visible(False); ax.grid(axis="y",visible=False)
plt.tight_layout(); plt.savefig("out/fig4_cox_forest.png",bbox_inches="tight"); plt.close()
print("figures rewritten in English")
