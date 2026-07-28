import pandas as pd, numpy as np, json
import matplotlib as mpl; mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from lifelines import KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts

BLU="#0072B2"; VER="#D55E00"; GREY="#7A7A7A"
plt.rcParams.update({"font.size":11,"font.family":"DejaVu Sans","axes.spines.top":False,
    "axes.spines.right":False,"axes.grid":True,"grid.color":"#ECECEC","grid.linewidth":0.7,
    "axes.axisbelow":True,"figure.dpi":150})

d = pd.read_csv("out/merged_raw.csv")
d = d[d["delay_months"].notna() & d["OS_MONTHS"].notna()]; d=d[d["delay_months"]>=0].copy()
d["T_dx"]=d["delay_months"]+d["OS_MONTHS"]; d["entry"]=d["delay_months"]

types=d["cancer_type"].value_counts().head(5).index.tolist()
lab={"Non-Small Cell Lung Cancer":"NSCLC","Colorectal Cancer":"Colorectal",
     "Breast Cancer":"Breast","Prostate Cancer":"Prostate","Pancreatic Cancer":"Pancreatic",
     "Bladder Cancer":"Bladder","Glioma":"Glioma","Renal Cell Carcinoma":"Renal"}

def fit(sub,lt):
    k=KaplanMeierFitter()
    if lt: k.fit(sub["T_dx"],sub["deceased"],entry=sub["entry"])
    else:  k.fit(sub["T_dx"],sub["deceased"])
    return k
def S(k,t):
    try: return float(k.survival_function_at_times(t).iloc[0])
    except Exception: return np.nan
def med(k):
    m=k.median_survival_time_; return m if np.isfinite(m) else np.nan

rows=[]
for ct in types:
    sub=d[d["cancer_type"]==ct]
    kn=fit(sub,False); kc=fit(sub,True)
    rows.append(dict(type=lab.get(ct,ct),n=len(sub),deaths=int(sub["deceased"].sum()),
        delay=round(sub["delay_months"].median(),1),
        mnaive=med(kn),mcorr=med(kc),
        s24n=S(kn,24)*100,s24c=S(kc,24)*100,
        s60n=S(kn,60)*100,s60c=S(kc,60)*100,
        s120n=S(kn,120)*100,s120c=S(kc,120)*100))
T=pd.DataFrame(rows)
T["over60"]=T["s60n"]-T["s60c"]; T["over24"]=T["s24n"]-T["s24c"]; T["over120"]=T["s120n"]-T["s120c"]
print(T[["type","n","deaths","delay","mnaive","mcorr","s60n","s60c","over60","over120"]].round(1).to_string(index=False))
T.to_csv("out/per_type.csv",index=False)

# ---------------- FIG 5: KM small multiples per cancer type ----------------
from matplotlib.lines import Line2D
fig,axes=plt.subplots(2,3,figsize=(12,8.4),sharex=True,sharey=True)
axl=axes.ravel()
for ax,ct in zip(axl,types):
    sub=d[d["cancer_type"]==ct]
    kn=fit(sub,False); kc=fit(sub,True)
    kn._label="Naive"; kc._label="L-trunc"
    xn=kn.survival_function_.index.values; yn=kn.survival_function_.iloc[:,0].values
    xc=kc.survival_function_.index.values; yc=kc.survival_function_.iloc[:,0].values
    ax.step(xn,yn,where="post",color=VER,lw=2.1)
    ax.step(xc,yc,where="post",color=BLU,lw=2.1)
    ax.axhline(0.5,color="#CCCCCC",lw=0.9,ls=":")
    ax.set_xlim(0,150); ax.set_ylim(0,1)
    ax.set_title(lab.get(ct,ct)+f"  (n={len(sub):,})",loc="left",fontsize=10.5,fontweight="bold")
    sub60n=S(kn,60)*100; sub60c=S(kc,60)*100
    ax.text(0.96,0.94,f"5-yr OS {sub60n:.0f}% vs {sub60c:.0f}%",transform=ax.transAxes,
            ha="right",va="top",fontsize=8.8,color="#333333")
    ax.set_xlabel("Months since diagnosis")
    add_at_risk_counts(kn,kc,ax=ax,labels=["Naive","L-trunc"],rows_to_show=["At risk"],xticks=[0,50,100,150])
# legend panel in the empty 6th cell
lg=axl[5]; lg.axis("off")
lg.legend([Line2D([0],[0],color=VER,lw=3),Line2D([0],[0],color=BLU,lw=3)],
          ["From diagnosis, naive","From diagnosis,\nleft-truncated"],
          loc="center",frameon=False,fontsize=11,handlelength=1.6)
lg.text(0.5,0.20,"Dotted line = median (50%)",transform=lg.transAxes,ha="center",fontsize=8.5,color="#777777")
for ax in [axl[0],axl[3]]: ax.set_ylabel("Overall survival")
plt.tight_layout(); plt.savefig("out/fig5_km_by_type.png",bbox_inches="tight"); plt.close()

# ---------------- FIG 6: heatmap of survival overestimation (pp) ----------------
cols=["24 mo","60 mo","120 mo"]
M=T.set_index("type")[["over24","over60","over120"]].values
fig,ax=plt.subplots(figsize=(6.6,4.2))
im=ax.imshow(M,cmap="Blues",aspect="auto",vmin=0,vmax=max(20,np.nanmax(M)))
ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols)
ax.set_yticks(range(len(T))); ax.set_yticklabels(T["type"])
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        v=M[i,j]
        ax.text(j,i,f"+{v:.0f}",ha="center",va="center",
                color="white" if v>np.nanmax(M)*0.55 else "#222222",fontsize=10,fontweight="bold")
ax.set_xlabel("Time since diagnosis"); ax.set_title("Overestimation of survival by naive analysis\n(percentage points, naive minus corrected)",
             loc="left",fontsize=10.5,fontweight="bold")
cb=fig.colorbar(im,ax=ax,fraction=0.046,pad=0.04); cb.set_label("Overestimation (pp)")
ax.set_xticks(np.arange(-.5,3,1),minor=True); ax.set_yticks(np.arange(-.5,len(T),1),minor=True)
ax.grid(which="minor",color="white",lw=2); ax.tick_params(which="minor",length=0)
plt.tight_layout(); plt.savefig("out/fig6_heatmap_bias.png",bbox_inches="tight"); plt.close()
print("\nfig5 and fig6 written")
