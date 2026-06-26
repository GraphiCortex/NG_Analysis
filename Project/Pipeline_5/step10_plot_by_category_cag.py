#!/usr/bin/env python3
import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import pearsonr

IN_MERGED = os.path.join("Project","New_Data_5","binding_vs_Oshikawa_CAG_expression.xlsx")
IN_GOXLS  = os.path.join("Project","New_Data","GO_results.xlsx")
OUT_DIR   = os.path.join("Project","New_Data_5","figures","cag")
os.makedirs(OUT_DIR, exist_ok=True)

FIGSIZE=(5.5,4.0)

merged = pd.read_excel(IN_MERGED, sheet_name="All_genes")
go = pd.read_excel(IN_GOXLS, sheet_name="Sh_2")

def geneset(cat):
    df = go[go["Category"]==cat]
    if df.empty: return set()
    return {g.strip().upper() for sub in df["geneID"].dropna().astype(str).str.split("/") for g in sub}

cats = {
    "Cell cycle":  os.path.join(OUT_DIR,"Scatter_Cell_cycle_cag.png"),
    "DNA repair":  os.path.join(OUT_DIR,"Scatter_DNA_repair_cag.png"),
    "Apoptosis":   os.path.join(OUT_DIR,"Scatter_Apoptosis_cag.png"),
    "Autophagy":   os.path.join(OUT_DIR,"Scatter_Autophagy_cag.png"),
    "Necroptosis": os.path.join(OUT_DIR,"Scatter_Necroptosis_cag.png"),
}

# common axis
lims=[]
subs={}
for c in cats:
    gs=geneset(c)
    s=merged[merged["Gene"].isin(gs)].copy()
    subs[c]=s
    if not s.empty:
        lims.append([s["mean_logFC_bind"].min(), s["mean_logFC_bind"].max(),
                     s["mean_logFC_expr"].min(), s["mean_logFC_expr"].max()])
if lims:
    arr=np.array(lims); X_LIM=(arr[:,0].min(), arr[:,1].max()); Y_LIM=(arr[:,2].min(), arr[:,3].max())
else:
    X_LIM=(-1,1); Y_LIM=(-1,1)

for c,out in cats.items():
    s=subs[c]
    x=s["mean_logFC_bind"] if not s.empty else pd.Series(dtype=float)
    y=s["mean_logFC_expr"] if not s.empty else pd.Series(dtype=float)
    fig,ax=plt.subplots(figsize=FIGSIZE); ax.scatter(x,y,s=18)
    n=len(s)
    if n>=2:
        m,b=np.polyfit(x,y,1); xs=np.linspace(*X_LIM,100); ax.plot(xs,m*xs+b,linewidth=1)
        r,p=pearsonr(x,y); txt=f"n = {n}\nr = {r:.2f}\np = {p:.2e}"
    elif n==1: txt="n = 1\n(no correlation)"
    else: txt="n = 0"
    ax.text(0.03,0.97,txt,transform=ax.transAxes,va="top",fontsize=9)
    ax.set_title(f"Binding vs. Expression (CAG): {c}",pad=8,fontsize=11)
    ax.set_xlabel("Binding log$_2$FC"); ax.set_ylabel("Expression log$_2$FC")
    ax.xaxis.grid(True,linestyle="--",linewidth=0.5,alpha=0.7)
    ax.yaxis.grid(True,linestyle="--",linewidth=0.5,alpha=0.7)
    ax.set_xlim(*X_LIM); ax.set_ylim(*Y_LIM)
    plt.tight_layout(); fig.savefig(out,dpi=300); plt.close(fig)
    print(f"→ Wrote {out}")
