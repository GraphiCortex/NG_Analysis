#!/usr/bin/env python3
import os, numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

FONG_XLSX = "Project/New_Data/binding_vs_TKO_expression.xlsx"
CAG_XLSX  = "Project/New_Data_5/binding_vs_Oshikawa_CAG_expression.xlsx"
GO_XLSX   = "Project/New_Data/GO_results.xlsx"
OUT_DIR   = "Project/New_Data_5/compare_fong_cag"
os.makedirs(OUT_DIR, exist_ok=True)

FAMS = ["Cell cycle","DNA repair","Apoptosis","Autophagy","Necroptosis"]

def need(df, cols, name):
    m=[c for c in cols if c not in df.columns]
    if m: raise ValueError(f"{name} missing: {m}")

def geneset(go, cat):
    df=go[go["Category"]==cat]
    if df.empty: return set()
    return {g.strip().upper() for sub in df["geneID"].dropna().astype(str).str.split("/") for g in sub}

def safe_stats(x,y):
    if len(x)<2: return (np.nan,np.nan,np.nan,np.nan)
    r,p=pearsonr(x,y); rs,ps=spearmanr(x,y); return (r,p,rs,ps)

fong = pd.read_excel(FONG_XLSX, sheet_name="All_genes")
cag  = pd.read_excel(CAG_XLSX,  sheet_name="All_genes")
go   = pd.read_excel(GO_XLSX,   sheet_name="Sh_2")

for df in (fong,cag):
    df["Gene"]=df["Gene"].astype(str).str.upper().str.strip()

need(fong,["Gene","mean_logFC_bind","mean_logFC_expr","adjP_expr"],"Fong")
need(cag, ["Gene","mean_logFC_bind","mean_logFC_expr","adjP_expr"],"CAG")

fong = fong.rename(columns={"mean_logFC_expr":"expr_fong","adjP_expr":"adjP_fong","mean_logFC_bind":"bind_logFC"})
cag  = cag.rename(columns={"mean_logFC_expr":"expr_cag","adjP_expr":"adjP_cag","mean_logFC_bind":"bind_logFC"})

rows=[]
for fam in FAMS:
    gs = geneset(go, fam)
    fsub=fong[fong["Gene"].isin(gs)].copy()
    csub=cag[cag["Gene"].isin(gs)].copy()
    F=set(fsub["Gene"]); C=set(csub["Gene"])
    ov=sorted(F & C); f_only=sorted(F - C); c_only=sorted(C - F)

    tag=fam.replace(" ","_")
    pd.DataFrame({"Gene":f_only}).to_csv(os.path.join(OUT_DIR,f"fong_only_{tag}.csv"),index=False)
    pd.DataFrame({"Gene":c_only}).to_csv(os.path.join(OUT_DIR,f"cag_only_{tag}.csv"),index=False)

    ovl=(fsub.set_index("Gene").loc[ov,["bind_logFC","expr_fong","adjP_fong"]]
         .join(csub.set_index("Gene").loc[ov,["expr_cag","adjP_cag"]]).reset_index())
    if not ovl.empty:
        ovl["sign_concordant"]=np.sign(ovl["expr_fong"])==np.sign(ovl["expr_cag"])
        ovl["delta_abs_expr"]=(ovl["expr_fong"].abs()-ovl["expr_cag"].abs()).abs()
        ovl=ovl[["Gene","bind_logFC","expr_fong","adjP_fong","expr_cag","adjP_cag","sign_concordant","delta_abs_expr"]]
        ovl.to_csv(os.path.join(OUT_DIR,f"overlap_{tag}.csv"),index=False)

        fig,ax=plt.subplots(figsize=(5.2,4.0))
        ax.scatter(ovl["expr_fong"],ovl["expr_cag"],s=20)
        if len(ovl)>=2:
            m,b=np.polyfit(ovl["expr_fong"],ovl["expr_cag"],1)
            xs=np.linspace(ovl["expr_fong"].min(),ovl["expr_fong"].max(),100); ax.plot(xs,m*xs+b,linewidth=1)
        r,p,rs,ps=safe_stats(ovl["expr_fong"],ovl["expr_cag"])
        ax.text(0.03,0.97,f"n={len(ovl)}\nPearson r={r:.3f} (p={p:.2e})\nSpearman ρ={rs:.3f} (p={ps:.2e})",
                transform=ax.transAxes,va="top",fontsize=9)
        ax.set_title(f"Fong (Rb-TKO) vs CAG expression • {fam}",pad=8)
        ax.set_xlabel("Fong log2FC (Rb-TKO)"); ax.set_ylabel("CAG log2FC")
        ax.xaxis.grid(True,linestyle="--",linewidth=0.5,alpha=0.7)
        ax.yaxis.grid(True,linestyle="--",linewidth=0.5,alpha=0.7)
        plt.tight_layout(); fig.savefig(os.path.join(OUT_DIR,f"scatter_{tag}.png"),dpi=300); plt.close(fig)
        concord=int(ovl["sign_concordant"].sum()); mean_delta=float(ovl["delta_abs_expr"].mean())
    else:
        r=p=rs=ps=np.nan; concord=0; mean_delta=np.nan

    rows.append({"family":fam,
                 "fong_only":len(f_only),"cag_only":len(c_only),"overlap":len(ov),
                 "overlap_concordant_signs":concord,"mean_delta_abs_expr_overlap":mean_delta,
                 "pearson_r_overlap":r,"pearson_p_overlap":p,
                 "spearman_rho_overlap":rs,"spearman_p_overlap":ps})

summary=pd.DataFrame(rows)
summary_path=os.path.join(OUT_DIR,"summary_per_family.csv")
summary.to_csv(summary_path,index=False)
print("✓ Wrote:", summary_path)
for fam in FAMS:
    tag=fam.replace(" ","_")
    print("  - overlap:",   os.path.join(OUT_DIR,f"overlap_{tag}.csv"))
    print("  - fong_only:", os.path.join(OUT_DIR,f"fong_only_{tag}.csv"))
    print("  - cag_only:",  os.path.join(OUT_DIR,f"cag_only_{tag}.csv"))
    print("  - scatter:",   os.path.join(OUT_DIR,f"scatter_{tag}.png"))
