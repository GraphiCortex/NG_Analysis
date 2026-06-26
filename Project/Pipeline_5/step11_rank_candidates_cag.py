#!/usr/bin/env python3
import os, pandas as pd
MERGED_XLSX = "Project/New_Data_5/binding_vs_Oshikawa_CAG_expression.xlsx"
GO_XLSX     = "Project/New_Data/GO_results.xlsx"
OUT_DIR     = "Project/New_Data_5/new_data"
os.makedirs(OUT_DIR, exist_ok=True)

families = ["Apoptosis","Autophagy","Necroptosis"]
outs = {
    "Apoptosis":   os.path.join(OUT_DIR,"apoptosis_candidates_ranked_cag.csv"),
    "Autophagy":   os.path.join(OUT_DIR,"autophagy_candidates_ranked_cag.csv"),
    "Necroptosis": os.path.join(OUT_DIR,"necroptosis_candidates_ranked_cag.csv"),
}

merged = pd.read_excel(MERGED_XLSX, sheet_name="All_genes")
go = pd.read_excel(GOXLSX:=GO_XLSX, sheet_name="Sh_2")

def geneset(cat):
    df = go[go["Category"]==cat]
    if df.empty: return set()
    return {g.strip().upper() for sub in df["geneID"].dropna().astype(str).str.split("/") for g in sub}

for fam in families:
    gs = geneset(fam)
    df = merged[merged["Gene"].isin(gs)].copy()
    df["abs_expr"]=df["mean_logFC_expr"].abs()
    df["abs_bind"]=df["mean_logFC_bind"].abs()
    ranked=df.sort_values(by=["abs_expr","abs_bind"],ascending=[False,False])
    ranked.to_csv(outs[fam],index=False,
                  columns=["Gene","mean_logFC_expr","adjP_expr",
                           "mean_logFC_bind","P_Value_bind","abs_expr","abs_bind"])
    print(f"→ {fam}: {outs[fam]}")
