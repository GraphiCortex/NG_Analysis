#!/usr/bin/env python3
"""
Rank candidate genes for three families in the Andrusiak legacy-filtered set:
  - Apoptosis
  - Autophagy
  - Necroptosis

Ranking:
  1) |expression log2FC| (descending)
  2) |binding   log2FC| (descending)
"""

import os
import pandas as pd

MERGED_XLSX = "Project/New_Data_6/binding_vs_Andrusiak_expression.xlsx"
GO_XLSX     = "Project/New_Data/GO_results.xlsx"
OUT_DIR     = "Project/New_Data_6/new_data"
os.makedirs(OUT_DIR, exist_ok=True)

FAMILIES = ["Apoptosis", "Autophagy", "Necroptosis"]
OUTFILES = {
    "Apoptosis":   os.path.join(OUT_DIR, "apoptosis_candidates_ranked_andrusiak.csv"),
    "Autophagy":   os.path.join(OUT_DIR, "autophagy_candidates_ranked_andrusiak.csv"),
    "Necroptosis": os.path.join(OUT_DIR, "necroptosis_candidates_ranked_andrusiak.csv"),
}

# 1) Load merged binding+expression table
merged = pd.read_excel(MERGED_XLSX, sheet_name="All_genes")
merged["Gene"] = merged["Gene"].astype(str).str.upper().str.strip()

# 2) Load GO results to build gene sets
go_full = pd.read_excel(GO_XLSX, sheet_name="Sh_2")

def build_geneset(category: str) -> set[str]:
    df = go_full[go_full["Category"] == category]
    if df.empty:
        return set()
    lists = df["geneID"].dropna().astype(str).str.split("/", expand=False)
    return {g.strip().upper() for sub in lists for g in sub}

# 3) Rank and write per family
for fam in FAMILIES:
    genes = build_geneset(fam)
    fam_df = merged[merged["Gene"].isin(genes)].copy()

    fam_df["abs_expr"] = fam_df["mean_logFC_expr"].abs()
    fam_df["abs_bind"] = fam_df["mean_logFC_bind"].abs()

    ranked = fam_df.sort_values(
        by=["abs_expr", "abs_bind"],
        ascending=[False, False]
    )

    ranked.to_csv(
        OUTFILES[fam],
        index=False,
        columns=[
            "Gene",
            "mean_logFC_expr",
            "mean_logFC_bind",
            "P_Value_bind",
            "adjP_bind",
            "abs_expr",
            "abs_bind"
        ]
    )
    print(f"→ {fam}: wrote {len(ranked)} candidates to {OUTFILES[fam]}")