#!/usr/bin/env python3
"""
Pipeline/step11_rank_celldeath_families.py

Rank candidate genes for three families:
  - Apoptosis
  - Autophagy
  - Necroptosis

Ranking per family:
  1) |expression log2FC| (descending)
  2) |binding   log2FC| (descending)

Writes one CSV per family.
"""

import os
import pandas as pd

MERGED_XLSX = "Project/New_Data/binding_vs_TKO_expression.xlsx"
GO_XLSX     = "Project/New_Data/GO_results.xlsx"
OUT_DIR     = "Project/New_Data/new_data"

FAMILIES = ["Apoptosis", "Autophagy", "Necroptosis"]
OUTFILES = {
    "Apoptosis":   os.path.join(OUT_DIR, "apoptosis_candidates_ranked.csv"),
    "Autophagy":   os.path.join(OUT_DIR, "autophagy_candidates_ranked.csv"),
    "Necroptosis": os.path.join(OUT_DIR, "necroptosis_candidates_ranked.csv"),
}

# 1) Load merged binding+expression table
merged = pd.read_excel(MERGED_XLSX, sheet_name="All_genes")

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
    fam_df = merged[ merged["Gene"].isin(genes) ].copy()

    # add ranking columns
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
            "mean_logFC_expr", "adjP_expr",
            "mean_logFC_bind", "P_Value_bind",
            "abs_expr", "abs_bind"
        ]
    )
    print(f"→ {fam}: wrote {len(ranked)} candidates to {OUTFILES[fam]}")
