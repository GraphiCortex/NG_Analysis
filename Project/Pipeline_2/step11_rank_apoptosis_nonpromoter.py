#!/usr/bin/env python3
"""
Pipeline_2 / step11_rank_apoptosis_nonpromoter.py

Rank apoptosis-related candidates from the NON-PROMOTER pipeline by:
1) |mean_logFC_expr| (desc), then
2) |mean_logFC_bind| (desc).

Inputs:
  - Project/New_Data_2/binding_vs_TKO_expression_nonpromoter.xlsx  (sheet 'All_genes')
  - Project/New_Data_2/GO_results_nonpromoter.xlsx                 (sheet 'Sh_2')

Output:
  - Project/New_Data_2/apoptosis_candidates_ranked_nonpromoter.csv
"""

import pandas as pd
from pathlib import Path

IN_DIR = Path("Project/New_Data_2")
merged_xlsx = IN_DIR / "binding_vs_TKO_expression_nonpromoter.xlsx"
go_xlsx     = IN_DIR / "GO_results_nonpromoter.xlsx"
out_csv     = IN_DIR / "apoptosis_candidates_ranked_nonpromoter.csv"

# 1) Load merged binding+expression table
merged = pd.read_excel(merged_xlsx, sheet_name="All_genes")

# 2) Load GO results and pull the *true* Apoptosis gene set (from Sh_2)
go_full = pd.read_excel(go_xlsx, sheet_name="Sh_2")
ap_df   = go_full[go_full["Category"] == "Apoptosis"]

# 3) Build uppercase gene set
if ap_df.empty or ap_df["geneID"].dropna().empty:
    print("No 'Apoptosis' terms found in GO workbook; nothing to rank.")
    pd.DataFrame(columns=[
        "Gene","mean_logFC_expr","adjP_expr","mean_logFC_bind","P_Value_bind",
        "abs_expr","abs_bind"
    ]).to_csv(out_csv, index=False)
    raise SystemExit(0)

ap_genes = {
    g.strip().upper()
    for sub in ap_df["geneID"].dropna().astype(str).str.split("/")
    for g in sub
}

# 4) Filter merged table to apoptosis genes (intersection)
ap_merged = merged[ merged["Gene"].str.upper().isin(ap_genes) ].copy()

# 5) Add ranking helpers
ap_merged["abs_expr"] = ap_merged["mean_logFC_expr"].abs()
ap_merged["abs_bind"] = ap_merged["mean_logFC_bind"].abs()

# 6) Sort by abs_expr desc, then abs_bind desc
ranked = ap_merged.sort_values(
    by=["abs_expr", "abs_bind"],
    ascending=[False, False]
)

# 7) Write full ranked list
cols = [
    "Gene",
    "mean_logFC_expr", "adjP_expr",
    "mean_logFC_bind", "P_Value_bind",
    "abs_expr", "abs_bind",
]
ranked.to_csv(out_csv, index=False, columns=[c for c in cols if c in ranked.columns])

print(f"→ Wrote {len(ranked)} apoptosis candidates to {out_csv}")
