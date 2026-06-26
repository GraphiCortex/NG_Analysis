#!/usr/bin/env python3
"""
Compare Julien-promoter × Oshikawa (MAP2) vs Julien-promoter × Oshikawa (CAG).

Inputs:
- Project/New_Data_4/binding_vs_Oshikawa_MAP2_expression.xlsx   (sheet 'All_genes')
- Project/New_Data_5/binding_vs_Oshikawa_CAG_expression.xlsx    (sheet 'All_genes')
- Project/New_Data/GO_results.xlsx                               (sheet 'Sh_2')

Outputs:
- Project/New_Data_5/compare_map2_cag/
    summary_per_family.csv
    overlap_<family>.csv
    map2_only_<family>.csv
    cag_only_<family>.csv
    scatter_<family>.png
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt

# ---------- paths ----------
MAP2_XLSX = "Project/New_Data_4/binding_vs_Oshikawa_MAP2_expression.xlsx"
CAG_XLSX  = "Project/New_Data_5/binding_vs_Oshikawa_CAG_expression.xlsx"
GO_XLSX   = "Project/New_Data/GO_results.xlsx"

OUT_DIR   = "Project/New_Data_5/compare_map2_cag"
os.makedirs(OUT_DIR, exist_ok=True)

FAMILIES = ["Cell cycle","DNA repair","Apoptosis","Autophagy","Necroptosis"]

# ---------- helpers ----------
def need_cols(df: pd.DataFrame, cols: list[str], name: str):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing columns: {missing}")

def build_geneset(go_full: pd.DataFrame, category: str) -> set[str]:
    df = go_full[go_full["Category"] == category]
    if df.empty:
        return set()
    lists = df["geneID"].dropna().astype(str).str.split("/", expand=False)
    return {g.strip().upper() for sub in lists for g in sub}

def safe_stats(x: pd.Series, y: pd.Series):
    if len(x) < 2:
        return (np.nan, np.nan, np.nan, np.nan)
    r, p = pearsonr(x, y)
    rs, ps = spearmanr(x, y)
    return (r, p, rs, ps)

# ---------- load data ----------
map2 = pd.read_excel(MAP2_XLSX, sheet_name="All_genes")
cag  = pd.read_excel(CAG_XLSX,  sheet_name="All_genes")
go_full = pd.read_excel(GO_XLSX, sheet_name="Sh_2")

for df in (map2, cag):
    df["Gene"] = df["Gene"].astype(str).str.upper().str.strip()

need_cols(map2, ["Gene","mean_logFC_bind","mean_logFC_expr","adjP_expr"], "MAP2 All_genes")
need_cols(cag,  ["Gene","mean_logFC_bind","mean_logFC_expr","adjP_expr"], "CAG All_genes")

# rename for clarity (binding is from Julien in both; keep one name)
map2 = map2.rename(columns={
    "mean_logFC_expr": "expr_map2",
    "adjP_expr":       "adjP_map2",
    "mean_logFC_bind": "bind_logFC"
})
cag = cag.rename(columns={
    "mean_logFC_expr": "expr_cag",
    "adjP_expr":       "adjP_cag",
    "mean_logFC_bind": "bind_logFC"
})

# ---------- per-family comparisons ----------
summary_rows = []
for fam in FAMILIES:
    genes = build_geneset(go_full, fam)
    map2_sub = map2[map2["Gene"].isin(genes)].copy()
    cag_sub  = cag[cag["Gene"].isin(genes)].copy()

    mset = set(map2_sub["Gene"])
    cset = set(cag_sub["Gene"])
    overlap   = sorted(mset & cset)
    map2_only = sorted(mset - cset)
    cag_only  = sorted(cset - mset)

    tag = fam.replace(" ","_")

    # save one-column lists
    pd.DataFrame({"Gene": map2_only}).to_csv(
        os.path.join(OUT_DIR, f"map2_only_{tag}.csv"), index=False
    )
    pd.DataFrame({"Gene": cag_only}).to_csv(
        os.path.join(OUT_DIR, f"cag_only_{tag}.csv"), index=False
    )

    # overlap table (keep one binding column)
    ovl = (map2_sub.set_index("Gene")
                    .loc[overlap, ["bind_logFC","expr_map2","adjP_map2"]]
                    .join(cag_sub.set_index("Gene")
                                  .loc[overlap, ["expr_cag","adjP_cag"]])
                    .reset_index())

    if not ovl.empty:
        ovl["sign_concordant"] = np.sign(ovl["expr_map2"]) == np.sign(ovl["expr_cag"])
        ovl["delta_abs_expr"]  = (ovl["expr_map2"].abs() - ovl["expr_cag"].abs()).abs()
        ovl = ovl[["Gene","bind_logFC","expr_map2","adjP_map2",
                   "expr_cag","adjP_cag","sign_concordant","delta_abs_expr"]]
        ovl.to_csv(os.path.join(OUT_DIR, f"overlap_{tag}.csv"), index=False)

        # scatter (MAP2 vs CAG)
        fig, ax = plt.subplots(figsize=(5.2, 4.0))
        ax.scatter(ovl["expr_map2"], ovl["expr_cag"], s=20)
        # regression line (guard tiny n)
        if len(ovl) >= 2:
            m, b = np.polyfit(ovl["expr_map2"], ovl["expr_cag"], 1)
            xs = np.linspace(ovl["expr_map2"].min(), ovl["expr_map2"].max(), 100)
            ax.plot(xs, m*xs + b, linewidth=1)
        r, p, rs, ps = safe_stats(ovl["expr_map2"], ovl["expr_cag"])
        txt = f"n={len(ovl)}\nPearson r={r:.3f} (p={p:.2e})\nSpearman ρ={rs:.3f} (p={ps:.2e})"
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=9)
        ax.set_title(f"MAP2 vs CAG expression • {fam}", pad=8)
        ax.set_xlabel("MAP2 log2FC")
        ax.set_ylabel("CAG log2FC")
        ax.xaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
        ax.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
        plt.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, f"scatter_{tag}.png"), dpi=300)
        plt.close(fig)

        concord = int(ovl["sign_concordant"].sum())
        mean_delta = float(ovl["delta_abs_expr"].mean())
    else:
        r = p = rs = ps = np.nan
        concord = 0
        mean_delta = np.nan

    summary_rows.append({
        "family": fam,
        "map2_only": len(map2_only),
        "cag_only": len(cag_only),
        "overlap": len(overlap),
        "overlap_concordant_signs": concord,
        "mean_delta_abs_expr_overlap": mean_delta,
        "pearson_r_overlap": r,
        "pearson_p_overlap": p,
        "spearman_rho_overlap": rs,
        "spearman_p_overlap": ps
    })

# write summary
summary = pd.DataFrame(summary_rows)
summary_path = os.path.join(OUT_DIR, "summary_per_family.csv")
summary.to_csv(summary_path, index=False)
print("✓ Wrote:", summary_path)
for fam in FAMILIES:
    tag = fam.replace(" ","_")
    print("  - overlap:",   os.path.join(OUT_DIR, f"overlap_{tag}.csv"))
    print("  - map2_only:", os.path.join(OUT_DIR, f"map2_only_{tag}.csv"))
    print("  - cag_only:",  os.path.join(OUT_DIR, f"cag_only_{tag}.csv"))
    print("  - scatter:",   os.path.join(OUT_DIR, f"scatter_{tag}.png"))
