#!/usr/bin/env python3
"""
Compare Julien-promoter × Fong vs Julien-promoter × Andrusiak.

Inputs:
- Project/New_Data/binding_vs_TKO_expression.xlsx          (sheet 'All_genes')
- Project/New_Data_6/binding_vs_Andrusiak_expression.xlsx  (sheet 'All_genes')
- Project/New_Data/GO_results.xlsx                         (sheet 'Sh_2')

Outputs:
- Project/New_Data_6/compare_fong_andrusiak/
    summary_per_family.csv
    overlap_<family>.csv
    fong_only_<family>.csv
    andrusiak_only_<family>.csv
    scatter_<family>.png
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import matplotlib.pyplot as plt

# ---------- paths ----------
FONG_XLSX = "Project/New_Data/binding_vs_TKO_expression.xlsx"
ANDR_XLSX = "Project/New_Data_6/binding_vs_Andrusiak_expression.xlsx"
GO_XLSX   = "Project/New_Data/GO_results.xlsx"
OUT_DIR   = "Project/New_Data_6/compare_fong_andrusiak"
os.makedirs(OUT_DIR, exist_ok=True)

FAMILIES = ["Cell cycle", "DNA repair", "Apoptosis", "Autophagy", "Necroptosis"]

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
fong = pd.read_excel(FONG_XLSX, sheet_name="All_genes")
andr = pd.read_excel(ANDR_XLSX, sheet_name="All_genes")
go_full = pd.read_excel(GO_XLSX, sheet_name="Sh_2")

for df in (fong, andr):
    df["Gene"] = df["Gene"].astype(str).str.upper().str.strip()

need_cols(fong, ["Gene", "mean_logFC_bind", "mean_logFC_expr"], "Fong All_genes")
need_cols(andr, ["Gene", "mean_logFC_bind", "mean_logFC_expr"], "Andrusiak All_genes")

# rename for clarity
fong = fong.rename(columns={
    "mean_logFC_expr": "expr_fong",
    "adjP_expr": "adjP_fong",
    "mean_logFC_bind": "bind_logFC"
})

andr = andr.rename(columns={
    "mean_logFC_expr": "expr_andr",
    "mean_logFC_bind": "bind_logFC"
})

# ---------- per-family comparisons ----------
summary_rows = []
for fam in FAMILIES:
    genes = build_geneset(go_full, fam)

    fong_sub = fong[fong["Gene"].isin(genes)].copy()
    andr_sub = andr[andr["Gene"].isin(genes)].copy()

    fset = set(fong_sub["Gene"])
    aset = set(andr_sub["Gene"])

    overlap = sorted(fset & aset)
    fong_only = sorted(fset - aset)
    andr_only = sorted(aset - fset)

    tag = fam.replace(" ", "_")

    pd.DataFrame({"Gene": fong_only}).to_csv(
        os.path.join(OUT_DIR, f"fong_only_{tag}.csv"),
        index=False
    )
    pd.DataFrame({"Gene": andr_only}).to_csv(
        os.path.join(OUT_DIR, f"andrusiak_only_{tag}.csv"),
        index=False
    )

    ovl = (
        fong_sub.set_index("Gene")
                .loc[overlap, ["bind_logFC", "expr_fong"]]
                .join(
                    andr_sub.set_index("Gene")
                            .loc[overlap, ["expr_andr"]]
                )
                .reset_index()
    )

    if not ovl.empty:
        ovl["sign_concordant"] = np.sign(ovl["expr_fong"]) == np.sign(ovl["expr_andr"])
        ovl["delta_abs_expr"] = (ovl["expr_fong"].abs() - ovl["expr_andr"].abs()).abs()

        ovl = ovl[[
            "Gene",
            "bind_logFC",
            "expr_fong",
            "expr_andr",
            "sign_concordant",
            "delta_abs_expr"
        ]]

        ovl.to_csv(os.path.join(OUT_DIR, f"overlap_{tag}.csv"), index=False)

        # scatter
        fig, ax = plt.subplots(figsize=(5.2, 4.0))
        ax.scatter(ovl["expr_fong"], ovl["expr_andr"], s=20)

        if len(ovl) >= 2:
            try:
                m, b = np.polyfit(ovl["expr_fong"], ovl["expr_andr"], 1)
                xs = np.linspace(ovl["expr_fong"].min(), ovl["expr_fong"].max(), 100)
                ax.plot(xs, m * xs + b, linewidth=1)
            except Exception:
                pass

        r, p, rs, ps = safe_stats(ovl["expr_fong"], ovl["expr_andr"])
        txt = f"n={len(ovl)}\nPearson r={r:.3f} (p={p:.2e})\nSpearman ρ={rs:.3f} (p={ps:.2e})"
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=9)
        ax.set_title(f"Fong vs Andrusiak expression • {fam}", pad=8)
        ax.set_xlabel("Fong log2FC (Rb-TKO)")
        ax.set_ylabel("Andrusiak log2FC")
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
        "fong_only": len(fong_only),
        "andrusiak_only": len(andr_only),
        "overlap": len(overlap),
        "overlap_concordant_signs": concord,
        "mean_delta_abs_expr_overlap": mean_delta,
        "pearson_r_overlap": r,
        "pearson_p_overlap": p,
        "spearman_rho_overlap": rs,
        "spearman_p_overlap": ps
    })

summary = pd.DataFrame(summary_rows)
summary_path = os.path.join(OUT_DIR, "summary_per_family.csv")
summary.to_csv(summary_path, index=False)

print("✓ Wrote:", summary_path)
for fam in FAMILIES:
    tag = fam.replace(" ", "_")
    print("  - overlap:", os.path.join(OUT_DIR, f"overlap_{tag}.csv"))
    print("  - fong_only:", os.path.join(OUT_DIR, f"fong_only_{tag}.csv"))
    print("  - andrusiak_only:", os.path.join(OUT_DIR, f"andrusiak_only_{tag}.csv"))
    print("  - scatter:", os.path.join(OUT_DIR, f"scatter_{tag}.png"))