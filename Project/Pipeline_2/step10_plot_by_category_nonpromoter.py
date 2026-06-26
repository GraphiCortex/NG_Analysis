#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

IN_DIR  = os.path.join("Project", "New_Data_2")
FIG_DIR = os.path.join(IN_DIR, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

merged_xlsx = os.path.join(IN_DIR, "binding_vs_TKO_expression_nonpromoter.xlsx")
go_xlsx     = os.path.join(IN_DIR, "GO_results_nonpromoter.xlsx")

# --- 1) Load merged binding+expression (All_genes sheet) ---
merged = pd.read_excel(merged_xlsx, sheet_name="All_genes")

# --- 2) Build gene set from GO workbook (Sh_2) ---
go_full = pd.read_excel(go_xlsx, sheet_name="Sh_2")

def build_geneset(category: str) -> set[str]:
    df = go_full[go_full["Category"] == category]
    if df.empty:
        return set()
    lists = df["geneID"].dropna().astype(str).str.split("/", expand=False)
    return {g.strip().upper() for sub in lists for g in sub}

# --- 3) Categories to plot (non-promoter) ---
cats = {
    "Cell cycle": os.path.join(FIG_DIR, "Scatter2_Cell_cycle_nonpromoter.png"),
    "DNA repair": os.path.join(FIG_DIR, "Scatter2_DNA_repair_nonpromoter.png"),
    "Apoptosis":  os.path.join(FIG_DIR, "Scatter2_Apoptosis_nonpromoter.png"),
}

# --- 4) Plot loop ---
for cat, out_png in cats.items():
    geneset = build_geneset(cat)
    if not geneset:
        print(f"Skipping {cat}: no genes in GO workbook.")
        continue

    sub = merged[merged["Gene"].isin(geneset)]
    x, y = sub["mean_logFC_bind"], sub["mean_logFC_expr"]
    n = len(sub)

    # Fixed small size so all figures align in the paper
    fig, ax = plt.subplots(figsize=(7, 3.2))
    ax.scatter(x, y, s=12, alpha=0.85)

    label_txt = f"n = {n}\n(n < 2)"
    if n >= 2:
        try:
            m, b = np.polyfit(x, y, 1)
            xs = np.linspace(x.min(), x.max(), 100)
            ax.plot(xs, m*xs + b, linewidth=1)
        except Exception:
            pass
        r, p = pearsonr(x, y)
        label_txt = f"n = {n}\nr = {r:.2f}\np = {p:.2e}"

    ax.text(0.03, 0.97, label_txt, transform=ax.transAxes, va="top", fontsize=9)
    ax.set_title(f"ChIP Binding vs. Expression: {cat} (non-promoter)", pad=8, fontsize=11)
    ax.set_xlabel("Binding log₂FC (non-promoter)", fontsize=10)
    ax.set_ylabel("Expression log₂FC (Rb-TKO)", fontsize=10)
    ax.xaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    plt.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"→ Wrote {out_png}")
