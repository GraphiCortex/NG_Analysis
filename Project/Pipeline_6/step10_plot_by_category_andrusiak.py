#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

IN_MERGED = os.path.join("Project", "New_Data_6", "binding_vs_Andrusiak_expression.xlsx")
IN_GOXLS  = os.path.join("Project", "New_Data", "GO_results.xlsx")
OUT_DIR   = os.path.join("Project", "New_Data_6", "figures", "andrusiak")
os.makedirs(OUT_DIR, exist_ok=True)

FIGSIZE = (5.5, 4.0)

# 1) Load merged binding+expression table
merged = pd.read_excel(IN_MERGED, sheet_name="All_genes")
merged["Gene"] = merged["Gene"].astype(str).str.upper().str.strip()

# 2) Load GO results (Sh_2) to build gene sets
go_full = pd.read_excel(IN_GOXLS, sheet_name="Sh_2")

def build_geneset(category: str) -> set[str]:
    df = go_full[go_full["Category"] == category]
    if df.empty:
        return set()
    lists = df["geneID"].dropna().astype(str).str.split("/", expand=False)
    return {g.strip().upper() for sub in lists for g in sub}

cats = {
    "Cell cycle":  os.path.join(OUT_DIR, "Scatter_Cell_cycle_andrusiak.png"),
    "DNA repair":  os.path.join(OUT_DIR, "Scatter_DNA_repair_andrusiak.png"),
    "Apoptosis":   os.path.join(OUT_DIR, "Scatter_Apoptosis_andrusiak.png"),
    "Autophagy":   os.path.join(OUT_DIR, "Scatter_Autophagy_andrusiak.png"),
    "Necroptosis": os.path.join(OUT_DIR, "Scatter_Necroptosis_andrusiak.png"),
}

# 3) Shared axis limits across all categories
xmins, xmaxs, ymins, ymaxs, subs = [], [], [], [], {}
for cat in cats:
    gs = build_geneset(cat)
    sub = merged[merged["Gene"].isin(gs)].copy()
    subs[cat] = sub
    if not sub.empty:
        xmins.append(sub["mean_logFC_bind"].min())
        xmaxs.append(sub["mean_logFC_bind"].max())
        ymins.append(sub["mean_logFC_expr"].min())
        ymaxs.append(sub["mean_logFC_expr"].max())

if xmins:
    x_lo, x_hi = min(xmins), max(xmaxs)
    y_lo, y_hi = min(ymins), max(ymaxs)
    x_pad = 0.05 * max(1e-6, abs(x_hi - x_lo))
    y_pad = 0.05 * max(1e-6, abs(y_hi - y_lo))
    X_LIM = (x_lo - x_pad, x_hi + x_pad)
    Y_LIM = (y_lo - y_pad, y_hi + y_pad)
else:
    X_LIM = (-1, 1)
    Y_LIM = (-1, 1)

# 4) Plot each category
for cat, out_png in cats.items():
    sub = subs[cat]
    x = sub["mean_logFC_bind"] if not sub.empty else pd.Series(dtype=float)
    y = sub["mean_logFC_expr"] if not sub.empty else pd.Series(dtype=float)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.scatter(x, y, s=18)

    n = len(sub)
    if n >= 2:
        m, b = np.polyfit(x, y, 1)
        xs = np.linspace(X_LIM[0], X_LIM[1], 100)
        ax.plot(xs, m * xs + b, linewidth=1)
        r, p = pearsonr(x, y)
        txt = f"n = {n}\nr = {r:.2f}\np = {p:.2e}"
    elif n == 1:
        txt = "n = 1\n(no correlation)"
    else:
        txt = "n = 0"

    ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=9)
    ax.set_title(f"Binding vs. Expression (Andrusiak): {cat}", pad=8, fontsize=11)
    ax.set_xlabel("Binding log$_2$FC", fontsize=10)
    ax.set_ylabel("Expression log$_2$FC", fontsize=10)
    ax.xaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    ax.yaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    ax.set_xlim(*X_LIM)
    ax.set_ylim(*Y_LIM)

    plt.tight_layout()
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print(f"→ Wrote {out_png}")