#!/usr/bin/env python3
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# 0. Settings
sns.set_theme(context="notebook", style="whitegrid", font_scale=1.1)
BASE = os.path.abspath(os.path.join(__file__, "..", ".."))
MASTER_F = os.path.join(BASE, "Data", "Processed", "integration",
                        "binding_expression_master.csv")

# 1. Load master table
df = pd.read_csv(MASTER_F)

# Drop any rows missing binding or expression
df = df.dropna(subset=["M_avg", "LFC"])

# 2. Overall correlation
r_all, p_all = pearsonr(df["M_avg"], df["LFC"])
print(f"Overall Pearson r = {r_all:.3f}, p = {p_all:.2e}")

# 3. Per-category correlation
cats = df["category"].unique()
corrs = {}
for cat in cats:
    sub = df[df["category"] == cat]
    r, p = pearsonr(sub["M_avg"], sub["LFC"])
    corrs[cat] = (r, p)
    print(f"{cat:<12}  r = {r:.3f}, p = {p:.2e}")

# 4. Faceted scatter with regression lines
g = sns.lmplot(
    data=df,
    x="M_avg", y="LFC",
    col="category",
    hue="category",
    palette="tab10",
    markers=["o","X","s"],
    height=4, aspect=1,
    scatter_kws={"alpha":0.6, "s":50},
    line_kws={"linewidth":2, "linestyle":"--"}
)
g.set_axis_labels("Average ChIP M-value", "Expression log₂-FC")
g.set_titles("{col_name}")

# Annotate each panel with its r and p
for ax, cat in zip(g.axes.flatten(), cats):
    r, p = corrs[cat]
    ax.text(
        0.05, 0.90,
        f"r = {r:.2f}\np = {p:.1e}",
        transform=ax.transAxes,
        ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.8)
    )
    # Label top 3 by |LFC| * M_avg
    subset = df[df["category"] == cat].copy()
    subset["score"] = subset["M_avg"].abs() * subset["LFC"].abs()
    top3 = subset.nlargest(3, "score")
    for _, row in top3.iterrows():
        ax.text(
            row["M_avg"], row["LFC"],
            row["Gene"],
            fontsize=8,
            weight="bold",
            va="bottom", ha="right"
        )

# 5. Overall plot title & save
plt.subplots_adjust(top=0.85)
g.fig.suptitle("ChIP Binding vs. Expression Change, by Category", fontsize=16)

outdir = os.path.join(BASE, "Data", "Processed", "figures")
os.makedirs(outdir, exist_ok=True)
outf = os.path.join(outdir, "binding_vs_expression_by_category.png")
g.savefig(outf, dpi=300)
print(f"✅ Figure saved to {outf}")
