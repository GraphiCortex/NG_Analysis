#!/usr/bin/env python3
# plot_cellcycle_top20.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap
import os

sns.set_theme(context="talk", style="whitegrid")

# 1. Load ORA results
infile = os.path.join(
    "Project", "Data", "Processed", "limma_results",
    "E2F4_cell_cycle_ORA.csv"
)
df = pd.read_csv(infile)

# 2. Take top 20 by Count
topn = 20
df_top = df.nlargest(topn, "Count").copy()

# 3. Wrap long descriptions
wrap_width = 25
df_top["WrappedDesc"] = (
    df_top["Description"]
      .astype(str)
      .apply(lambda s: "\n".join(textwrap.wrap(s, width=wrap_width)))
)

# 4. Sort so biggest bar is at top
df_top = df_top.sort_values("Count")

# 5. Plot
fig, ax = plt.subplots(figsize=(12, 12))
plt.subplots_adjust(left=0.35, right=0.95, top=0.88, bottom=0.05)

sns.barplot(
    x="Count", y="WrappedDesc",
    data=df_top,
    palette="muted",
    ax=ax
)

# 6. Annotate
max_count = df_top["Count"].max()
for i, v in enumerate(df_top["Count"]):
    ax.text(
        v + max_count*0.02,
        i,
        str(v),
        va="center",
        fontsize=12,
        color="black"
    )

# 7. Labels & title
ax.set_title(
    f"Top {topn} cell-cycle GO terms enriched among E2F4 targets",
    fontsize=18, pad=12
)
ax.set_xlabel("Number of E2F4 target promoters", fontsize=14)
ax.set_ylabel("")
ax.tick_params(axis="y", labelsize=12)
ax.tick_params(axis="x", labelsize=12)

ax.invert_yaxis()

# 8. Save
outdir = os.path.join("Project","Data","Processed","figures")
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, f"cellcycle_top{topn}_go_bar.png")
fig.savefig(outfile, dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"✅ Top {topn} cell-cycle GO bar-plot saved to {outfile}")
