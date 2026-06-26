#!/usr/bin/env python3
import pandas as pd, matplotlib.pyplot as plt, textwrap, os, seaborn as sns

sns.set_theme(context="talk", style="whitegrid")
ora_files = {
  "Apoptosis": "E2F4_cell_death_ORA.csv",
  "Cell-cycle": "E2F4_cell_cycle_ORA.csv",
  "DNA-repair": "E2F4_DNA_repair_ORA.csv"
}

for category, fname in ora_files.items():
    df = pd.read_csv(os.path.join("Project","Data","Processed","limma_results", fname))
    df["WrappedDesc"] = df["Description"].apply(lambda s: "\n".join(textwrap.wrap(s, width=30)))
    df = df.sort_values("Count")
    fig, ax = plt.subplots(figsize=(14,13))
    plt.subplots_adjust(left=0.4, right=0.95, top=0.9)
    sns.barplot(x="Count", y="WrappedDesc", data=df, palette="deep", ax=ax)
    maxc = df["Count"].max()
    for i,v in enumerate(df["Count"]):
        ax.text(v + maxc*0.01, i, str(v), va="center", fontsize=12)
    ax.set_title(f"{category}-related GO terms enriched among E2F4 targets", fontsize=18, pad=15)
    ax.set_xlabel("Number of E2F4 target promoters", fontsize=14)
    ax.set_ylabel("")
    ax.tick_params(axis="both", labelsize=12)
    ax.invert_yaxis()
    outdir = os.path.join("Project","Data","Processed","figures")
    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f"{category.lower().replace('-','_')}_go_bar_clean.png")
    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"✅ {category} plot saved to {outfile}")
