import pandas as pd

# 1. Load expression set
expr = pd.read_csv('Project/Data/Processed/expression_integration/all_three_categories.csv')
expr_genes = set(expr['Gene'])

# 2. Load merged master table
master = pd.read_csv('Project/Data/Processed/integration/binding_expression_master.csv')

# 3. Unique genes in each
expr_n = len(expr_genes)
master_n = master['Gene'].nunique()

# 4. Genes with any binding data (non‐NA M_avg)
binding_genes = set(master.loc[master['M_avg'].notna(), 'Gene'])
binding_n = len(binding_genes & expr_genes)

# 5. Report
print(f"Distinct genes in expression set:       {expr_n}")
print(f"Distinct genes in master merge:         {master_n}")
print(f"Distinct genes with binding (M_avg):    {binding_n}")
print(f"Genes missing binding entirely:         {expr_n - binding_n}")

#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# 1. Load the master table
master = pd.read_csv("Project/Data/Processed/integration/binding_expression_master.csv")

# 2. Scatter plot M_avg vs LFC, colored by significance
sns.set_theme(style="whitegrid", context="talk")
plt.figure(figsize=(7, 6))
ax = sns.scatterplot(
    data=master,
    x="M_avg",
    y="LFC",
    hue="sig",
    style="category",
    palette={True: "tab:red", False: "tab:blue"},
    alpha=0.7,
    edgecolor=None
)

# 3. A regression line (optional)
sns.regplot(
    data=master, x="M_avg", y="LFC",
    scatter=False, ax=ax,
    line_kws={"color":"grey", "linewidth":1, "linestyle":"--"}
)

# 4. Labels & title
ax.set_title("ChIP binding (M_avg) vs. Expression change (LFC)", pad=15)
ax.set_xlabel("Average ChIP M-value")
ax.set_ylabel("log₂ fold‐change (Fong2022)")

# 5. Save
outdir = os.path.join("Project","Data","Processed","figures")
os.makedirs(outdir, exist_ok=True)
outfile = os.path.join(outdir, "binding_vs_expression_scatter.png")
plt.tight_layout()
plt.savefig(outfile, dpi=300)
plt.show()
