#!/usr/bin/env python3
"""
Prepare the legacy filtered Andrusiak dataset for downstream integration.

Input:
- Project/Andrusiak Rbf Corticals 1.5 Fold FDR5 RV-NG microarray data.xlsx

Output:
- Project/New_Data_6/Andrusiak_filtered_full.csv

Notes
-----
This is NOT a fresh limma/TREAT re-analysis.
The workbook is already a legacy filtered hit list (paper-level 1.5-fold / FDR5% style output),
so this step only reformats it into a gene-level table usable by the Julien-integration pipeline.
"""

import os
import math
import numpy as np
import pandas as pd
import glob

candidates = glob.glob(os.path.join("Project", "Andrusiak*.xlsx"))
if not candidates:
    raise FileNotFoundError("Could not find any Andrusiak .xlsx file under Project/")
IN_XLSX = candidates[0]
print(f"Using input file: {IN_XLSX}")
OUT_DIR = os.path.join("Project", "New_Data_6")
OUT_CSV = os.path.join(OUT_DIR, "Andrusiak_filtered_full.csv")
os.makedirs(OUT_DIR, exist_ok=True)

# ---------- load ----------
# Default to first sheet
raw = pd.read_excel(IN_XLSX, sheet_name=0, header=None)

# The file has a two-row header structure.
# Row 0 = broad labels / annotation names
# Row 1 = sample labels (Control / Rb f/f +Cre) or blanks
#
# We'll build robust column names:
# - keep annotation names from row 0
# - for sample columns, append row 1 label
row0 = raw.iloc[0].astype(str).fillna("")
row1 = raw.iloc[1].astype(str).fillna("")

cols = []
seen = {}
for a, b in zip(row0, row1):
    a = a.strip()
    b = b.strip()

    if a.lower() == "nan":
        a = ""
    if b.lower() == "nan":
        b = ""

    if a and b:
        name = f"{a}__{b}"
    elif a:
        name = a
    elif b:
        name = b
    else:
        name = "unnamed"

    # make unique
    seen[name] = seen.get(name, 0) + 1
    if seen[name] > 1:
        name = f"{name}_{seen[name]}"
    cols.append(name)

df = raw.iloc[2:].copy()
df.columns = cols
df = df.reset_index(drop=True)

# ---------- identify key columns ----------
# expected annotation columns
gene_col = None
fc_col = None
p_col = None
padj_col = None

for c in df.columns:
    cl = c.lower()
    if gene_col is None and "gene_symbol" in cl:
        gene_col = c
    if fc_col is None and "fold change" in cl:
        fc_col = c
    if p_col is None and ("observed score" in cl or cl == "observed score(d)"):
        p_col = c
    if padj_col is None and ("expected score" in cl or "dexp" in cl):
        padj_col = c

if gene_col is None:
    raise ValueError(f"Could not find gene symbol column. Columns were:\n{list(df.columns)}")
if fc_col is None:
    raise ValueError(f"Could not find unlogged fold-change column. Columns were:\n{list(df.columns)}")

# sample columns: anything whose second header included Control or Cre
sample_cols = [
    c for c in df.columns
    if ("control" in c.lower()) or ("cre" in c.lower())
]

# ---------- clean core fields ----------
df[gene_col] = df[gene_col].astype(str).str.strip()
df = df[df[gene_col].notna()]
df = df[df[gene_col] != ""]
df = df[df[gene_col].str.lower() != "nan"]

# numeric conversion
df[fc_col] = pd.to_numeric(df[fc_col], errors="coerce")

for c in sample_cols:
    df[c] = pd.to_numeric(df[c], errors="coerce")

if p_col is not None:
    df[p_col] = pd.to_numeric(df[p_col], errors="coerce")
if padj_col is not None:
    df[padj_col] = pd.to_numeric(df[padj_col], errors="coerce")

# ---------- convert unlogged fold change -> signed log2FC ----------
# We assume the column is ratio-like:
#   >1  means up in Cre
#   <1  means down in Cre
#
# Then signed log2FC is:
#   log2(fc) if fc >= 1
#   -log2(1/fc) if 0 < fc < 1
# which is simply log2(fc) for positive fc.
def safe_log2_fc(x):
    if pd.isna(x) or x <= 0:
        return np.nan
    return math.log2(x)

df["log2FoldChange"] = df[fc_col].apply(safe_log2_fc)

# baseMean from available sample columns
if sample_cols:
    df["baseMean"] = df[sample_cols].mean(axis=1, skipna=True)
else:
    df["baseMean"] = np.nan

# ---------- build provisional pvalue / padj ----------
# IMPORTANT:
# The legacy workbook does not provide standard p-value / adjusted p-value
# columns in the same format as our limma/TREAT outputs.
# The detected "Observed score(d)" and "Expected score (dExp)" are SAM-style
# scores, not p-values/FDR, so we do NOT treat them as pvalue/padj.
df["pvalue"] = np.nan
df["padj"] = np.nan

# ---------- collapse to gene level ----------
# Multiple probes per gene can exist; keep the row with largest |log2FC|
# and use lowest available p/padj as summaries.
def first_valid(series):
    s = series.dropna()
    return s.iloc[0] if len(s) else np.nan

df["Gene"] = df[gene_col].astype(str).str.upper().str.strip()
df["abs_log2FC"] = df["log2FoldChange"].abs()

# representative row per gene = largest |log2FC|
rep_idx = df.groupby("Gene")["abs_log2FC"].idxmax()
rep = df.loc[rep_idx, ["Gene", "log2FoldChange", "baseMean"]].copy()

# summary stats per gene
summ = (
    df.groupby("Gene", as_index=False)
      .agg(
          pvalue=("pvalue", "min"),
          padj=("padj", "min")
      )
)

out = rep.merge(summ, on="Gene", how="left")

# reorder
out = out[["Gene", "log2FoldChange", "pvalue", "padj", "baseMean"]]
out = out.sort_values("log2FoldChange", ascending=False).reset_index(drop=True)

# ---------- write ----------
out.to_csv(OUT_CSV, index=False)

# ---------- report ----------
print(f"→ Wrote {OUT_CSV}")
print(f"Rows written: {len(out)}")
print(f"Genes with |log2FC| >= 0.25: {(out['log2FoldChange'].abs() >= 0.25).sum()}")
print(f"Genes with non-null pvalue: {out['pvalue'].notna().sum()}")
print(f"Genes with non-null padj:   {out['padj'].notna().sum()}")

print("\nTop 10 rows preview:")
print(out.head(10).to_string(index=False))