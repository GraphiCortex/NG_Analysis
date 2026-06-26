#!/usr/bin/env python3
"""
integrate_binding_expression.py

1) Load the three replicate M‐value tables from Chip–chip (Julian2016).  
2) For each, keep only probe ID (ID_REF) and binding (M), rename to M1/M2/M3.  
3) Join them on probe ID → one table with M1, M2, M3 per probe.  
4) Read the GPL_probe2gene mapping (ID_REF → SYMBOL), drop unmapped probes.  
5) Collapse to gene level by averaging the three M’s.  
6) Load our expression‐integration master (Gene, LFC, adjP, sig, category).  
7) Merge binding + expression on “Gene” → one master table.  
8) Write out to CSV.
"""

import os
import pandas as pd

# ——— 0) Base directories —————————————————————————————
# assume script lives in Project/Scripts/
BASE = os.path.abspath(os.path.join(__file__, os.pardir, os.pardir))
PROC = os.path.join(BASE, "Data", "Processed")
META = os.path.join(BASE, "MetaData")

# ——— 1) Files to read ————————————————————————————————
expr_file   = os.path.join(PROC, "expression_integration", "all_three_categories.csv")
mapping_f   = os.path.join(META, "GPL_probe2gene.csv")
m1_file     = os.path.join(PROC, "Julian2016_rep1_Mvalues.csv")
m2_file     = os.path.join(PROC, "Julian2016_rep2_Mvalues.csv")
m3_file     = os.path.join(PROC, "Julian2016_rep3_Mvalues.csv")

# ——— 2) Helper to load each M‐value file ————————————————————
def load_M(path, colname):
    """Read only ID_REF + M, rename M → colname."""
    df = pd.read_csv(path, usecols=["ID_REF","M"])
    df = df.rename(columns={"M": colname})
    return df

# load all three replicates
m1 = load_M(m1_file, "M_rep1")
m2 = load_M(m2_file, "M_rep2")
m3 = load_M(m3_file, "M_rep3")

# ——— 3) Merge probes across replicates ————————————————————
m = m1.merge(m2, on="ID_REF", how="inner") \
      .merge(m3, on="ID_REF", how="inner")

# ——— 4) Map probes → gene symbols ——————————————————————
# the GPL_probe2gene.csv should have ID_REF and SYMBOL columns
map_df = pd.read_csv(mapping_f, usecols=["ID_REF","SYMBOL"])
map_df = map_df.rename(columns={"SYMBOL":"Gene"})

# join, drop probes without a Gene
m = m.merge(map_df, on="ID_REF", how="left")
m = m.dropna(subset=["Gene"])

# ——— 5) Collapse to gene level by averaging the three M’s ——————
gene_bind = (
    m
    .groupby("Gene")[["M_rep1","M_rep2","M_rep3"]]
    .mean()
    .reset_index()
)
# also add an explicit column for the average across the three
gene_bind["M_avg"] = gene_bind[["M_rep1","M_rep2","M_rep3"]].mean(axis=1)

# ——— 6) Load expression‐integration results ——————————————
expr = pd.read_csv(expr_file)
# expect columns: Gene, LFC, adjP, sig, category

# ——— 7) Merge binding ↔ expression —————————————————————
master = expr.merge(gene_bind, on="Gene", how="left")

# ——— 8) Write out final table ————————————————————————
outdir = os.path.join(PROC, "integration")
os.makedirs(outdir, exist_ok=True)
outf   = os.path.join(outdir, "binding_expression_master.csv")
master.to_csv(outf, index=False)

print(f"✅ Wrote {len(master)} rows to {outf}")
