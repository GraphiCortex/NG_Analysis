import pandas as pd
import os

# 1) Load your expression data
expr = pd.read_csv(
    os.path.join("Project", "Data", "Raw", "microarray", "Fong2022",
                 "Fong 2022 L2FC1 microarray data.csv"),
    usecols=["Row.names", "log2FoldChange", "padj"]
).rename(columns={
    "Row.names": "Gene",
    "log2FoldChange": "LFC",
    "padj": "adjP"
})

# 2) Load the three gene-set CSVs
base = os.path.join("Project", "Data", "Processed", "gene_sets")
mitotic    = pd.read_csv(os.path.join(base, "E2F4_mitotic_genes.csv"))
dna_repair = pd.read_csv(os.path.join(base, "E2F4_DNArepair_genes.csv"))
apoptosis  = pd.read_csv(os.path.join(base, "E2F4_apoptosis_genes.csv"))

# 3) Merge & flag significant expression
def merge_and_flag(df_genes, category):
    df = df_genes.merge(expr, on="Gene", how="left")
    df["sig"] = (df["LFC"].abs() >= 0.5) & (df["adjP"] < 0.10)
    df["category"] = category
    return df

mitotic_df   = merge_and_flag(mitotic,   "Cell-cycle")
repair_df    = merge_and_flag(dna_repair,"DNA-repair")
apoptosis_df = merge_and_flag(apoptosis, "Apoptosis")

# 4) Save outputs
outdir = os.path.join("Project", "Data", "Processed", "expression_integration")
os.makedirs(outdir, exist_ok=True)

mitotic_df   .to_csv(os.path.join(outdir, "mitotic_expression.csv"),   index=False)
repair_df    .to_csv(os.path.join(outdir, "DNArepair_expression.csv"), index=False)
apoptosis_df .to_csv(os.path.join(outdir, "apoptosis_expression.csv"), index=False)
pd.concat([mitotic_df, repair_df, apoptosis_df], ignore_index=True) \
  .to_csv(os.path.join(outdir, "all_three_categories.csv"), index=False)

print("Merged and flagged expression results saved to", outdir)
