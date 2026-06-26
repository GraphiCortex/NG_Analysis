import pandas as pd
# replace path with your full (unfiltered) promoter probe→gene map
m = pd.read_csv("Project/New_Data/E2F4_vs_E2F3_significant_genes.csv")
m["GENE_SYMBOL"] = m["GENE_SYMBOL"].str.upper().str.strip()
g = ["BAX","BAK1","BID","CASP9","CASP3","CASP7","CASP8","PARP1","PARP2","PARP3","PARP4","BCL2","BCL2L1","BIRC5","PIDD1"]
print(m[m["GENE_SYMBOL"].isin(g)].groupby("GENE_SYMBOL").size())
