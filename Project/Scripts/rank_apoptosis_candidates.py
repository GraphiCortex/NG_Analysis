import pandas as pd

# 1) Load the merged binding+expression table
master = pd.read_csv("Project/Data/Processed/integration/binding_expression_master.csv")

# 2) Subset to apoptosis category
apop = master[master.category == "Apoptosis"].copy()

# 3) Compute absolute LFC
apop["absLFC"] = apop.LFC.abs()

# 4) Rank by absLFC then by M_avg
apop = apop.sort_values(["absLFC", "M_avg"], ascending=[False, False])

# 5) Show top 5 candidates
print(apop[["Gene","LFC","adjP","M_avg"]].head(5))
