import pandas as pd

# point this at whichever file you want to check first
fn = "Project/Data/Processed/expression_integration/apoptosis_expression.csv"
df = pd.read_csv(fn)

# 1. count of missing values in LFC & adjP
print("Missing LFC:",    df["LFC"].isna().sum())
print("Missing adjP:",   df["adjP"].isna().sum())

# 2. ensure 'sig' is boolean and count both states
print("\n'sig' dtype:", df["sig"].dtype)
print(df["sig"].value_counts(dropna=False))

# 3. spot-check ranges
print("\nLFC  min/max:",  df["LFC"].min(), "/", df["LFC"].max())
print("adjP min/max:",   df["adjP"].min(), "/", df["adjP"].max())
