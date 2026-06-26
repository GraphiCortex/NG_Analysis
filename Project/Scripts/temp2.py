import pandas as pd

# Point to one of the integrated files, e.g. the mitotic one:
df = pd.read_csv("Project/Data/Processed/expression_integration/mitotic_expression.csv")
print(df.columns)
print(df.head())
