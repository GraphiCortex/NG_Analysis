import pandas as pd

# 1) Load just the first few rows to inspect columns
df = pd.read_csv(
    "../Data/Raw/microarray/Fong2022/Fong 2022 L2FC1 microarray data.csv",
    nrows=5
)
print("PATH USED → ../Data/microarray/Fong2022/Fong 2022 L2FC1 microarray data.csv")
print("COLUMNS:", list(df.columns))
