#!/usr/bin/env python3
import pandas as pd
import numpy as np
import glob
import os
from scipy.stats import pearsonr
import matplotlib.pyplot as plt

def quantile_norm(df: pd.DataFrame) -> pd.DataFrame:
    """
    Perform quantile normalization across the columns of df.
    Each column is sorted, averaged across columns, then
    values are reassigned to original ranks.
    """
    # 1. Extract the raw matrix
    mat = df.values.astype(float)
    # 2. Sort each column
    sorted_mat = np.sort(mat, axis=0)
    # 3. Compute the mean of each row across all columns
    mean_sorted = np.mean(sorted_mat, axis=1)
    # 4. Get the rank (1..n) of each value within its column
    ranks = df.rank(method='min').astype(int).values
    # 5. Build the normalized matrix by mapping rank -> mean_sorted
    #    subtract 1 because ranks are 1-based
    mat_qn = mean_sorted[ranks - 1]
    # 6. Return as DataFrame with original index/columns
    return pd.DataFrame(mat_qn, index=df.index, columns=df.columns)

# — adjust this path if your processed CSVs live elsewhere —
processed_dir = 'Project/Data/Processed'

# 1. Find all your *_Mvalues.csv files
csv_paths = sorted(glob.glob(os.path.join(processed_dir, '*_Mvalues.csv')))
if not csv_paths:
    raise FileNotFoundError(f"No files matching '*_Mvalues.csv' in {processed_dir}")

# 2. Load each into a DataFrame, renaming "M" → sample name
dfs = []
for path in csv_paths:
    sample = os.path.basename(path).replace('.csv','')
    df = pd.read_csv(path, usecols=['ID_REF','M'])
    df = df.rename(columns={'M': sample})
    dfs.append(df)

# 3. Merge all DataFrames on ID_REF
merged = dfs[0]
for df in dfs[1:]:
    merged = merged.merge(df, on='ID_REF', how='inner')
merged = merged.set_index('ID_REF')
print(f"\nLoaded M-values matrix: {merged.shape[0]} probes × {merged.shape[1]} arrays")

# 4. Pearson correlations BEFORE normalization
corr_before = merged.corr()
print("\nPearson correlations BEFORE quantile-norm:")
print(corr_before)

# 5. Quantile-normalize
merged_qn = quantile_norm(merged)

# 6. Pearson correlations AFTER normalization
corr_after = merged_qn.corr()
print("\nPearson correlations AFTER quantile-norm:")
print(corr_after)

# 7. Scatterplot of the first two arrays, before & after
if merged.shape[1] >= 2:
    c1, c2 = merged.columns[:2]
    # Before
    plt.figure(figsize=(5,5))
    plt.scatter(merged[c1], merged[c2], alpha=0.4)
    plt.plot([-3,3],[-3,3],'--',color='grey')
    plt.title(f'Before QN: r={corr_before.loc[c1,c2]:.2f}')
    plt.xlabel(c1); plt.ylabel(c2)
    plt.xlim(-3,3); plt.ylim(-3,3)
    plt.savefig(os.path.join(processed_dir, 'QC_scatter_before.png'))
    plt.close()
    # After
    plt.figure(figsize=(5,5))
    plt.scatter(merged_qn[c1], merged_qn[c2], alpha=0.4)
    plt.plot([-3,3],[-3,3],'--',color='grey')
    plt.title(f'After QN: r={corr_after.loc[c1,c2]:.2f}')
    plt.xlabel(c1); plt.ylabel(c2)
    plt.xlim(-3,3); plt.ylim(-3,3)
    plt.savefig(os.path.join(processed_dir, 'QC_scatter_after.png'))
    plt.close()

# 8. Save correlation matrices
corr_before.to_csv(os.path.join(processed_dir, 'chipchip_corr_before_QN.csv'))
corr_after.to_csv(os.path.join(processed_dir, 'chipchip_corr_after_QN.csv'))

print("\nQC complete.")
print(f" • Before-QN correlations → {processed_dir}/chipchip_corr_before_QN.csv")
print(f" • After-QN correlations  → {processed_dir}/chipchip_corr_after_QN.csv")
print(f" • Scatter plots          → {processed_dir}/QC_scatter_before.png, QC_scatter_after.png\n")
