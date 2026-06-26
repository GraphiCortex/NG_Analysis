import pandas as pd
import itertools
import matplotlib.pyplot as plt

# Paths to your three M-values CSVs (they include the raw IP/INPUT columns)
paths = {
    'rep1': 'Project/Data/Processed/Julian2016_rep1_Mvalues.csv',
    'rep2': 'Project/Data/Processed/Julian2016_rep2_Mvalues.csv',
    'rep3': 'Project/Data/Processed/Julian2016_rep3_Mvalues.csv',
}

# Load each into a DataFrame of ID_REF, IP, INPUT
dfs = {}
for name, p in paths.items():
    df = pd.read_csv(p, usecols=['ID_REF','IP','INPUT'])
    df = df.rename(columns={'IP':f'IP_{name}', 'INPUT':f'IN_{name}'})
    dfs[name] = df

# Merge all three on ID_REF
merged = dfs['rep1']
for name in ['rep2','rep3']:
    merged = merged.merge(dfs[name], on='ID_REF')

# Compute and print all pairwise correlations
for a, b in itertools.combinations(['rep1','rep2','rep3'], 2):
    ip_r = merged[f'IP_{a}'].corr(merged[f'IP_{b}'])
    in_r = merged[f'IN_{a}'].corr(merged[f'IN_{b}'])
    print(f"{a.upper()} vs {b.upper()}  IP–IP r = {ip_r:.3f},  INPUT–INPUT r = {in_r:.3f}")

# Plot IP scatter for each pair
fig, axes = plt.subplots(1, 3, figsize=(15,5))
for ax,(a,b) in zip(axes, itertools.combinations(['rep1','rep2','rep3'], 2)):
    ax.scatter(merged[f'IP_{a}'], merged[f'IP_{b}'], s=2, alpha=0.3)
    mn = min(merged[f'IP_{a}'].min(), merged[f'IP_{b}'].min())
    mx = max(merged[f'IP_{a}'].max(), merged[f'IP_{b}'].max())
    ax.plot([mn,mx],[mn,mx], '--', color='grey')
    ax.set_xlabel(f'IP_{a}')
    ax.set_ylabel(f'IP_{b}')
    ax.set_title(f'{a.upper()} vs {b.upper()} IP r={merged[f"IP_{a}"].corr(merged[f"IP_{b}"]):.2f}')

plt.tight_layout()
plt.show()
