#!/usr/bin/env python3
# save as tools/compare_apoptosis_rank_lists.py
import pandas as pd

old = pd.read_csv("Project/New_Data/apoptosis_candidates_ranked.csv")
new = pd.read_csv("Project/New_Data/new_data/apoptosis_candidates_ranked.csv")

# Normalize symbols just in case
for df in (old,new):
    df["Gene"] = df["Gene"].astype(str).str.upper().str.strip()

# Minimal columns to keep
keep = ["Gene","mean_logFC_expr","adjP_expr","mean_logFC_bind","P_Value_bind"]
old = old[keep].drop_duplicates(subset=["Gene"])
new = new[keep].drop_duplicates(subset=["Gene"])

old["old_rank"] = old["mean_logFC_expr"].abs().rank(ascending=False, method="first")
old["old_tie2"] = old["mean_logFC_bind"].abs().rank(ascending=False, method="first")
new["new_rank"] = new["mean_logFC_expr"].abs().rank(ascending=False, method="first")
new["new_tie2"] = new["mean_logFC_bind"].abs().rank(ascending=False, method="first")

# Set diffs
old_only = sorted(set(old.Gene) - set(new.Gene))
new_only = sorted(set(new.Gene) - set(old.Gene))
both     = sorted(set(old.Gene) & set(new.Gene))

print(f"Old count: {len(old)} | New count: {len(new)}")
print(f"Shared: {len(both)} | Old-only: {len(old_only)} | New-only: {len(new_only)}")

pd.Series(old_only, name="old_only").to_csv("Project/New_Data/diag_old_only.csv", index=False)
pd.Series(new_only, name="new_only").to_csv("Project/New_Data/diag_new_only.csv", index=False)

# Rank deltas for shared
merged = (old.set_index("Gene")[["old_rank","old_tie2","mean_logFC_expr","mean_logFC_bind"]]
            .join(new.set_index("Gene")[["new_rank","new_tie2","mean_logFC_expr","mean_logFC_bind"]],
                  how="inner", lsuffix="_old", rsuffix="_new")
            .reset_index())
merged["rank_delta"] = merged["new_rank"] - merged["old_rank"]
merged.sort_values("rank_delta").to_csv("Project/New_Data/diag_rank_deltas.csv", index=False)
print("→ Wrote diagnostics under Project/New_Data/: diag_old_only.csv, diag_new_only.csv, diag_rank_deltas.csv")
