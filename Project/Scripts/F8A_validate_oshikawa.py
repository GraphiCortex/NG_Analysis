from pathlib import Path
import numpy as np
import pandas as pd


# ============================================================
# Paths
# ============================================================

ROOT = Path(".")

OLD = (
    ROOT
    / "Data"
    / "Raw"
    / "microarray"
    / "oshikawa2013"
    / "Oshikawa 2013 L2FC1.csv"
)

NEW_DIR = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F8A_Oshikawa"
    / "DE"
)

OUT_DIR = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F8A_Oshikawa"
    / "Validation"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Helpers
# ============================================================

def clean_symbol(x):
    if pd.isna(x):
        return np.nan

    x = str(x).strip().upper()

    if x in {
        "",
        "NA",
        "NAN",
        "---",
        "NULL",
    }:
        return np.nan

    return x


def collapse_old(df):
    """
    Collapse old probe-level Oshikawa table to gene level.

    For genes represented by multiple probes:
      - expression effect = median logFC
      - adjusted P = minimum adjusted P
      - raw P = minimum raw P
      - probe count retained
    """

    x = df.copy()

    x["Gene"] = x["GENE_SYMBOL"].map(clean_symbol)

    x["logFC"] = pd.to_numeric(
        x["logFC"],
        errors="coerce",
    )

    x["P.Value"] = pd.to_numeric(
        x["P.Value"],
        errors="coerce",
    )

    x["adj.P.Val"] = pd.to_numeric(
        x["adj.P.Val"],
        errors="coerce",
    )

    x = x.dropna(
        subset=[
            "Gene",
            "logFC",
        ]
    )

    collapsed = (
        x.groupby("Gene", as_index=False)
        .agg(
            old_logFC=("logFC", "median"),
            old_P=("P.Value", "min"),
            old_adjP=("adj.P.Val", "min"),
            old_probe_count=("ID", "count"),
        )
    )

    return collapsed


def load_new(path, prefix):
    x = pd.read_csv(path)

    x["Gene"] = (
        x["Gene"]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    x["logFC"] = pd.to_numeric(
        x["logFC"],
        errors="coerce",
    )

    x["P.Value"] = pd.to_numeric(
        x["P.Value"],
        errors="coerce",
    )

    x["adj.P.Val"] = pd.to_numeric(
        x["adj.P.Val"],
        errors="coerce",
    )

    x = x[
        [
            "Gene",
            "logFC",
            "P.Value",
            "adj.P.Val",
        ]
    ].copy()

    x = x.rename(
        columns={
            "logFC": f"{prefix}_logFC",
            "P.Value": f"{prefix}_P",
            "adj.P.Val": f"{prefix}_adjP",
        }
    )

    return x


def compare(old, new, prefix):
    m = old.merge(
        new,
        on="Gene",
        how="inner",
    )

    new_lfc = f"{prefix}_logFC"

    valid = (
        m["old_logFC"].notna()
        &
        m[new_lfc].notna()
    )

    m = m.loc[valid].copy()

    pearson = m[
        ["old_logFC", new_lfc]
    ].corr(
        method="pearson"
    ).iloc[0, 1]

    spearman = m[
        ["old_logFC", new_lfc]
    ].corr(
        method="spearman"
    ).iloc[0, 1]

    sign_agree = (
        np.sign(m["old_logFC"])
        ==
        np.sign(m[new_lfc])
    ).mean()

    mae = (
        m["old_logFC"]
        -
        m[new_lfc]
    ).abs().mean()

    rmse = np.sqrt(
        np.mean(
            (
                m["old_logFC"]
                -
                m[new_lfc]
            ) ** 2
        )
    )

    return {
        "comparison": prefix,
        "genes_compared": len(m),
        "pearson": pearson,
        "spearman": spearman,
        "sign_agreement": sign_agree,
        "mean_abs_error": mae,
        "rmse": rmse,
    }, m


# ============================================================
# Load data
# ============================================================

old_probe = pd.read_csv(OLD)

print("\nOLD OSHIKAWA FILE")
print("=" * 80)

print(
    "Rows:",
    len(old_probe)
)

print(
    "Unique probe IDs:",
    old_probe["ID"].nunique()
)

print(
    "Nonmissing symbols:",
    old_probe["GENE_SYMBOL"].notna().sum()
)

print(
    "Unique nonmissing symbols:",
    old_probe["GENE_SYMBOL"].dropna().nunique()
)


old = collapse_old(old_probe)

print(
    "Collapsed gene-level symbols:",
    len(old)
)


cag = load_new(
    NEW_DIR
    / "CAG_TKO_vs_Control_full.csv",
    "CAG",
)

map2 = load_new(
    NEW_DIR
    / "MAP2_TKO_vs_Control_full.csv",
    "MAP2",
)


# ============================================================
# Compare old file to both possible contrasts
# ============================================================

stats_cag, merged_cag = compare(
    old,
    cag,
    "CAG",
)

stats_map2, merged_map2 = compare(
    old,
    map2,
    "MAP2",
)

stats = pd.DataFrame(
    [
        stats_cag,
        stats_map2,
    ]
)

stats.to_csv(
    OUT_DIR
    / "old_vs_new_comparison_summary.csv",
    index=False,
)


# ============================================================
# Determine which new contrast old file resembles more closely
# ============================================================

winner = (
    stats
    .sort_values(
        [
            "pearson",
            "spearman",
        ],
        ascending=False,
    )
    .iloc[0]
)


# ============================================================
# Export merged comparisons
# ============================================================

merged_cag.to_csv(
    OUT_DIR
    / "old_vs_new_CAG.csv",
    index=False,
)

merged_map2.to_csv(
    OUT_DIR
    / "old_vs_new_MAP2.csv",
    index=False,
)


# ============================================================
# Largest disagreements
# ============================================================

for prefix, merged in [
    ("CAG", merged_cag),
    ("MAP2", merged_map2),
]:

    x = merged.copy()

    x["difference"] = (
        x["old_logFC"]
        -
        x[f"{prefix}_logFC"]
    )

    x["abs_difference"] = (
        x["difference"].abs()
    )

    x["sign_match"] = (
        np.sign(x["old_logFC"])
        ==
        np.sign(x[f"{prefix}_logFC"])
    )

    x.sort_values(
        "abs_difference",
        ascending=False,
    ).head(100).to_csv(
        OUT_DIR
        / f"largest_disagreements_{prefix}.csv",
        index=False,
    )


# ============================================================
# Console output
# ============================================================

print("\n")
print("=" * 80)
print("OLD vs NEW OSHIKAWA VALIDATION")
print("=" * 80)

display = stats.copy()

for col in [
    "pearson",
    "spearman",
    "sign_agreement",
    "mean_abs_error",
    "rmse",
]:
    display[col] = display[col].round(4)

print(
    display.to_string(
        index=False
    )
)

print("\nBest match for historical L2FC1 file:")
print(
    f"    {winner['comparison']}"
)

print(
    f"    Pearson = {winner['pearson']:.4f}"
)

print(
    f"    Spearman = {winner['spearman']:.4f}"
)

print(
    f"    Sign agreement = "
    f"{winner['sign_agreement']:.1%}"
)


print("\nOutputs:")
print(OUT_DIR)

print(
    "\n[NEXT] If the historical file agrees with one of the "
    "fresh contrasts, use the fresh raw-derived results as "
    "the canonical Oshikawa dataset for F8B."
)