from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


# ============================================================
# F10A — Adult NSC residual decomposition
#
# Goal:
#
# Learn the conserved Rb-loss response from genes that:
#
#   1. are strict Fong Rb-TKO DEGs,
#   2. respond significantly and reciprocally to E2F1/3 loss,
#   3. move in the same direction in Fong, CAG and MAP2.
#
# Then estimate:
#
#   adult residual =
#       observed Fong Rb-TKO response
#       -
#       response expected from the conserved backbone.
#
# Large residuals nominate genes/programs whose adult behavior
# is poorly explained by the conserved Rb/E2F backbone.
#
# IMPORTANT:
# Residual != proof of adult-specific causality.
# ============================================================


ROOT = Path(".")

INFILE = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F8B_ContextMatrix"
    / "F8B_fong_strict_complete_contexts.csv"
)

OUT = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F10A_AdultResiduals"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Helpers
# ============================================================

def as_bool(series):

    if series.dtype == bool:
        return series

    return (
        series.astype(str)
        .str.lower()
        .map({
            "true": True,
            "false": False,
        })
        .fillna(False)
    )


# ============================================================
# Load
# ============================================================

df = pd.read_csv(INFILE)


for col in [
    "Adult_E2F_reciprocal_strong",
    "Rb_loss_sign_conserved",
    "CAG_FDR05",
    "MAP2_FDR05",
]:

    if col in df.columns:
        df[col] = as_bool(df[col])


# ============================================================
# Shared cortical Rb-loss response
#
# CAG and MAP2 are on the same reprocessed platform and are
# strongly correlated, so their mean gives a compact shared
# cortical response.
# ============================================================

df["Osh_shared_LFC"] = (
    df[
        [
            "CAG_LFC",
            "MAP2_LFC",
        ]
    ]
    .mean(axis=1)
)


df["Osh_CAG_MAP2_disagreement"] = (
    df["MAP2_LFC"]
    -
    df["CAG_LFC"]
)


# ============================================================
# Define high-confidence conserved backbone training set
# ============================================================

train_mask = (
    df["Adult_E2F_reciprocal_strong"]
    &
    df["Rb_loss_sign_conserved"]
    &
    (
        df["CAG_FDR05"]
        |
        df["MAP2_FDR05"]
    )
)


train = df.loc[
    train_mask
].copy()


print()
print("=" * 105)
print("F10A — ADULT NSC RESIDUAL DECOMPOSITION")
print("=" * 105)

print()
print(
    "Complete-context Fong genes:",
    len(df)
)

print(
    "High-confidence backbone training genes:",
    len(train)
)


if len(train) < 20:

    raise RuntimeError(
        "Backbone training set is unexpectedly small. "
        "Inspect F8B flags before proceeding."
    )


# ============================================================
# Robust regression
#
# Adult Fong response:
#
#   y = intercept + slope * shared Oshikawa response
#
# Theil-Sen is used so that strongly divergent genes do not
# dictate the backbone mapping.
# ============================================================

slope, intercept, slope_lo, slope_hi = (
    stats.theilslopes(
        train["RbTKO_LFC"].astype(float),
        train["Osh_shared_LFC"].astype(float),
        alpha=0.95,
    )
)


print()
print("ROBUST BACKBONE MAPPING")
print("-" * 105)

print(
    f"Fong_expected = "
    f"{intercept:.4f} "
    f"+ {slope:.4f} * Osh_shared"
)

print(
    f"Theil-Sen slope 95% CI: "
    f"[{slope_lo:.4f}, {slope_hi:.4f}]"
)


# ============================================================
# Predict adult response + residuals
# ============================================================

df["RbTKO_expected_from_backbone"] = (
    intercept
    +
    slope
    *
    df["Osh_shared_LFC"]
)


df["Adult_residual"] = (
    df["RbTKO_LFC"]
    -
    df["RbTKO_expected_from_backbone"]
)


# ============================================================
# Robust residual scale from backbone training genes
# ============================================================

train_residual = (
    train["RbTKO_LFC"]
    -
    (
        intercept
        +
        slope
        *
        train["Osh_shared_LFC"]
    )
)


center = np.median(
    train_residual
)

mad = np.median(
    np.abs(
        train_residual
        -
        center
    )
)

robust_sigma = (
    1.4826
    *
    mad
)


if robust_sigma <= 1e-8:
    raise RuntimeError(
        "Residual MAD is zero or nearly zero."
    )


df["Adult_residual_z"] = (
    (
        df["Adult_residual"]
        -
        center
    )
    /
    robust_sigma
)


# ============================================================
# Descriptive residual classes
#
# These are ranking aids, NOT formal biological categories.
# ============================================================

df["Residual_class"] = "Backbone_compatible"

df.loc[
    df["Adult_residual_z"] >= 2,
    "Residual_class",
] = "Adult_amplified"

df.loc[
    df["Adult_residual_z"] <= -2,
    "Residual_class",
] = "Adult_attenuated"


# Actual sign reversal relative to shared cortical response

df["Adult_vs_Osh_sign_reversal"] = (
    np.sign(
        df["RbTKO_LFC"]
    )
    !=
    np.sign(
        df["Osh_shared_LFC"]
    )
)


# Require a nontrivial Oshikawa response so tiny values do not
# become meaningless "reversals".

df["Adult_vs_Osh_strong_reversal"] = (
    df["Adult_vs_Osh_sign_reversal"]
    &
    (
        df["Osh_shared_LFC"].abs()
        >= 0.5
    )
)


# ============================================================
# E2F1/3 relationship
# ============================================================

def e2f_class(row):

    if pd.isna(
        row["E2F13_LFC"]
    ):
        return "Untestable"

    if (
        pd.isna(
            row["E2F13_padj"]
        )
        or
        row["E2F13_padj"] >= 0.05
    ):
        return "E2F13_not_significant"

    if (
        np.sign(
            row["RbTKO_LFC"]
        )
        ==
        -np.sign(
            row["E2F13_LFC"]
        )
    ):
        return "E2F_reciprocal"

    if (
        np.sign(
            row["RbTKO_LFC"]
        )
        ==
        np.sign(
            row["E2F13_LFC"]
        )
    ):
        return "E2F_same_direction"

    return "Other"


df["E2F_relationship"] = (
    df.apply(
        e2f_class,
        axis=1,
    )
)


# ============================================================
# Rank divergence
# ============================================================

df["Abs_adult_residual_z"] = (
    df["Adult_residual_z"]
    .abs()
)


ranked = (
    df.sort_values(
        "Abs_adult_residual_z",
        ascending=False,
    )
    .copy()
)


ranked.to_csv(
    OUT
    / "F10A_all_adult_residuals_ranked.csv",
    index=False,
)


# ============================================================
# Strong residuals
# ============================================================

strong = ranked[
    ranked[
        "Abs_adult_residual_z"
    ]
    >= 2
].copy()


strong.to_csv(
    OUT
    / "F10A_strong_adult_divergence.csv",
    index=False,
)


reversals = ranked[
    ranked[
        "Adult_vs_Osh_strong_reversal"
    ]
].copy()


reversals.to_csv(
    OUT
    / "F10A_adult_vs_cortical_reversals.csv",
    index=False,
)


# ============================================================
# Candidate genes from F6/F7
# ============================================================

CANDIDATES = [
    "E2F1",
    "E2F2",
    "MYC",
    "DDIT3",

    "PCNA",
    "CLSPN",
    "TOPBP1",
    "CHEK1",
    "CHEK2",
    "H2AFX",
    "BRCA1",
    "CDK1",
    "CDK2",

    "FAS",
    "TNFRSF10B",

    "PMAIP1",
    "BBC3",

    "ZBP1",
    "PYCARD",
    "MLKL",
    "EIF2AK2",

    "BDNF",
    "KDR",
    "FGF2",
    "IRS2",
    "ANGPT1",
    "MET",
    "CXCL12",

    "AURKB",
    "BIRC5",
    "BUB1",
    "MELK",
    "KIF14",
    "MYBL2",

    "NR4A1",
    "NR4A2",
    "NR4A3",
]


candidate_table = ranked[
    ranked["Gene"].isin(
        CANDIDATES
    )
].copy()


candidate_table.to_csv(
    OUT
    / "F10A_candidate_residuals.csv",
    index=False,
)


# ============================================================
# Backbone fit statistics
# ============================================================

train_prediction = (
    intercept
    +
    slope
    *
    train["Osh_shared_LFC"]
)


pearson = stats.pearsonr(
    train["Osh_shared_LFC"],
    train["RbTKO_LFC"],
).statistic

spearman = stats.spearmanr(
    train["Osh_shared_LFC"],
    train["RbTKO_LFC"],
).statistic


# ============================================================
# Console
# ============================================================

print()
print(
    f"Training Pearson r  = {pearson:.4f}"
)

print(
    f"Training Spearman r = {spearman:.4f}"
)

print(
    f"Robust residual sigma = {robust_sigma:.4f}"
)


print()
print("RESIDUAL COUNTS")
print("-" * 105)

print(
    df[
        "Residual_class"
    ]
    .value_counts()
    .to_string()
)

print(
    "\nStrong adult-vs-cortical sign reversals:",
    int(
        df[
            "Adult_vs_Osh_strong_reversal"
        ].sum()
    )
)


# ============================================================
# Top residuals
# ============================================================

show = [
    "Gene",
    "RbTKO_LFC",
    "E2F13_LFC",
    "E2F13_padj",
    "CAG_LFC",
    "MAP2_LFC",
    "Osh_shared_LFC",
    "RbTKO_expected_from_backbone",
    "Adult_residual",
    "Adult_residual_z",
    "Residual_class",
    "E2F_relationship",
]


print()
print("TOP 40 ADULT-DIVERGENT GENES")
print("-" * 105)

print(
    ranked[
        show
    ]
    .head(40)
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.3g}",
    )
)


print()
print("KEY CANDIDATES")
print("-" * 105)

print(
    candidate_table
    .sort_values(
        "Abs_adult_residual_z",
        ascending=False,
    )[
        show
    ]
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.3g}",
    )
)


# ============================================================
# Summary
# ============================================================

with open(
    OUT / "F10A_summary.txt",
    "w",
    encoding="utf-8",
) as f:

    f.write(
        "F10A — ADULT NSC RESIDUAL DECOMPOSITION\n"
    )

    f.write(
        "======================================\n\n"
    )

    f.write(
        f"Complete-context genes: {len(df)}\n"
    )

    f.write(
        f"Backbone training genes: {len(train)}\n"
    )

    f.write(
        f"Theil-Sen intercept: {intercept:.6f}\n"
    )

    f.write(
        f"Theil-Sen slope: {slope:.6f}\n"
    )

    f.write(
        f"Training Pearson: {pearson:.6f}\n"
    )

    f.write(
        f"Training Spearman: {spearman:.6f}\n"
    )

    f.write(
        f"Residual robust sigma: {robust_sigma:.6f}\n"
    )

    f.write(
        f"Strong residual genes |z|>=2: {len(strong)}\n"
    )

    f.write(
        "Strong adult/cortical sign reversals: "
        f"{int(df['Adult_vs_Osh_strong_reversal'].sum())}\n\n"
    )

    f.write(
        "Interpretation:\n"
    )

    f.write(
        "- Residuals quantify divergence from the conserved "
        "Rb-loss transcriptional backbone.\n"
    )

    f.write(
        "- They do not by themselves prove adult-specific "
        "causal regulation.\n"
    )

    f.write(
        "- High residual genes are candidates for the adult "
        "state/fate-routing layer to be tested in F10B/F11.\n"
    )


print()
print("Outputs written to:")
print(OUT)

print()
print(
    "[NEXT] Test whether adult-divergent residuals concentrate "
    "in coherent death/survival programs and causal network "
    "positions rather than appearing as isolated genes."
)