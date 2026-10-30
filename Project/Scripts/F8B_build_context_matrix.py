from pathlib import Path
import numpy as np
import pandas as pd


# ============================================================
# F8B — Multi-context perturbation response matrix
#
# Integrates:
#   1. Adult Rb-TKO (Fong 2022)
#   2. Adult E2F1/3-DKO (Fong 2022)
#   3. Oshikawa CAG-TKO
#   4. Oshikawa MAP2-TKO
#
# IMPORTANT:
#   Absence from the strict Rb-TKO workbook is encoded as NA,
#   NOT zero.
# ============================================================


ROOT = Path(".")

RAW_FONG = (
    ROOT
    / "Data"
    / "Raw"
    / "microarray"
    / "Fong2022"
)

OSH_DIR = (
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
    / "F8B_ContextMatrix"
)

OUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Helpers
# ============================================================

def clean_gene(x):
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


def sign_class(x, threshold=0):
    if pd.isna(x):
        return "NA"

    if x > threshold:
        return "Up"

    if x < -threshold:
        return "Down"

    return "Near_zero"


# ============================================================
# 1. Fong Rb-TKO
# ============================================================

rb_path = (
    RAW_FONG
    / "Fong 2022 L2FC1 microarray data.xlsx"
)

rb = pd.read_excel(
    rb_path,
    sheet_name="Sig L2FC1",
)

# In this workbook:
#   Unnamed: 2 = gene symbol
#   Row.names  = Ensembl ID

rb = rb.rename(
    columns={
        "Unnamed: 2": "Gene",
        "Row.names": "Ensembl_RbTKO",
        "log2FoldChange": "RbTKO_LFC",
        "pvalue": "RbTKO_P",
        "padj": "RbTKO_padj",
    }
)

rb["Gene"] = rb["Gene"].map(clean_gene)

rb = rb.dropna(
    subset=["Gene", "RbTKO_LFC"]
)

rb["RbTKO_LFC"] = pd.to_numeric(
    rb["RbTKO_LFC"],
    errors="coerce",
)

rb["RbTKO_padj"] = pd.to_numeric(
    rb["RbTKO_padj"],
    errors="coerce",
)

rb = (
    rb.sort_values(
        "RbTKO_padj",
        na_position="last",
    )
    .drop_duplicates(
        "Gene",
        keep="first",
    )
)

rb = rb[
    [
        "Gene",
        "Ensembl_RbTKO",
        "RbTKO_LFC",
        "RbTKO_P",
        "RbTKO_padj",
    ]
].copy()


# ============================================================
# 2. Fong E2F1/3-DKO
# ============================================================

e2f_path = (
    RAW_FONG
    / "Fong 2022 E2F1.3 microarray data.xlsx"
)

e2f = pd.read_excel(
    e2f_path,
    sheet_name="GSE190764_FDSD_FC.IHW",
)

e2f = e2f.rename(
    columns={
        "name": "Gene",
        "Row.names": "Ensembl_E2F13",
        "log2FoldChange": "E2F13_LFC",
        "pvalue": "E2F13_P",
        "padj": "E2F13_padj",
    }
)

e2f["Gene"] = e2f["Gene"].map(clean_gene)

e2f["E2F13_LFC"] = pd.to_numeric(
    e2f["E2F13_LFC"],
    errors="coerce",
)

e2f["E2F13_P"] = pd.to_numeric(
    e2f["E2F13_P"],
    errors="coerce",
)

e2f["E2F13_padj"] = pd.to_numeric(
    e2f["E2F13_padj"],
    errors="coerce",
)

e2f = e2f.dropna(
    subset=[
        "Gene",
        "E2F13_LFC",
    ]
)

e2f = (
    e2f.sort_values(
        "E2F13_padj",
        na_position="last",
    )
    .drop_duplicates(
        "Gene",
        keep="first",
    )
)

e2f = e2f[
    [
        "Gene",
        "Ensembl_E2F13",
        "E2F13_LFC",
        "E2F13_P",
        "E2F13_padj",
    ]
].copy()


# ============================================================
# 3. Oshikawa CAG
# ============================================================

cag = pd.read_csv(
    OSH_DIR
    / "CAG_TKO_vs_Control_full.csv"
)

cag["Gene"] = cag["Gene"].map(clean_gene)

cag = cag.rename(
    columns={
        "logFC": "CAG_LFC",
        "P.Value": "CAG_P",
        "adj.P.Val": "CAG_padj",
    }
)

cag = cag[
    [
        "Gene",
        "CAG_LFC",
        "CAG_P",
        "CAG_padj",
    ]
].copy()


# ============================================================
# 4. Oshikawa MAP2
# ============================================================

map2 = pd.read_csv(
    OSH_DIR
    / "MAP2_TKO_vs_Control_full.csv"
)

map2["Gene"] = map2["Gene"].map(clean_gene)

map2 = map2.rename(
    columns={
        "logFC": "MAP2_LFC",
        "P.Value": "MAP2_P",
        "adj.P.Val": "MAP2_padj",
    }
)

map2 = map2[
    [
        "Gene",
        "MAP2_LFC",
        "MAP2_P",
        "MAP2_padj",
    ]
].copy()


# ============================================================
# 5. Formal CAG-vs-MAP2 interaction
# ============================================================

interaction = pd.read_csv(
    OSH_DIR
    / "MAP2_vs_CAG_interaction_full.csv"
)

interaction["Gene"] = interaction["Gene"].map(clean_gene)

interaction = interaction.rename(
    columns={
        "logFC": "Osh_ContextInteraction_LFC",
        "P.Value": "Osh_ContextInteraction_P",
        "adj.P.Val": "Osh_ContextInteraction_padj",
    }
)

interaction = interaction[
    [
        "Gene",
        "Osh_ContextInteraction_LFC",
        "Osh_ContextInteraction_P",
        "Osh_ContextInteraction_padj",
    ]
].copy()


# ============================================================
# 6. Full union matrix
# ============================================================

frames = [
    rb,
    e2f,
    cag,
    map2,
    interaction,
]

master = frames[0]

for x in frames[1:]:
    master = master.merge(
        x,
        on="Gene",
        how="outer",
    )


# ============================================================
# 7. Basic significance flags
# ============================================================

master["RbTKO_strict"] = (
    master["RbTKO_LFC"].notna()
)

master["E2F13_FDR05"] = (
    master["E2F13_padj"] < 0.05
)

master["E2F13_strict"] = (
    (master["E2F13_padj"] < 0.05)
    &
    (master["E2F13_LFC"].abs() >= 1)
)

master["CAG_FDR05"] = (
    master["CAG_padj"] < 0.05
)

master["CAG_strict"] = (
    (master["CAG_padj"] < 0.05)
    &
    (master["CAG_LFC"].abs() >= 1)
)

master["MAP2_FDR05"] = (
    master["MAP2_padj"] < 0.05
)

master["MAP2_strict"] = (
    (master["MAP2_padj"] < 0.05)
    &
    (master["MAP2_LFC"].abs() >= 1)
)


# ============================================================
# 8. Signs
# ============================================================

for prefix in [
    "RbTKO",
    "E2F13",
    "CAG",
    "MAP2",
]:

    master[f"{prefix}_Direction"] = (
        master[f"{prefix}_LFC"]
        .apply(sign_class)
    )


# ============================================================
# 9. Rb-TKO vs E2F1/3 reciprocity
# ============================================================

has_adult_pair = (
    master["RbTKO_LFC"].notna()
    &
    master["E2F13_LFC"].notna()
)

master["Adult_sign_relation"] = "Untestable"

same = (
    has_adult_pair
    &
    (
        np.sign(master["RbTKO_LFC"])
        ==
        np.sign(master["E2F13_LFC"])
    )
)

opposite = (
    has_adult_pair
    &
    (
        np.sign(master["RbTKO_LFC"])
        ==
        -np.sign(master["E2F13_LFC"])
    )
)

master.loc[
    same,
    "Adult_sign_relation",
] = "Same_direction"

master.loc[
    opposite,
    "Adult_sign_relation",
] = "E2F_reciprocal"


# Require significant E2F1/3 perturbation for strong reciprocity

master["Adult_E2F_reciprocal_strong"] = (
    opposite
    &
    master["E2F13_FDR05"]
)


# ============================================================
# 10. Conservation of Rb loss across Fong / CAG / MAP2
# ============================================================

three_available = (
    master["RbTKO_LFC"].notna()
    &
    master["CAG_LFC"].notna()
    &
    master["MAP2_LFC"].notna()
)

same_three_sign = (
    three_available
    &
    (
        np.sign(master["RbTKO_LFC"])
        ==
        np.sign(master["CAG_LFC"])
    )
    &
    (
        np.sign(master["RbTKO_LFC"])
        ==
        np.sign(master["MAP2_LFC"])
    )
)

master["Rb_loss_sign_conserved"] = (
    same_three_sign
)


# ============================================================
# 11. CAG-vs-MAP2 descriptive categories
#
# IMPORTANT:
# "threshold_only" deliberately avoids claiming a formally
# significant interaction.
# ============================================================

master["Oshikawa_pattern"] = "No_strict_pattern"

both_strict = (
    master["CAG_strict"]
    &
    master["MAP2_strict"]
)

same_osh_sign = (
    both_strict
    &
    (
        np.sign(master["CAG_LFC"])
        ==
        np.sign(master["MAP2_LFC"])
    )
)

opposite_osh_sign = (
    both_strict
    &
    (
        np.sign(master["CAG_LFC"])
        !=
        np.sign(master["MAP2_LFC"])
    )
)

master.loc[
    same_osh_sign,
    "Oshikawa_pattern",
] = "Conserved_strict"

master.loc[
    opposite_osh_sign,
    "Oshikawa_pattern",
] = "Strict_direction_reversal"

master.loc[
    master["CAG_strict"]
    &
    ~master["MAP2_strict"],
    "Oshikawa_pattern",
] = "CAG_threshold_only"

master.loc[
    ~master["CAG_strict"]
    &
    master["MAP2_strict"],
    "Oshikawa_pattern",
] = "MAP2_threshold_only"


# Formal interaction remains a separate statistical variable

master["Oshikawa_interaction_FDR05"] = (
    master["Osh_ContextInteraction_padj"]
    < 0.05
)

master["Oshikawa_interaction_rawP05"] = (
    master["Osh_ContextInteraction_P"]
    < 0.05
)

master["Oshikawa_large_nominal_interaction"] = (
    (master["Osh_ContextInteraction_P"] < 0.05)
    &
    (
        master[
            "Osh_ContextInteraction_LFC"
        ].abs()
        >= 1
    )
)


# ============================================================
# 12. Useful response magnitudes
# ============================================================

master["RbTKO_vs_E2F13_opposition_score"] = (
    master["RbTKO_LFC"]
    -
    master["E2F13_LFC"]
)

master["CAG_MAP2_difference"] = (
    master["MAP2_LFC"]
    -
    master["CAG_LFC"]
)


# ============================================================
# 13. Sort + export master
# ============================================================

master = master.sort_values(
    "Gene"
).reset_index(
    drop=True
)

master.to_csv(
    OUT_DIR
    / "F8B_all_context_union.csv",
    index=False,
)


# ============================================================
# 14. Fong-strict anchored matrix
#
# This is our primary adult Rb-loss discovery universe.
# ============================================================

fong_anchor = master[
    master["RbTKO_strict"]
].copy()

fong_anchor.to_csv(
    OUT_DIR
    / "F8B_fong_strict_context_matrix.csv",
    index=False,
)


# ============================================================
# 15. Fully observed Fong-strict genes
# ============================================================

complete = fong_anchor[
    fong_anchor[
        [
            "E2F13_LFC",
            "CAG_LFC",
            "MAP2_LFC",
        ]
    ]
    .notna()
    .all(axis=1)
].copy()

complete.to_csv(
    OUT_DIR
    / "F8B_fong_strict_complete_contexts.csv",
    index=False,
)


# ============================================================
# 16. Summaries
# ============================================================

print()
print("=" * 92)
print("F8B — FOUR-CONTEXT EMPIRICAL RESPONSE MATRIX")
print("=" * 92)

print(
    "\nUnion genes:",
    len(master),
)

print(
    "Fong strict Rb-TKO genes:",
    len(fong_anchor),
)

print(
    "Fong-strict genes observed in E2F13 + CAG + MAP2:",
    len(complete),
)


# Adult E2F reciprocity among Fong genes

adult_testable = fong_anchor[
    fong_anchor["E2F13_LFC"].notna()
].copy()

adult_sig = fong_anchor[
    fong_anchor["E2F13_FDR05"]
].copy()

print()
print("ADULT Rb-TKO vs E2F1/3-DKO")
print("-" * 92)

print(
    "Fong strict genes with E2F1/3 measurement:",
    len(adult_testable),
)

print(
    "Fong strict genes with significant E2F1/3 response:",
    len(adult_sig),
)

if len(adult_sig) > 0:

    reciprocal_n = (
        adult_sig[
            "Adult_E2F_reciprocal_strong"
        ].sum()
    )

    print(
        "Strong reciprocal genes:",
        reciprocal_n,
        f"({reciprocal_n / len(adult_sig):.1%})",
    )

    adult_corr = adult_sig[
        [
            "RbTKO_LFC",
            "E2F13_LFC",
        ]
    ].corr(
        method="pearson"
    ).iloc[0, 1]

    print(
        "Pearson correlation among significant paired genes:",
        round(adult_corr, 4),
    )


# Conserved Rb-loss response

conserved = fong_anchor[
    fong_anchor[
        "Rb_loss_sign_conserved"
    ]
].copy()

print()
print("RB-LOSS CONSERVATION")
print("-" * 92)

print(
    "Fong strict genes with same sign in Fong + CAG + MAP2:",
    len(conserved),
)


# Oshikawa pattern counts

print()
print("OSHIKAWA DESCRIPTIVE PATTERNS")
print("-" * 92)

print(
    master[
        "Oshikawa_pattern"
    ]
    .value_counts()
    .to_string()
)


print()
print("Formal interaction FDR < 0.05:")
print(
    master[
        "Oshikawa_interaction_FDR05"
    ].sum()
)

print(
    "Nominal interaction P < 0.05:",
    master[
        "Oshikawa_interaction_rawP05"
    ].sum()
)

print(
    "Nominal interaction P < 0.05 & |interaction LFC| >= 1:",
    master[
        "Oshikawa_large_nominal_interaction"
    ].sum()
)


# ============================================================
# 17. Candidate genes of immediate interest
# ============================================================

candidates = [
    "E2F1",
    "E2F2",
    "MYC",
    "CDK1",
    "CDK2",
    "PCNA",
    "CLSPN",
    "TOPBP1",
    "CHEK1",
    "CHEK2",
    "H2AFX",
    "BRCA1",
    "DDIT3",
    "FAS",
    "TNFRSF10B",
    "PMAIP1",
    "BBC3",
    "ZBP1",
    "PYCARD",
    "MLKL",
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
    "NR4A3",
]

candidate_table = master[
    master["Gene"].isin(
        candidates
    )
].copy()

candidate_table.to_csv(
    OUT_DIR
    / "F8B_key_candidate_contexts.csv",
    index=False,
)


show_cols = [
    "Gene",
    "RbTKO_LFC",
    "E2F13_LFC",
    "E2F13_padj",
    "CAG_LFC",
    "CAG_padj",
    "MAP2_LFC",
    "MAP2_padj",
    "Osh_ContextInteraction_LFC",
    "Osh_ContextInteraction_padj",
    "Adult_sign_relation",
    "Rb_loss_sign_conserved",
]

print()
print("KEY CANDIDATES")
print("-" * 92)

print(
    candidate_table[
        show_cols
    ]
    .sort_values("Gene")
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.3g}",
    )
)


# ============================================================
# 18. Summary file
# ============================================================

summary_path = (
    OUT_DIR
    / "F8B_summary.txt"
)

with open(
    summary_path,
    "w",
    encoding="utf-8",
) as f:

    f.write(
        "F8B — FOUR-CONTEXT EMPIRICAL RESPONSE MATRIX\n"
    )

    f.write(
        "============================================\n\n"
    )

    f.write(
        f"Union genes: {len(master)}\n"
    )

    f.write(
        f"Fong strict genes: {len(fong_anchor)}\n"
    )

    f.write(
        "Fong strict genes with all three additional "
        f"context measurements: {len(complete)}\n\n"
    )

    if len(adult_sig) > 0:

        reciprocal_n = (
            adult_sig[
                "Adult_E2F_reciprocal_strong"
            ].sum()
        )

        f.write(
            "Adult RbTKO/E2F1-3 reciprocity:\n"
        )

        f.write(
            f"  significant paired genes = {len(adult_sig)}\n"
        )

        f.write(
            f"  reciprocal = {reciprocal_n}\n"
        )

        f.write(
            f"  reciprocal fraction = "
            f"{reciprocal_n / len(adult_sig):.4f}\n"
        )

        f.write(
            f"  Pearson r = {adult_corr:.4f}\n\n"
        )

    f.write(
        "Important interpretation:\n"
    )

    f.write(
        "- RbTKO absence is NA, never assumed to be zero.\n"
    )

    f.write(
        "- CAG_threshold_only and MAP2_threshold_only are "
        "descriptive threshold classes, not formal interactions.\n"
    )

    f.write(
        "- Formal developmental-context inference remains the "
        "MAP2-vs-CAG interaction contrast.\n"
    )


print()
print("Outputs written to:")
print(OUT_DIR)

print()
print(
    "[NEXT] F8C: test network modules and biological state "
    "programs for coherent context interactions, then derive "
    "the empirical state variables for the dynamical model."
)