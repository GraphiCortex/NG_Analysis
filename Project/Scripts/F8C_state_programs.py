from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats


# ============================================================
# F8C — Empirical state-program inference
#
# Goals:
#
# 1. Collapse mechanistically related genes into biological
#    state variables.
#
# 2. Orient each state relative to the adult Fong Rb-TKO
#    response.
#
# 3. Test CAG-vs-MAP2 context interaction at the MODULE level
#    using the actual 12 Oshikawa biological replicates.
#
# 4. Describe each module across:
#
#       Rb-TKO
#       E2F1/3-DKO
#       CAG-TKO
#       MAP2-TKO
#
# ============================================================


ROOT = Path(".")

F8A = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F8A_Oshikawa"
)

F8B = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F8B_ContextMatrix"
)

OUT = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F8C_StatePrograms"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 1. Mechanistically defined state programs
# ============================================================
#
# These are deliberately compact.
#
# They are based on the F2/F4/F5 architecture already derived
# in the project, rather than on genes selected from Oshikawa.
#
# Overlap between programs is allowed because biological state
# variables need not be disjoint.
# ============================================================

PROGRAMS = {

    "E2F_CELL_CYCLE": [
        "E2F1",
        "E2F2",
        "CDK1",
        "CDK2",
        "CCND1",
        "CDKN1A",
        "CDKN2D",
        "MYBL2",
    ],

    "MITOTIC_EXECUTION": [
        "AURKB",
        "BIRC5",
        "BUB1",
        "CDK1",
        "CIT",
        "ECT2",
        "KIF14",
        "MAD2L1",
        "MELK",
        "PLK1",
        "RAD21",
        "TOP2A",
        "TPX2",
    ],

    "REPLICATION_ATR": [
        "PCNA",
        "CLSPN",
        "TOPBP1",
        "CHEK1",
        "H2AFX",
        "MCM2",
        "MCM10",
        "POLD1",
        "POLE",
        "RPA2",
        "TIMELESS",
        "TIPIN",
        "CHAF1A",
        "CHAF1B",
        "DONSON",
        "LIG1",
    ],

    "DNA_REPAIR_CHECKPOINT": [
        "BRCA1",
        "BRCA2",
        "CHEK1",
        "CHEK2",
        "BLM",
        "PALB2",
        "XRCC2",
        "ATAD5",
        "BARD1",
        "FIGNL1",
        "EXO1",
    ],

    "DEATH_RECEPTOR_STRESS": [
        "DDIT3",
        "FAS",
        "TNFRSF10B",
        "EIF2AK2",
        "MAP3K5",
        "TRIB3",
        "ATF5",
        "TNFRSF22",
    ],

    "INTRINSIC_APOPTOSIS": [
        "PMAIP1",
        "BBC3",
    ],

    "REGULATED_DEATH": [
        "ZBP1",
        "PYCARD",
        "MLKL",
        "EIF2AK2",
        "TNF",
    ],

    "RTK_TROPHIC_STATE": [
        "ANGPT1",
        "BDNF",
        "FGF2",
        "IRS2",
        "KDR",
        "MET",
        "CXCL12",
        "PTK2B",
        "RASGRF1",
        "RASGRF2",
        "THBS1",
        "UNC5B",
    ],

    "INNATE_STRESS": [
        "DTX3L",
        "EIF2AK2",
        "IFIT2",
        "IFIT3",
        "IGTP",
        "IIGP1",
        "IRGM1",
        "PARP9",
        "TRIM25",
    ],

    "NR4A_STATE": [
        "NR4A1",
        "NR4A2",
        "NR4A3",
    ],
}


# ============================================================
# 2. Load expression and metadata
# ============================================================

expr = pd.read_csv(
    F8A
    / "oshikawa_gene_level_expression.csv"
)

meta = pd.read_csv(
    F8A
    / "sample_metadata.csv"
)

context = pd.read_csv(
    F8B
    / "F8B_all_context_union.csv"
)


expr["Gene"] = (
    expr["Gene"]
    .astype(str)
    .str.upper()
    .str.strip()
)

context["Gene"] = (
    context["Gene"]
    .astype(str)
    .str.upper()
    .str.strip()
)

expr = expr.set_index("Gene")

context = context.set_index("Gene")


sample_order = meta["Sample"].tolist()

missing_samples = [
    s
    for s in sample_order
    if s not in expr.columns
]

if missing_samples:
    raise RuntimeError(
        f"Expression matrix missing samples: {missing_samples}"
    )


expr = expr[
    sample_order
].copy()


# ============================================================
# 3. Metadata coding
# ============================================================

meta = meta.set_index(
    "Sample"
).loc[
    sample_order
].reset_index()


meta["TKO"] = (
    meta["Group"]
    .str.endswith("TKO")
    .astype(int)
)

meta["MAP2"] = (
    meta["Group"]
    .str.startswith("MAP2")
    .astype(int)
)

meta["Interaction"] = (
    meta["TKO"]
    *
    meta["MAP2"]
)


# Design:
#
# intercept = CAG control
# beta_TKO = CAG TKO effect
# beta_MAP2 = baseline MAP2-vs-CAG control difference
# beta_INT = MAP2 TKO effect - CAG TKO effect

X = np.column_stack([
    np.ones(len(meta)),
    meta["TKO"].values,
    meta["MAP2"].values,
    meta["Interaction"].values,
])


# ============================================================
# 4. BH correction
# ============================================================

def bh_adjust(pvalues):

    p = np.asarray(
        pvalues,
        dtype=float,
    )

    out = np.full(
        len(p),
        np.nan,
    )

    valid = np.isfinite(p)

    pv = p[valid]

    if len(pv) == 0:
        return out

    order = np.argsort(pv)

    ranked = pv[order]

    m = len(ranked)

    adj = (
        ranked
        *
        m
        /
        np.arange(
            1,
            m + 1,
        )
    )

    adj = np.minimum.accumulate(
        adj[::-1]
    )[::-1]

    adj = np.minimum(
        adj,
        1,
    )

    temp = np.empty(m)

    temp[order] = adj

    out[valid] = temp

    return out


# ============================================================
# 5. Linear model helper
# ============================================================

def fit_module_model(y):

    y = np.asarray(
        y,
        dtype=float,
    )

    beta = np.linalg.lstsq(
        X,
        y,
        rcond=None,
    )[0]

    residuals = (
        y
        -
        X @ beta
    )

    df_resid = (
        len(y)
        -
        X.shape[1]
    )

    sigma2 = (
        np.sum(
            residuals ** 2
        )
        /
        df_resid
    )

    XtX_inv = np.linalg.inv(
        X.T @ X
    )

    cov_beta = (
        sigma2
        *
        XtX_inv
    )


    def contrast(c):

        c = np.asarray(
            c,
            dtype=float,
        )

        estimate = float(
            c @ beta
        )

        se = float(
            np.sqrt(
                c
                @ cov_beta
                @ c
            )
        )

        if se == 0:
            tval = np.nan
            pval = np.nan

        else:

            tval = (
                estimate
                /
                se
            )

            pval = (
                2
                *
                stats.t.sf(
                    abs(tval),
                    df=df_resid,
                )
            )

        return (
            estimate,
            se,
            tval,
            pval,
        )


    # CAG TKO - CAG Control
    cag = contrast(
        [0, 1, 0, 0]
    )

    # MAP2 TKO - MAP2 Control
    map2 = contrast(
        [0, 1, 0, 1]
    )

    # Difference-of-differences
    interaction = contrast(
        [0, 0, 0, 1]
    )


    return {
        "beta": beta,
        "sigma2": sigma2,
        "df_resid": df_resid,
        "CAG": cag,
        "MAP2": map2,
        "Interaction": interaction,
    }


# ============================================================
# 6. Build Fong-oriented Oshikawa module scores
# ============================================================

score_table = meta[
    [
        "Sample",
        "Group",
    ]
].copy()

module_results = []

coverage_rows = []


for module_name, genes in PROGRAMS.items():

    genes = [
        g.upper()
        for g in genes
    ]

    # --------------------------------------------
    # Must have:
    #   - Oshikawa expression
    #   - Fong Rb-TKO direction
    # --------------------------------------------

    available = []

    for gene in genes:

        if (
            gene in expr.index
            and
            gene in context.index
            and
            pd.notna(
                context.loc[
                    gene,
                    "RbTKO_LFC"
                ]
            )
        ):
            available.append(
                gene
            )


    coverage_rows.append({
        "Program": module_name,
        "Requested_genes": len(genes),
        "Genes_used": len(available),
        "Coverage_fraction": (
            len(available)
            /
            len(genes)
        ),
        "Genes": ";".join(
            available
        ),
    })


    if len(available) < 2:

        print(
            f"[WARNING] {module_name}: "
            f"only {len(available)} usable genes."
        )

        continue


    # --------------------------------------------
    # Expression matrix:
    # genes x samples
    # --------------------------------------------

    M = (
        expr.loc[
            available,
            sample_order,
        ]
        .astype(float)
        .copy()
    )


    # --------------------------------------------
    # Z-score every gene across the 12 samples
    # --------------------------------------------

    means = M.mean(
        axis=1
    )

    sds = M.std(
        axis=1,
        ddof=1,
    )


    valid_variance = (
        sds > 0
    )

    M = M.loc[
        valid_variance
    ]

    means = means.loc[
        valid_variance
    ]

    sds = sds.loc[
        valid_variance
    ]


    Z = M.sub(
        means,
        axis=0,
    ).div(
        sds,
        axis=0,
    )


    used = Z.index.tolist()


    # --------------------------------------------
    # Orient every gene according to adult Rb-TKO
    #
    # Positive module score =
    # more "adult Rb-TKO-like" for that program
    # --------------------------------------------

    fong_sign = np.sign(
        context.loc[
            used,
            "RbTKO_LFC"
        ].astype(float)
    )

    oriented = Z.mul(
        fong_sign,
        axis=0,
    )


    score = oriented.mean(
        axis=0
    )


    score_table[
        module_name
    ] = (
        score.loc[
            sample_order
        ].values
    )


    # --------------------------------------------
    # Fit module interaction model
    # --------------------------------------------

    model = fit_module_model(
        score.loc[
            sample_order
        ].values
    )


    cag_est, cag_se, cag_t, cag_p = (
        model["CAG"]
    )

    map2_est, map2_se, map2_t, map2_p = (
        model["MAP2"]
    )

    int_est, int_se, int_t, int_p = (
        model["Interaction"]
    )


    # --------------------------------------------
    # Group means
    # --------------------------------------------

    temp = score_table[
        [
            "Group",
            module_name,
        ]
    ]

    group_means = (
        temp.groupby(
            "Group"
        )[module_name]
        .mean()
    )


    module_results.append({

        "Program": module_name,

        "Genes_used": len(used),

        "CAG_Control_mean":
            group_means.get(
                "CAG_Control",
                np.nan,
            ),

        "CAG_TKO_mean":
            group_means.get(
                "CAG_TKO",
                np.nan,
            ),

        "MAP2_Control_mean":
            group_means.get(
                "MAP2_Control",
                np.nan,
            ),

        "MAP2_TKO_mean":
            group_means.get(
                "MAP2_TKO",
                np.nan,
            ),

        "CAG_TKO_effect":
            cag_est,

        "CAG_effect_P":
            cag_p,

        "MAP2_TKO_effect":
            map2_est,

        "MAP2_effect_P":
            map2_p,

        "MAP2_minus_CAG_interaction":
            int_est,

        "Interaction_P":
            int_p,
    })


module_results = pd.DataFrame(
    module_results
)


# ============================================================
# 7. Multiple-testing correction
# ============================================================

module_results[
    "Interaction_FDR"
] = bh_adjust(
    module_results[
        "Interaction_P"
    ].values
)

module_results[
    "CAG_effect_FDR"
] = bh_adjust(
    module_results[
        "CAG_effect_P"
    ].values
)

module_results[
    "MAP2_effect_FDR"
] = bh_adjust(
    module_results[
        "MAP2_effect_P"
    ].values
)


# ============================================================
# 8. Four-context module response geometry
# ============================================================
#
# For each gene:
#
# oriented response =
#     sign(Fong Rb-TKO LFC) * context LFC
#
# Therefore:
#
#   > 0 = same direction as adult Rb-TKO
#   < 0 = opposite direction
#
# For E2F1/3-DKO, a strongly NEGATIVE value is what we expect
# from an E2F-dependent Rb-loss program.
# ============================================================

contexts = {
    "RbTKO": (
        "RbTKO_LFC",
        None,
    ),

    "E2F13": (
        "E2F13_LFC",
        "E2F13_padj",
    ),

    "CAG": (
        "CAG_LFC",
        "CAG_padj",
    ),

    "MAP2": (
        "MAP2_LFC",
        "MAP2_padj",
    ),
}


response_rows = []


for module_name, genes in PROGRAMS.items():

    genes = [
        g.upper()
        for g in genes
    ]

    genes = [
        g
        for g in genes
        if (
            g in context.index
            and
            pd.notna(
                context.loc[
                    g,
                    "RbTKO_LFC"
                ]
            )
        )
    ]


    for ctx_name, (
        lfc_col,
        padj_col,
    ) in contexts.items():

        vals = []

        vals_sig = []

        same_sign = []

        same_sign_sig = []


        for gene in genes:

            lfc = context.loc[
                gene,
                lfc_col,
            ]

            fong_lfc = context.loc[
                gene,
                "RbTKO_LFC",
            ]

            if pd.isna(lfc):
                continue

            oriented_value = (
                np.sign(
                    fong_lfc
                )
                *
                float(lfc)
            )

            vals.append(
                oriented_value
            )

            same_sign.append(
                oriented_value > 0
            )


            # RbTKO file is already strict by construction

            if ctx_name == "RbTKO":

                vals_sig.append(
                    oriented_value
                )

                same_sign_sig.append(
                    oriented_value > 0
                )

            elif (
                padj_col is not None
                and
                pd.notna(
                    context.loc[
                        gene,
                        padj_col,
                    ]
                )
                and
                context.loc[
                    gene,
                    padj_col,
                ] < 0.05
            ):

                vals_sig.append(
                    oriented_value
                )

                same_sign_sig.append(
                    oriented_value > 0
                )


        response_rows.append({

            "Program":
                module_name,

            "Context":
                ctx_name,

            "N":
                len(vals),

            "Oriented_mean":
                np.mean(vals)
                if vals
                else np.nan,

            "Oriented_median":
                np.median(vals)
                if vals
                else np.nan,

            "Adult_sign_agreement":
                np.mean(same_sign)
                if same_sign
                else np.nan,

            "N_significant":
                len(vals_sig),

            "Oriented_mean_significant":
                np.mean(vals_sig)
                if vals_sig
                else np.nan,

            "Adult_sign_agreement_significant":
                np.mean(
                    same_sign_sig
                )
                if same_sign_sig
                else np.nan,
        })


response = pd.DataFrame(
    response_rows
)


# ============================================================
# 9. Export
# ============================================================

score_table.to_csv(
    OUT
    / "F8C_oshikawa_sample_state_scores.csv",
    index=False,
)

module_results.to_csv(
    OUT
    / "F8C_module_interaction_results.csv",
    index=False,
)

pd.DataFrame(
    coverage_rows
).to_csv(
    OUT
    / "F8C_program_gene_coverage.csv",
    index=False,
)

response.to_csv(
    OUT
    / "F8C_four_context_program_responses.csv",
    index=False,
)


# ============================================================
# 10. Console summary
# ============================================================

print()
print("=" * 105)
print("F8C — EMPIRICAL STATE-PROGRAM ANALYSIS")
print("=" * 105)


print()
print("A. OSHIKAWA MODULE-LEVEL RESPONSE")
print("-" * 105)

display_cols = [
    "Program",
    "Genes_used",
    "CAG_TKO_effect",
    "CAG_effect_P",
    "MAP2_TKO_effect",
    "MAP2_effect_P",
    "MAP2_minus_CAG_interaction",
    "Interaction_P",
    "Interaction_FDR",
]

display = (
    module_results[
        display_cols
    ]
    .sort_values(
        "Interaction_P"
    )
    .copy()
)

print(
    display.to_string(
        index=False,
        float_format=lambda x: f"{x:.4g}",
    )
)


print()
print("Interpretation:")
print(
    "  Positive CAG/MAP2 effect = movement toward the "
    "adult Rb-TKO direction for that state."
)

print(
    "  Negative effect = movement opposite to the "
    "adult Rb-TKO state."
)

print(
    "  Interaction = MAP2 TKO effect minus CAG TKO effect."
)


# ============================================================
# Four-context table
# ============================================================

print()
print()
print("B. FOUR-CONTEXT RESPONSE GEOMETRY")
print("-" * 105)


pivot_mean = response.pivot(
    index="Program",
    columns="Context",
    values="Oriented_mean",
)

pivot_sig = response.pivot(
    index="Program",
    columns="Context",
    values="Oriented_mean_significant",
)

for c in [
    "RbTKO",
    "E2F13",
    "CAG",
    "MAP2",
]:

    if c not in pivot_mean.columns:
        pivot_mean[c] = np.nan

    if c not in pivot_sig.columns:
        pivot_sig[c] = np.nan


pivot_mean = pivot_mean[
    [
        "RbTKO",
        "E2F13",
        "CAG",
        "MAP2",
    ]
]

pivot_sig = pivot_sig[
    [
        "RbTKO",
        "E2F13",
        "CAG",
        "MAP2",
    ]
]


print()
print("All measured genes")
print()

print(
    pivot_mean.to_string(
        float_format=lambda x: f"{x:.3f}",
    )
)


print()
print("Significant genes only")
print()

print(
    pivot_sig.to_string(
        float_format=lambda x: f"{x:.3f}",
    )
)


# ============================================================
# 11. Candidate empirical state variables
# ============================================================

print()
print()
print("C. CANDIDATE STATE VARIABLES")
print("-" * 105)


for _, row in (
    module_results
    .sort_values(
        "Interaction_P"
    )
    .iterrows()
):

    module = row["Program"]

    print(
        f"{module:28s} "
        f"CAG={row['CAG_TKO_effect']:+.3f}  "
        f"MAP2={row['MAP2_TKO_effect']:+.3f}  "
        f"INT={row['MAP2_minus_CAG_interaction']:+.3f}  "
        f"Pint={row['Interaction_P']:.4g}  "
        f"FDR={row['Interaction_FDR']:.4g}"
    )


print()
print("Outputs written to:")
print(OUT)

print()
print(
    "[NEXT] Use these program-level responses to decide which "
    "state variables are empirically distinguishable enough "
    "to enter the first dynamical model."
)