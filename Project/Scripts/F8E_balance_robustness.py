from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd


# ============================================================
# F8E — Robustness analysis for the mitotic/replication
#       imbalance discovered in F8D
#
# Tests:
#   1. Mean-score reproduction
#   2. Median-score reproduction
#   3. Exact within-context permutation test
#   4. Leave-one-gene-out robustness
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
    / "F8E_BalanceRobustness"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Programs
# ============================================================

MITOTIC = [
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
]

REPLICATION = [
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
]

E2F_CYCLE = [
    "E2F1",
    "E2F2",
    "CDK1",
    "CDK2",
    "CCND1",
    "CDKN1A",
    "CDKN2D",
    "MYBL2",
]


# ============================================================
# Load data
# ============================================================

expr = pd.read_csv(
    F8A / "oshikawa_gene_level_expression.csv"
)

meta = pd.read_csv(
    F8A / "sample_metadata.csv"
)

ctx = pd.read_csv(
    F8B / "F8B_all_context_union.csv"
)


expr["Gene"] = (
    expr["Gene"]
    .astype(str)
    .str.upper()
    .str.strip()
)

ctx["Gene"] = (
    ctx["Gene"]
    .astype(str)
    .str.upper()
    .str.strip()
)

expr = expr.set_index("Gene")
ctx = ctx.set_index("Gene")


samples = meta["Sample"].tolist()

expr = expr[samples]


# ============================================================
# Build oriented gene-level z-scores
# ============================================================

def oriented_zscores(genes):

    usable = [
        g for g in genes
        if (
            g in expr.index
            and
            g in ctx.index
            and
            pd.notna(
                ctx.loc[g, "RbTKO_LFC"]
            )
        )
    ]

    M = (
        expr.loc[
            usable,
            samples
        ]
        .astype(float)
        .copy()
    )

    sd = M.std(
        axis=1,
        ddof=1
    )

    M = M.loc[
        sd > 0
    ]

    Z = M.sub(
        M.mean(axis=1),
        axis=0
    ).div(
        M.std(axis=1, ddof=1),
        axis=0
    )

    signs = np.sign(
        ctx.loc[
            Z.index,
            "RbTKO_LFC"
        ].astype(float)
    )

    return Z.mul(
        signs,
        axis=0
    )


mitotic_z = oriented_zscores(
    MITOTIC
)

replication_z = oriented_zscores(
    REPLICATION
)

cycle_z = oriented_zscores(
    E2F_CYCLE
)


print()
print("=" * 100)
print("F8E — BALANCE ROBUSTNESS")
print("=" * 100)

print()
print(
    "Usable mitotic genes:",
    len(mitotic_z)
)

print(
    "Usable replication genes:",
    len(replication_z)
)

print(
    "Usable E2F/cycle genes:",
    len(cycle_z)
)


# ============================================================
# Module score helper
# ============================================================

def module_score(Z, method="mean"):

    if method == "mean":
        return Z.mean(axis=0)

    if method == "median":
        return Z.median(axis=0)

    raise ValueError(method)


# ============================================================
# Difference-of-differences helper
# ============================================================

def interaction_effect(
    values,
    tko_labels=None
):

    values = pd.Series(
        values,
        index=samples,
        dtype=float,
    )

    m = meta.copy()

    if tko_labels is None:

        m["TKO_test"] = (
            m["Group"]
            .str.endswith("TKO")
            .astype(int)
        )

    else:

        m["TKO_test"] = [
            int(tko_labels[s])
            for s in samples
        ]


    cag = m[
        m["Group"]
        .str.startswith("CAG")
    ].copy()

    map2 = m[
        m["Group"]
        .str.startswith("MAP2")
    ].copy()


    cag_tko_samples = cag.loc[
        cag["TKO_test"] == 1,
        "Sample"
    ]

    cag_ctrl_samples = cag.loc[
        cag["TKO_test"] == 0,
        "Sample"
    ]

    map2_tko_samples = map2.loc[
        map2["TKO_test"] == 1,
        "Sample"
    ]

    map2_ctrl_samples = map2.loc[
        map2["TKO_test"] == 0,
        "Sample"
    ]


    cag_effect = (
        values.loc[
            cag_tko_samples
        ].mean()
        -
        values.loc[
            cag_ctrl_samples
        ].mean()
    )

    map2_effect = (
        values.loc[
            map2_tko_samples
        ].mean()
        -
        values.loc[
            map2_ctrl_samples
        ].mean()
    )


    return (
        cag_effect,
        map2_effect,
        map2_effect - cag_effect,
    )


# ============================================================
# 1. Mean and median aggregation
# ============================================================

aggregation_rows = []


for method in [
    "mean",
    "median",
]:

    mit = module_score(
        mitotic_z,
        method
    )

    rep = module_score(
        replication_z,
        method
    )

    cyc = module_score(
        cycle_z,
        method
    )


    balances = {

        "MITOTIC_MINUS_REPLICATION":
            mit - rep,

        "MITOTIC_MINUS_E2F_CYCLE":
            mit - cyc,
    }


    for name, values in balances.items():

        cag, map2, interaction = (
            interaction_effect(
                values
            )
        )

        aggregation_rows.append({
            "Aggregation":
                method,

            "Balance":
                name,

            "CAG_effect":
                cag,

            "MAP2_effect":
                map2,

            "Interaction":
                interaction,
        })


aggregation = pd.DataFrame(
    aggregation_rows
)


# ============================================================
# 2. Exact permutation test
#
# Preserve:
#   3 TKO / 3 control within CAG
#   3 TKO / 3 control within MAP2
#
# Number permutations:
#   C(6,3) * C(6,3) = 400
# ============================================================

cag_samples = (
    meta.loc[
        meta["Group"]
        .str.startswith("CAG"),
        "Sample"
    ]
    .tolist()
)

map2_samples = (
    meta.loc[
        meta["Group"]
        .str.startswith("MAP2"),
        "Sample"
    ]
    .tolist()
)


cag_assignments = list(
    combinations(
        cag_samples,
        3
    )
)

map2_assignments = list(
    combinations(
        map2_samples,
        3
    )
)


perm_rows = []


for method in [
    "mean",
    "median",
]:

    mit = module_score(
        mitotic_z,
        method
    )

    rep = module_score(
        replication_z,
        method
    )

    cyc = module_score(
        cycle_z,
        method
    )


    balance_map = {

        "MITOTIC_MINUS_REPLICATION":
            mit - rep,

        "MITOTIC_MINUS_E2F_CYCLE":
            mit - cyc,
    }


    for balance_name, values in balance_map.items():

        _, _, observed = (
            interaction_effect(
                values
            )
        )

        permuted = []


        for cag_tko in cag_assignments:

            cag_tko = set(
                cag_tko
            )

            for map2_tko in map2_assignments:

                map2_tko = set(
                    map2_tko
                )

                labels = {}

                for s in cag_samples:

                    labels[s] = (
                        s in cag_tko
                    )

                for s in map2_samples:

                    labels[s] = (
                        s in map2_tko
                    )


                _, _, val = (
                    interaction_effect(
                        values,
                        labels
                    )
                )

                permuted.append(
                    val
                )


        permuted = np.asarray(
            permuted
        )

        # Exact two-sided randomization P-value
        empirical_p = (
            np.sum(
                np.abs(permuted)
                >=
                abs(observed)
            )
            /
            len(permuted)
        )


        perm_rows.append({

            "Aggregation":
                method,

            "Balance":
                balance_name,

            "Observed_interaction":
                observed,

            "Permutations":
                len(permuted),

            "Exact_permutation_P":
                empirical_p,

            "Null_mean":
                permuted.mean(),

            "Null_SD":
                permuted.std(
                    ddof=1
                ),
        })


permutation = pd.DataFrame(
    perm_rows
)


# ============================================================
# 3. Leave-one-gene-out:
#    mitotic minus replication
# ============================================================

jackknife_rows = []


# Full reference

full_mit = (
    mitotic_z.mean(
        axis=0
    )
)

full_rep = (
    replication_z.mean(
        axis=0
    )
)

full_balance = (
    full_mit
    -
    full_rep
)

_, _, full_interaction = (
    interaction_effect(
        full_balance
    )
)


# Drop one mitotic gene at a time

for gene in mitotic_z.index:

    score_mit = (
        mitotic_z
        .drop(
            index=gene
        )
        .mean(
            axis=0
        )
    )

    balance = (
        score_mit
        -
        full_rep
    )

    cag, map2, inter = (
        interaction_effect(
            balance
        )
    )

    jackknife_rows.append({

        "Dropped_from":
            "MITOTIC",

        "Dropped_gene":
            gene,

        "CAG_effect":
            cag,

        "MAP2_effect":
            map2,

        "Interaction":
            inter,

        "Same_sign_as_full":
            np.sign(inter)
            ==
            np.sign(full_interaction),
    })


# Drop one replication gene at a time

for gene in replication_z.index:

    score_rep = (
        replication_z
        .drop(
            index=gene
        )
        .mean(
            axis=0
        )
    )

    balance = (
        full_mit
        -
        score_rep
    )

    cag, map2, inter = (
        interaction_effect(
            balance
        )
    )

    jackknife_rows.append({

        "Dropped_from":
            "REPLICATION",

        "Dropped_gene":
            gene,

        "CAG_effect":
            cag,

        "MAP2_effect":
            map2,

        "Interaction":
            inter,

        "Same_sign_as_full":
            np.sign(inter)
            ==
            np.sign(full_interaction),
    })


jackknife = pd.DataFrame(
    jackknife_rows
)


# ============================================================
# 4. Summarize jackknife
# ============================================================

jack_summary = pd.DataFrame([{

    "Full_interaction":
        full_interaction,

    "N_leave_one_out":
        len(jackknife),

    "Min_interaction":
        jackknife[
            "Interaction"
        ].min(),

    "Max_interaction":
        jackknife[
            "Interaction"
        ].max(),

    "Mean_interaction":
        jackknife[
            "Interaction"
        ].mean(),

    "SD_interaction":
        jackknife[
            "Interaction"
        ].std(
            ddof=1
        ),

    "Fraction_same_sign":
        jackknife[
            "Same_sign_as_full"
        ].mean(),
}])


# ============================================================
# Export
# ============================================================

aggregation.to_csv(
    OUT
    / "F8E_mean_median_robustness.csv",
    index=False,
)

permutation.to_csv(
    OUT
    / "F8E_exact_permutation.csv",
    index=False,
)

jackknife.to_csv(
    OUT
    / "F8E_leave_one_gene_out.csv",
    index=False,
)

jack_summary.to_csv(
    OUT
    / "F8E_jackknife_summary.csv",
    index=False,
)


# ============================================================
# Console
# ============================================================

print()
print("A. MEAN / MEDIAN AGGREGATION")
print("-" * 100)

print(
    aggregation.to_string(
        index=False,
        float_format=lambda x: f"{x:+.4f}",
    )
)


print()
print("B. EXACT WITHIN-CONTEXT PERMUTATION")
print("-" * 100)

print(
    permutation.to_string(
        index=False,
        float_format=lambda x: f"{x:.6g}",
    )
)


print()
print("C. LEAVE-ONE-GENE-OUT ROBUSTNESS")
print("-" * 100)

print(
    jack_summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.5f}",
    )
)


print()
print("Most influential deletions:")

jackknife["Shift_from_full"] = (
    jackknife["Interaction"]
    -
    full_interaction
).abs()

print(
    jackknife.sort_values(
        "Shift_from_full",
        ascending=False
    )[
        [
            "Dropped_from",
            "Dropped_gene",
            "Interaction",
            "Shift_from_full",
            "Same_sign_as_full",
        ]
    ]
    .head(15)
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.5f}",
    )
)


print()
print("Outputs written to:")
print(OUT)

print()
print(
    "[NEXT] If the imbalance survives aggregation, exact "
    "permutation, and gene jackknifing, promote it to a "
    "constraint in the first dynamical model."
)