from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd


# ============================================================
# F11G — Minimal molecular control sets
#
# Goal
# ----
# Identify molecular nodes that cover the largest fraction of
# evidence supporting the dynamically important state edges.
#
# This is NOT a quantitative prediction of knockout efficacy.
#
# A molecular node "covers" an edge when it is either the
# source or target of that edge.
#
# We focus primarily on:
#
#   U -> D
#   D -> U
#   D -> D
#
# because these were the state-level interactions implicated
# most directly in commitment/rescue.
#
# U -> U is retained as a secondary conserved-backbone role.
# ============================================================


ROOT = Path(".")

INFILE = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F11F_MolecularLeverage"
    / "F11F_switch_molecular_edges.csv"
)

OUT = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F11G_MolecularControlSets"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


df = pd.read_csv(INFILE)


# ============================================================
# 1. Roles
# ============================================================

CORE_ROLES = [
    "U_TO_D_GATE",
    "D_TO_U_RETURN",
    "D_SELF_REINFORCEMENT",
]

ALL_ROLES = [
    "U_TO_D_GATE",
    "D_TO_U_RETURN",
    "D_SELF_REINFORCEMENT",
    "U_SELF_REINFORCEMENT",
]


df = df[
    df["Switch_role"].isin(
        ALL_ROLES
    )
].copy()


# ============================================================
# 2. Edge weights
#
# Use F11F exploratory priority score directly.
#
# This is still only a ranking heuristic.
# ============================================================

df["Edge_weight"] = pd.to_numeric(
    df["Exploratory_priority_score"],
    errors="coerce",
).fillna(0.0)


# Unique molecular relation within role

df = (
    df.sort_values(
        "Edge_weight",
        ascending=False,
    )
    .drop_duplicates(
        [
            "Switch_role",
            "Source",
            "Target",
        ]
    )
    .reset_index(drop=True)
)


# ============================================================
# 3. Total evidence per role
# ============================================================

role_totals = (
    df.groupby(
        "Switch_role"
    )[
        "Edge_weight"
    ]
    .sum()
    .to_dict()
)


core = df[
    df[
        "Switch_role"
    ].isin(
        CORE_ROLES
    )
].copy()


core_total = (
    core["Edge_weight"]
    .sum()
)


all_total = (
    df["Edge_weight"]
    .sum()
)


# ============================================================
# 4. Candidate molecular nodes
# ============================================================

genes = sorted(
    set(
        df["Source"]
    )
    |
    set(
        df["Target"]
    )
)


print()
print("=" * 105)
print("F11G — MINIMAL MOLECULAR CONTROL SETS")
print("=" * 105)

print()
print(
    "Switch genes:",
    len(genes),
)

print(
    "Core switch edges:",
    len(core),
)

print(
    "All switch edges:",
    len(df),
)


# ============================================================
# 5. Coverage function
# ============================================================

def evaluate_gene_set(
    selected,
):

    selected = set(
        selected
    )


    hit = (
        df["Source"].isin(
            selected
        )
        |
        df["Target"].isin(
            selected
        )
    )


    covered = df[
        hit
    ]


    core_covered = covered[
        covered[
            "Switch_role"
        ].isin(
            CORE_ROLES
        )
    ]


    row = {

        "Genes":
            ";".join(
                sorted(selected)
            ),

        "Set_size":
            len(selected),

        "Edges_covered":
            len(covered),

        "Core_edges_covered":
            len(core_covered),

        "Weighted_core_coverage":
            (
                core_covered[
                    "Edge_weight"
                ].sum()
                /
                core_total
                if core_total > 0
                else np.nan
            ),

        "Weighted_all_coverage":
            (
                covered[
                    "Edge_weight"
                ].sum()
                /
                all_total
                if all_total > 0
                else np.nan
            ),
    }


    role_count = 0


    for role in ALL_ROLES:

        role_df = df[
            df[
                "Switch_role"
            ]
            ==
            role
        ]


        role_hit = covered[
            covered[
                "Switch_role"
            ]
            ==
            role
        ]


        denom = role_totals.get(
            role,
            0,
        )


        coverage = (
            role_hit[
                "Edge_weight"
            ].sum()
            /
            denom
            if denom > 0
            else np.nan
        )


        row[
            f"Coverage_{role}"
        ] = coverage


        if (
            pd.notna(
                coverage
            )
            and
            coverage > 0
        ):

            role_count += 1


    row[
        "Roles_touched"
    ] = role_count


    return row


# ============================================================
# 6. Single genes
# ============================================================

single_rows = [
    evaluate_gene_set(
        [gene]
    )
    for gene in genes
]


singles = pd.DataFrame(
    single_rows
).sort_values(
    [
        "Weighted_core_coverage",
        "Roles_touched",
    ],
    ascending=False,
)


singles.to_csv(
    OUT
    / "F11G_single_gene_control.csv",
    index=False,
)


# ============================================================
# 7. Restrict combination search
#
# Use top 15 single nodes to keep combinations interpretable.
# ============================================================

TOP_N = min(
    15,
    len(singles),
)


top_genes = (
    singles
    .head(TOP_N)[
        "Genes"
    ]
    .tolist()
)


combo_rows = []


for size in [
    2,
    3,
]:

    for combo in combinations(
        top_genes,
        size,
    ):

        combo_rows.append(
            evaluate_gene_set(
                combo
            )
        )


combos = pd.DataFrame(
    combo_rows
)


combos = combos.sort_values(
    [
        "Weighted_core_coverage",
        "Weighted_all_coverage",
        "Roles_touched",
    ],
    ascending=False,
)


combos.to_csv(
    OUT
    / "F11G_combination_control.csv",
    index=False,
)


# ============================================================
# 8. Pareto-like minimal sets
#
# For each set size, report best core coverage.
# ============================================================

best_rows = []


for size in [
    1,
    2,
    3,
]:

    source = (
        singles
        if size == 1
        else combos[
            combos[
                "Set_size"
            ]
            ==
            size
        ]
    )


    if source.empty:
        continue


    best = source.iloc[
        source[
            "Weighted_core_coverage"
        ]
        .argmax()
    ]


    best_rows.append(
        best
    )


best_sets = pd.DataFrame(
    best_rows
)


best_sets.to_csv(
    OUT
    / "F11G_best_control_sets.csv",
    index=False,
)


# ============================================================
# 9. Find sets covering all three core roles
# ============================================================

core_coverage_cols = [
    f"Coverage_{role}"
    for role in CORE_ROLES
]


def covers_all_core(row):

    return all(
        (
            pd.notna(
                row[col]
            )
            and
            row[col] > 0
        )
        for col in core_coverage_cols
    )


singles[
    "Covers_all_core_roles"
] = singles.apply(
    covers_all_core,
    axis=1,
)


combos[
    "Covers_all_core_roles"
] = combos.apply(
    covers_all_core,
    axis=1,
)


all_core_sets = pd.concat(
    [
        singles[
            singles[
                "Covers_all_core_roles"
            ]
        ],
        combos[
            combos[
                "Covers_all_core_roles"
            ]
        ],
    ],
    ignore_index=True,
)


all_core_sets = (
    all_core_sets.sort_values(
        [
            "Set_size",
            "Weighted_core_coverage",
        ],
        ascending=[
            True,
            False,
        ],
    )
)


all_core_sets.to_csv(
    OUT
    / "F11G_all_core_role_sets.csv",
    index=False,
)


# ============================================================
# 10. Console output
# ============================================================

SHOW = [
    "Genes",
    "Set_size",
    "Edges_covered",
    "Core_edges_covered",
    "Roles_touched",
    "Weighted_core_coverage",
    "Weighted_all_coverage",
    "Coverage_U_TO_D_GATE",
    "Coverage_D_TO_U_RETURN",
    "Coverage_D_SELF_REINFORCEMENT",
    "Coverage_U_SELF_REINFORCEMENT",
]


print()
print("TOP SINGLE-GENE CONTROL NODES")
print("-" * 105)

print(
    singles[
        SHOW
    ]
    .head(15)
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.3f}",
    )
)


print()
print("TOP TWO-GENE CONTROL SETS")
print("-" * 105)

print(
    combos[
        combos[
            "Set_size"
        ]
        ==
        2
    ][
        SHOW
    ]
    .head(15)
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.3f}",
    )
)


print()
print("TOP THREE-GENE CONTROL SETS")
print("-" * 105)

print(
    combos[
        combos[
            "Set_size"
        ]
        ==
        3
    ][
        SHOW
    ]
    .head(15)
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.3f}",
    )
)


print()
print("SMALLEST SETS TOUCHING ALL THREE CORE SWITCH ROLES")
print("-" * 105)

if all_core_sets.empty:

    print(
        "No tested set touched all three core roles."
    )

else:

    print(
        all_core_sets[
            SHOW
        ]
        .head(20)
        .to_string(
            index=False,
            float_format=lambda x: f"{x:.3f}",
        )
    )


print()
print("BEST SET AT EACH SIZE")
print("-" * 105)

print(
    best_sets[
        SHOW
    ]
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.3f}",
    )
)


print()
print("Outputs written to:")
print(OUT)

print()
print(
    "[NEXT] Convert the strongest molecular control sets "
    "into explicit experimental perturbation hypotheses, "
    "with predicted readouts for U, D, T and cell fate."
)