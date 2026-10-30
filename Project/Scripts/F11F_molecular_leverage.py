from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# F11F — Molecular translation of dynamical leverage points
#
# Goal
# ----
# Translate the state-level rescue predictions:
#
#   U -> D
#   D -> U
#   D -> D
#   U -> U
#
# back into the actual signed molecular relations present in
# the corrected F5 network.
#
# We integrate:
#
#   - direct signed F5 molecular edges
#   - supporting records
#   - representative path usage
#   - representative pair usage
#   - number of evidence layers
#   - F10A adult residuals
#
# IMPORTANT
# ---------
# The resulting ranking is exploratory.
#
# It is NOT:
#   - a drug ranking
#   - causal proof
#   - a calibrated intervention strength
#
# It identifies concrete molecular relations that could
# instantiate the dynamical feedback architecture.
# ============================================================


ROOT = Path(".")

F11A = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F11A_StateCausalGraph"
)

F10A = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F10A_AdultResiduals"
)

F5 = (
    ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F5_Mechanistic"
    / "Corrected"
)

OUT = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F11F_MolecularLeverage"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 1. Load
# ============================================================

edges = pd.read_csv(
    F11A
    / "F11A_state_molecular_edges.csv"
)

membership = pd.read_csv(
    F11A
    / "F11A_state_membership.csv"
)

residuals = pd.read_csv(
    F10A
    / "F10A_all_adult_residuals_ranked.csv"
)

pairs = pd.read_csv(
    F5
    / "pair_robustness.csv"
)


# ============================================================
# 2. Normalize symbols
# ============================================================

for frame, cols in [

    (
        edges,
        [
            "Source",
            "Target",
        ],
    ),

    (
        membership,
        [
            "Gene",
        ],
    ),

    (
        residuals,
        [
            "Gene",
        ],
    ),

    (
        pairs,
        [
            "Source",
            "Target",
        ],
    ),
]:

    for col in cols:

        frame[col] = (
            frame[col]
            .astype(str)
            .str.upper()
            .str.strip()
        )


# ============================================================
# 3. State membership lookup
# ============================================================

state_lookup = (
    membership
    .drop_duplicates(
        "Gene"
    )
    .set_index(
        "Gene"
    )[
        "State"
    ]
    .to_dict()
)


# ============================================================
# 4. Residual lookup
# ============================================================

RESIDUAL_COLUMNS = [

    "Adult_residual_z",

    "Residual_class",

    "E2F_relationship",

    "RbTKO_LFC",

    "E2F13_LFC",

    "CAG_LFC",

    "MAP2_LFC",
]


available_residual_columns = [

    col
    for col in RESIDUAL_COLUMNS
    if col in residuals.columns
]


res = (
    residuals[
        [
            "Gene",
            *available_residual_columns,
        ]
    ]
    .drop_duplicates(
        "Gene"
    )
    .set_index(
        "Gene"
    )
)


# ============================================================
# 5. Switch relations supported by F11E
# ============================================================

SWITCH_EDGES = {

    "U_TO_D_GATE":
        (
            "U_BACKBONE",
            "D_DEATH",
            1,
        ),

    "D_TO_U_RETURN":
        (
            "D_DEATH",
            "U_BACKBONE",
            1,
        ),

    "D_SELF_REINFORCEMENT":
        (
            "D_DEATH",
            "D_DEATH",
            1,
        ),

    "U_SELF_REINFORCEMENT":
        (
            "U_BACKBONE",
            "U_BACKBONE",
            1,
        ),
}


candidate_frames = []


for role, (
    source_state,
    target_state,
    required_effect,
) in SWITCH_EDGES.items():

    x = edges[
        (
            edges[
                "Source_state"
            ]
            ==
            source_state
        )
        &
        (
            edges[
                "Target_state"
            ]
            ==
            target_state
        )
        &
        (
            pd.to_numeric(
                edges[
                    "Effect"
                ],
                errors="coerce",
            )
            ==
            required_effect
        )
    ].copy()


    x[
        "Switch_role"
    ] = role


    candidate_frames.append(
        x
    )


candidates = pd.concat(
    candidate_frames,
    ignore_index=True,
)


if candidates.empty:

    raise RuntimeError(
        "No signed molecular edges matched the "
        "state-level switch architecture."
    )


# ============================================================
# 6. Evidence layer count
# ============================================================

def count_layers(value):

    if pd.isna(value):
        return 0

    parts = {

        item.strip()

        for item in str(
            value
        ).split(";")

        if item.strip()
    }

    return len(
        parts
    )


candidates[
    "Evidence_layer_count"
] = candidates[
    "Layers"
].apply(
    count_layers
)


# ============================================================
# 7. Add F10A evidence for both endpoints
# ============================================================

for prefix, gene_col in [

    (
        "Source",
        "Source",
    ),

    (
        "Target",
        "Target",
    ),
]:

    temp = (
        res.copy()
        .add_prefix(
            f"{prefix}_"
        )
    )


    candidates = candidates.merge(

        temp,

        left_on=gene_col,

        right_index=True,

        how="left",
    )


# ============================================================
# 8. Exact pair robustness where available
#
# pair_robustness was designed around selected F5
# source-target families, so many feedback edges will not have
# an exact entry. Absence here is NOT absence of support.
# ============================================================

pair_small = (
    pairs[
        [
            "Source",
            "Target",
            "Testable",
            "Support_tier",
            "Minimum_external_any_setting",
            "Shortest_edges_at_min_external",
        ]
    ]
    .drop_duplicates(
        [
            "Source",
            "Target",
        ]
    )
)


candidates = candidates.merge(

    pair_small,

    on=[
        "Source",
        "Target",
    ],

    how="left",
)


# ============================================================
# 9. Adult-context component
#
# Different switch roles should use residual information
# differently.
#
# U -> D:
#     Does the death target show adult-specific amplification?
#
# D -> U:
#     Does the death source show adult-specific amplification?
#
# D -> D:
#     Are both death nodes adult-amplified?
#
# U -> U:
#     This is conserved backbone circuitry, so residual
#     amplification is NOT required.
# ============================================================

def adult_context_score(row):

    role = row[
        "Switch_role"
    ]


    source_z = row.get(
        "Source_Adult_residual_z",
        np.nan,
    )

    target_z = row.get(
        "Target_Adult_residual_z",
        np.nan,
    )


    source_z = (
        float(source_z)
        if pd.notna(source_z)
        else 0.0
    )

    target_z = (
        float(target_z)
        if pd.notna(target_z)
        else 0.0
    )


    if role == "U_TO_D_GATE":

        return max(
            target_z,
            0.0,
        )


    if role == "D_TO_U_RETURN":

        return max(
            source_z,
            0.0,
        )


    if role == "D_SELF_REINFORCEMENT":

        return (
            max(
                source_z,
                0.0,
            )
            +
            max(
                target_z,
                0.0,
            )
        ) / 2


    if role == "U_SELF_REINFORCEMENT":

        return 0.0


    return 0.0


candidates[
    "Adult_context_score"
] = candidates.apply(
    adult_context_score,
    axis=1,
)


# ============================================================
# 10. Transparent exploratory evidence score
#
# Network term:
#
#   log(1 + records)
# + 0.5 log(1 + path usage)
# + 0.5 log(1 + pair usage)
# + 0.25 * evidence layers
#
# Adult-context term:
#
# + 0.75 * relevant positive residual z
#
# This score is only a ranking heuristic.
# ============================================================

for col in [

    "Supporting_records",

    "Representative_path_usage",

    "Representative_pair_usage",
]:

    candidates[col] = pd.to_numeric(
        candidates[col],
        errors="coerce",
    ).fillna(0)


candidates[
    "Network_support_score"
] = (

    np.log1p(
        candidates[
            "Supporting_records"
        ]
    )

    +

    0.5
    *
    np.log1p(
        candidates[
            "Representative_path_usage"
        ]
    )

    +

    0.5
    *
    np.log1p(
        candidates[
            "Representative_pair_usage"
        ]
    )

    +

    0.25
    *
    candidates[
        "Evidence_layer_count"
    ]
)


candidates[
    "Exploratory_priority_score"
] = (

    candidates[
        "Network_support_score"
    ]

    +

    0.75
    *
    candidates[
        "Adult_context_score"
    ]
)


# ============================================================
# 11. Human-readable relation
# ============================================================

candidates[
    "Relation"
] = (

    candidates[
        "Source"
    ]

    +
    " -> "

    +
    candidates[
        "Target"
    ]
)


# ============================================================
# 12. Priority tier
#
# This is intentionally descriptive rather than statistical.
# ============================================================

def priority_tier(row):

    strong_network = (

        row[
            "Supporting_records"
        ]
        >=
        2

        or

        row[
            "Evidence_layer_count"
        ]
        >=
        2
    )


    adult_supported = (

        row[
            "Adult_context_score"
        ]
        >=
        2
    )


    if (
        strong_network
        and
        adult_supported
    ):

        return "A"


    if (
        strong_network
        or
        adult_supported
    ):

        return "B"


    return "C"


candidates[
    "Exploratory_priority_tier"
] = candidates.apply(
    priority_tier,
    axis=1,
)


# ============================================================
# 13. Sort + export
# ============================================================

candidates = (
    candidates.sort_values(

        [
            "Switch_role",
            "Exploratory_priority_score",
        ],

        ascending=[
            True,
            False,
        ],
    )
    .reset_index(
        drop=True
    )
)


candidates.to_csv(

    OUT
    / "F11F_switch_molecular_edges.csv",

    index=False,
)


# ============================================================
# 14. Node-level leverage
#
# Count how often each molecular node participates in one or
# more of the switch relations.
# ============================================================

node_rows = []


for _, row in candidates.iterrows():

    for endpoint in [

        "Source",
        "Target",
    ]:

        gene = row[
            endpoint
        ]


        node_rows.append({

            "Gene":
                gene,

            "State":
                state_lookup.get(
                    gene,
                    np.nan,
                ),

            "Endpoint":
                endpoint,

            "Switch_role":
                row[
                    "Switch_role"
                ],

            "Partner":
                (
                    row[
                        "Target"
                    ]
                    if endpoint
                    ==
                    "Source"
                    else
                    row[
                        "Source"
                    ]
                ),

            "Edge_priority_score":
                row[
                    "Exploratory_priority_score"
                ],

            "Supporting_records":
                row[
                    "Supporting_records"
                ],

            "Representative_pair_usage":
                row[
                    "Representative_pair_usage"
                ],
        })


node_long = pd.DataFrame(
    node_rows
)


node_summary_rows = []


for gene, sub in node_long.groupby(
    "Gene"
):

    state = (
        sub[
            "State"
        ]
        .dropna()
    )


    state = (
        state.iloc[0]
        if len(state)
        else np.nan
    )


    residual_z = np.nan

    residual_class = np.nan


    if gene in res.index:

        if (
            "Adult_residual_z"
            in res.columns
        ):

            residual_z = res.loc[
                gene,
                "Adult_residual_z"
            ]


        if (
            "Residual_class"
            in res.columns
        ):

            residual_class = res.loc[
                gene,
                "Residual_class"
            ]


    roles = sorted(
        sub[
            "Switch_role"
        ]
        .unique()
    )


    partners = sorted(
        sub[
            "Partner"
        ]
        .unique()
    )


    node_summary_rows.append({

        "Gene":
            gene,

        "State":
            state,

        "Switch_role_count":
            len(
                roles
            ),

        "Switch_roles":
            ";".join(
                roles
            ),

        "Distinct_partners":
            len(
                partners
            ),

        "Partners":
            ";".join(
                partners
            ),

        "Incident_switch_edges":
            len(
                sub
            ),

        "Total_edge_priority":
            sub[
                "Edge_priority_score"
            ].sum(),

        "Maximum_edge_priority":
            sub[
                "Edge_priority_score"
            ].max(),

        "Total_support_records":
            sub[
                "Supporting_records"
            ].sum(),

        "Total_pair_usage":
            sub[
                "Representative_pair_usage"
            ].sum(),

        "Adult_residual_z":
            residual_z,

        "Residual_class":
            residual_class,
    })


node_summary = pd.DataFrame(
    node_summary_rows
)


node_summary = (
    node_summary.sort_values(

        [
            "Switch_role_count",
            "Total_edge_priority",
            "Distinct_partners",
        ],

        ascending=False,
    )
    .reset_index(
        drop=True
    )
)


node_summary.to_csv(

    OUT
    / "F11F_switch_node_priorities.csv",

    index=False,
)


# ============================================================
# 15. Role summaries
# ============================================================

role_summary = (

    candidates.groupby(
        "Switch_role"
    )

    .agg(

        Molecular_edges=(
            "Relation",
            "count",
        ),

        Distinct_sources=(
            "Source",
            "nunique",
        ),

        Distinct_targets=(
            "Target",
            "nunique",
        ),

        Median_network_support=(
            "Network_support_score",
            "median",
        ),

        Median_adult_context=(
            "Adult_context_score",
            "median",
        ),

        Maximum_priority=(
            "Exploratory_priority_score",
            "max",
        ),
    )

    .reset_index()
)


role_summary.to_csv(

    OUT
    / "F11F_switch_role_summary.csv",

    index=False,
)


# ============================================================
# 16. Console
# ============================================================

print()
print("=" * 110)
print("F11F — MOLECULAR LEVERAGE MAP")
print("=" * 110)


print()
print("SWITCH-ROLE SUMMARY")
print("-" * 110)


print(
    role_summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.3f}",
    )
)


# ============================================================
# 17. Top molecular edges by dynamical role
# ============================================================

SHOW = [

    "Relation",

    "Exploratory_priority_tier",

    "Supporting_records",

    "Evidence_layer_count",

    "Representative_path_usage",

    "Representative_pair_usage",

    "Source_Adult_residual_z",

    "Target_Adult_residual_z",

    "Adult_context_score",

    "Exploratory_priority_score",
]


for role in SWITCH_EDGES:

    print()
    print()
    print(role)
    print("-" * 110)


    x = candidates[
        candidates[
            "Switch_role"
        ]
        ==
        role
    ]


    if x.empty:

        print(
            "No matching signed molecular edges."
        )

        continue


    show_cols = [

        col
        for col in SHOW
        if col in x.columns
    ]


    print(
        x[
            show_cols
        ]
        .head(20)
        .to_string(
            index=False,
            float_format=lambda y: f"{y:.3f}",
        )
    )


# ============================================================
# 18. Node priorities
# ============================================================

print()
print()
print("TOP SWITCH NODES")
print("-" * 110)


print(
    node_summary
    .head(30)
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.3f}",
    )
)


# ============================================================
# 19. Adult-amplified death nodes participating in switch
# ============================================================

adult_death = node_summary[
    (
        node_summary[
            "State"
        ]
        ==
        "D_DEATH"
    )
    &
    (
        pd.to_numeric(
            node_summary[
                "Adult_residual_z"
            ],
            errors="coerce",
        )
        >=
        2
    )
]


print()
print()
print("ADULT-AMPLIFIED DEATH NODES IN THE SWITCH")
print("-" * 110)


if adult_death.empty:

    print(
        "None."
    )

else:

    print(
        adult_death.to_string(
            index=False,
            float_format=lambda x: f"{x:.3f}",
        )
    )


# ============================================================
# 20. Final note
# ============================================================

print()
print("Outputs written to:")
print(OUT)

print()
print(
    "[NEXT] Inspect which concrete molecules instantiate "
    "U->D, D->U and D->D. Then build perturbation hypotheses "
    "around molecular targets rather than abstract state edges."
)