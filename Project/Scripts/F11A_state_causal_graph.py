from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx


# ============================================================
# F11A — State-level causal graph reconstruction
#
# Goal:
#   Collapse the corrected F5 molecular network into a small
#   signed state graph suitable for dynamical modeling.
#
# States:
#
#   U = conserved Rb/E2F-cycle/replication backbone
#   D = apoptotic / death-stress program
#   I = innate / regulated-death reinforcement
#   C = adult-amplified compensatory trophic program
#   T = adult-suppressed trophic/quiescence-support program
#
# We then search for:
#
#   - positive/negative state interactions
#   - self-reinforcement
#   - feedback cycles
#   - candidate positive feedback capable of multistability
#
# IMPORTANT:
#   Absence of a feedback loop in this graph is NOT proof that
#   the biological system lacks it. F5 itself was constructed
#   from a restricted source->target mechanistic search.
# ============================================================


ROOT = Path(".")

F5 = (
    ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F5_Mechanistic"
    / "Corrected"
)

F10A = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F10A_AdultResiduals"
)

OUT = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F11A_StateCausalGraph"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 1. Load network + residual data
# ============================================================

edges = pd.read_csv(
    F5 / "representative_union_edges.csv"
)

pairs = pd.read_csv(
    F5 / "pair_robustness.csv"
)

residuals = pd.read_csv(
    F10A / "F10A_all_adult_residuals_ranked.csv"
)


# Standardize gene symbols

for col in [
    "Source",
    "Target",
]:

    edges[col] = (
        edges[col]
        .astype(str)
        .str.upper()
        .str.strip()
    )


for col in [
    "Source",
    "Target",
]:

    pairs[col] = (
        pairs[col]
        .astype(str)
        .str.upper()
        .str.strip()
    )


residuals["Gene"] = (
    residuals["Gene"]
    .astype(str)
    .str.upper()
    .str.strip()
)

residuals = (
    residuals
    .drop_duplicates("Gene")
    .set_index("Gene")
)


# Numeric conversion

for col in [
    "Effect",
    "Supporting_records",
    "Representative_path_usage",
    "Representative_pair_usage",
]:

    edges[col] = pd.to_numeric(
        edges[col],
        errors="coerce",
    )


# ============================================================
# 2. Define empirically motivated states
# ============================================================


# ------------------------------------------------------------
# U — conserved Rb/E2F-cycle/replication backbone
# ------------------------------------------------------------

BACKBONE = [
    "E2F1",
    "E2F2",
    "CDK1",
    "CDK2",
    "PCNA",
    "CLSPN",
    "TOPBP1",
    "CHEK1",
    "CHEK2",
    "H2AFX",
    "BRCA1",
    "AURKB",
    "BIRC5",
    "BUB1",
    "MELK",
    "KIF14",
    "MYBL2",
]


# ------------------------------------------------------------
# D — apoptosis / death-stress
#
# Kept separate from innate/regulated death below.
# ------------------------------------------------------------

DEATH = [
    "MYC",
    "DDIT3",
    "FAS",
    "TNFRSF10B",
    "PMAIP1",
    "BBC3",
    "ATF5",
    "TRIB3",
    "MAP3K5",
    "TNFRSF22",
]


# ------------------------------------------------------------
# I — innate / regulated-death reinforcement
# ------------------------------------------------------------

INNATE_REGULATED = [
    "ZBP1",
    "PYCARD",
    "MLKL",
    "EIF2AK2",
    "TNF",
    "DTX3L",
    "IFIT2",
    "IFIT3",
    "IGTP",
    "IIGP1",
    "IRGM1",
    "PARP9",
    "TRIM25",
]


# ------------------------------------------------------------
# Trophic genes
#
# F10B showed that the old "trophic" module is NOT coherent.
#
# We therefore split it empirically:
#
#   C = adult-amplified compensatory trophic response
#   T = adult-suppressed trophic/quiescence-support response
#
# Only |residual z| >= 1 is used.
# ------------------------------------------------------------

TROPHIC_CANDIDATES = [
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
]


COMPENSATORY = []

SUPPRESSED_TROPHIC = []


for gene in TROPHIC_CANDIDATES:

    if gene not in residuals.index:
        continue

    z = residuals.loc[
        gene,
        "Adult_residual_z"
    ]

    if pd.isna(z):
        continue

    if z >= 1:

        COMPENSATORY.append(
            gene
        )

    elif z <= -1:

        SUPPRESSED_TROPHIC.append(
            gene
        )


STATES = {

    "U_BACKBONE":
        BACKBONE,

    "D_DEATH":
        DEATH,

    "I_INNATE_REGULATED":
        INNATE_REGULATED,

    "C_COMPENSATORY":
        COMPENSATORY,

    "T_TROPHIC_SUPPORT":
        SUPPRESSED_TROPHIC,
}


# ============================================================
# 3. Verify exclusive state membership
# ============================================================

gene_to_state = {}

overlaps = []


for state, genes in STATES.items():

    for gene in genes:

        gene = gene.upper()

        if gene in gene_to_state:

            overlaps.append(
                (
                    gene,
                    gene_to_state[gene],
                    state,
                )
            )

        else:

            gene_to_state[
                gene
            ] = state


if overlaps:

    print()
    print(
        "[WARNING] Overlapping state membership detected:"
    )

    for item in overlaps:
        print(item)

    raise RuntimeError(
        "States must be exclusive for F11A."
    )


# ============================================================
# 4. Export state membership
# ============================================================

membership_rows = []


for state, genes in STATES.items():

    for gene in genes:

        residual_z = np.nan

        if gene in residuals.index:

            residual_z = (
                residuals.loc[
                    gene,
                    "Adult_residual_z"
                ]
            )


        membership_rows.append({

            "State":
                state,

            "Gene":
                gene,

            "Adult_residual_z":
                residual_z,
        })


membership = pd.DataFrame(
    membership_rows
)


membership.to_csv(
    OUT
    / "F11A_state_membership.csv",
    index=False,
)


# ============================================================
# 5. Deduplicate molecular edges across F5 settings
#
# Same molecular relation may appear under E4_X0, E4_X1, etc.
# We do not want to count that as independent biological
# evidence.
# ============================================================

def join_unique(values):

    out = []

    for value in values:

        if pd.isna(value):
            continue

        for item in str(value).split(";"):

            item = item.strip()

            if (
                item
                and
                item not in out
            ):
                out.append(item)

    return ";".join(
        sorted(out)
    )


edge_unique = (
    edges.groupby(
        [
            "Source",
            "Target",
            "Effect",
        ],
        as_index=False,
    )
    .agg({

        "Setting":
            join_unique,

        "Layers":
            join_unique,

        "Supporting_records":
            "max",

        "Representative_path_usage":
            "max",

        "Representative_pair_usage":
            "max",
    })
)


# ============================================================
# 6. Map molecular edges onto states
# ============================================================

edge_unique[
    "Source_state"
] = edge_unique[
    "Source"
].map(
    gene_to_state
)


edge_unique[
    "Target_state"
] = edge_unique[
    "Target"
].map(
    gene_to_state
)


state_molecular = edge_unique[
    edge_unique[
        "Source_state"
    ].notna()
    &
    edge_unique[
        "Target_state"
    ].notna()
].copy()


state_molecular[
    "Effect_label_F11"
] = np.select(

    [
        state_molecular[
            "Effect"
        ] > 0,

        state_molecular[
            "Effect"
        ] < 0,
    ],

    [
        "positive",
        "negative",
    ],

    default="unsigned",
)


state_molecular.to_csv(
    OUT
    / "F11A_state_molecular_edges.csv",
    index=False,
)


# ============================================================
# 7. Aggregate molecular evidence by state -> state
# ============================================================

state_rows = []


for (
    source_state,
    target_state
), sub in state_molecular.groupby(
    [
        "Source_state",
        "Target_state",
    ]
):

    positive = sub[
        sub["Effect"] > 0
    ]

    negative = sub[
        sub["Effect"] < 0
    ]

    unsigned = sub[
        sub["Effect"] == 0
    ]


    pos_support = positive[
        "Supporting_records"
    ].fillna(0).sum()

    neg_support = negative[
        "Supporting_records"
    ].fillna(0).sum()

    unsigned_support = unsigned[
        "Supporting_records"
    ].fillna(0).sum()


    signed_support = (
        pos_support
        +
        neg_support
    )


    if signed_support == 0:

        dominant_sign = 0

        dominant_label = (
            "unsigned_only"
        )

        sign_purity = np.nan

    else:

        if pos_support > neg_support:

            dominant_sign = 1

            dominant_label = (
                "positive"
            )

        elif neg_support > pos_support:

            dominant_sign = -1

            dominant_label = (
                "negative"
            )

        else:

            dominant_sign = 0

            dominant_label = (
                "conflicting"
            )


        sign_purity = (
            max(
                pos_support,
                neg_support
            )
            /
            signed_support
        )


    molecular_pairs = (
        sub[
            [
                "Source",
                "Target",
            ]
        ]
        .drop_duplicates()
    )


    state_rows.append({

        "Source_state":
            source_state,

        "Target_state":
            target_state,

        "Molecular_edges":
            len(
                molecular_pairs
            ),

        "Positive_edges":
            len(
                positive
            ),

        "Negative_edges":
            len(
                negative
            ),

        "Unsigned_edges":
            len(
                unsigned
            ),

        "Positive_support_records":
            pos_support,

        "Negative_support_records":
            neg_support,

        "Unsigned_support_records":
            unsigned_support,

        "Representative_path_usage":
            sub[
                "Representative_path_usage"
            ].fillna(0).sum(),

        "Representative_pair_usage":
            sub[
                "Representative_pair_usage"
            ].fillna(0).sum(),

        "Distinct_source_genes":
            sub[
                "Source"
            ].nunique(),

        "Distinct_target_genes":
            sub[
                "Target"
            ].nunique(),

        "Dominant_sign_numeric":
            dominant_sign,

        "Dominant_sign":
            dominant_label,

        "Sign_purity":
            sign_purity,

        "Molecular_relations":
            "; ".join(
                sorted(
                    (
                        sub["Source"]
                        +
                        "->"
                        +
                        sub["Target"]
                    )
                    .drop_duplicates()
                )
            ),
    })


state_summary = pd.DataFrame(
    state_rows
)


if not state_summary.empty:

    state_summary[
        "Sign_resolved"
    ] = (
        state_summary[
            "Dominant_sign_numeric"
        ]
        .isin(
            [-1, 1]
        )
        &
        (
            state_summary[
                "Sign_purity"
            ]
            >= 0.67
        )
    )


    state_summary = (
        state_summary.sort_values(
            [
                "Source_state",
                "Target_state",
            ]
        )
        .reset_index(
            drop=True
        )
    )


state_summary.to_csv(
    OUT
    / "F11A_state_edge_summary.csv",
    index=False,
)


# ============================================================
# 8. Map pair-robustness evidence onto states
#
# Pair robustness is SOURCE->TARGET reachability evidence.
# It is complementary to the direct edge aggregation above.
# ============================================================

pairs[
    "Source_state"
] = pairs[
    "Source"
].map(
    gene_to_state
)


pairs[
    "Target_state"
] = pairs[
    "Target"
].map(
    gene_to_state
)


state_pairs = pairs[
    pairs[
        "Source_state"
    ].notna()
    &
    pairs[
        "Target_state"
    ].notna()
].copy()


pair_rows = []


for (
    source_state,
    target_state
), sub in state_pairs.groupby(
    [
        "Source_state",
        "Target_state",
    ]
):

    tiers = (
        sub[
            "Support_tier"
        ]
        .value_counts()
        .to_dict()
    )


    pair_rows.append({

        "Source_state":
            source_state,

        "Target_state":
            target_state,

        "Pairs":
            len(sub),

        "Tier_A":
            tiers.get(
                "A",
                0
            ),

        "Tier_B":
            tiers.get(
                "B",
                0
            ),

        "Tier_C":
            tiers.get(
                "C",
                0
            ),

        "Tier_D":
            tiers.get(
                "D",
                0
            ),

        "Unreachable":
            tiers.get(
                "UNREACHABLE",
                0
            ),

        "Direct_pairs":
            (
                pd.to_numeric(
                    sub[
                        "Shortest_edges_at_min_external"
                    ],
                    errors="coerce",
                )
                ==
                1
            ).sum(),

        "Median_shortest_edges":
            pd.to_numeric(
                sub[
                    "Shortest_edges_at_min_external"
                ],
                errors="coerce",
            ).median(),

        "Source_genes":
            ";".join(
                sorted(
                    sub[
                        "Source"
                    ]
                    .unique()
                )
            ),

        "Target_genes":
            ";".join(
                sorted(
                    sub[
                        "Target"
                    ]
                    .unique()
                )
            ),
    })


pair_summary = pd.DataFrame(
    pair_rows
)


pair_summary.to_csv(
    OUT
    / "F11A_state_pair_robustness.csv",
    index=False,
)


# ============================================================
# 9. Build signed state graph
# ============================================================

G = nx.DiGraph()


for state in STATES:

    G.add_node(
        state
    )


if not state_summary.empty:

    for _, row in state_summary.iterrows():

        G.add_edge(

            row[
                "Source_state"
            ],

            row[
                "Target_state"
            ],

            sign=int(
                row[
                    "Dominant_sign_numeric"
                ]
            ),

            sign_label=row[
                "Dominant_sign"
            ],

            sign_purity=row[
                "Sign_purity"
            ],

            molecular_edges=int(
                row[
                    "Molecular_edges"
                ]
            ),
        )


# ============================================================
# 10. Search state-level feedback cycles
#
# Cycles of length <= 4.
#
# Cycle sign:
#
#   product(edge signs) = +1  -> positive feedback
#   product(edge signs) = -1  -> negative feedback
#   any unresolved edge       -> unknown
#
# Positive feedback is a necessary structural ingredient for
# multistability in many regulatory systems, although its
# presence alone does not prove bistability.
# ============================================================

cycle_rows = []


cycles = list(
    nx.simple_cycles(
        G
    )
)


# Canonicalize cycles to avoid representation issues

seen_cycles = set()


for cycle in cycles:

    if len(cycle) > 4:
        continue


    # Rotate so lexicographically smallest node is first.

    rotations = []

    for i in range(len(cycle)):

        rotations.append(
            tuple(
                cycle[i:]
                +
                cycle[:i]
            )
        )


    canonical = min(
        rotations
    )


    if canonical in seen_cycles:
        continue


    seen_cycles.add(
        canonical
    )


    signs = []

    purities = []

    edge_labels = []


    for i in range(
        len(cycle)
    ):

        source = cycle[i]

        target = cycle[
            (i + 1)
            %
            len(cycle)
        ]


        data = G.get_edge_data(
            source,
            target
        )


        sign = data[
            "sign"
        ]

        signs.append(
            sign
        )

        purities.append(
            data[
                "sign_purity"
            ]
        )

        edge_labels.append(
            f"{source}->{target}"
        )


    if all(
        s in [-1, 1]
        for s in signs
    ):

        cycle_sign_numeric = int(
            np.prod(
                signs
            )
        )


        cycle_sign = (
            "positive_feedback"
            if cycle_sign_numeric > 0
            else
            "negative_feedback"
        )

    else:

        cycle_sign_numeric = 0

        cycle_sign = (
            "sign_unresolved"
        )


    cycle_rows.append({

        "Cycle":
            " -> ".join(
                list(cycle)
                +
                [cycle[0]]
            ),

        "Length":
            len(cycle),

        "Cycle_sign":
            cycle_sign,

        "Cycle_sign_numeric":
            cycle_sign_numeric,

        "Minimum_edge_sign_purity":
            np.nanmin(
                purities
            )
            if any(
                pd.notna(x)
                for x in purities
            )
            else np.nan,

        "Edges":
            ";".join(
                edge_labels
            ),
    })


cycle_df = pd.DataFrame(
    cycle_rows
)


if not cycle_df.empty:

    cycle_df = (
        cycle_df.sort_values(
            [
                "Cycle_sign",
                "Length",
            ]
        )
    )


cycle_df.to_csv(
    OUT
    / "F11A_feedback_cycles.csv",
    index=False,
)


# ============================================================
# 11. Signed adjacency matrix
# ============================================================

state_names = list(
    STATES.keys()
)


adj = pd.DataFrame(
    0,
    index=state_names,
    columns=state_names,
    dtype=int,
)


if not state_summary.empty:

    for _, row in state_summary.iterrows():

        adj.loc[
            row[
                "Source_state"
            ],
            row[
                "Target_state"
            ],
        ] = int(
            row[
                "Dominant_sign_numeric"
            ]
        )


adj.to_csv(
    OUT
    / "F11A_signed_state_adjacency.csv"
)


# ============================================================
# 12. Console output
# ============================================================

print()
print("=" * 110)
print("F11A — STATE-LEVEL CAUSAL GRAPH")
print("=" * 110)


print()
print("STATE MEMBERSHIP")
print("-" * 110)


for state, genes in STATES.items():

    print(
        f"{state:24s} "
        f"n={len(genes):2d} : "
        +
        ", ".join(
            genes
        )
    )


print()
print()
print("DIRECT MOLECULAR STATE EDGES")
print("-" * 110)


if state_summary.empty:

    print(
        "No direct molecular edges connect the defined states."
    )

else:

    cols = [
        "Source_state",
        "Target_state",
        "Molecular_edges",
        "Positive_edges",
        "Negative_edges",
        "Unsigned_edges",
        "Dominant_sign",
        "Sign_purity",
        "Representative_pair_usage",
    ]


    print(
        state_summary[
            cols
        ]
        .to_string(
            index=False,
            float_format=lambda x: f"{x:.3f}",
        )
    )


print()
print()
print("PAIR-ROBUSTNESS STATE CONNECTIONS")
print("-" * 110)


if pair_summary.empty:

    print(
        "No pair-robustness rows connect the defined states."
    )

else:

    print(
        pair_summary.to_string(
            index=False,
            float_format=lambda x: f"{x:.3f}",
        )
    )


print()
print()
print("FEEDBACK CYCLES")
print("-" * 110)


if cycle_df.empty:

    print(
        "No state-level cycles were recovered from the "
        "current F5 representative-edge graph."
    )

else:

    print(
        cycle_df.to_string(
            index=False,
            float_format=lambda x: f"{x:.3f}",
        )
    )


# ============================================================
# 13. Dynamical interpretation
# ============================================================

positive_cycles = pd.DataFrame()

negative_cycles = pd.DataFrame()


if not cycle_df.empty:

    positive_cycles = cycle_df[
        cycle_df[
            "Cycle_sign"
        ]
        ==
        "positive_feedback"
    ]

    negative_cycles = cycle_df[
        cycle_df[
            "Cycle_sign"
        ]
        ==
        "negative_feedback"
    ]


print()
print()
print("DYNAMICAL INTERPRETATION")
print("-" * 110)


print(
    "Positive-feedback cycles:",
    len(
        positive_cycles
    )
)

print(
    "Negative-feedback cycles:",
    len(
        negative_cycles
    )
)


if len(
    positive_cycles
) > 0:

    print()
    print(
        "The current molecular graph contains signed "
        "positive-feedback architecture."
    )

    print(
        "This makes multistability structurally plausible, "
        "although it does NOT prove multiple attractors."
    )

else:

    print()
    print(
        "The current F5 path-restricted graph does not "
        "establish a signed positive-feedback cycle."
    )

    print(
        "Therefore bistability should NOT be inserted into "
        "the ODE model by assumption."
    )


if len(
    negative_cycles
) > 0:

    print(
        "Negative-feedback architecture is also present and "
        "may support stabilization/adaptation depending on "
        "nonlinearity and time scales."
    )


print()
print(
    "IMPORTANT: F5 was constructed around selected "
    "source-to-target paths."
)

print(
    "Failure to recover reverse feedback here is therefore "
    "a limitation of current network coverage, not evidence "
    "that the biological feedback does not exist."
)


# ============================================================
# 14. Summary file
# ============================================================

with open(
    OUT
    / "F11A_summary.txt",
    "w",
    encoding="utf-8",
) as f:

    f.write(
        "F11A — STATE-LEVEL CAUSAL GRAPH\n"
    )

    f.write(
        "================================\n\n"
    )

    f.write(
        f"State-mapped molecular edges: "
        f"{len(state_molecular)}\n"
    )

    f.write(
        f"State-level directed relations: "
        f"{len(state_summary)}\n"
    )

    f.write(
        f"Pair-robustness state relations: "
        f"{len(pair_summary)}\n"
    )

    f.write(
        f"Feedback cycles <=4 states: "
        f"{len(cycle_df)}\n"
    )

    f.write(
        f"Positive-feedback cycles: "
        f"{len(positive_cycles)}\n"
    )

    f.write(
        f"Negative-feedback cycles: "
        f"{len(negative_cycles)}\n\n"
    )

    f.write(
        "Interpretation:\n"
    )

    f.write(
        "- Positive feedback is structurally compatible with "
        "multistability but does not itself prove bistability.\n"
    )

    f.write(
        "- No feedback should be invented merely to create "
        "attractors.\n"
    )

    f.write(
        "- The F5 network is path-restricted, so missing "
        "feedback may reflect incomplete network coverage.\n"
    )


print()
print("Outputs written to:")
print(OUT)

print()
print(
    "[NEXT] Use the recovered signed state topology to define "
    "the minimal adult-fate ODE system. Then solve fixed points, "
    "classify Jacobian eigenvalues, and draw phase portraits."
)