from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# F10B — Adult residual program enrichment
#
# Question:
#
# Are adult-vs-conserved residuals coherently concentrated
# within mechanistic fate programs?
#
# This is an exploratory internal validation step.
#
# IMPORTANT:
# Several programs were originally motivated by the Fong
# analysis itself, so enrichment here is NOT independent
# external validation. The purpose is to determine whether
# candidate state variables are internally coherent enough
# to justify dynamical modeling.
# ============================================================


ROOT = Path(".")

INFILE = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F10A_AdultResiduals"
    / "F10A_all_adult_residuals_ranked.csv"
)

OUT = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F10B_ResidualPrograms"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 1. Programs
#
# Kept aligned with programs/candidates already defined before
# seeing the F10A residual ranking.
# ============================================================

PROGRAMS = {

    # --------------------------------------------------------
    # Negative-control-ish conserved backbone
    # --------------------------------------------------------

    "CONSERVED_BACKBONE": [
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
    ],


    # --------------------------------------------------------
    # Fate / stress programs
    # --------------------------------------------------------

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

    "NR4A_STATE": [
        "NR4A1",
        "NR4A2",
        "NR4A3",
    ],


    # --------------------------------------------------------
    # Pre-existing F7-style composite pro-death hypothesis
    # --------------------------------------------------------

    "FATE_DEATH_COMPOSITE": [
        "MYC",
        "DDIT3",
        "FAS",
        "TNFRSF10B",
        "PMAIP1",
        "BBC3",
        "ZBP1",
        "PYCARD",
        "MLKL",
        "EIF2AK2",
    ],
}


# ============================================================
# 2. Load residual universe
# ============================================================

df = pd.read_csv(
    INFILE
)

df["Gene"] = (
    df["Gene"]
    .astype(str)
    .str.upper()
    .str.strip()
)

df = (
    df.dropna(
        subset=[
            "Adult_residual_z",
        ]
    )
    .drop_duplicates(
        "Gene"
    )
    .copy()
)


universe_genes = set(
    df["Gene"]
)

residual_lookup = (
    df.set_index(
        "Gene"
    )[
        "Adult_residual_z"
    ]
)


print()
print("=" * 105)
print("F10B — ADULT RESIDUAL PROGRAM ANALYSIS")
print("=" * 105)

print()
print(
    "Residual universe:",
    len(df),
    "genes"
)


# ============================================================
# 3. Background properties
# ============================================================

background_mean = (
    df["Adult_residual_z"]
    .mean()
)

background_median = (
    df["Adult_residual_z"]
    .median()
)

background_positive = (
    df["Adult_residual_z"]
    > 0
).mean()


print()
print(
    f"Background mean residual z   = "
    f"{background_mean:+.4f}"
)

print(
    f"Background median residual z = "
    f"{background_median:+.4f}"
)

print(
    f"Background positive fraction = "
    f"{background_positive:.3f}"
)


# ============================================================
# 4. Permutation setup
# ============================================================

rng = np.random.default_rng(
    20261022
)

N_PERM = 50000

background_values = (
    df["Adult_residual_z"]
    .astype(float)
    .values
)


# ============================================================
# 5. Program analysis
# ============================================================

results = []

gene_rows = []

jackknife_rows = []


for program, requested_genes in PROGRAMS.items():

    requested_genes = [
        g.upper()
        for g in requested_genes
    ]

    genes = [
        g
        for g in requested_genes
        if g in universe_genes
    ]


    if len(genes) < 2:

        print(
            f"[WARNING] {program}: "
            f"only {len(genes)} genes available."
        )

        continue


    values = (
        residual_lookup.loc[
            genes
        ]
        .astype(float)
        .values
    )


    observed_mean = (
        np.mean(values)
    )

    observed_median = (
        np.median(values)
    )

    positive_fraction = (
        np.mean(
            values > 0
        )
    )

    strong_positive_fraction = (
        np.mean(
            values >= 2
        )
    )

    strong_negative_fraction = (
        np.mean(
            values <= -2
        )
    )


    # ========================================================
    # Random gene-set null
    #
    # Sample equal-sized sets from the same 1206-gene universe.
    # ========================================================

    null_means = np.empty(
        N_PERM,
        dtype=float,
    )


    n = len(values)


    for i in range(N_PERM):

        sample = rng.choice(
            background_values,
            size=n,
            replace=False,
        )

        null_means[i] = (
            sample.mean()
        )


    # Two-sided empirical P

    centered_observed = (
        observed_mean
        -
        background_mean
    )

    centered_null = (
        null_means
        -
        background_mean
    )


    permutation_p_two_sided = (
        1
        +
        np.sum(
            np.abs(
                centered_null
            )
            >=
            abs(
                centered_observed
            )
        )
    ) / (
        N_PERM
        +
        1
    )


    # Directional empirical P

    if (
        observed_mean
        >=
        background_mean
    ):

        permutation_p_directional = (
            1
            +
            np.sum(
                null_means
                >=
                observed_mean
            )
        ) / (
            N_PERM
            +
            1
        )

        direction = (
            "Adult_amplified"
        )

    else:

        permutation_p_directional = (
            1
            +
            np.sum(
                null_means
                <=
                observed_mean
            )
        ) / (
            N_PERM
            +
            1
        )

        direction = (
            "Adult_attenuated"
        )


    # Standardized program effect relative to null distribution

    null_sd = np.std(
        null_means,
        ddof=1,
    )

    program_z = (
        (
            observed_mean
            -
            np.mean(
                null_means
            )
        )
        /
        null_sd
    )


    results.append({

        "Program":
            program,

        "Requested_genes":
            len(
                requested_genes
            ),

        "Genes_available":
            len(
                genes
            ),

        "Genes":
            ";".join(
                genes
            ),

        "Mean_residual_z":
            observed_mean,

        "Median_residual_z":
            observed_median,

        "Positive_fraction":
            positive_fraction,

        "Strong_positive_fraction":
            strong_positive_fraction,

        "Strong_negative_fraction":
            strong_negative_fraction,

        "Direction":
            direction,

        "Null_mean":
            np.mean(
                null_means
            ),

        "Null_SD":
            null_sd,

        "Program_enrichment_z":
            program_z,

        "Permutation_P_directional":
            permutation_p_directional,

        "Permutation_P_two_sided":
            permutation_p_two_sided,
    })


    # ========================================================
    # Gene-level contributions
    # ========================================================

    for gene in genes:

        gene_rows.append({

            "Program":
                program,

            "Gene":
                gene,

            "Adult_residual_z":
                residual_lookup.loc[
                    gene
                ],
        })


    # ========================================================
    # Leave-one-gene-out robustness
    # ========================================================

    if len(genes) >= 3:

        for gene in genes:

            remaining = [
                g
                for g in genes
                if g != gene
            ]

            loo_mean = (
                residual_lookup.loc[
                    remaining
                ]
                .mean()
            )

            jackknife_rows.append({

                "Program":
                    program,

                "Dropped_gene":
                    gene,

                "LOO_mean_residual_z":
                    loo_mean,

                "Same_direction":
                    (
                        np.sign(
                            loo_mean
                            -
                            background_mean
                        )
                        ==
                        np.sign(
                            observed_mean
                            -
                            background_mean
                        )
                    ),
            })


# ============================================================
# 6. Results dataframe
# ============================================================

results = pd.DataFrame(
    results
)


# ============================================================
# 7. BH adjustment
# ============================================================

def bh_adjust(values):

    p = np.asarray(
        values,
        dtype=float,
    )

    order = np.argsort(
        p
    )

    ranked = p[
        order
    ]

    n = len(
        ranked
    )

    adjusted = (
        ranked
        *
        n
        /
        np.arange(
            1,
            n + 1
        )
    )

    adjusted = np.minimum.accumulate(
        adjusted[::-1]
    )[::-1]

    adjusted = np.minimum(
        adjusted,
        1
    )

    result = np.empty(
        n
    )

    result[
        order
    ] = adjusted

    return result


results[
    "Permutation_FDR_directional"
] = bh_adjust(
    results[
        "Permutation_P_directional"
    ].values
)


results[
    "Permutation_FDR_two_sided"
] = bh_adjust(
    results[
        "Permutation_P_two_sided"
    ].values
)


# ============================================================
# 8. Jackknife summaries
# ============================================================

jackknife = pd.DataFrame(
    jackknife_rows
)


jack_summaries = []


for program in (
    jackknife[
        "Program"
    ]
    .unique()
):

    x = jackknife[
        jackknife[
            "Program"
        ]
        ==
        program
    ]

    jack_summaries.append({

        "Program":
            program,

        "N_leave_one_out":
            len(x),

        "LOO_min_mean":
            x[
                "LOO_mean_residual_z"
            ].min(),

        "LOO_max_mean":
            x[
                "LOO_mean_residual_z"
            ].max(),

        "LOO_mean":
            x[
                "LOO_mean_residual_z"
            ].mean(),

        "Fraction_same_direction":
            x[
                "Same_direction"
            ].mean(),
    })


jack_summary = pd.DataFrame(
    jack_summaries
)


results = results.merge(
    jack_summary,
    on="Program",
    how="left",
)


# ============================================================
# 9. Export
# ============================================================

results = results.sort_values(
    "Permutation_P_two_sided"
)


results.to_csv(
    OUT
    / "F10B_program_residual_enrichment.csv",
    index=False,
)


pd.DataFrame(
    gene_rows
).to_csv(
    OUT
    / "F10B_program_gene_contributions.csv",
    index=False,
)


jackknife.to_csv(
    OUT
    / "F10B_program_leave_one_out.csv",
    index=False,
)


# ============================================================
# 10. Console output
# ============================================================

print()
print()
print("PROGRAM RESIDUAL ENRICHMENT")
print("-" * 105)


display_cols = [
    "Program",
    "Genes_available",
    "Mean_residual_z",
    "Median_residual_z",
    "Program_enrichment_z",
    "Direction",
    "Permutation_P_two_sided",
    "Permutation_FDR_two_sided",
    "Fraction_same_direction",
]


print(
    results[
        display_cols
    ]
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.4g}",
    )
)


# ============================================================
# 11. Gene contributions
# ============================================================

print()
print()
print("GENE CONTRIBUTIONS")
print("-" * 105)


gene_table = pd.DataFrame(
    gene_rows
)


for program in results[
    "Program"
]:

    x = (
        gene_table[
            gene_table[
                "Program"
            ]
            ==
            program
        ]
        .sort_values(
            "Adult_residual_z",
            ascending=False,
        )
    )

    print()
    print(
        program
    )

    print(
        x[
            [
                "Gene",
                "Adult_residual_z",
            ]
        ]
        .to_string(
            index=False,
            float_format=lambda y: f"{y:+.3f}",
        )
    )


# ============================================================
# 12. Dynamical-state recommendation
# ============================================================

print()
print()
print("CANDIDATE DYNAMICAL STATE VARIABLES")
print("-" * 105)


supported = results[
    (
        results[
            "Permutation_FDR_two_sided"
        ]
        <
        0.10
    )
    &
    (
        results[
            "Fraction_same_direction"
        ]
        .fillna(1)
        >=
        0.75
    )
].copy()


if supported.empty:

    print(
        "No program passes the provisional "
        "FDR/robustness criterion."
    )

else:

    print(
        supported[
            [
                "Program",
                "Mean_residual_z",
                "Direction",
                "Permutation_FDR_two_sided",
                "Fraction_same_direction",
            ]
        ]
        .to_string(
            index=False,
            float_format=lambda x: f"{x:.4g}",
        )
    )


print()
print("Outputs written to:")
print(OUT)

print()
print(
    "[NEXT] Promote only coherent residual programs—not "
    "individual dramatic genes—into the first adult fate "
    "dynamical system."
)