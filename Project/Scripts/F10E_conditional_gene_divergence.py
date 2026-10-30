#!/usr/bin/env python

"""
F10E_conditional_gene_divergence.py
===================================

Gene-level conditional divergence analysis.

PURPOSE
-------
F10D showed that broad program-level adult-specific enrichment weakens
substantially after controlling for:

    1. adult Rb-TKO response direction,
    2. |adult Rb-TKO log2FC|,
    3. adult baseline expression,
    4. CAG-control baseline expression,
    5. MAP2-control baseline expression.

F10E now asks:

    Which INDIVIDUAL genes remain unusually adult-divergent relative
    to genes with similar response magnitude and baseline detectability?

This is primarily a ROBUSTNESS / RANKING analysis.

It should NOT be interpreted as establishing causal adult-specific genes.

METHOD
------
For each gene:

    - restrict potential comparators to genes with the same adult
      Rb-TKO direction;

    - find the nearest neighbors in the four-dimensional matching space:

          |adult Rb-TKO LFC|
          adult baseline expression
          CAG-control expression percentile
          MAP2-control expression percentile

    - compare the gene's original F10A adult residual z-score with
      the residual distribution among those matched neighbors.

The procedure is repeated at:

    k = 25
    k = 50
    k = 100

A gene is considered "robustly high" only if its conditional excess
remains positive across all matching neighborhoods.

IMPORTANT
---------
The local empirical P values are exploratory because:

    - neighborhoods overlap,
    - genes are not independent,
    - matching neighborhoods are data-dependent.

The main output is therefore the stability of the conditional ranking,
not a claim that any individual gene is statistically proven causal.

INPUT
-----
Data/Processed/Dynamics/F10D_CrossContextMatchedNull/
    F10D_analysis_universe.csv

Data/Processed/Dynamics/F10B_ResidualPrograms/
    F10B_program_gene_contributions.csv

OUTPUT
------
Data/Processed/Dynamics/F10E_ConditionalGeneDivergence/

    F10E_all_genes.csv
    F10E_top_positive.csv
    F10E_top_negative.csv
    F10E_candidate_audit.csv
    F10E_summary.txt
"""


from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# DEFAULT PATHS
# ============================================================

DEFAULT_UNIVERSE = (
    "Data/Processed/Dynamics/"
    "F10D_CrossContextMatchedNull/"
    "F10D_analysis_universe.csv"
)

DEFAULT_PROGRAMS = (
    "Data/Processed/Dynamics/"
    "F10B_ResidualPrograms/"
    "F10B_program_gene_contributions.csv"
)

DEFAULT_OUTDIR = (
    "Data/Processed/Dynamics/"
    "F10E_ConditionalGeneDivergence"
)


# ============================================================
# MATCHING VARIABLES
# ============================================================

MATCH_COLUMNS = [
    "Match_z_abs_lfc",
    "Match_z_adult_expr",
    "Match_z_cag_rank",
    "Match_z_map2_rank",
]


# ============================================================
# CANDIDATE AUDIT PANEL
# ============================================================

CANDIDATE_GENES = [
    # Current MYC-centered hypothesis
    "MYC",
    "DDIT3",
    "MAP3K5",

    # Death receptor / apoptosis
    "FAS",
    "TNFRSF10B",
    "PMAIP1",
    "BBC3",

    # Innate / regulated death
    "ZBP1",
    "PYCARD",
    "MLKL",
    "EIF2AK2",

    # NR4A
    "NR4A1",
    "NR4A2",
    "NR4A3",

    # Classical p53 axis
    "TP53",
    "CDKN2A",
    "MDM2",
    "BAX",
    "BCL2",
    "CASP3",

    # Backbone / E2F
    "E2F1",
    "E2F2",
    "CDK1",
    "CDK2",
    "CHEK1",
    "CHEK2",
    "CLSPN",
    "TOPBP1",
    "BRCA1",
    "H2AFX",
]


# ============================================================
# UTILITIES
# ============================================================

def require_columns(
    df: pd.DataFrame,
    columns: list[str],
    label: str,
) -> None:

    missing = [
        column
        for column in columns
        if column not in df.columns
    ]

    if missing:

        raise ValueError(
            f"\n{label} is missing required columns:\n"
            +
            "\n".join(
                f"  - {column}"
                for column in missing
            )
        )


def bh_fdr(
    pvalues: np.ndarray,
) -> np.ndarray:

    pvalues = np.asarray(
        pvalues,
        dtype=float,
    )

    n = len(
        pvalues
    )

    order = np.argsort(
        pvalues
    )

    ranked = (
        pvalues[
            order
        ]
    )

    adjusted = (
        ranked
        *
        n
        /
        np.arange(
            1,
            n + 1,
        )
    )

    adjusted = (
        np.minimum.accumulate(
            adjusted[::-1]
        )[::-1]
    )

    adjusted = np.clip(
        adjusted,
        0.0,
        1.0,
    )

    output = np.empty(
        n,
        dtype=float,
    )

    output[
        order
    ] = adjusted

    return output


def robust_scale(
    values: np.ndarray,
) -> float:

    values = np.asarray(
        values,
        dtype=float,
    )

    median = np.median(
        values
    )

    mad = (
        1.4826
        *
        np.median(
            np.abs(
                values
                -
                median
            )
        )
    )

    if (
        not np.isfinite(
            mad
        )
        or mad <= 1e-12
    ):

        mad = np.std(
            values,
            ddof=1,
        )

    if (
        not np.isfinite(
            mad
        )
        or mad <= 1e-12
    ):

        mad = 1.0

    return float(
        mad
    )


# ============================================================
# PROGRAM MEMBERSHIP
# ============================================================

def build_program_membership(
    programs: pd.DataFrame,
) -> dict[str, str]:

    require_columns(
        programs,
        [
            "Program",
            "Gene",
        ],
        "F10B program table",
    )

    membership = (
        programs
        .dropna(
            subset=[
                "Program",
                "Gene",
            ]
        )
        .groupby(
            "Gene"
        )[
            "Program"
        ]
        .apply(
            lambda x:
            ";".join(
                sorted(
                    set(
                        x.astype(str)
                    )
                )
            )
        )
        .to_dict()
    )

    return membership


# ============================================================
# CONDITIONAL NEIGHBOR ANALYSIS
# ============================================================

def analyze_one_k(
    universe: pd.DataFrame,
    k: int,
) -> pd.DataFrame:

    genes = (
        universe[
            "Gene"
        ]
        .astype(str)
        .to_numpy()
    )

    directions = (
        universe[
            "RbTKO_Direction"
        ]
        .astype(str)
        .to_numpy()
    )

    residuals = (
        universe[
            "Adult_residual_z"
        ]
        .to_numpy(
            dtype=float
        )
    )

    X = (
        universe[
            MATCH_COLUMNS
        ]
        .to_numpy(
            dtype=float
        )
    )

    rows = []

    for index in range(
        len(
            universe
        )
    ):

        eligible = (
            directions
            ==
            directions[
                index
            ]
        )

        eligible[
            index
        ] = False

        candidate_indices = np.where(
            eligible
        )[0]

        if len(
            candidate_indices
        ) == 0:

            continue

        delta = (
            X[
                candidate_indices
            ]
            -
            X[
                index
            ]
        )

        distance_squared = (
            np.sum(
                delta * delta,
                axis=1,
            )
        )

        kk = min(
            k,
            len(
                candidate_indices
            ),
        )

        nearest_local = (
            np.argpartition(
                distance_squared,
                kk - 1,
            )[:kk]
        )

        nearest_local = (
            nearest_local[
                np.argsort(
                    distance_squared[
                        nearest_local
                    ]
                )
            ]
        )

        neighbors = (
            candidate_indices[
                nearest_local
            ]
        )

        neighbor_residuals = (
            residuals[
                neighbors
            ]
        )

        observed = float(
            residuals[
                index
            ]
        )

        local_median = float(
            np.median(
                neighbor_residuals
            )
        )

        local_mean = float(
            np.mean(
                neighbor_residuals
            )
        )

        local_scale = robust_scale(
            neighbor_residuals
        )

        conditional_excess = (
            observed
            -
            local_median
        )

        conditional_z = (
            conditional_excess
            /
            local_scale
        )

        # ----------------------------------------------------
        # Percentile inside matched neighborhood
        # ----------------------------------------------------

        percentile = (
            100.0
            *
            (
                0.5
                +
                np.count_nonzero(
                    neighbor_residuals
                    <
                    observed
                )
            )
            /
            (
                len(
                    neighbor_residuals
                )
                +
                1
            )
        )

        # ----------------------------------------------------
        # Exploratory local P values
        # ----------------------------------------------------

        p_two_sided = (
            1
            +
            np.count_nonzero(
                np.abs(
                    neighbor_residuals
                    -
                    local_median
                )
                >=
                abs(
                    observed
                    -
                    local_median
                )
            )
        ) / (
            len(
                neighbor_residuals
            )
            +
            1
        )

        p_high = (
            1
            +
            np.count_nonzero(
                neighbor_residuals
                >=
                observed
            )
        ) / (
            len(
                neighbor_residuals
            )
            +
            1
        )

        p_low = (
            1
            +
            np.count_nonzero(
                neighbor_residuals
                <=
                observed
            )
        ) / (
            len(
                neighbor_residuals
            )
            +
            1
        )

        mean_neighbor_distance = float(
            np.mean(
                np.sqrt(
                    distance_squared[
                        nearest_local
                    ]
                )
            )
        )

        max_neighbor_distance = float(
            np.max(
                np.sqrt(
                    distance_squared[
                        nearest_local
                    ]
                )
            )
        )

        rows.append(
            {
                "Gene":
                    genes[
                        index
                    ],

                "K":
                    k,

                "Adult_residual_z":
                    observed,

                "Local_null_median":
                    local_median,

                "Local_null_mean":
                    local_mean,

                "Local_null_scale":
                    local_scale,

                "Conditional_excess":
                    conditional_excess,

                "Conditional_z":
                    conditional_z,

                "Local_percentile":
                    percentile,

                "Local_P_two_sided":
                    p_two_sided,

                "Local_P_high":
                    p_high,

                "Local_P_low":
                    p_low,

                "Mean_neighbor_distance":
                    mean_neighbor_distance,

                "Max_neighbor_distance":
                    max_neighbor_distance,
            }
        )

    results = pd.DataFrame(
        rows
    )

    results[
        "Local_FDR_two_sided"
    ] = bh_fdr(
        results[
            "Local_P_two_sided"
        ]
        .to_numpy()
    )

    return results


# ============================================================
# COMBINE K VALUES
# ============================================================

def combine_k_results(
    universe: pd.DataFrame,
    all_k_results: pd.DataFrame,
    k_values: list[int],
    membership: dict[str, str],
) -> pd.DataFrame:

    base_columns = [
        "Gene",
        "RbTKO_Direction",
        "RbTKO_LFC",
        "Adult_residual",
        "Adult_residual_z",
        "Adult_control_mean",
        "CAG_control_mean",
        "MAP2_control_mean",
        "CAG_control_rank",
        "MAP2_control_rank",
    ]

    available_base_columns = [
        column
        for column in base_columns
        if column in universe.columns
    ]

    combined = (
        universe[
            available_base_columns
        ]
        .copy()
    )

    combined[
        "Programs"
    ] = (
        combined[
            "Gene"
        ]
        .astype(str)
        .map(
            membership
        )
        .fillna(
            ""
        )
    )

    for k in k_values:

        subset = (
            all_k_results[
                all_k_results[
                    "K"
                ]
                ==
                k
            ]
            .copy()
        )

        keep = [
            "Gene",
            "Local_null_median",
            "Conditional_excess",
            "Conditional_z",
            "Local_percentile",
            "Local_P_two_sided",
            "Local_FDR_two_sided",
            "Mean_neighbor_distance",
        ]

        subset = subset[
            keep
        ]

        subset = subset.rename(
            columns={
                "Local_null_median":
                    f"Local_null_median_k{k}",

                "Conditional_excess":
                    f"Conditional_excess_k{k}",

                "Conditional_z":
                    f"Conditional_z_k{k}",

                "Local_percentile":
                    f"Local_percentile_k{k}",

                "Local_P_two_sided":
                    f"Local_P_two_sided_k{k}",

                "Local_FDR_two_sided":
                    f"Local_FDR_two_sided_k{k}",

                "Mean_neighbor_distance":
                    f"Mean_neighbor_distance_k{k}",
            }
        )

        combined = combined.merge(
            subset,
            on="Gene",
            how="left",
            validate="one_to_one",
        )

    conditional_z_columns = [
        f"Conditional_z_k{k}"
        for k in k_values
    ]

    percentile_columns = [
        f"Local_percentile_k{k}"
        for k in k_values
    ]

    excess_columns = [
        f"Conditional_excess_k{k}"
        for k in k_values
    ]

    combined[
        "Median_conditional_z"
    ] = (
        combined[
            conditional_z_columns
        ]
        .median(
            axis=1
        )
    )

    combined[
        "Min_conditional_z"
    ] = (
        combined[
            conditional_z_columns
        ]
        .min(
            axis=1
        )
    )

    combined[
        "Max_conditional_z"
    ] = (
        combined[
            conditional_z_columns
        ]
        .max(
            axis=1
        )
    )

    combined[
        "Median_conditional_excess"
    ] = (
        combined[
            excess_columns
        ]
        .median(
            axis=1
        )
    )

    combined[
        "Min_local_percentile"
    ] = (
        combined[
            percentile_columns
        ]
        .min(
            axis=1
        )
    )

    combined[
        "Max_local_percentile"
    ] = (
        combined[
            percentile_columns
        ]
        .max(
            axis=1
        )
    )

    combined[
        "Positive_across_all_k"
    ] = (
        combined[
            conditional_z_columns
        ]
        .gt(
            0
        )
        .all(
            axis=1
        )
    )

    combined[
        "Negative_across_all_k"
    ] = (
        combined[
            conditional_z_columns
        ]
        .lt(
            0
        )
        .all(
            axis=1
        )
    )

    # --------------------------------------------------------
    # Exploratory robustness flags.
    #
    # These are descriptive thresholds, NOT formal
    # significance definitions.
    # --------------------------------------------------------

    combined[
        "Robust_high_90pct"
    ] = (
        combined[
            "Positive_across_all_k"
        ]
        &
        (
            combined[
                "Min_local_percentile"
            ]
            >=
            90.0
        )
    )

    combined[
        "Robust_high_z1"
    ] = (
        combined[
            "Positive_across_all_k"
        ]
        &
        (
            combined[
                "Min_conditional_z"
            ]
            >=
            1.0
        )
    )

    combined[
        "Robust_low_10pct"
    ] = (
        combined[
            "Negative_across_all_k"
        ]
        &
        (
            combined[
                "Max_local_percentile"
            ]
            <=
            10.0
        )
    )

    combined[
        "Absolute_robust_score"
    ] = (
        combined[
            "Median_conditional_z"
        ]
        .abs()
    )

    combined = (
        combined
        .sort_values(
            "Absolute_robust_score",
            ascending=False,
        )
        .reset_index(
            drop=True
        )
    )

    combined[
        "Robust_rank"
    ] = (
        np.arange(
            1,
            len(
                combined
            )
            +
            1,
        )
    )

    return combined


# ============================================================
# CANDIDATE AUDIT
# ============================================================

def build_candidate_audit(
    combined: pd.DataFrame,
) -> pd.DataFrame:

    rows = []

    indexed = (
        combined
        .set_index(
            "Gene"
        )
    )

    for gene in CANDIDATE_GENES:

        if gene in indexed.index:

            row = (
                indexed
                .loc[
                    gene
                ]
                .to_dict()
            )

            row[
                "Gene"
            ] = gene

            row[
                "Present_in_F10_universe"
            ] = True

        else:

            row = {
                "Gene":
                    gene,

                "Present_in_F10_universe":
                    False,
            }

        rows.append(
            row
        )

    audit = pd.DataFrame(
        rows
    )

    first_columns = [
        "Gene",
        "Present_in_F10_universe",
    ]

    remaining_columns = [
        column
        for column in audit.columns
        if column
        not in first_columns
    ]

    audit = audit[
        first_columns
        +
        remaining_columns
    ]

    return audit


# ============================================================
# SUMMARY
# ============================================================

def write_summary(
    combined: pd.DataFrame,
    candidate_audit: pd.DataFrame,
    k_values: list[int],
    output_path: Path,
) -> None:

    lines = []

    lines.append(
        "F10E — CONDITIONAL GENE-LEVEL DIVERGENCE"
    )

    lines.append(
        "=" * 48
    )

    lines.append("")

    lines.append(
        "Purpose:"
    )

    lines.append(
        "Identify genes whose adult-vs-cortical residual "
        "remains unusual relative to genes matched on adult "
        "effect size and baseline expression/detectability."
    )

    lines.append("")

    lines.append(
        "Matching neighborhoods:"
    )

    for k in k_values:

        lines.append(
            f"  - k = {k}"
        )

    lines.append("")

    lines.append(
        "IMPORTANT:"
    )

    lines.append(
        "This is a conditional ranking / robustness analysis."
    )

    lines.append(
        "It does not establish causality or biological "
        "state membership."
    )

    lines.append("")

    lines.append(
        "ROBUST POSITIVE GENES"
    )

    lines.append(
        "-" * 48
    )

    positive = (
        combined[
            combined[
                "Robust_high_90pct"
            ]
        ]
        .sort_values(
            "Median_conditional_z",
            ascending=False,
        )
        .head(
            30
        )
    )

    if len(
        positive
    ) == 0:

        lines.append(
            "None met the exploratory robust-high criterion."
        )

    else:

        for _, row in (
            positive
            .iterrows()
        ):

            lines.append(
                f"{row['Gene']}: "
                f"original_z="
                f"{row['Adult_residual_z']:.3f}, "
                f"conditional_z_median="
                f"{row['Median_conditional_z']:.3f}, "
                f"min_percentile="
                f"{row['Min_local_percentile']:.1f}"
            )

    lines.append("")

    lines.append(
        "ROBUST NEGATIVE GENES"
    )

    lines.append(
        "-" * 48
    )

    negative = (
        combined[
            combined[
                "Robust_low_10pct"
            ]
        ]
        .sort_values(
            "Median_conditional_z"
        )
        .head(
            30
        )
    )

    if len(
        negative
    ) == 0:

        lines.append(
            "None met the exploratory robust-low criterion."
        )

    else:

        for _, row in (
            negative
            .iterrows()
        ):

            lines.append(
                f"{row['Gene']}: "
                f"original_z="
                f"{row['Adult_residual_z']:.3f}, "
                f"conditional_z_median="
                f"{row['Median_conditional_z']:.3f}, "
                f"max_percentile="
                f"{row['Max_local_percentile']:.1f}"
            )

    lines.append("")

    lines.append(
        "PRE-SPECIFIED CANDIDATE AUDIT"
    )

    lines.append(
        "-" * 48
    )

    for _, row in (
        candidate_audit
        .iterrows()
    ):

        gene = row[
            "Gene"
        ]

        present = row[
            "Present_in_F10_universe"
        ]

        if not present:

            lines.append(
                f"{gene}: not present in complete-context universe"
            )

        else:

            lines.append(
                f"{gene}: "
                f"original_z="
                f"{row['Adult_residual_z']:.3f}, "
                f"conditional_z_median="
                f"{row['Median_conditional_z']:.3f}, "
                f"rank="
                f"{int(row['Robust_rank'])}"
            )

    lines.append("")

    lines.append(
        "INTERPRETATION GUARDRAILS"
    )

    lines.append(
        "-" * 48
    )

    lines.append(
        "- A high conditional rank means a gene remains unusual "
        "relative to phenotypically/expression-matched genes."
    )

    lines.append(
        "- It does NOT prove the gene causes the adult phenotype."
    )

    lines.append(
        "- Failure of MYC or another prior candidate to remain "
        "highly ranked should weaken, not be ignored by, the "
        "subsequent mechanistic model."
    )

    lines.append(
        "- TP53 absence from this transcriptomic ranking would not "
        "exclude p53 pathway activity, because p53 is extensively "
        "regulated post-translationally."
    )

    output_path.write_text(
        "\n".join(
            lines
        ),
        encoding="utf-8",
    )


# ============================================================
# ARGUMENTS
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--universe",
        default=DEFAULT_UNIVERSE,
    )

    parser.add_argument(
        "--programs",
        default=DEFAULT_PROGRAMS,
    )

    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
    )

    parser.add_argument(
        "--k",
        default="25,50,100",
        help=(
            "Comma-separated nearest-neighbor "
            "matching sizes."
        ),
    )

    return parser.parse_args()


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_args()

    universe_path = Path(
        args.universe
    )

    programs_path = Path(
        args.programs
    )

    output_directory = Path(
        args.outdir
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    universe = pd.read_csv(
        universe_path
    )

    programs = pd.read_csv(
        programs_path
    )

    required_columns = [
        "Gene",
        "RbTKO_Direction",
        "RbTKO_LFC",
        "Adult_residual_z",
        *MATCH_COLUMNS,
    ]

    require_columns(
        universe,
        required_columns,
        "F10D analysis universe",
    )

    k_values = [
        int(
            value.strip()
        )

        for value
        in args.k.split(
            ","
        )

        if value.strip()
    ]

    membership = (
        build_program_membership(
            programs
        )
    )

    # ========================================================
    # RUN EACH NEIGHBORHOOD SIZE
    # ========================================================

    all_results = []

    for k in k_values:

        print(
            f"Running conditional divergence: "
            f"k = {k}"
        )

        result = analyze_one_k(
            universe,
            k,
        )

        all_results.append(
            result
        )

    all_k_results = pd.concat(
        all_results,
        ignore_index=True,
    )

    # ========================================================
    # COMBINE ACROSS K
    # ========================================================

    combined = combine_k_results(
        universe=universe,
        all_k_results=all_k_results,
        k_values=k_values,
        membership=membership,
    )

    combined.to_csv(
        output_directory
        /
        "F10E_all_genes.csv",
        index=False,
    )

    # ========================================================
    # TOP POSITIVE / NEGATIVE
    # ========================================================

    top_positive = (
        combined
        .sort_values(
            "Median_conditional_z",
            ascending=False,
        )
        .head(
            50
        )
    )

    top_positive.to_csv(
        output_directory
        /
        "F10E_top_positive.csv",
        index=False,
    )

    top_negative = (
        combined
        .sort_values(
            "Median_conditional_z",
            ascending=True,
        )
        .head(
            50
        )
    )

    top_negative.to_csv(
        output_directory
        /
        "F10E_top_negative.csv",
        index=False,
    )

    # ========================================================
    # CANDIDATE AUDIT
    # ========================================================

    candidate_audit = (
        build_candidate_audit(
            combined
        )
    )

    candidate_audit.to_csv(
        output_directory
        /
        "F10E_candidate_audit.csv",
        index=False,
    )

    # ========================================================
    # SUMMARY
    # ========================================================

    write_summary(
        combined=combined,
        candidate_audit=candidate_audit,
        k_values=k_values,
        output_path=(
            output_directory
            /
            "F10E_summary.txt"
        ),
    )

    # ========================================================
    # CONSOLE OUTPUT
    # ========================================================

    print()
    print(
        "Top conditionally HIGH genes"
    )

    print(
        "----------------------------"
    )

    display_columns = [
        "Gene",
        "RbTKO_LFC",
        "Adult_residual_z",
        "Median_conditional_z",
        "Min_local_percentile",
        "Robust_high_90pct",
        "Programs",
    ]

    print(
        top_positive[
            display_columns
        ]
        .head(
            25
        )
        .to_string(
            index=False
        )
    )

    print()
    print(
        "Pre-specified candidate audit"
    )

    print(
        "-----------------------------"
    )

    candidate_display_columns = [
        column
        for column in [
            "Gene",
            "Present_in_F10_universe",
            "RbTKO_LFC",
            "Adult_residual_z",
            "Median_conditional_z",
            "Min_local_percentile",
            "Robust_rank",
            "Programs",
        ]
        if column
        in candidate_audit.columns
    ]

    print(
        candidate_audit[
            candidate_display_columns
        ]
        .to_string(
            index=False
        )
    )

    print()

    print(
        "Outputs written to:"
    )

    print(
        output_directory
    )

    print()

    print(
        "IMPORTANT:"
    )

    print(
        "Do not interpret Robust_high_90pct as a "
        "formal significance call. It is an exploratory "
        "stability flag across matched neighborhoods."
    )


if __name__ == "__main__":
    main()