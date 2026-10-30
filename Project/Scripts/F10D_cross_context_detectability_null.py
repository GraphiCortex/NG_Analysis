#!/usr/bin/env python

"""
F10D_cross_context_detectability_null.py
========================================

Cross-context detectability-aware robustness test for F10C.

PURPOSE
-------
F10C matched program genes against null genes with:

    1. the same adult Rb-TKO direction,
    2. similar |adult Rb-TKO log2FC|,
    3. similar adult baseline expression.

F10D adds two further matching variables:

    4. baseline expression in CAG-control cortical cells,
    5. baseline expression in MAP2-control cortical cells.

This tests the alternative explanation that an apparently
"adult-specific" program might simply consist of genes that are
poorly expressed / poorly detectable in the cortical datasets.

IMPORTANT
---------
Cortical expression is represented using within-Oshikawa expression
percentiles rather than raw intensity values. This avoids pretending
that the adult and cortical expression platforms are directly
comparable on the same numerical scale.

PRIMARY INFERENCE
-----------------
Two-sided matched-null P values with Benjamini-Hochberg FDR.

Directional P values are also reported, but are SECONDARY.

OUTPUTS
-------
Data/Processed/Dynamics/F10D_CrossContextMatchedNull/

    F10D_analysis_universe.csv
    F10D_cortical_sample_mapping.csv
    F10D_matching_diagnostics.csv
    F10D_primary.csv
    F10D_sensitivity.csv
    F10D_program_detectability_audit.csv
    F10D_null_distributions.csv.gz
    F10D_summary.txt

    figures/
        F10D_matched_null_summary.png
        F10D_INNATE_STRESS_null.png
        F10D_FATE_DEATH_COMPOSITE_null.png
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# DEFAULT PATHS
# ============================================================

DEFAULT_RESIDUALS = (
    "Data/Processed/Dynamics/F10A_AdultResiduals/"
    "F10A_all_adult_residuals_ranked.csv"
)

DEFAULT_PROGRAMS = (
    "Data/Processed/Dynamics/F10B_ResidualPrograms/"
    "F10B_program_gene_contributions.csv"
)

DEFAULT_FONG = (
    "Data/Raw/microarray/Fong2022/"
    "Fong 2022 L2FC1 microarray data.xlsx"
)

DEFAULT_OSH_EXPR = (
    "Data/Processed/Dynamics/F8A_Oshikawa/"
    "oshikawa_gene_level_expression.csv"
)

DEFAULT_DESIGN = (
    "Data/Processed/Dynamics/F8A_Oshikawa/"
    "design_matrix.csv"
)

DEFAULT_OUTDIR = (
    "Data/Processed/Dynamics/"
    "F10D_CrossContextMatchedNull"
)


# ============================================================
# GENERAL UTILITIES
# ============================================================

def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction.
    """

    p = np.asarray(pvalues, dtype=float)

    n = len(p)

    order = np.argsort(p)

    ranked = p[order]

    q = (
        ranked
        * n
        / np.arange(1, n + 1)
    )

    q = np.minimum.accumulate(
        q[::-1]
    )[::-1]

    q = np.clip(
        q,
        0.0,
        1.0,
    )

    output = np.empty(
        n,
        dtype=float,
    )

    output[order] = q

    return output


def robust_z(
    x: pd.Series,
) -> pd.Series:
    """
    Robust z-score using median and MAD.
    """

    arr = pd.to_numeric(
        x,
        errors="coerce",
    ).to_numpy(
        dtype=float
    )

    median = np.nanmedian(
        arr
    )

    mad = (
        1.4826
        *
        np.nanmedian(
            np.abs(
                arr - median
            )
        )
    )

    if (
        not np.isfinite(mad)
        or mad <= 1e-12
    ):

        mad = np.nanstd(
            arr,
            ddof=1,
        )

    if (
        not np.isfinite(mad)
        or mad <= 1e-12
    ):

        mad = 1.0

    return pd.Series(
        (
            arr - median
        )
        / mad,
        index=x.index,
    )


def require_columns(
    df: pd.DataFrame,
    columns: list[str],
    label: str,
) -> None:
    """
    Stop with a useful message if required columns are missing.
    """

    missing = [
        col
        for col in columns
        if col not in df.columns
    ]

    if missing:

        raise ValueError(
            f"\n{label} is missing required columns:\n"
            +
            "\n".join(
                f"  - {col}"
                for col in missing
            )
        )


def normalize_name(
    value: str,
) -> str:
    """
    Normalize a string for loose sample-name matching.
    """

    return re.sub(
        r"[^a-z0-9]+",
        "",
        str(value).lower(),
    )


# ============================================================
# READ / INTERPRET OSHIKAWA DESIGN MATRIX
# ============================================================

def read_design_matrix(
    path: Path,
) -> pd.DataFrame:

    design = pd.read_csv(
        path
    )

    first_col = design.columns[0]

    # R write.csv() often creates an unnamed row-name column.
    if (
        str(first_col).startswith(
            "Unnamed"
        )
        or normalize_name(
            first_col
        )
        in {
            "x",
            "sample",
            "sampleid",
        }
    ):

        design = design.rename(
            columns={
                first_col:
                "Sample"
            }
        )

    elif "Sample" not in design.columns:

        # Fallback:
        # if all columns after the first look numeric,
        # assume first column contains sample names.

        other = (
            design.iloc[:, 1:]
            .apply(
                pd.to_numeric,
                errors="coerce",
            )
        )

        numeric_fraction = (
            other
            .notna()
            .mean()
            .mean()
        )

        if numeric_fraction > 0.90:

            design = design.rename(
                columns={
                    first_col:
                    "Sample"
                }
            )

    if "Sample" not in design.columns:

        raise ValueError(
            "\nCould not identify the sample-name "
            "column in design_matrix.csv.\n\n"
            f"Columns found:\n{list(design.columns)}"
        )

    design["Sample"] = (
        design["Sample"]
        .astype(str)
    )

    return design


def infer_control_columns(
    design: pd.DataFrame,
) -> tuple[str, str]:
    """
    Infer CAG-control and MAP2-control design columns.

    Supports names such as:

        cCont
        CAG_control
        groupcCont

    and:

        mCont
        MAP2_control
        groupmCont
    """

    candidates = [
        col
        for col in design.columns
        if col != "Sample"
    ]

    normalized = {
        col:
        normalize_name(col)
        for col in candidates
    }

    def score_cag(
        col: str,
    ) -> int:

        name = normalized[col]

        score = 0

        if "ccont" in name:
            score += 20

        if (
            "cag" in name
            and (
                "cont" in name
                or "ctrl" in name
                or "control" in name
            )
        ):
            score += 20

        if (
            name.startswith("c")
            and (
                "cont" in name
                or "ctrl" in name
            )
        ):
            score += 10

        if (
            "tko" in name
            or name.endswith("ko")
        ):
            score -= 20

        return score


    def score_map2(
        col: str,
    ) -> int:

        name = normalized[col]

        score = 0

        if "mcont" in name:
            score += 20

        if (
            "map2" in name
            and (
                "cont" in name
                or "ctrl" in name
                or "control" in name
            )
        ):
            score += 20

        if (
            name.startswith("m")
            and (
                "cont" in name
                or "ctrl" in name
            )
        ):
            score += 10

        if (
            "tko" in name
            or name.endswith("ko")
        ):
            score -= 20

        return score


    cag_control = max(
        candidates,
        key=score_cag,
    )

    map2_control = max(
        candidates,
        key=score_map2,
    )

    if score_cag(
        cag_control
    ) <= 0:

        raise ValueError(
            "\nCould not infer the CAG-control "
            "column from design_matrix.csv.\n\n"
            f"Columns found:\n{candidates}"
        )

    if score_map2(
        map2_control
    ) <= 0:

        raise ValueError(
            "\nCould not infer the MAP2-control "
            "column from design_matrix.csv.\n\n"
            f"Columns found:\n{candidates}"
        )

    return (
        cag_control,
        map2_control,
    )


def samples_for_group(
    design: pd.DataFrame,
    group_column: str,
) -> list[str]:

    values = pd.to_numeric(
        design[group_column],
        errors="coerce",
    )

    samples = (
        design.loc[
            values > 0.5,
            "Sample",
        ]
        .astype(str)
        .tolist()
    )

    if len(samples) == 0:

        raise ValueError(
            f"\nNo samples were assigned to "
            f"{group_column!r}."
        )

    return samples


# ============================================================
# READ OSHIKAWA GENE-LEVEL EXPRESSION
# ============================================================

def infer_gene_column(
    expression: pd.DataFrame,
) -> str:

    preferred = [
        "Gene",
        "gene",
        "Symbol",
        "SYMBOL",
        "GeneSymbol",
        "gene_symbol",
        "Gene_Symbol",
    ]

    for col in preferred:

        if col in expression.columns:

            return col

    # Fallback:
    # first string/object column.

    for col in expression.columns:

        if expression[col].dtype == object:

            return col

    raise ValueError(
        "\nCould not identify a gene-symbol "
        "column in oshikawa_gene_level_expression.csv."
    )


def map_design_samples_to_expression(
    design_samples: list[str],
    expression_columns: list[str],
) -> dict[str, str]:

    mapping = {}

    normalized_expression = {
        col:
        normalize_name(col)
        for col in expression_columns
    }

    for sample in design_samples:

        # Exact match first.
        if sample in expression_columns:

            mapping[sample] = sample

            continue

        normalized_sample = (
            normalize_name(sample)
        )

        # Exact normalized match.
        normalized_matches = [
            col
            for col, norm_col
            in normalized_expression.items()
            if norm_col == normalized_sample
        ]

        if len(
            normalized_matches
        ) == 1:

            mapping[sample] = (
                normalized_matches[0]
            )

            continue

        # Containment fallback.
        containment_matches = [
            col
            for col, norm_col
            in normalized_expression.items()
            if (
                normalized_sample in norm_col
                or norm_col in normalized_sample
            )
        ]

        if len(
            containment_matches
        ) == 1:

            mapping[sample] = (
                containment_matches[0]
            )

            continue

        raise ValueError(
            "\nCould not map design sample "
            f"{sample!r} to exactly one expression column.\n\n"
            f"Potential matches:\n"
            f"{containment_matches}"
        )

    return mapping


def load_oshikawa_baseline_expression(
    expression_path: Path,
    design_path: Path,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
]:

    expression = pd.read_csv(
        expression_path
    )

    design = read_design_matrix(
        design_path
    )

    gene_column = infer_gene_column(
        expression
    )

    (
        cag_control_column,
        map2_control_column,
    ) = infer_control_columns(
        design
    )

    cag_control_samples = (
        samples_for_group(
            design,
            cag_control_column,
        )
    )

    map2_control_samples = (
        samples_for_group(
            design,
            map2_control_column,
        )
    )

    expression_columns = [
        col
        for col in expression.columns
        if col != gene_column
    ]

    all_control_samples = (
        cag_control_samples
        +
        map2_control_samples
    )

    mapping = (
        map_design_samples_to_expression(
            all_control_samples,
            expression_columns,
        )
    )

    cag_expression_columns = [
        mapping[sample]
        for sample
        in cag_control_samples
    ]

    map2_expression_columns = [
        mapping[sample]
        for sample
        in map2_control_samples
    ]

    # Force expression values to numeric.
    for col in (
        cag_expression_columns
        +
        map2_expression_columns
    ):

        expression[col] = (
            pd.to_numeric(
                expression[col],
                errors="coerce",
            )
        )

    baseline = pd.DataFrame(
        {
            "Gene":
                expression[
                    gene_column
                ].astype(str),

            "CAG_control_mean":
                expression[
                    cag_expression_columns
                ].mean(
                    axis=1
                ),

            "MAP2_control_mean":
                expression[
                    map2_expression_columns
                ].mean(
                    axis=1
                ),
        }
    )

    # If a gene appears multiple times,
    # collapse it conservatively by mean.
    baseline = (
        baseline
        .dropna(
            subset=[
                "Gene",
                "CAG_control_mean",
                "MAP2_control_mean",
            ]
        )
        .groupby(
            "Gene",
            as_index=False,
        )
        .mean(
            numeric_only=True
        )
    )

    # --------------------------------------------------------
    # IMPORTANT:
    # Use within-platform percentiles instead of raw intensity.
    # --------------------------------------------------------

    baseline[
        "CAG_control_rank"
    ] = (
        baseline[
            "CAG_control_mean"
        ]
        .rank(
            method="average",
            pct=True,
        )
    )

    baseline[
        "MAP2_control_rank"
    ] = (
        baseline[
            "MAP2_control_mean"
        ]
        .rank(
            method="average",
            pct=True,
        )
    )

    # Save explicit sample mapping for manual verification.
    sample_map_rows = []

    for sample in cag_control_samples:

        sample_map_rows.append(
            {
                "Design_sample":
                    sample,

                "Expression_column":
                    mapping[sample],

                "Group":
                    "CAG_control",

                "Design_group_column":
                    cag_control_column,
            }
        )

    for sample in map2_control_samples:

        sample_map_rows.append(
            {
                "Design_sample":
                    sample,

                "Expression_column":
                    mapping[sample],

                "Group":
                    "MAP2_control",

                "Design_group_column":
                    map2_control_column,
            }
        )

    sample_mapping = pd.DataFrame(
        sample_map_rows
    )

    return (
        baseline,
        sample_mapping,
    )


# ============================================================
# BUILD ANALYSIS UNIVERSE
# ============================================================

def load_analysis_universe(
    residual_path: Path,
    program_path: Path,
    fong_path: Path,
    fong_sheet: str,
    oshikawa_expression_path: Path,
    design_path: Path,
):

    residuals = pd.read_csv(
        residual_path
    )

    programs = pd.read_csv(
        program_path
    )

    require_columns(
        residuals,
        [
            "Gene",
            "Ensembl_RbTKO",
            "RbTKO_LFC",
            "RbTKO_Direction",
            "CAG_LFC",
            "MAP2_LFC",
            "Adult_residual_z",
        ],
        "F10A residual table",
    )

    require_columns(
        programs,
        [
            "Program",
            "Gene",
            "Adult_residual_z",
        ],
        "F10B program contribution table",
    )

    # ========================================================
    # Adult baseline expression
    # ========================================================

    fong = pd.read_excel(
        fong_path,
        sheet_name=fong_sheet,
    )

    adult_control_columns = [
        "RB_THC_1",
        "RB_THC_2",
        "RB_THC_3",
    ]

    require_columns(
        fong,
        [
            "Row.names",
            *adult_control_columns,
        ],
        "Fong expression table",
    )

    adult = fong[
        [
            "Row.names",
            *adult_control_columns,
        ]
    ].copy()

    adult = adult.rename(
        columns={
            "Row.names":
            "Ensembl_RbTKO"
        }
    )

    adult = (
        adult
        .dropna(
            subset=[
                "Ensembl_RbTKO"
            ]
        )
        .drop_duplicates(
            "Ensembl_RbTKO",
            keep="first",
        )
    )

    merged = residuals.merge(
        adult,
        on="Ensembl_RbTKO",
        how="left",
        validate="many_to_one",
    )

    merged[
        "Adult_control_mean"
    ] = (
        merged[
            adult_control_columns
        ]
        .mean(
            axis=1
        )
    )

    merged[
        "Adult_control_log1p"
    ] = np.log1p(
        merged[
            "Adult_control_mean"
        ]
    )

    merged[
        "Abs_RbTKO_LFC"
    ] = (
        merged[
            "RbTKO_LFC"
        ]
        .abs()
    )

    # ========================================================
    # Cortical baseline expression
    # ========================================================

    (
        cortical_baseline,
        sample_mapping,
    ) = load_oshikawa_baseline_expression(
        oshikawa_expression_path,
        design_path,
    )

    merged["Gene"] = (
        merged["Gene"]
        .astype(str)
    )

    cortical_baseline["Gene"] = (
        cortical_baseline["Gene"]
        .astype(str)
    )

    merged = merged.merge(
        cortical_baseline,
        on="Gene",
        how="left",
        validate="many_to_one",
    )

    required_complete = [
        "Gene",
        "RbTKO_LFC",
        "RbTKO_Direction",
        "CAG_LFC",
        "MAP2_LFC",
        "Adult_residual_z",
        "Adult_control_log1p",
        "CAG_control_rank",
        "MAP2_control_rank",
    ]

    before = len(
        merged
    )

    merged = (
        merged
        .dropna(
            subset=required_complete
        )
        .copy()
    )

    after = len(
        merged
    )

    print(
        f"Complete cross-context matching universe: "
        f"{after}/{before} F10 genes"
    )

    if after < before:

        print(
            f"WARNING: {before - after} genes were removed "
            "because one or more baseline-expression "
            "variables were unavailable."
        )

    # ========================================================
    # Standardize matching variables
    # ========================================================

    merged[
        "Match_z_abs_lfc"
    ] = robust_z(
        merged[
            "Abs_RbTKO_LFC"
        ]
    )

    merged[
        "Match_z_adult_expr"
    ] = robust_z(
        merged[
            "Adult_control_log1p"
        ]
    )

    merged[
        "Match_z_cag_rank"
    ] = robust_z(
        merged[
            "CAG_control_rank"
        ]
    )

    merged[
        "Match_z_map2_rank"
    ] = robust_z(
        merged[
            "MAP2_control_rank"
        ]
    )

    merged = merged.reset_index(
        drop=True
    )

    return (
        merged,
        programs,
        sample_mapping,
    )


# ============================================================
# MATCHING
# ============================================================

MATCH_COLUMNS = [
    "Match_z_abs_lfc",
    "Match_z_adult_expr",
    "Match_z_cag_rank",
    "Match_z_map2_rank",
]


def get_program_indices(
    universe: pd.DataFrame,
    programs: pd.DataFrame,
    program_name: str,
) -> np.ndarray:

    genes = (
        programs.loc[
            programs["Program"]
            ==
            program_name,
            "Gene",
        ]
        .dropna()
        .astype(str)
        .unique()
    )

    index_by_gene = {
        str(gene):
        i

        for i, gene
        in enumerate(
            universe[
                "Gene"
            ].astype(str)
        )
    }

    indices = [
        index_by_gene[gene]

        for gene in genes

        if gene
        in index_by_gene
    ]

    return np.asarray(
        indices,
        dtype=int,
    )


def build_neighbor_lists(
    universe: pd.DataFrame,
    program_indices: np.ndarray,
    k: int,
    exclude_all_program_genes: bool,
    all_program_genes: set[str],
) -> list[np.ndarray]:

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

    X = (
        universe[
            MATCH_COLUMNS
        ]
        .to_numpy(
            dtype=float
        )
    )

    current_program_genes = set(
        genes[
            program_indices
        ].tolist()
    )

    if exclude_all_program_genes:

        excluded_genes = (
            all_program_genes
        )

    else:

        excluded_genes = (
            current_program_genes
        )

    exclusion_mask = np.isin(
        genes,
        list(
            excluded_genes
        ),
    )

    neighbor_lists = []

    for idx in program_indices:

        eligible = (
            (
                directions
                ==
                directions[idx]
            )
            &
            (
                ~exclusion_mask
            )
        )

        candidates = np.where(
            eligible
        )[0]

        if len(
            candidates
        ) == 0:

            raise RuntimeError(
                f"No eligible matched-null candidates "
                f"for gene {genes[idx]}."
            )

        delta = (
            X[candidates]
            -
            X[idx]
        )

        distance_squared = np.sum(
            delta * delta,
            axis=1,
        )

        kk = min(
            k,
            len(
                candidates
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

        neighbor_lists.append(
            candidates[
                nearest_local
            ]
        )

    return neighbor_lists


# ============================================================
# MONTE CARLO MATCHED NULL
# ============================================================

def simulate_matched_null(
    universe: pd.DataFrame,
    program_indices: np.ndarray,
    neighbor_lists: list[np.ndarray],
    iterations: int,
    seed: int,
):

    rng = np.random.default_rng(
        seed
    )

    residual = (
        universe[
            "Adult_residual_z"
        ]
        .to_numpy(
            dtype=float
        )
    )

    abs_lfc = (
        universe[
            "Abs_RbTKO_LFC"
        ]
        .to_numpy(
            dtype=float
        )
    )

    adult_expression = (
        universe[
            "Adult_control_log1p"
        ]
        .to_numpy(
            dtype=float
        )
    )

    cag_rank = (
        universe[
            "CAG_control_rank"
        ]
        .to_numpy(
            dtype=float
        )
    )

    map2_rank = (
        universe[
            "MAP2_control_rank"
        ]
        .to_numpy(
            dtype=float
        )
    )

    n_program_genes = len(
        program_indices
    )

    null_residual = np.empty(
        iterations,
        dtype=float,
    )

    null_abs_lfc = np.empty(
        iterations,
        dtype=float,
    )

    null_adult_expression = np.empty(
        iterations,
        dtype=float,
    )

    null_cag_rank = np.empty(
        iterations,
        dtype=float,
    )

    null_map2_rank = np.empty(
        iterations,
        dtype=float,
    )

    for iteration in range(
        iterations
    ):

        used = set()

        selected = []

        # Randomize gene order to avoid systematic
        # first-choice advantage when candidate pools overlap.
        order = rng.permutation(
            n_program_genes
        )

        for local_index in order:

            candidates = (
                neighbor_lists[
                    local_index
                ]
            )

            available = [
                int(candidate)

                for candidate
                in candidates

                if int(candidate)
                not in used
            ]

            if not available:

                raise RuntimeError(
                    "\nUnique matching failed because "
                    "no unused candidate remained.\n"
                    "Increase --k and rerun."
                )

            chosen = int(
                available[
                    rng.integers(
                        0,
                        len(
                            available
                        ),
                    )
                ]
            )

            used.add(
                chosen
            )

            selected.append(
                chosen
            )

        selected = np.asarray(
            selected,
            dtype=int,
        )

        null_residual[
            iteration
        ] = (
            residual[
                selected
            ]
            .mean()
        )

        null_abs_lfc[
            iteration
        ] = (
            abs_lfc[
                selected
            ]
            .mean()
        )

        null_adult_expression[
            iteration
        ] = (
            adult_expression[
                selected
            ]
            .mean()
        )

        null_cag_rank[
            iteration
        ] = (
            cag_rank[
                selected
            ]
            .mean()
        )

        null_map2_rank[
            iteration
        ] = (
            map2_rank[
                selected
            ]
            .mean()
        )

    return {
        "residual":
            null_residual,

        "abs_lfc":
            null_abs_lfc,

        "adult_expression":
            null_adult_expression,

        "cag_rank":
            null_cag_rank,

        "map2_rank":
            null_map2_rank,
    }


# ============================================================
# SUMMARIZE ONE PROGRAM
# ============================================================

def summarize_program(
    universe: pd.DataFrame,
    program_name: str,
    program_indices: np.ndarray,
    null: dict[str, np.ndarray],
):

    observed_residual = float(
        universe.loc[
            program_indices,
            "Adult_residual_z",
        ]
        .mean()
    )

    null_residual = (
        null[
            "residual"
        ]
    )

    null_mean = float(
        null_residual.mean()
    )

    null_sd = float(
        null_residual.std(
            ddof=1
        )
    )

    if null_sd > 0:

        matched_effect_z = (
            observed_residual
            -
            null_mean
        ) / null_sd

    else:

        matched_effect_z = (
            np.nan
        )

    # --------------------------------------------------------
    # PRIMARY TEST:
    # two-sided relative to matched-null center
    # --------------------------------------------------------

    primary_p = (
        1
        +
        np.count_nonzero(
            np.abs(
                null_residual
                -
                null_mean
            )
            >=
            abs(
                observed_residual
                -
                null_mean
            )
        )
    ) / (
        len(
            null_residual
        )
        +
        1
    )

    # --------------------------------------------------------
    # SECONDARY directional test
    # --------------------------------------------------------

    if observed_residual >= null_mean:

        directional_p = (
            1
            +
            np.count_nonzero(
                null_residual
                >=
                observed_residual
            )
        ) / (
            len(
                null_residual
            )
            +
            1
        )

        direction = (
            "Above_matched_null"
        )

    else:

        directional_p = (
            1
            +
            np.count_nonzero(
                null_residual
                <=
                observed_residual
            )
        ) / (
            len(
                null_residual
            )
            +
            1
        )

        direction = (
            "Below_matched_null"
        )

    def observed_mean(
        column: str,
    ) -> float:

        return float(
            universe.loc[
                program_indices,
                column,
            ]
            .mean()
        )

    result = {
        "Program":
            program_name,

        "Genes_available":
            int(
                len(
                    program_indices
                )
            ),

        "Observed_mean_residual_z":
            observed_residual,

        "Matched_null_mean":
            null_mean,

        "Matched_null_SD":
            null_sd,

        "Matched_effect_z":
            matched_effect_z,

        "Matched_q025":
            float(
                np.quantile(
                    null_residual,
                    0.025,
                )
            ),

        "Matched_q975":
            float(
                np.quantile(
                    null_residual,
                    0.975,
                )
            ),

        "Matched_percentile":
            float(
                100.0
                *
                np.mean(
                    null_residual
                    <
                    observed_residual
                )
            ),

        "Primary_P_two_sided":
            float(
                primary_p
            ),

        "Secondary_P_directional":
            float(
                directional_p
            ),

        "Direction_vs_null":
            direction,

        # ----------------------------------------------------
        # MATCHING DIAGNOSTICS
        # ----------------------------------------------------

        "Program_mean_abs_LFC":
            observed_mean(
                "Abs_RbTKO_LFC"
            ),

        "Null_mean_abs_LFC":
            float(
                null[
                    "abs_lfc"
                ]
                .mean()
            ),

        "Delta_abs_LFC":
            (
                observed_mean(
                    "Abs_RbTKO_LFC"
                )
                -
                float(
                    null[
                        "abs_lfc"
                    ]
                    .mean()
                )
            ),

        "Program_mean_adult_logexpr":
            observed_mean(
                "Adult_control_log1p"
            ),

        "Null_mean_adult_logexpr":
            float(
                null[
                    "adult_expression"
                ]
                .mean()
            ),

        "Delta_adult_logexpr":
            (
                observed_mean(
                    "Adult_control_log1p"
                )
                -
                float(
                    null[
                        "adult_expression"
                    ]
                    .mean()
                )
            ),

        "Program_mean_CAG_rank":
            observed_mean(
                "CAG_control_rank"
            ),

        "Null_mean_CAG_rank":
            float(
                null[
                    "cag_rank"
                ]
                .mean()
            ),

        "Delta_CAG_rank":
            (
                observed_mean(
                    "CAG_control_rank"
                )
                -
                float(
                    null[
                        "cag_rank"
                    ]
                    .mean()
                )
            ),

        "Program_mean_MAP2_rank":
            observed_mean(
                "MAP2_control_rank"
            ),

        "Null_mean_MAP2_rank":
            float(
                null[
                    "map2_rank"
                ]
                .mean()
            ),

        "Delta_MAP2_rank":
            (
                observed_mean(
                    "MAP2_control_rank"
                )
                -
                float(
                    null[
                        "map2_rank"
                    ]
                    .mean()
                )
            ),
    }

    return result


# ============================================================
# RUN ONE MATCHING SETTING
# ============================================================

def run_setting(
    universe: pd.DataFrame,
    programs: pd.DataFrame,
    k: int,
    iterations: int,
    seed: int,
    exclude_all_program_genes: bool,
):

    all_program_genes = set(
        programs[
            "Gene"
        ]
        .dropna()
        .astype(str)
    )

    program_names = list(
        programs[
            "Program"
        ]
        .dropna()
        .astype(str)
        .unique()
    )

    rows = []

    null_store = {}

    for i, program_name in enumerate(
        program_names
    ):

        indices = get_program_indices(
            universe,
            programs,
            program_name,
        )

        if len(
            indices
        ) == 0:

            continue

        neighbors = build_neighbor_lists(
            universe=universe,
            program_indices=indices,
            k=k,
            exclude_all_program_genes=(
                exclude_all_program_genes
            ),
            all_program_genes=(
                all_program_genes
            ),
        )

        null = simulate_matched_null(
            universe=universe,
            program_indices=indices,
            neighbor_lists=neighbors,
            iterations=iterations,
            seed=(
                seed
                +
                i * 1009
            ),
        )

        row = summarize_program(
            universe=universe,
            program_name=program_name,
            program_indices=indices,
            null=null,
        )

        row["K"] = k

        row[
            "Exclude_all_program_genes"
        ] = (
            exclude_all_program_genes
        )

        rows.append(
            row
        )

        null_store[
            program_name
        ] = null[
            "residual"
        ]

    results = pd.DataFrame(
        rows
    )

    results[
        "Primary_FDR_two_sided"
    ] = bh_fdr(
        results[
            "Primary_P_two_sided"
        ]
        .to_numpy()
    )

    results[
        "Secondary_FDR_directional"
    ] = bh_fdr(
        results[
            "Secondary_P_directional"
        ]
        .to_numpy()
    )

    return (
        results,
        null_store,
    )


# ============================================================
# DETECTABILITY AUDIT
# ============================================================

def build_detectability_audit(
    universe: pd.DataFrame,
    programs: pd.DataFrame,
) -> pd.DataFrame:

    program_genes = (
        programs[
            [
                "Program",
                "Gene",
            ]
        ]
        .copy()
    )

    relevant = (
        universe[
            [
                "Gene",
                "RbTKO_LFC",
                "Adult_residual_z",
                "Adult_control_mean",
                "CAG_control_mean",
                "MAP2_control_mean",
                "CAG_control_rank",
                "MAP2_control_rank",
            ]
        ]
        .copy()
    )

    audit = (
        program_genes
        .merge(
            relevant,
            on="Gene",
            how="left",
            validate="many_to_one",
        )
    )

    audit[
        "CAG_bottom_10pct"
    ] = (
        audit[
            "CAG_control_rank"
        ]
        <=
        0.10
    )

    audit[
        "MAP2_bottom_10pct"
    ] = (
        audit[
            "MAP2_control_rank"
        ]
        <=
        0.10
    )

    audit[
        "Either_cortical_bottom_10pct"
    ] = (
        audit[
            "CAG_bottom_10pct"
        ]
        |
        audit[
            "MAP2_bottom_10pct"
        ]
    )

    return audit


# ============================================================
# PLOTTING
# ============================================================

def plot_summary(
    results: pd.DataFrame,
    output_path: Path,
) -> None:

    data = (
        results
        .sort_values(
            "Observed_mean_residual_z"
        )
        .reset_index(
            drop=True
        )
    )

    y = np.arange(
        len(
            data
        )
    )

    lower_error = (
        data[
            "Matched_null_mean"
        ]
        -
        data[
            "Matched_q025"
        ]
    )

    upper_error = (
        data[
            "Matched_q975"
        ]
        -
        data[
            "Matched_null_mean"
        ]
    )

    fig, ax = plt.subplots(
        figsize=(
            9.2,
            6.2,
        )
    )

    ax.errorbar(
        data[
            "Matched_null_mean"
        ],
        y,
        xerr=np.vstack(
            [
                lower_error,
                upper_error,
            ]
        ),
        fmt="o",
        capsize=4,
        label=(
            "Matched-null mean "
            "and 95% interval"
        ),
    )

    ax.scatter(
        data[
            "Observed_mean_residual_z"
        ],
        y,
        marker="x",
        s=70,
        label=(
            "Observed program mean"
        ),
    )

    ax.axvline(
        0,
        linewidth=1,
    )

    ax.set_yticks(
        y
    )

    ax.set_yticklabels(
        data[
            "Program"
        ]
    )

    ax.set_xlabel(
        "Mean adult residual z"
    )

    ax.set_title(
        "F10D: cross-context detectability-matched null"
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=300,
    )

    plt.close(
        fig
    )


def plot_null_histogram(
    program_name: str,
    observed: float,
    null_values: np.ndarray,
    output_path: Path,
) -> None:

    fig, ax = plt.subplots(
        figsize=(
            7.6,
            5.2,
        )
    )

    ax.hist(
        null_values,
        bins=45,
        density=True,
        alpha=0.75,
    )

    ax.axvline(
        observed,
        linewidth=2,
        label=(
            f"Observed = "
            f"{observed:.3f}"
        ),
    )

    ax.axvline(
        null_values.mean(),
        linestyle="--",
        linewidth=2,
        label=(
            "Matched-null mean = "
            f"{null_values.mean():.3f}"
        ),
    )

    ax.set_xlabel(
        "Mean adult residual z"
    )

    ax.set_ylabel(
        "Density"
    )

    ax.set_title(
        f"F10D matched null: "
        f"{program_name}"
    )

    ax.legend()

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=300,
    )

    plt.close(
        fig
    )


# ============================================================
# WRITE NULL DISTRIBUTIONS
# ============================================================

def write_null_distributions(
    null_store: dict[
        str,
        np.ndarray,
    ],
    output_path: Path,
) -> None:

    with gzip.open(
        output_path,
        "wt",
    ) as handle:

        handle.write(
            "Program,"
            "Iteration,"
            "Matched_null_mean_residual_z\n"
        )

        for (
            program_name,
            values,
        ) in null_store.items():

            safe_name = (
                program_name
                .replace(
                    ",",
                    "_",
                )
            )

            for (
                iteration,
                value,
            ) in enumerate(
                values
            ):

                handle.write(
                    f"{safe_name},"
                    f"{iteration},"
                    f"{value:.10g}\n"
                )


# ============================================================
# COMMAND-LINE ARGUMENTS
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--residuals",
        default=DEFAULT_RESIDUALS,
    )

    parser.add_argument(
        "--programs",
        default=DEFAULT_PROGRAMS,
    )

    parser.add_argument(
        "--fong",
        default=DEFAULT_FONG,
    )

    parser.add_argument(
        "--fong-sheet",
        default="Sig L2FC1",
    )

    parser.add_argument(
        "--osh-expression",
        default=DEFAULT_OSH_EXPR,
    )

    parser.add_argument(
        "--design",
        default=DEFAULT_DESIGN,
    )

    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
    )

    parser.add_argument(
        "--iterations",
        type=int,
        default=50_000,
    )

    parser.add_argument(
        "--sensitivity-iterations",
        type=int,
        default=10_000,
    )

    parser.add_argument(
        "--k",
        type=int,
        default=25,
    )

    parser.add_argument(
        "--sensitivity-k",
        default="10,25,50,100",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20261022,
    )

    return parser.parse_args()


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_args()

    output_directory = Path(
        args.outdir
    )

    figure_directory = (
        output_directory
        /
        "figures"
    )

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    figure_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ========================================================
    # LOAD DATA
    # ========================================================

    (
        universe,
        programs,
        sample_mapping,
    ) = load_analysis_universe(
        residual_path=Path(
            args.residuals
        ),
        program_path=Path(
            args.programs
        ),
        fong_path=Path(
            args.fong
        ),
        fong_sheet=(
            args.fong_sheet
        ),
        oshikawa_expression_path=Path(
            args.osh_expression
        ),
        design_path=Path(
            args.design
        ),
    )

    # ========================================================
    # SAVE SAMPLE MAPPING
    # ========================================================

    sample_mapping.to_csv(
        output_directory
        /
        "F10D_cortical_sample_mapping.csv",
        index=False,
    )

    print()
    print(
        "Inferred cortical control sample mapping"
    )
    print(
        "----------------------------------------"
    )
    print(
        sample_mapping.to_string(
            index=False
        )
    )
    print()

    # ========================================================
    # SAVE ANALYSIS UNIVERSE
    # ========================================================

    universe.to_csv(
        output_directory
        /
        "F10D_analysis_universe.csv",
        index=False,
    )

    # ========================================================
    # DETECTABILITY AUDIT
    # ========================================================

    detectability_audit = (
        build_detectability_audit(
            universe,
            programs,
        )
    )

    detectability_audit.to_csv(
        output_directory
        /
        "F10D_program_detectability_audit.csv",
        index=False,
    )

    # ========================================================
    # PRIMARY ANALYSIS
    # ========================================================

    (
        primary,
        null_store,
    ) = run_setting(
        universe=universe,
        programs=programs,
        k=args.k,
        iterations=args.iterations,
        seed=args.seed,
        exclude_all_program_genes=False,
    )

    primary.to_csv(
        output_directory
        /
        "F10D_primary.csv",
        index=False,
    )

    # ========================================================
    # MATCHING DIAGNOSTICS
    # ========================================================

    diagnostic_columns = [
        "Program",
        "Genes_available",
        "K",
        "Delta_abs_LFC",
        "Delta_adult_logexpr",
        "Delta_CAG_rank",
        "Delta_MAP2_rank",
    ]

    primary[
        diagnostic_columns
    ].to_csv(
        output_directory
        /
        "F10D_matching_diagnostics.csv",
        index=False,
    )

    # ========================================================
    # SENSITIVITY ANALYSIS
    # ========================================================

    k_values = [
        int(
            value.strip()
        )

        for value
        in args.sensitivity_k.split(
            ","
        )

        if value.strip()
    ]

    sensitivity_frames = []

    setting_counter = 0

    for k in k_values:

        for exclude_all in [
            False,
            True,
        ]:

            (
                sensitivity_result,
                _,
            ) = run_setting(
                universe=universe,
                programs=programs,
                k=k,
                iterations=(
                    args.sensitivity_iterations
                ),
                seed=(
                    args.seed
                    +
                    1_000_000
                    +
                    setting_counter
                    *
                    100_003
                ),
                exclude_all_program_genes=(
                    exclude_all
                ),
            )

            sensitivity_frames.append(
                sensitivity_result
            )

            setting_counter += 1

    sensitivity = pd.concat(
        sensitivity_frames,
        ignore_index=True,
    )

    sensitivity.to_csv(
        output_directory
        /
        "F10D_sensitivity.csv",
        index=False,
    )

    # ========================================================
    # SAVE NULL DISTRIBUTIONS
    # ========================================================

    write_null_distributions(
        null_store,
        output_directory
        /
        "F10D_null_distributions.csv.gz",
    )

    # ========================================================
    # FIGURES
    # ========================================================

    plot_summary(
        primary,
        figure_directory
        /
        "F10D_matched_null_summary.png",
    )

    lookup = (
        primary
        .set_index(
            "Program"
        )
    )

    for program_name in [
        "INNATE_STRESS",
        "FATE_DEATH_COMPOSITE",
    ]:

        if (
            program_name
            in null_store
            and
            program_name
            in lookup.index
        ):

            plot_null_histogram(
                program_name=(
                    program_name
                ),
                observed=float(
                    lookup.loc[
                        program_name,
                        "Observed_mean_residual_z",
                    ]
                ),
                null_values=(
                    null_store[
                        program_name
                    ]
                ),
                output_path=(
                    figure_directory
                    /
                    (
                        f"F10D_"
                        f"{program_name}"
                        f"_null.png"
                    )
                ),
            )

    # ========================================================
    # TEXT SUMMARY
    # ========================================================

    lines = []

    lines.append(
        "F10D — CROSS-CONTEXT DETECTABILITY-MATCHED NULL"
    )

    lines.append(
        "=" * 55
    )

    lines.append("")

    lines.append(
        f"Analysis universe: "
        f"{len(universe)} genes"
    )

    lines.append("")

    lines.append(
        "Matching variables:"
    )

    lines.append(
        "  - exact adult Rb-TKO direction"
    )

    lines.append(
        "  - |adult Rb-TKO log2FC|"
    )

    lines.append(
        "  - adult control expression "
        "(log1p mean RB_THC_1..3)"
    )

    lines.append(
        "  - CAG-control expression percentile"
    )

    lines.append(
        "  - MAP2-control expression percentile"
    )

    lines.append("")

    lines.append(
        "Primary inference:"
    )

    lines.append(
        "  two-sided matched-null P "
        "+ Benjamini-Hochberg FDR"
    )

    lines.append("")

    lines.append(
        "Secondary inference:"
    )

    lines.append(
        "  directional matched-null P "
        "+ Benjamini-Hochberg FDR"
    )

    lines.append("")

    lines.append(
        "PRIMARY RESULTS"
    )

    lines.append(
        "-" * 55
    )

    ordered_primary = (
        primary
        .sort_values(
            "Primary_P_two_sided"
        )
    )

    for _, row in (
        ordered_primary
        .iterrows()
    ):

        lines.append(
            f"{row['Program']}: "
            f"obs="
            f"{row['Observed_mean_residual_z']:.3f}, "
            f"null="
            f"{row['Matched_null_mean']:.3f}, "
            f"P2="
            f"{row['Primary_P_two_sided']:.5g}, "
            f"FDR2="
            f"{row['Primary_FDR_two_sided']:.5g}, "
            f"Pdir="
            f"{row['Secondary_P_directional']:.5g}, "
            f"FDRdir="
            f"{row['Secondary_FDR_directional']:.5g}"
        )

    lines.append("")

    lines.append(
        "SENSITIVITY RANGES"
    )

    lines.append(
        "-" * 55
    )

    for program_name in (
        primary[
            "Program"
        ]
    ):

        subset = sensitivity[
            sensitivity[
                "Program"
            ]
            ==
            program_name
        ]

        lines.append(
            f"{program_name}: "
            f"P2 range "
            f"["
            f"{subset['Primary_P_two_sided'].min():.5g}, "
            f"{subset['Primary_P_two_sided'].max():.5g}"
            f"], "
            f"FDR2 range "
            f"["
            f"{subset['Primary_FDR_two_sided'].min():.5g}, "
            f"{subset['Primary_FDR_two_sided'].max():.5g}"
            f"]"
        )

    lines.append("")

    lines.append(
        "INTERPRETATION GUARDRAILS"
    )

    lines.append(
        "-" * 55
    )

    lines.append(
        "- This analysis controls for baseline "
        "detectability in adult, CAG-control, "
        "and MAP2-control contexts."
    )

    lines.append(
        "- It does not make the adult and "
        "Agilent cortical effect sizes directly "
        "equivalent."
    )

    lines.append(
        "- A robust signal should remain "
        "directionally stable across matching "
        "neighborhoods and exclusion rules."
    )

    lines.append(
        "- Overlapping biological programs are "
        "correlated tests rather than independent "
        "replications."
    )

    lines.append(
        "- Residual enrichment remains associative "
        "and does not establish causality."
    )

    (
        output_directory
        /
        "F10D_summary.txt"
    ).write_text(
        "\n".join(
            lines
        ),
        encoding="utf-8",
    )

    # ========================================================
    # CONSOLE OUTPUT
    # ========================================================

    display_columns = [
        "Program",
        "Genes_available",
        "Observed_mean_residual_z",
        "Matched_null_mean",
        "Matched_effect_z",
        "Primary_P_two_sided",
        "Primary_FDR_two_sided",
        "Secondary_P_directional",
        "Secondary_FDR_directional",
    ]

    print(
        primary[
            display_columns
        ]
        .sort_values(
            "Primary_P_two_sided"
        )
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
        "Check F10D_cortical_sample_mapping.csv "
        "or the mapping printed above before "
        "interpreting the results."
    )


if __name__ == "__main__":
    main()