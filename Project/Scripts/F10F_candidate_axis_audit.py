#!/usr/bin/env python

"""
F10F_candidate_axis_audit.py
============================

Candidate-axis audit across the adult Rb-TKO, adult E2F1/3-DKO,
Oshikawa CAG/MAP2, and F10E conditional-divergence analyses.

PURPOSE
-------
The previous analyses weakened the original broad adult-death / MYC-switch
interpretation.

F10F therefore performs a deliberately descriptive audit of candidate axes
WITHOUT constructing another mechanistic model.

Questions:

1. Which candidate genes are present in the strict adult Rb-TKO list?
2. What happens to those genes after E2F1/3 deletion?
3. Are their cortical CAG/MAP2 responses conserved or divergent?
4. Are they conditionally unusual after the F10E matched-neighborhood test?
5. What happened to the p53 axis?
6. Does the evidence actually support a MYC-centered, NR4A-centered,
   p53-centered, or other interpretation?

IMPORTANT DATA LIMITATION
-------------------------
The available adult Rb-TKO workbook contains only:

    "Sig L2FC1"

Therefore:

    RbTKO_strict_present = False

means ONLY:

    "gene is not present in the supplied strict adult Rb-TKO table"

It does NOT mean:

    - the gene was not measured,
    - the gene was not expressed,
    - or the gene had no biological role.

The E2F1/3-DKO workbook contains a full sheet:

    "GSE190764_FDSD_FC.IHW"

and a raw-p<0.05 target sheet:

    "Significant targets p<0.05"

The Oshikawa response matrix is the full gene-level cortical response.

OUTPUT
------
Data/Processed/Dynamics/F10F_CandidateAxisAudit/

    F10F_candidate_axis_audit.csv
    F10F_candidate_axis_long.csv
    F10F_candidate_axis_summary.txt
"""


from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# DEFAULT PATHS
# ============================================================

DEFAULT_RB = (
    "Data/Raw/microarray/Fong2022/"
    "Fong 2022 L2FC1 microarray data.xlsx"
)

DEFAULT_E2F = (
    "Data/Raw/microarray/Fong2022/"
    "Fong 2022 E2F1.3 microarray data.xlsx"
)

DEFAULT_OSH_RESPONSE = (
    "Data/Processed/Dynamics/F8A_Oshikawa/"
    "oshikawa_gene_response_matrix.csv"
)

DEFAULT_OSH_EXPRESSION = (
    "Data/Processed/Dynamics/F8A_Oshikawa/"
    "oshikawa_gene_level_expression.csv"
)

DEFAULT_F10E = (
    "Data/Processed/Dynamics/F10E_ConditionalGeneDivergence/"
    "F10E_all_genes.csv"
)

DEFAULT_OUTDIR = (
    "Data/Processed/Dynamics/"
    "F10F_CandidateAxisAudit"
)


# ============================================================
# PRE-SPECIFIED CANDIDATE AXES
# ============================================================

CANDIDATE_AXES = {

    "MYC_DEATH_STRESS": [
        "MYC",
        "DDIT3",
        "MAP3K5",
        "FAS",
        "TNFRSF10B",
        "PMAIP1",
        "BBC3",
    ],

    "NR4A": [
        "NR4A1",
        "NR4A2",
        "NR4A3",
    ],

    "P53_AXIS": [
        "TRP53",
        "CDKN2A",
        "MDM2",
        "CDKN1A",
        "BAX",
        "BBC3",
        "PMAIP1",
    ],

    "INNATE_REGULATED_DEATH": [
        "ZBP1",
        "PYCARD",
        "MLKL",
        "EIF2AK2",
        "DTX3L",
        "IFIT2",
        "IFIT3",
        "IGTP",
        "IIGP1",
        "IRGM1",
        "PARP9",
        "TRIM25",
    ],

    "E2F_BACKBONE": [
        "E2F1",
        "E2F2",
        "CDK1",
        "CDK2",
        "PCNA",
        "CLSPN",
        "TOPBP1",
        "CHEK1",
        "CHEK2",
        "BRCA1",
        "H2AFX",
        "AURKB",
        "BIRC5",
        "BUB1",
        "MELK",
        "KIF14",
        "MYBL2",
    ],

    "TROPHIC_CONTEXT": [
        "BDNF",
        "FGF2",
        "KDR",
        "CXCL12",
        "ANGPT1",
        "IRS2",
        "MET",
        "THBS1",
        "PTK2B",
    ],
}


# ============================================================
# GENERAL UTILITIES
# ============================================================

def normalize_gene(
    value,
):

    if pd.isna(
        value
    ):
        return np.nan

    return str(
        value
    ).strip().upper()


def direction(
    value,
    tolerance=1e-12,
):

    if pd.isna(
        value
    ):
        return "NA"

    value = float(
        value
    )

    if value > tolerance:
        return "Up"

    if value < -tolerance:
        return "Down"

    return "Flat"


def same_sign(
    a,
    b,
):

    if (
        pd.isna(a)
        or pd.isna(b)
    ):
        return np.nan

    a = float(a)
    b = float(b)

    if (
        a == 0
        or b == 0
    ):
        return False

    return bool(
        np.sign(a)
        ==
        np.sign(b)
    )


def opposite_sign(
    a,
    b,
):

    if (
        pd.isna(a)
        or pd.isna(b)
    ):
        return np.nan

    a = float(a)
    b = float(b)

    if (
        a == 0
        or b == 0
    ):
        return False

    return bool(
        np.sign(a)
        ==
        -np.sign(b)
    )


def safe_mean(
    row,
    columns,
):

    values = pd.to_numeric(
        row[
            columns
        ],
        errors="coerce",
    )

    return values.mean()


def collapse_by_gene(
    dataframe,
    gene_column,
    priority_column=None,
    ascending=True,
):

    df = dataframe.copy()

    df[
        "Gene_norm"
    ] = (
        df[
            gene_column
        ]
        .map(
            normalize_gene
        )
    )

    df = (
        df
        .dropna(
            subset=[
                "Gene_norm"
            ]
        )
        .copy()
    )

    df[
        "_rows_for_symbol"
    ] = (
        df.groupby(
            "Gene_norm"
        )[
            "Gene_norm"
        ]
        .transform(
            "size"
        )
    )

    if (
        priority_column is not None
        and
        priority_column in df.columns
    ):

        df[
            "_priority"
        ] = pd.to_numeric(
            df[
                priority_column
            ],
            errors="coerce",
        )

        df = (
            df
            .sort_values(
                [
                    "Gene_norm",
                    "_priority",
                ],
                ascending=[
                    True,
                    ascending,
                ],
                na_position="last",
            )
            .drop_duplicates(
                "Gene_norm",
                keep="first",
            )
        )

    else:

        df = (
            df
            .drop_duplicates(
                "Gene_norm",
                keep="first",
            )
        )

    return df


# ============================================================
# ADULT RB-TKO
# ============================================================

def load_rb_tko(
    path,
):

    rb = pd.read_excel(
        path,
        sheet_name="Sig L2FC1",
    )

    required = [
        "Row.names",
        "Unnamed: 2",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj",
        "RB_THC_1",
        "RB_THC_2",
        "RB_THC_3",
        "RB_TKO_1",
        "RB_TKO_2",
        "RB_TKO_3",
    ]

    missing = [
        col
        for col in required
        if col not in rb.columns
    ]

    if missing:

        raise ValueError(
            "Rb-TKO workbook missing columns:\n"
            +
            "\n".join(
                missing
            )
        )

    rb = rb.rename(
        columns={
            "Unnamed: 2":
            "Gene_symbol",

            "Row.names":
            "Ensembl",
        }
    )

    # If duplicate symbols occur, retain the row with
    # the largest absolute adult Rb-TKO effect.
    rb[
        "_abs_lfc"
    ] = (
        pd.to_numeric(
            rb[
                "log2FoldChange"
            ],
            errors="coerce",
        )
        .abs()
    )

    rb = collapse_by_gene(
        rb,
        gene_column="Gene_symbol",
        priority_column="_abs_lfc",
        ascending=False,
    )

    rb[
        "RbTKO_control_mean"
    ] = (
        rb[
            [
                "RB_THC_1",
                "RB_THC_2",
                "RB_THC_3",
            ]
        ]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
        .mean(
            axis=1
        )
    )

    rb[
        "RbTKO_perturbed_mean"
    ] = (
        rb[
            [
                "RB_TKO_1",
                "RB_TKO_2",
                "RB_TKO_3",
            ]
        ]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
        .mean(
            axis=1
        )
    )

    keep = [
        "Gene_norm",
        "Ensembl",
        "_rows_for_symbol",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj",
        "RbTKO_control_mean",
        "RbTKO_perturbed_mean",
    ]

    rb = rb[
        keep
    ].rename(
        columns={
            "_rows_for_symbol":
            "RbTKO_symbol_rows",

            "baseMean":
            "RbTKO_baseMean",

            "log2FoldChange":
            "RbTKO_LFC",

            "pvalue":
            "RbTKO_pvalue",

            "padj":
            "RbTKO_padj",
        }
    )

    rb[
        "RbTKO_strict_present"
    ] = True

    return rb


# ============================================================
# E2F1/3-DKO
# ============================================================

def load_e2f_dko(
    path,
):

    full = pd.read_excel(
        path,
        sheet_name="GSE190764_FDSD_FC.IHW",
    )

    sig = pd.read_excel(
        path,
        sheet_name="Significant targets p<0.05",
    )

    required = [
        "Row.names",
        "name",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj",
        "EHT_A",
        "EHT_B",
        "EHT_C",
        "EKO_A",
        "EKO_B",
        "EKO_C",
    ]

    missing = [
        col
        for col in required
        if col not in full.columns
    ]

    if missing:

        raise ValueError(
            "E2F1/3-DKO workbook missing columns:\n"
            +
            "\n".join(
                missing
            )
        )

    full = full.rename(
        columns={
            "name":
            "Gene_symbol",

            "Row.names":
            "Ensembl_E2F",
        }
    )

    # Prefer the row with the smallest finite raw p-value.
    full = collapse_by_gene(
        full,
        gene_column="Gene_symbol",
        priority_column="pvalue",
        ascending=True,
    )

    full[
        "E2F_control_mean"
    ] = (
        full[
            [
                "EHT_A",
                "EHT_B",
                "EHT_C",
            ]
        ]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
        .mean(
            axis=1
        )
    )

    full[
        "E2F_DKO_mean"
    ] = (
        full[
            [
                "EKO_A",
                "EKO_B",
                "EKO_C",
            ]
        ]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
        .mean(
            axis=1
        )
    )

    sig_genes = set(
        sig[
            "name"
        ]
        .map(
            normalize_gene
        )
        .dropna()
    )

    full[
        "E2F_raw_p05_listed"
    ] = (
        full[
            "Gene_norm"
        ]
        .isin(
            sig_genes
        )
    )

    keep = [
        "Gene_norm",
        "Ensembl_E2F",
        "_rows_for_symbol",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj",
        "E2F_control_mean",
        "E2F_DKO_mean",
        "E2F_raw_p05_listed",
    ]

    full = full[
        keep
    ].rename(
        columns={
            "_rows_for_symbol":
            "E2F_symbol_rows",

            "baseMean":
            "E2F_baseMean",

            "log2FoldChange":
            "E2F_DKO_LFC",

            "pvalue":
            "E2F_DKO_pvalue",

            "padj":
            "E2F_DKO_padj",
        }
    )

    full[
        "E2F_measured"
    ] = True

    return full


# ============================================================
# OSHIKAWA RESPONSE
# ============================================================

def load_oshikawa_response(
    path,
):

    osh = pd.read_csv(
        path
    )

    required = [
        "Gene",
        "logFC_CAG",
        "P.Value_CAG",
        "adj.P.Val_CAG",
        "logFC_MAP2",
        "P.Value_MAP2",
        "adj.P.Val_MAP2",
        "logFC_Interaction",
        "P.Value_Interaction",
        "adj.P.Val_Interaction",
    ]

    missing = [
        col
        for col in required
        if col not in osh.columns
    ]

    if missing:

        raise ValueError(
            "Oshikawa response matrix missing columns:\n"
            +
            "\n".join(
                missing
            )
        )

    osh[
        "Gene_norm"
    ] = (
        osh[
            "Gene"
        ]
        .map(
            normalize_gene
        )
    )

    osh = (
        osh
        .dropna(
            subset=[
                "Gene_norm"
            ]
        )
        .drop_duplicates(
            "Gene_norm",
            keep="first",
        )
    )

    osh = osh.rename(
        columns={
            "logFC_CAG":
            "CAG_LFC",

            "P.Value_CAG":
            "CAG_pvalue",

            "adj.P.Val_CAG":
            "CAG_padj",

            "logFC_MAP2":
            "MAP2_LFC",

            "P.Value_MAP2":
            "MAP2_pvalue",

            "adj.P.Val_MAP2":
            "MAP2_padj",

            "logFC_Interaction":
            "Interaction_LFC",

            "P.Value_Interaction":
            "Interaction_pvalue",

            "adj.P.Val_Interaction":
            "Interaction_padj",
        }
    )

    keep = [
        "Gene_norm",
        "CAG_LFC",
        "CAG_pvalue",
        "CAG_padj",
        "MAP2_LFC",
        "MAP2_pvalue",
        "MAP2_padj",
        "Interaction_LFC",
        "Interaction_pvalue",
        "Interaction_padj",
    ]

    osh = osh[
        keep
    ]

    osh[
        "Osh_measured"
    ] = True

    return osh


# ============================================================
# OSHIKAWA BASELINE EXPRESSION
# ============================================================

def load_oshikawa_expression(
    path,
):

    expr = pd.read_csv(
        path
    )

    if "Gene" not in expr.columns:

        raise ValueError(
            "Expected 'Gene' column in "
            "oshikawa_gene_level_expression.csv"
        )

    cag_controls = [
        col
        for col in expr.columns
        if "cCont" in col
    ]

    map2_controls = [
        col
        for col in expr.columns
        if "mCont" in col
    ]

    if len(
        cag_controls
    ) != 3:

        raise ValueError(
            "Expected exactly 3 CAG-control columns. "
            f"Found: {cag_controls}"
        )

    if len(
        map2_controls
    ) != 3:

        raise ValueError(
            "Expected exactly 3 MAP2-control columns. "
            f"Found: {map2_controls}"
        )

    expr[
        "Gene_norm"
    ] = (
        expr[
            "Gene"
        ]
        .map(
            normalize_gene
        )
    )

    expr[
        "CAG_control_mean"
    ] = (
        expr[
            cag_controls
        ]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
        .mean(
            axis=1
        )
    )

    expr[
        "MAP2_control_mean"
    ] = (
        expr[
            map2_controls
        ]
        .apply(
            pd.to_numeric,
            errors="coerce",
        )
        .mean(
            axis=1
        )
    )

    expr = (
        expr[
            [
                "Gene_norm",
                "CAG_control_mean",
                "MAP2_control_mean",
            ]
        ]
        .dropna(
            subset=[
                "Gene_norm"
            ]
        )
        .groupby(
            "Gene_norm",
            as_index=False,
        )
        .mean(
            numeric_only=True
        )
    )

    expr[
        "CAG_control_percentile"
    ] = (
        expr[
            "CAG_control_mean"
        ]
        .rank(
            pct=True,
            method="average",
        )
    )

    expr[
        "MAP2_control_percentile"
    ] = (
        expr[
            "MAP2_control_mean"
        ]
        .rank(
            pct=True,
            method="average",
        )
    )

    return expr


# ============================================================
# F10E CONDITIONAL DIVERGENCE
# ============================================================

def load_f10e(
    path,
):

    path = Path(
        path
    )

    if not path.exists():

        print(
            "WARNING: F10E file not found. "
            "Conditional-divergence fields will be empty."
        )

        return pd.DataFrame(
            columns=[
                "Gene_norm"
            ]
        )

    f10e = pd.read_csv(
        path
    )

    f10e[
        "Gene_norm"
    ] = (
        f10e[
            "Gene"
        ]
        .map(
            normalize_gene
        )
    )

    desired = [
        "Gene_norm",
        "Adult_residual_z",
        "Median_conditional_z",
        "Min_conditional_z",
        "Max_conditional_z",
        "Min_local_percentile",
        "Max_local_percentile",
        "Positive_across_all_k",
        "Negative_across_all_k",
        "Robust_high_90pct",
        "Robust_high_z1",
        "Robust_rank",
    ]

    desired = [
        col
        for col in desired
        if col in f10e.columns
    ]

    return (
        f10e[
            desired
        ]
        .drop_duplicates(
            "Gene_norm",
            keep="first",
        )
    )


# ============================================================
# BUILD CANDIDATE MASTER
# ============================================================

def build_candidate_master():

    rows = []

    for (
        axis,
        genes,
    ) in CANDIDATE_AXES.items():

        for gene in genes:

            rows.append(
                {
                    "Axis":
                        axis,

                    "Gene":
                        gene,

                    "Gene_norm":
                        normalize_gene(
                            gene
                        ),
                }
            )

    candidates = pd.DataFrame(
        rows
    )

    # Preserve axis membership for genes that occur in more than one axis.
    candidates = (
        candidates
        .groupby(
            [
                "Gene",
                "Gene_norm",
            ],
            as_index=False,
        )[
            "Axis"
        ]
        .apply(
            lambda x:
            ";".join(
                sorted(
                    set(
                        x
                    )
                )
            )
        )
        .rename(
            columns={
                "Axis":
                "Axes"
            }
        )
    )

    return candidates


# ============================================================
# MERGE ALL CONTEXTS
# ============================================================

def build_audit(
    candidates,
    rb,
    e2f,
    osh,
    osh_expr,
    f10e,
):

    audit = candidates.copy()

    for table in [
        rb,
        e2f,
        osh,
        osh_expr,
        f10e,
    ]:

        audit = audit.merge(
            table,
            on="Gene_norm",
            how="left",
            validate="one_to_one",
        )

    # --------------------------------------------------------
    # Explicit presence flags
    # --------------------------------------------------------

    audit[
        "RbTKO_strict_present"
    ] = (
        audit[
            "RbTKO_strict_present"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
    )

    audit[
        "E2F_measured"
    ] = (
        audit[
            "E2F_measured"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
    )

    audit[
        "Osh_measured"
    ] = (
        audit[
            "Osh_measured"
        ]
        .fillna(
            False
        )
        .astype(
            bool
        )
    )

    if (
        "E2F_raw_p05_listed"
        in audit.columns
    ):

        audit[
            "E2F_raw_p05_listed"
        ] = (
            audit[
                "E2F_raw_p05_listed"
            ]
            .fillna(
                False
            )
            .astype(
                bool
            )
        )

    # --------------------------------------------------------
    # Directions
    # --------------------------------------------------------

    audit[
        "RbTKO_direction"
    ] = (
        audit[
            "RbTKO_LFC"
        ]
        .apply(
            direction
        )
    )

    audit[
        "E2F_DKO_direction"
    ] = (
        audit[
            "E2F_DKO_LFC"
        ]
        .apply(
            direction
        )
    )

    audit[
        "CAG_direction"
    ] = (
        audit[
            "CAG_LFC"
        ]
        .apply(
            direction
        )
    )

    audit[
        "MAP2_direction"
    ] = (
        audit[
            "MAP2_LFC"
        ]
        .apply(
            direction
        )
    )

    # --------------------------------------------------------
    # E2F reciprocity
    # --------------------------------------------------------

    audit[
        "E2F_reciprocal_to_RbTKO"
    ] = [
        opposite_sign(
            rb_lfc,
            e2f_lfc,
        )

        for (
            rb_lfc,
            e2f_lfc,
        ) in zip(
            audit[
                "RbTKO_LFC"
            ],
            audit[
                "E2F_DKO_LFC"
            ],
        )
    ]

    # --------------------------------------------------------
    # Cortical conservation
    # --------------------------------------------------------

    audit[
        "CAG_same_sign_as_RbTKO"
    ] = [
        same_sign(
            rb_lfc,
            cag_lfc,
        )

        for (
            rb_lfc,
            cag_lfc,
        ) in zip(
            audit[
                "RbTKO_LFC"
            ],
            audit[
                "CAG_LFC"
            ],
        )
    ]

    audit[
        "MAP2_same_sign_as_RbTKO"
    ] = [
        same_sign(
            rb_lfc,
            map_lfc,
        )

        for (
            rb_lfc,
            map_lfc,
        ) in zip(
            audit[
                "RbTKO_LFC"
            ],
            audit[
                "MAP2_LFC"
            ],
        )
    ]

    # --------------------------------------------------------
    # Statistical flags
    # --------------------------------------------------------

    audit[
        "E2F_raw_p05"
    ] = (
        pd.to_numeric(
            audit[
                "E2F_DKO_pvalue"
            ],
            errors="coerce",
        )
        <
        0.05
    )

    audit[
        "E2F_FDR05"
    ] = (
        pd.to_numeric(
            audit[
                "E2F_DKO_padj"
            ],
            errors="coerce",
        )
        <
        0.05
    )

    audit[
        "CAG_FDR05"
    ] = (
        pd.to_numeric(
            audit[
                "CAG_padj"
            ],
            errors="coerce",
        )
        <
        0.05
    )

    audit[
        "MAP2_FDR05"
    ] = (
        pd.to_numeric(
            audit[
                "MAP2_padj"
            ],
            errors="coerce",
        )
        <
        0.05
    )

    # --------------------------------------------------------
    # Explicit interpretation field for missing adult genes
    # --------------------------------------------------------

    audit[
        "Adult_Rb_status"
    ] = np.where(
        audit[
            "RbTKO_strict_present"
        ],
        "Present_in_strict_RbTKO_LFC1_list",
        "Absent_from_supplied_strict_RbTKO_list",
    )

    return audit


# ============================================================
# LONG FORMAT
# ============================================================

def build_long_table(
    audit,
):

    rows = []

    for _, row in audit.iterrows():

        gene = row[
            "Gene"
        ]

        axes = row[
            "Axes"
        ]

        # ----------------------------------------------------
        # Adult Rb-TKO
        # ----------------------------------------------------

        rows.append(
            {
                "Gene":
                    gene,

                "Axes":
                    axes,

                "Context":
                    "Adult_RbTKO",

                "Available":
                    row[
                        "RbTKO_strict_present"
                    ],

                "LFC":
                    row.get(
                        "RbTKO_LFC",
                        np.nan,
                    ),

                "P_value":
                    row.get(
                        "RbTKO_pvalue",
                        np.nan,
                    ),

                "FDR":
                    row.get(
                        "RbTKO_padj",
                        np.nan,
                    ),

                "Baseline_expression":
                    row.get(
                        "RbTKO_control_mean",
                        np.nan,
                    ),

                "Perturbed_expression":
                    row.get(
                        "RbTKO_perturbed_mean",
                        np.nan,
                    ),

                "Availability_note":
                    row[
                        "Adult_Rb_status"
                    ],
            }
        )

        # ----------------------------------------------------
        # Adult E2F1/3-DKO
        # ----------------------------------------------------

        rows.append(
            {
                "Gene":
                    gene,

                "Axes":
                    axes,

                "Context":
                    "Adult_E2F1_3_DKO",

                "Available":
                    row[
                        "E2F_measured"
                    ],

                "LFC":
                    row.get(
                        "E2F_DKO_LFC",
                        np.nan,
                    ),

                "P_value":
                    row.get(
                        "E2F_DKO_pvalue",
                        np.nan,
                    ),

                "FDR":
                    row.get(
                        "E2F_DKO_padj",
                        np.nan,
                    ),

                "Baseline_expression":
                    row.get(
                        "E2F_control_mean",
                        np.nan,
                    ),

                "Perturbed_expression":
                    row.get(
                        "E2F_DKO_mean",
                        np.nan,
                    ),

                "Availability_note":
                    "Full_E2F_DKO_table",
            }
        )

        # ----------------------------------------------------
        # CAG-TKO
        # ----------------------------------------------------

        rows.append(
            {
                "Gene":
                    gene,

                "Axes":
                    axes,

                "Context":
                    "CAG_TKO",

                "Available":
                    row[
                        "Osh_measured"
                    ],

                "LFC":
                    row.get(
                        "CAG_LFC",
                        np.nan,
                    ),

                "P_value":
                    row.get(
                        "CAG_pvalue",
                        np.nan,
                    ),

                "FDR":
                    row.get(
                        "CAG_padj",
                        np.nan,
                    ),

                "Baseline_expression":
                    row.get(
                        "CAG_control_mean",
                        np.nan,
                    ),

                "Perturbed_expression":
                    np.nan,

                "Availability_note":
                    "Oshikawa_gene_level_matrix",
            }
        )

        # ----------------------------------------------------
        # MAP2-TKO
        # ----------------------------------------------------

        rows.append(
            {
                "Gene":
                    gene,

                "Axes":
                    axes,

                "Context":
                    "MAP2_TKO",

                "Available":
                    row[
                        "Osh_measured"
                    ],

                "LFC":
                    row.get(
                        "MAP2_LFC",
                        np.nan,
                    ),

                "P_value":
                    row.get(
                        "MAP2_pvalue",
                        np.nan,
                    ),

                "FDR":
                    row.get(
                        "MAP2_padj",
                        np.nan,
                    ),

                "Baseline_expression":
                    row.get(
                        "MAP2_control_mean",
                        np.nan,
                    ),

                "Perturbed_expression":
                    np.nan,

                "Availability_note":
                    "Oshikawa_gene_level_matrix",
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# SUMMARY TEXT
# ============================================================

def format_number(
    value,
    digits=3,
):

    if pd.isna(
        value
    ):
        return "NA"

    return (
        f"{float(value):.{digits}f}"
    )


def write_summary(
    audit,
    output_path,
):

    indexed = (
        audit
        .set_index(
            "Gene"
        )
    )

    lines = []

    lines.append(
        "F10F — CANDIDATE-AXIS AUDIT"
    )

    lines.append(
        "=" * 50
    )

    lines.append("")

    lines.append(
        "IMPORTANT DATA LIMITATION"
    )

    lines.append(
        "-" * 50
    )

    lines.append(
        "The supplied adult Rb-TKO workbook contains only the "
        "'Sig L2FC1' sheet."
    )

    lines.append(
        "Therefore absence from that workbook means absence from "
        "the supplied strict adult DEG list, not absence from the "
        "original experiment."
    )

    lines.append("")

    # ========================================================
    # MYC
    # ========================================================

    lines.append(
        "MYC AUDIT"
    )

    lines.append(
        "-" * 50
    )

    if "MYC" in indexed.index:

        row = indexed.loc[
            "MYC"
        ]

        lines.append(
            "Adult Rb-TKO LFC: "
            +
            format_number(
                row.get(
                    "RbTKO_LFC"
                )
            )
        )

        lines.append(
            "E2F1/3-DKO LFC: "
            +
            format_number(
                row.get(
                    "E2F_DKO_LFC"
                )
            )
            +
            ", p="
            +
            format_number(
                row.get(
                    "E2F_DKO_pvalue"
                ),
                digits=4,
            )
        )

        lines.append(
            "CAG LFC: "
            +
            format_number(
                row.get(
                    "CAG_LFC"
                )
            )
        )

        lines.append(
            "MAP2 LFC: "
            +
            format_number(
                row.get(
                    "MAP2_LFC"
                )
            )
        )

        lines.append(
            "F10E median conditional z: "
            +
            format_number(
                row.get(
                    "Median_conditional_z"
                )
            )
        )

        lines.append(
            "E2F reciprocal to adult Rb-TKO: "
            +
            str(
                row.get(
                    "E2F_reciprocal_to_RbTKO"
                )
            )
        )

    lines.append("")

    # ========================================================
    # NR4A
    # ========================================================

    lines.append(
        "NR4A AUDIT"
    )

    lines.append(
        "-" * 50
    )

    for gene in [
        "NR4A1",
        "NR4A2",
        "NR4A3",
    ]:

        if gene not in indexed.index:
            continue

        row = indexed.loc[
            gene
        ]

        lines.append(
            f"{gene}: "
            f"RbTKO="
            f"{format_number(row.get('RbTKO_LFC'))}, "
            f"E2F_DKO="
            f"{format_number(row.get('E2F_DKO_LFC'))}, "
            f"CAG="
            f"{format_number(row.get('CAG_LFC'))}, "
            f"MAP2="
            f"{format_number(row.get('MAP2_LFC'))}, "
            f"conditional_z="
            f"{format_number(row.get('Median_conditional_z'))}"
        )

    lines.append("")

    # ========================================================
    # P53
    # ========================================================

    lines.append(
        "P53-AXIS AVAILABILITY AUDIT"
    )

    lines.append(
        "-" * 50
    )

    for gene in [
        "TRP53",
        "CDKN2A",
        "MDM2",
        "CDKN1A",
        "BAX",
        "BBC3",
        "PMAIP1",
    ]:

        if gene not in indexed.index:
            continue

        row = indexed.loc[
            gene
        ]

        lines.append(
            f"{gene}: "
            f"adult_strict={row.get('RbTKO_strict_present')}, "
            f"E2F_measured={row.get('E2F_measured')}, "
            f"E2F_LFC={format_number(row.get('E2F_DKO_LFC'))}, "
            f"Osh_measured={row.get('Osh_measured')}, "
            f"CAG={format_number(row.get('CAG_LFC'))}, "
            f"MAP2={format_number(row.get('MAP2_LFC'))}"
        )

    lines.append("")

    # ========================================================
    # CONDITIONAL TOP CANDIDATES
    # ========================================================

    if (
        "Median_conditional_z"
        in audit.columns
    ):

        lines.append(
            "TOP PRE-SPECIFIED CANDIDATES BY F10E CONDITIONAL Z"
        )

        lines.append(
            "-" * 50
        )

        ranked = (
            audit
            .dropna(
                subset=[
                    "Median_conditional_z"
                ]
            )
            .sort_values(
                "Median_conditional_z",
                ascending=False,
            )
        )

        for _, row in (
            ranked
            .head(
                15
            )
            .iterrows()
        ):

            lines.append(
                f"{row['Gene']}: "
                f"conditional_z="
                f"{format_number(row['Median_conditional_z'])}, "
                f"adult_RbTKO="
                f"{format_number(row.get('RbTKO_LFC'))}"
            )

    lines.append("")

    lines.append(
        "INTERPRETATION GUARDRAILS"
    )

    lines.append(
        "-" * 50
    )

    lines.append(
        "- This audit is descriptive and does not establish causal "
        "regulatory edges."
    )

    lines.append(
        "- A gene can be mechanistically important without being "
        "transcriptionally differential."
    )

    lines.append(
        "- TRP53 pathway activity cannot be excluded by lack of TRP53 "
        "mRNA differential expression because p53 is strongly "
        "post-translationally regulated."
    )

    lines.append(
        "- Prior-knowledge regulatory edges should not be treated as "
        "active in adult NSCs unless perturbation data support them."
    )

    lines.append(
        "- F11 should only be rebuilt after this audit is interpreted."
    )

    Path(
        output_path
    ).write_text(
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
        "--rb",
        default=DEFAULT_RB,
    )

    parser.add_argument(
        "--e2f",
        default=DEFAULT_E2F,
    )

    parser.add_argument(
        "--osh-response",
        default=DEFAULT_OSH_RESPONSE,
    )

    parser.add_argument(
        "--osh-expression",
        default=DEFAULT_OSH_EXPRESSION,
    )

    parser.add_argument(
        "--f10e",
        default=DEFAULT_F10E,
    )

    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
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

    output_directory.mkdir(
        parents=True,
        exist_ok=True,
    )

    print(
        "Loading adult Rb-TKO strict table..."
    )

    rb = load_rb_tko(
        args.rb
    )

    print(
        "Loading full E2F1/3-DKO table..."
    )

    e2f = load_e2f_dko(
        args.e2f
    )

    print(
        "Loading Oshikawa response matrix..."
    )

    osh = load_oshikawa_response(
        args.osh_response
    )

    print(
        "Loading Oshikawa baseline expression..."
    )

    osh_expr = load_oshikawa_expression(
        args.osh_expression
    )

    print(
        "Loading F10E conditional ranking..."
    )

    f10e = load_f10e(
        args.f10e
    )

    candidates = (
        build_candidate_master()
    )

    audit = build_audit(
        candidates=candidates,
        rb=rb,
        e2f=e2f,
        osh=osh,
        osh_expr=osh_expr,
        f10e=f10e,
    )

    # ========================================================
    # SAVE WIDE TABLE
    # ========================================================

    audit.to_csv(
        output_directory
        /
        "F10F_candidate_axis_audit.csv",
        index=False,
    )

    # ========================================================
    # SAVE LONG TABLE
    # ========================================================

    long_table = build_long_table(
        audit
    )

    long_table.to_csv(
        output_directory
        /
        "F10F_candidate_axis_long.csv",
        index=False,
    )

    # ========================================================
    # SUMMARY
    # ========================================================

    write_summary(
        audit=audit,
        output_path=(
            output_directory
            /
            "F10F_candidate_axis_summary.txt"
        ),
    )

    # ========================================================
    # CONSOLE DISPLAY
    # ========================================================

    display_columns = [
        "Gene",
        "Axes",

        "RbTKO_strict_present",
        "RbTKO_LFC",

        "E2F_measured",
        "E2F_DKO_LFC",
        "E2F_DKO_pvalue",
        "E2F_DKO_padj",
        "E2F_reciprocal_to_RbTKO",

        "CAG_LFC",
        "CAG_padj",

        "MAP2_LFC",
        "MAP2_padj",

        "Median_conditional_z",
        "Min_local_percentile",
        "Robust_rank",
    ]

    display_columns = [
        column
        for column
        in display_columns
        if column
        in audit.columns
    ]

    print()
    print(
        "============================================================"
    )

    print(
        "CANDIDATE-AXIS AUDIT"
    )

    print(
        "============================================================"
    )

    print()

    print(
        audit[
            display_columns
        ]
        .sort_values(
            [
                "Axes",
                "Gene",
            ]
        )
        .to_string(
            index=False
        )
    )

    # ========================================================
    # SMALL FOCUSED PANELS
    # ========================================================

    print()
    print(
        "============================================================"
    )

    print(
        "MYC / DDIT3 / MAP3K5"
    )

    print(
        "============================================================"
    )

    focused = audit[
        audit[
            "Gene"
        ]
        .isin(
            [
                "MYC",
                "DDIT3",
                "MAP3K5",
            ]
        )
    ]

    print(
        focused[
            display_columns
        ]
        .to_string(
            index=False
        )
    )

    print()
    print(
        "============================================================"
    )

    print(
        "NR4A"
    )

    print(
        "============================================================"
    )

    focused = audit[
        audit[
            "Gene"
        ]
        .isin(
            [
                "NR4A1",
                "NR4A2",
                "NR4A3",
            ]
        )
    ]

    print(
        focused[
            display_columns
        ]
        .to_string(
            index=False
        )
    )

    print()
    print(
        "============================================================"
    )

    print(
        "P53 AXIS"
    )

    print(
        "============================================================"
    )

    focused = audit[
        audit[
            "Gene"
        ]
        .isin(
            [
                "TRP53",
                "CDKN2A",
                "MDM2",
                "CDKN1A",
                "BAX",
                "BBC3",
                "PMAIP1",
            ]
        )
    ]

    print(
        focused[
            display_columns
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
        "RbTKO_strict_present=False means absent from the supplied "
        "'Sig L2FC1' workbook only. Do NOT interpret it as "
        "'not measured' or 'not expressed'."
    )


if __name__ == "__main__":
    main()