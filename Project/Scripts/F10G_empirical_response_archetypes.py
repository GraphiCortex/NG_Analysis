#!/usr/bin/env python

"""
F10G_empirical_response_archetypes.py
=====================================

Data-first empirical response archetypes for complete-context Rb-loss genes.

RATIONALE
---------
F10C-F10F showed that predefined adult death/stress programs and the
MYC-centered switch interpretation were not robust to selection-aware
and detectability-aware controls.

F10G therefore removes biological state labels from the clustering step.

The analysis asks:

    Among genes already known to respond strongly in adult Rb-TKO cells,
    what empirical response patterns are observed across:

        1. adult E2F1/3-DKO,
        2. CAG-driven cortical Rb loss,
        3. MAP2-driven cortical Rb loss?

IMPORTANT DESIGN CHOICE
-----------------------
The adult Rb-TKO sheet available to us is already selected for strong
adult effects.

Therefore naive clustering across all four LFC columns would tend to
rediscover the trivial adult Up/Down split.

To reduce this circularity, F10G:

    - stratifies genes by adult Rb-TKO direction first;
    - clusters WITHIN adult-Up and adult-Down strata;
    - uses only E2F1/3-DKO, CAG, and MAP2 responses for clustering;
    - rank-Gaussian normalizes each response dimension within each stratum.

Adult Rb-TKO magnitude is retained only as an annotation.

Thus these are response archetypes CONDITIONAL ON adult Rb-TKO
responsiveness. They are not unbiased transcriptome-wide states.

MODEL SELECTION
---------------
For each adult-direction stratum, k = 2..max_k is evaluated using:

    - silhouette score;
    - perturbation stability using adjusted Rand index;
    - minimum cluster size.

The primary k is the k with the largest silhouette score.

The perturbation-stability metric measures algorithmic robustness only.
It is NOT biological replicate robustness.

A separate supported sign-pattern census is also generated using the
precomputed F10A strict-response flags. This provides a non-clustering
cross-check.

OUTPUT
------
Data/Processed/Dynamics/F10G_EmpiricalResponseArchetypes/

    F10G_model_selection.csv
    F10G_cluster_assignments.csv
    F10G_cluster_summary.csv
    F10G_supported_pattern_census.csv
    F10G_candidate_placement.csv
    F10G_summary.txt

    figures/
        F10G_model_selection_Up.png
        F10G_model_selection_Down.png
        F10G_centroids_Up.png
        F10G_centroids_Down.png
"""

from __future__ import annotations

import argparse
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# PATHS
# ============================================================

DEFAULT_INPUT = (
    "Data/Processed/Dynamics/F10A_AdultResiduals/"
    "F10A_all_adult_residuals_ranked.csv"
)

DEFAULT_OUTDIR = (
    "Data/Processed/Dynamics/"
    "F10G_EmpiricalResponseArchetypes"
)


# ============================================================
# FEATURES
# ============================================================

CLUSTER_FEATURES = [
    "E2F13_LFC",
    "CAG_LFC",
    "MAP2_LFC",
]


# ============================================================
# CANDIDATES TO TRACK
# ============================================================

CANDIDATE_GENES = [
    "MYC",
    "DDIT3",
    "MAP3K5",
    "FAS",
    "TNFRSF10B",
    "PMAIP1",
    "BBC3",
    "ZBP1",
    "PYCARD",
    "MLKL",
    "EIF2AK2",

    "NR4A1",
    "NR4A2",
    "NR4A3",

    "TRP53",
    "CDKN1A",
    "CDKN2A",
    "MDM2",
    "BAX",

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

    "FGF2",
    "KDR",
    "CXCL12",
    "BDNF",
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
            f"{label} is missing required columns:\n"
            +
            "\n".join(
                f"  - {column}"
                for column in missing
            )
        )


def rank_gaussian(
    series: pd.Series,
) -> pd.Series:
    """
    Convert values to normal scores using within-stratum ranks.

        p_i = (rank_i - 0.5) / n
        z_i = Phi^{-1}(p_i)

    This preserves response ordering while reducing scale dominance.
    """

    values = pd.to_numeric(
        series,
        errors="coerce",
    )

    ranks = values.rank(
        method="average",
    )

    n = int(
        values.notna().sum()
    )

    probabilities = (
        ranks - 0.5
    ) / n

    normal = NormalDist()

    transformed = probabilities.map(
        lambda p:
        normal.inv_cdf(
            float(p)
        )
        if pd.notna(p)
        else np.nan
    )

    return transformed


# ============================================================
# K-MEANS++
# ============================================================

def kmeans_pp_once(
    X: np.ndarray,
    k: int,
    rng: np.random.Generator,
    max_iter: int = 300,
    tol: float = 1e-8,
):
    """
    Minimal k-means++ implementation.

    No sklearn dependency required.
    """

    n, d = X.shape

    centers = np.empty(
        (k, d),
        dtype=float,
    )

    first = int(
        rng.integers(
            0,
            n,
        )
    )

    centers[0] = X[first]

    min_dist2 = np.sum(
        (
            X
            -
            centers[0]
        )
        ** 2,
        axis=1,
    )

    for center_index in range(
        1,
        k,
    ):

        total = float(
            min_dist2.sum()
        )

        if total <= 0:

            chosen = int(
                rng.integers(
                    0,
                    n,
                )
            )

        else:

            probabilities = (
                min_dist2
                /
                total
            )

            chosen = int(
                rng.choice(
                    n,
                    p=probabilities,
                )
            )

        centers[
            center_index
        ] = X[
            chosen
        ]

        new_dist2 = np.sum(
            (
                X
                -
                centers[
                    center_index
                ]
            )
            ** 2,
            axis=1,
        )

        min_dist2 = np.minimum(
            min_dist2,
            new_dist2,
        )

    previous_inertia = np.inf

    for _ in range(
        max_iter
    ):

        dist2 = np.sum(
            (
                X[:, None, :]
                -
                centers[None, :, :]
            )
            ** 2,
            axis=2,
        )

        labels = np.argmin(
            dist2,
            axis=1,
        )

        # Repair empty clusters if needed.
        for cluster in range(
            k
        ):

            if not np.any(
                labels == cluster
            ):

                nearest_center_dist = np.min(
                    dist2,
                    axis=1,
                )

                farthest = int(
                    np.argmax(
                        nearest_center_dist
                    )
                )

                labels[
                    farthest
                ] = cluster

        new_centers = np.vstack(
            [
                X[
                    labels == cluster
                ].mean(
                    axis=0
                )

                for cluster
                in range(
                    k
                )
            ]
        )

        inertia = float(
            np.sum(
                (
                    X
                    -
                    new_centers[
                        labels
                    ]
                )
                ** 2
            )
        )

        if np.isfinite(
            previous_inertia
        ):

            change = abs(
                previous_inertia
                -
                inertia
            )

            scale = max(
                1.0,
                previous_inertia,
            )

            if change <= (
                tol
                *
                scale
            ):

                centers = (
                    new_centers
                )

                break

        centers = (
            new_centers
        )

        previous_inertia = (
            inertia
        )

    final_dist2 = np.sum(
        (
            X[:, None, :]
            -
            centers[None, :, :]
        )
        ** 2,
        axis=2,
    )

    labels = np.argmin(
        final_dist2,
        axis=1,
    )

    inertia = float(
        final_dist2[
            np.arange(
                n
            ),
            labels,
        ].sum()
    )

    return (
        labels,
        centers,
        inertia,
    )


def best_kmeans(
    X: np.ndarray,
    k: int,
    n_init: int,
    seed: int,
):
    """
    Repeat k-means++ initialization and retain the minimum-inertia fit.
    """

    best = None

    for init in range(
        n_init
    ):

        rng = np.random.default_rng(
            seed
            +
            init
            *
            1009
        )

        result = kmeans_pp_once(
            X,
            k,
            rng,
        )

        if (
            best is None
            or
            result[2]
            <
            best[2]
        ):

            best = (
                result
            )

    return best


# ============================================================
# DISTANCE + SILHOUETTE
# ============================================================

def pairwise_euclidean(
    X: np.ndarray,
) -> np.ndarray:

    squared_norm = np.sum(
        X * X,
        axis=1,
        keepdims=True,
    )

    dist2 = (
        squared_norm
        +
        squared_norm.T
        -
        2.0
        *
        X
        @
        X.T
    )

    dist2 = np.maximum(
        dist2,
        0.0,
    )

    return np.sqrt(
        dist2
    )


def silhouette_from_distance(
    distance_matrix: np.ndarray,
    labels: np.ndarray,
) -> float:

    labels = np.asarray(
        labels
    )

    clusters = np.unique(
        labels
    )

    values = np.zeros(
        len(
            labels
        ),
        dtype=float,
    )

    for i in range(
        len(
            labels
        )
    ):

        same = (
            labels
            ==
            labels[i]
        )

        same_count = int(
            same.sum()
        )

        if same_count > 1:

            a = float(
                distance_matrix[
                    i,
                    same,
                ].sum()
                /
                (
                    same_count
                    -
                    1
                )
            )

        else:

            a = 0.0

        b = np.inf

        for cluster in clusters:

            if cluster == labels[i]:

                continue

            other = (
                labels
                ==
                cluster
            )

            candidate = float(
                distance_matrix[
                    i,
                    other,
                ].mean()
            )

            b = min(
                b,
                candidate,
            )

        denominator = max(
            a,
            b,
        )

        if denominator > 0:

            values[i] = (
                b - a
            ) / denominator

        else:

            values[i] = 0.0

    return float(
        values.mean()
    )


# ============================================================
# ADJUSTED RAND INDEX
# ============================================================

def adjusted_rand_index(
    labels_a: np.ndarray,
    labels_b: np.ndarray,
) -> float:

    labels_a = np.asarray(
        labels_a
    )

    labels_b = np.asarray(
        labels_b
    )

    unique_a, encoded_a = np.unique(
        labels_a,
        return_inverse=True,
    )

    unique_b, encoded_b = np.unique(
        labels_b,
        return_inverse=True,
    )

    contingency = np.zeros(
        (
            len(
                unique_a
            ),
            len(
                unique_b
            ),
        ),
        dtype=int,
    )

    np.add.at(
        contingency,
        (
            encoded_a,
            encoded_b,
        ),
        1,
    )

    def choose2(
        array,
    ):

        array = np.asarray(
            array,
            dtype=float,
        )

        return (
            array
            *
            (
                array - 1.0
            )
            /
            2.0
        )

    sum_comb = float(
        choose2(
            contingency
        ).sum()
    )

    row_sums = (
        contingency.sum(
            axis=1
        )
    )

    column_sums = (
        contingency.sum(
            axis=0
        )
    )

    n = len(
        labels_a
    )

    total_pairs = (
        n
        *
        (
            n - 1
        )
        /
        2.0
    )

    if total_pairs <= 0:

        return 1.0

    expected = (
        choose2(
            row_sums
        ).sum()
        *
        choose2(
            column_sums
        ).sum()
        /
        total_pairs
    )

    maximum = (
        choose2(
            row_sums
        ).sum()
        +
        choose2(
            column_sums
        ).sum()
    ) / 2.0

    denominator = (
        maximum
        -
        expected
    )

    if abs(
        denominator
    ) <= 1e-12:

        return 1.0

    return float(
        (
            sum_comb
            -
            expected
        )
        /
        denominator
    )


# ============================================================
# SMALL-NOISE STABILITY
# ============================================================

def algorithmic_stability(
    X: np.ndarray,
    reference_labels: np.ndarray,
    k: int,
    runs: int,
    noise_sd: float,
    n_init: int,
    seed: int,
):
    """
    Small-noise perturbation stability.

    This tests algorithmic robustness only.
    It is NOT biological replicate uncertainty.
    """

    scores = []

    for run in range(
        runs
    ):

        rng = np.random.default_rng(
            seed
            +
            run
            *
            100003
        )

        perturbed = (
            X
            +
            rng.normal(
                loc=0.0,
                scale=noise_sd,
                size=X.shape,
            )
        )

        labels, _, _ = best_kmeans(
            X=perturbed,
            k=k,
            n_init=max(
                3,
                n_init // 5,
            ),
            seed=(
                seed
                +
                7_000_000
                +
                run
                *
                1009
            ),
        )

        scores.append(
            adjusted_rand_index(
                reference_labels,
                labels,
            )
        )

    return (
        float(
            np.mean(
                scores
            )
        ),

        float(
            np.min(
                scores
            )
        ),

        float(
            np.max(
                scores
            )
        ),
    )


# ============================================================
# RESPONSE LABELS
# ============================================================

def response_word(
    value: float,
    threshold: float = 0.25,
) -> str:

    if pd.isna(
        value
    ):

        return "NA"

    if value >= threshold:

        return "Up"

    if value <= -threshold:

        return "Down"

    return "Near0"


def supported_state(
    strict_value,
    direction_value,
) -> str:

    if not bool(
        strict_value
    ):

        return "NS"

    if direction_value in [
        "Up",
        "Down",
    ]:

        return str(
            direction_value
        )

    return "NS"


# ============================================================
# PREPARE ONE ADULT DIRECTION
# ============================================================

def prepare_stratum(
    df: pd.DataFrame,
    adult_direction: str,
):

    subset = (
        df[
            df[
                "RbTKO_Direction"
            ]
            ==
            adult_direction
        ]
        .copy()
        .reset_index(
            drop=True
        )
    )

    for feature in CLUSTER_FEATURES:

        subset[
            f"RG_{feature}"
        ] = rank_gaussian(
            subset[
                feature
            ]
        )

    normalized_columns = [
        f"RG_{feature}"
        for feature
        in CLUSTER_FEATURES
    ]

    subset = (
        subset
        .dropna(
            subset=normalized_columns
        )
        .reset_index(
            drop=True
        )
    )

    X = (
        subset[
            normalized_columns
        ]
        .to_numpy(
            dtype=float
        )
    )

    return (
        subset,
        X,
    )


# ============================================================
# MODEL SELECTION
# ============================================================

def evaluate_model_selection(
    X: np.ndarray,
    adult_direction: str,
    min_k: int,
    max_k: int,
    n_init: int,
    stability_runs: int,
    noise_sd: float,
    seed: int,
):

    distance_matrix = (
        pairwise_euclidean(
            X
        )
    )

    rows = []

    fits = {}

    actual_max_k = min(
        max_k,
        len(
            X
        )
        -
        1,
    )

    for k in range(
        min_k,
        actual_max_k
        +
        1,
    ):

        (
            labels,
            centers,
            inertia,
        ) = best_kmeans(
            X=X,
            k=k,
            n_init=n_init,
            seed=(
                seed
                +
                k
                *
                1009
            ),
        )

        silhouette = (
            silhouette_from_distance(
                distance_matrix,
                labels,
            )
        )

        (
            mean_ari,
            min_ari,
            max_ari,
        ) = algorithmic_stability(
            X=X,
            reference_labels=labels,
            k=k,
            runs=stability_runs,
            noise_sd=noise_sd,
            n_init=n_init,
            seed=(
                seed
                +
                k
                *
                100003
            ),
        )

        cluster_sizes = (
            np.bincount(
                labels,
                minlength=k,
            )
        )

        rows.append(
            {
                "Adult_direction":
                    adult_direction,

                "K":
                    k,

                "Silhouette":
                    silhouette,

                "Perturbation_ARI_mean":
                    mean_ari,

                "Perturbation_ARI_min":
                    min_ari,

                "Perturbation_ARI_max":
                    max_ari,

                "Inertia":
                    inertia,

                "Smallest_cluster":
                    int(
                        cluster_sizes.min()
                    ),

                "Largest_cluster":
                    int(
                        cluster_sizes.max()
                    ),
            }
        )

        fits[
            k
        ] = (
            labels,
            centers,
            inertia,
        )

    selection = pd.DataFrame(
        rows
    )

    best_k = int(
        selection.loc[
            selection[
                "Silhouette"
            ].idxmax(),
            "K",
        ]
    )

    selection[
        "Selected_by_silhouette"
    ] = (
        selection[
            "K"
        ]
        ==
        best_k
    )

    return (
        selection,
        fits,
        best_k,
    )


# ============================================================
# DETERMINISTIC CLUSTER LABELING
# ============================================================

def relabel_clusters(
    subset: pd.DataFrame,
    labels: np.ndarray,
):

    temp = (
        subset[
            [
                "E2F13_LFC",
                "CAG_LFC",
                "MAP2_LFC",
            ]
        ]
        .copy()
    )

    temp[
        "_cluster"
    ] = labels

    means = (
        temp
        .groupby(
            "_cluster"
        )[
            [
                "E2F13_LFC",
                "CAG_LFC",
                "MAP2_LFC",
            ]
        ]
        .mean()
        .sort_values(
            [
                "E2F13_LFC",
                "CAG_LFC",
                "MAP2_LFC",
            ]
        )
    )

    ordered_original_ids = (
        means.index.tolist()
    )

    mapping = {
        original:
        new_index + 1

        for new_index, original
        in enumerate(
            ordered_original_ids
        )
    }

    relabeled = np.asarray(
        [
            mapping[
                label
            ]
            for label
            in labels
        ],
        dtype=int,
    )

    return (
        relabeled,
        mapping,
    )


# ============================================================
# CLUSTER SUMMARY
# ============================================================

def build_cluster_summary(
    subset: pd.DataFrame,
    labels: np.ndarray,
    adult_direction: str,
):

    data = subset.copy()

    data[
        "Archetype_ID"
    ] = labels

    rows = []

    for cluster in sorted(
        data[
            "Archetype_ID"
        ].unique()
    ):

        group = (
            data[
                data[
                    "Archetype_ID"
                ]
                ==
                cluster
            ]
            .copy()
        )

        mean_rb = float(
            group[
                "RbTKO_LFC"
            ].mean()
        )

        mean_e2f = float(
            group[
                "E2F13_LFC"
            ].mean()
        )

        mean_cag = float(
            group[
                "CAG_LFC"
            ].mean()
        )

        mean_map2 = float(
            group[
                "MAP2_LFC"
            ].mean()
        )

        descriptive_label = (
            f"Adult{adult_direction}"
            f"__E2F{response_word(mean_e2f)}"
            f"__CAG{response_word(mean_cag)}"
            f"__MAP2{response_word(mean_map2)}"
        )

        centroid = np.array(
            [
                group[
                    "RG_E2F13_LFC"
                ].mean(),

                group[
                    "RG_CAG_LFC"
                ].mean(),

                group[
                    "RG_MAP2_LFC"
                ].mean(),
            ]
        )

        X = (
            group[
                [
                    "RG_E2F13_LFC",
                    "RG_CAG_LFC",
                    "RG_MAP2_LFC",
                ]
            ]
            .to_numpy(
                dtype=float
            )
        )

        distance = np.sqrt(
            np.sum(
                (
                    X
                    -
                    centroid
                )
                ** 2,
                axis=1,
            )
        )

        representative_order = (
            np.argsort(
                distance
            )[:10]
        )

        representative_genes = ";".join(
            group.iloc[
                representative_order
            ][
                "Gene"
            ]
            .astype(str)
            .tolist()
        )

        strongest_adult = ";".join(
            group.assign(
                _abs_adult=(
                    group[
                        "RbTKO_LFC"
                    ].abs()
                )
            )
            .sort_values(
                "_abs_adult",
                ascending=False,
            )
            .head(
                10
            )[
                "Gene"
            ]
            .astype(str)
            .tolist()
        )

        rows.append(
            {
                "Adult_direction":
                    adult_direction,

                "Archetype_ID":
                    int(
                        cluster
                    ),

                "Archetype_label":
                    descriptive_label,

                "N_genes":
                    int(
                        len(
                            group
                        )
                    ),

                "Mean_RbTKO_LFC":
                    mean_rb,

                "Median_RbTKO_LFC":
                    float(
                        group[
                            "RbTKO_LFC"
                        ].median()
                    ),

                "Mean_E2F13_LFC":
                    mean_e2f,

                "Median_E2F13_LFC":
                    float(
                        group[
                            "E2F13_LFC"
                        ].median()
                    ),

                "Mean_CAG_LFC":
                    mean_cag,

                "Median_CAG_LFC":
                    float(
                        group[
                            "CAG_LFC"
                        ].median()
                    ),

                "Mean_MAP2_LFC":
                    mean_map2,

                "Median_MAP2_LFC":
                    float(
                        group[
                            "MAP2_LFC"
                        ].median()
                    ),

                "Frac_E2F13_FDR05":
                    float(
                        group[
                            "E2F13_FDR05"
                        ]
                        .astype(bool)
                        .mean()
                    ),

                "Frac_E2F13_strict":
                    float(
                        group[
                            "E2F13_strict"
                        ]
                        .astype(bool)
                        .mean()
                    ),

                "Frac_Adult_E2F_reciprocal_strong":
                    float(
                        group[
                            "Adult_E2F_reciprocal_strong"
                        ]
                        .astype(bool)
                        .mean()
                    ),

                "Frac_CAG_FDR05":
                    float(
                        group[
                            "CAG_FDR05"
                        ]
                        .astype(bool)
                        .mean()
                    ),

                "Frac_CAG_strict":
                    float(
                        group[
                            "CAG_strict"
                        ]
                        .astype(bool)
                        .mean()
                    ),

                "Frac_MAP2_FDR05":
                    float(
                        group[
                            "MAP2_FDR05"
                        ]
                        .astype(bool)
                        .mean()
                    ),

                "Frac_MAP2_strict":
                    float(
                        group[
                            "MAP2_strict"
                        ]
                        .astype(bool)
                        .mean()
                    ),

                "Frac_Rb_loss_sign_conserved":
                    float(
                        group[
                            "Rb_loss_sign_conserved"
                        ]
                        .astype(bool)
                        .mean()
                    ),

                "Median_Adult_residual_z":
                    float(
                        group[
                            "Adult_residual_z"
                        ].median()
                    ),

                "Representative_genes":
                    representative_genes,

                "Strongest_adult_effect_genes":
                    strongest_adult,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# SUPPORTED SIGN-PATTERN CENSUS
# ============================================================

def build_supported_pattern_census(
    df: pd.DataFrame,
):

    data = df.copy()

    data[
        "E2F_supported_state"
    ] = [
        supported_state(
            strict,
            direction,
        )

        for strict, direction
        in zip(
            data[
                "E2F13_strict"
            ],
            data[
                "E2F13_Direction"
            ],
        )
    ]

    data[
        "CAG_supported_state"
    ] = [
        supported_state(
            strict,
            direction,
        )

        for strict, direction
        in zip(
            data[
                "CAG_strict"
            ],
            data[
                "CAG_Direction"
            ],
        )
    ]

    data[
        "MAP2_supported_state"
    ] = [
        supported_state(
            strict,
            direction,
        )

        for strict, direction
        in zip(
            data[
                "MAP2_strict"
            ],
            data[
                "MAP2_Direction"
            ],
        )
    ]

    data[
        "Supported_pattern"
    ] = (
        "Adult_"
        +
        data[
            "RbTKO_Direction"
        ].astype(str)
        +
        "__E2F_"
        +
        data[
            "E2F_supported_state"
        ].astype(str)
        +
        "__CAG_"
        +
        data[
            "CAG_supported_state"
        ].astype(str)
        +
        "__MAP2_"
        +
        data[
            "MAP2_supported_state"
        ].astype(str)
    )

    census = (
        data
        .groupby(
            [
                "RbTKO_Direction",
                "Supported_pattern",
            ],
            as_index=False,
        )
        .agg(
            N_genes=(
                "Gene",
                "size",
            ),

            Mean_RbTKO_LFC=(
                "RbTKO_LFC",
                "mean",
            ),

            Mean_E2F13_LFC=(
                "E2F13_LFC",
                "mean",
            ),

            Mean_CAG_LFC=(
                "CAG_LFC",
                "mean",
            ),

            Mean_MAP2_LFC=(
                "MAP2_LFC",
                "mean",
            ),

            Median_Adult_residual_z=(
                "Adult_residual_z",
                "median",
            ),
        )
        .sort_values(
            [
                "RbTKO_Direction",
                "N_genes",
            ],
            ascending=[
                True,
                False,
            ],
        )
    )

    totals = (
        census
        .groupby(
            "RbTKO_Direction"
        )[
            "N_genes"
        ]
        .transform(
            "sum"
        )
    )

    census[
        "Fraction_within_adult_direction"
    ] = (
        census[
            "N_genes"
        ]
        /
        totals
    )

    return census


# ============================================================
# CANDIDATE PLACEMENT
# ============================================================

def build_candidate_placement(
    assignments: pd.DataFrame,
):

    indexed = (
        assignments
        .set_index(
            "Gene"
        )
    )

    rows = []

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
                "Present_in_F10G_universe"
            ] = True

        else:

            row = {
                "Gene":
                    gene,

                "Present_in_F10G_universe":
                    False,
            }

        rows.append(
            row
        )

    candidate = pd.DataFrame(
        rows
    )

    first = [
        "Gene",
        "Present_in_F10G_universe",
    ]

    rest = [
        column
        for column
        in candidate.columns
        if column not in first
    ]

    return candidate[
        first
        +
        rest
    ]


# ============================================================
# FIGURES
# ============================================================

def plot_model_selection(
    selection: pd.DataFrame,
    adult_direction: str,
    output_path: Path,
):

    fig, ax = plt.subplots(
        figsize=(
            7.5,
            5.2,
        )
    )

    ax.plot(
        selection[
            "K"
        ],
        selection[
            "Silhouette"
        ],
        marker="o",
        label="Silhouette",
    )

    ax.plot(
        selection[
            "K"
        ],
        selection[
            "Perturbation_ARI_mean"
        ],
        marker="o",
        label="Perturbation ARI",
    )

    ax.set_xlabel(
        "Number of clusters (k)"
    )

    ax.set_ylabel(
        "Score"
    )

    ax.set_title(
        f"F10G model selection: adult Rb-TKO {adult_direction}"
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


def plot_centroids(
    summary: pd.DataFrame,
    adult_direction: str,
    output_path: Path,
):

    data = (
        summary[
            summary[
                "Adult_direction"
            ]
            ==
            adult_direction
        ]
        .sort_values(
            "Archetype_ID"
        )
        .copy()
    )

    if len(
        data
    ) == 0:

        return

    matrix = (
        data[
            [
                "Mean_E2F13_LFC",
                "Mean_CAG_LFC",
                "Mean_MAP2_LFC",
            ]
        ]
        .to_numpy(
            dtype=float
        )
    )

    fig, ax = plt.subplots(
        figsize=(
            7.2,
            max(
                3.5,
                0.75
                *
                len(
                    data
                )
                +
                1.5,
            ),
        )
    )

    image = ax.imshow(
        matrix,
        aspect="auto",
    )

    ax.set_xticks(
        [
            0,
            1,
            2,
        ]
    )

    ax.set_xticklabels(
        [
            "E2F1/3-DKO",
            "CAG Rb-loss",
            "MAP2 Rb-loss",
        ]
    )

    ax.set_yticks(
        np.arange(
            len(
                data
            )
        )
    )

    ax.set_yticklabels(
        [
            f"A{int(row.Archetype_ID)} "
            f"(n={int(row.N_genes)})"

            for row
            in data.itertuples()
        ]
    )

    ax.set_title(
        f"F10G raw-LFC archetype centroids: adult {adult_direction}"
    )

    fig.colorbar(
        image,
        ax=ax,
        label="Mean log2 fold change",
    )

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=300,
    )

    plt.close(
        fig
    )


# ============================================================
# SUMMARY FILE
# ============================================================

def write_summary(
    model_selection: pd.DataFrame,
    cluster_summary: pd.DataFrame,
    pattern_census: pd.DataFrame,
    candidate_placement: pd.DataFrame,
    output_path: Path,
):

    lines = []

    lines.append(
        "F10G — EMPIRICAL RESPONSE ARCHETYPES"
    )

    lines.append(
        "=" * 52
    )

    lines.append("")

    lines.append(
        "DESIGN"
    )

    lines.append(
        "-" * 52
    )

    lines.append(
        "Genes were stratified by adult Rb-TKO direction before clustering."
    )

    lines.append(
        "Clustering used rank-Gaussian normalized E2F1/3-DKO, CAG, "
        "and MAP2 responses only."
    )

    lines.append(
        "Adult Rb-TKO magnitude was retained as annotation rather than "
        "as a clustering feature."
    )

    lines.append("")

    lines.append(
        "MODEL SELECTION"
    )

    lines.append(
        "-" * 52
    )

    for adult_direction in [
        "Up",
        "Down",
    ]:

        subset = (
            model_selection[
                model_selection[
                    "Adult_direction"
                ]
                ==
                adult_direction
            ]
        )

        if len(
            subset
        ) == 0:

            continue

        chosen = (
            subset[
                subset[
                    "Selected_by_silhouette"
                ]
            ]
            .iloc[
                0
            ]
        )

        lines.append(
            f"Adult {adult_direction}: "
            f"selected k={int(chosen['K'])}; "
            f"silhouette={chosen['Silhouette']:.3f}; "
            f"perturbation ARI mean="
            f"{chosen['Perturbation_ARI_mean']:.3f}; "
            f"smallest cluster="
            f"{int(chosen['Smallest_cluster'])}."
        )

    lines.append("")

    lines.append(
        "SELECTED ARCHETYPES"
    )

    lines.append(
        "-" * 52
    )

    for _, row in (
        cluster_summary
        .sort_values(
            [
                "Adult_direction",
                "Archetype_ID",
            ]
        )
        .iterrows()
    ):

        lines.append(
            f"Adult {row['Adult_direction']} / "
            f"A{int(row['Archetype_ID'])}: "
            f"n={int(row['N_genes'])}; "
            f"E2F={row['Mean_E2F13_LFC']:.3f}; "
            f"CAG={row['Mean_CAG_LFC']:.3f}; "
            f"MAP2={row['Mean_MAP2_LFC']:.3f}; "
            f"strong E2F reciprocity="
            f"{row['Frac_Adult_E2F_reciprocal_strong']:.3f}; "
            f"Rb-loss sign conserved="
            f"{row['Frac_Rb_loss_sign_conserved']:.3f}."
        )

    lines.append("")

    lines.append(
        "LARGEST SUPPORTED SIGN PATTERNS"
    )

    lines.append(
        "-" * 52
    )

    for adult_direction in [
        "Up",
        "Down",
    ]:

        subset = (
            pattern_census[
                pattern_census[
                    "RbTKO_Direction"
                ]
                ==
                adult_direction
            ]
            .head(
                10
            )
        )

        lines.append(
            f"Adult {adult_direction}:"
        )

        for _, row in (
            subset
            .iterrows()
        ):

            lines.append(
                f"  {row['Supported_pattern']}: "
                f"n={int(row['N_genes'])} "
                f"("
                f"{100.0 * row['Fraction_within_adult_direction']:.1f}%"
                f")"
            )

    lines.append("")

    lines.append(
        "PRE-SPECIFIED CANDIDATES"
    )

    lines.append(
        "-" * 52
    )

    for _, row in (
        candidate_placement
        .iterrows()
    ):

        if not bool(
            row[
                "Present_in_F10G_universe"
            ]
        ):

            lines.append(
                f"{row['Gene']}: not present in the complete-context "
                f"adult-selected F10G universe."
            )

        else:

            lines.append(
                f"{row['Gene']}: "
                f"adult={row['RbTKO_Direction']}, "
                f"archetype=A{int(row['Archetype_ID'])}, "
                f"RbTKO={row['RbTKO_LFC']:.3f}, "
                f"E2F={row['E2F13_LFC']:.3f}, "
                f"CAG={row['CAG_LFC']:.3f}, "
                f"MAP2={row['MAP2_LFC']:.3f}."
            )

    lines.append("")

    lines.append(
        "INTERPRETATION GUARDRAILS"
    )

    lines.append(
        "-" * 52
    )

    lines.append(
        "- These are empirical response archetypes, not biological cell states."
    )

    lines.append(
        "- The analysis universe is selected on strong adult Rb-TKO response; "
        "cluster frequencies are conditional on that selection."
    )

    lines.append(
        "- Rank normalization reduces scale differences but does not remove "
        "cell-type, platform, timing, recombination, or statistical-power differences."
    )

    lines.append(
        "- Perturbation ARI measures algorithmic stability under small numerical "
        "perturbations, not experimental replicate uncertainty."
    )

    lines.append(
        "- No cluster should receive a biological name until independent "
        "functional enrichment and network tests are performed."
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
        "--input",
        default=DEFAULT_INPUT,
    )

    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
    )

    parser.add_argument(
        "--min-k",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--max-k",
        type=int,
        default=6,
    )

    parser.add_argument(
        "--n-init",
        type=int,
        default=30,
    )

    parser.add_argument(
        "--stability-runs",
        type=int,
        default=20,
    )

    parser.add_argument(
        "--noise-sd",
        type=float,
        default=0.05,
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20261023,
    )

    return parser.parse_args()


# ============================================================
# MAIN
# ============================================================

def main():

    args = parse_args()

    input_path = Path(
        args.input
    )

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

    data = pd.read_csv(
        input_path
    )

    required = [
        "Gene",
        "RbTKO_LFC",
        "E2F13_LFC",
        "CAG_LFC",
        "MAP2_LFC",
        "RbTKO_Direction",
        "E2F13_Direction",
        "CAG_Direction",
        "MAP2_Direction",
        "E2F13_FDR05",
        "E2F13_strict",
        "CAG_FDR05",
        "CAG_strict",
        "MAP2_FDR05",
        "MAP2_strict",
        "Adult_E2F_reciprocal_strong",
        "Rb_loss_sign_conserved",
        "Adult_residual_z",
    ]

    require_columns(
        data,
        required,
        "F10A residual table",
    )

    data = (
        data
        .dropna(
            subset=[
                "Gene",
                "RbTKO_LFC",
                "E2F13_LFC",
                "CAG_LFC",
                "MAP2_LFC",
                "RbTKO_Direction",
            ]
        )
        .copy()
    )

    print()
    print(
        "F10G empirical response archetypes"
    )

    print(
        "---------------------------------"
    )

    print(
        f"Complete-context input genes: "
        f"{len(data)}"
    )

    print()

    model_selection_frames = []

    assignment_frames = []

    summary_frames = []

    for (
        direction_index,
        adult_direction,
    ) in enumerate(
        [
            "Up",
            "Down",
        ]
    ):

        (
            subset,
            X,
        ) = prepare_stratum(
            data,
            adult_direction,
        )

        print(
            f"Adult {adult_direction}: "
            f"{len(subset)} genes"
        )

        (
            selection,
            fits,
            best_k,
        ) = evaluate_model_selection(
            X=X,
            adult_direction=adult_direction,
            min_k=args.min_k,
            max_k=args.max_k,
            n_init=args.n_init,
            stability_runs=args.stability_runs,
            noise_sd=args.noise_sd,
            seed=(
                args.seed
                +
                direction_index
                *
                10_000_000
            ),
        )

        model_selection_frames.append(
            selection
        )

        print(
            selection[
                [
                    "K",
                    "Silhouette",
                    "Perturbation_ARI_mean",
                    "Smallest_cluster",
                    "Selected_by_silhouette",
                ]
            ]
            .to_string(
                index=False
            )
        )

        raw_labels = (
            fits[
                best_k
            ][0]
        )

        (
            labels,
            _,
        ) = relabel_clusters(
            subset,
            raw_labels,
        )

        subset[
            "Archetype_ID"
        ] = labels

        subset[
            "Archetype_key"
        ] = (
            adult_direction
            +
            "_A"
            +
            subset[
                "Archetype_ID"
            ]
            .astype(str)
        )

        assignment_frames.append(
            subset
        )

        summary = (
            build_cluster_summary(
                subset=subset,
                labels=labels,
                adult_direction=adult_direction,
            )
        )

        summary_frames.append(
            summary
        )

        plot_model_selection(
            selection=selection,
            adult_direction=adult_direction,
            output_path=(
                figure_directory
                /
                (
                    f"F10G_model_selection_"
                    f"{adult_direction}.png"
                )
            ),
        )

        print(
            f"Selected k for adult "
            f"{adult_direction}: "
            f"{best_k}"
        )

        print()

    model_selection = (
        pd.concat(
            model_selection_frames,
            ignore_index=True,
        )
    )

    assignments = (
        pd.concat(
            assignment_frames,
            ignore_index=True,
        )
    )

    cluster_summary = (
        pd.concat(
            summary_frames,
            ignore_index=True,
        )
    )

    pattern_census = (
        build_supported_pattern_census(
            data
        )
    )

    candidate_placement = (
        build_candidate_placement(
            assignments
        )
    )

    # ========================================================
    # SAVE TABLES
    # ========================================================

    model_selection.to_csv(
        output_directory
        /
        "F10G_model_selection.csv",
        index=False,
    )

    assignments.to_csv(
        output_directory
        /
        "F10G_cluster_assignments.csv",
        index=False,
    )

    cluster_summary.to_csv(
        output_directory
        /
        "F10G_cluster_summary.csv",
        index=False,
    )

    pattern_census.to_csv(
        output_directory
        /
        "F10G_supported_pattern_census.csv",
        index=False,
    )

    candidate_placement.to_csv(
        output_directory
        /
        "F10G_candidate_placement.csv",
        index=False,
    )

    # ========================================================
    # FIGURES
    # ========================================================

    for adult_direction in [
        "Up",
        "Down",
    ]:

        plot_centroids(
            summary=cluster_summary,
            adult_direction=adult_direction,
            output_path=(
                figure_directory
                /
                (
                    f"F10G_centroids_"
                    f"{adult_direction}.png"
                )
            ),
        )

    # ========================================================
    # TEXT SUMMARY
    # ========================================================

    write_summary(
        model_selection=model_selection,
        cluster_summary=cluster_summary,
        pattern_census=pattern_census,
        candidate_placement=candidate_placement,
        output_path=(
            output_directory
            /
            "F10G_summary.txt"
        ),
    )

    # ========================================================
    # CONSOLE OUTPUT
    # ========================================================

    print(
        "Selected archetypes"
    )

    print(
        "-------------------"
    )

    display_summary = (
        cluster_summary[
            [
                "Adult_direction",
                "Archetype_ID",
                "Archetype_label",
                "N_genes",
                "Mean_E2F13_LFC",
                "Mean_CAG_LFC",
                "Mean_MAP2_LFC",
                "Frac_Adult_E2F_reciprocal_strong",
                "Frac_Rb_loss_sign_conserved",
                "Representative_genes",
            ]
        ]
    )

    print(
        display_summary
        .to_string(
            index=False
        )
    )

    print()

    print(
        "Candidate placement"
    )

    print(
        "-------------------"
    )

    candidate_columns = [
        column
        for column
        in [
            "Gene",
            "Present_in_F10G_universe",
            "RbTKO_Direction",
            "Archetype_key",
            "RbTKO_LFC",
            "E2F13_LFC",
            "E2F13_padj",
            "CAG_LFC",
            "CAG_padj",
            "MAP2_LFC",
            "MAP2_padj",
            "Adult_E2F_reciprocal_strong",
            "Adult_residual_z",
        ]
        if column
        in candidate_placement.columns
    ]

    print(
        candidate_placement[
            candidate_columns
        ]
        .to_string(
            index=False
        )
    )

    print()

    print(
        "Largest supported response patterns"
    )

    print(
        "-----------------------------------"
    )

    print(
        pattern_census[
            [
                "RbTKO_Direction",
                "Supported_pattern",
                "N_genes",
                "Fraction_within_adult_direction",
            ]
        ]
        .groupby(
            "RbTKO_Direction",
            group_keys=False,
        )
        .head(
            10
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
        "These clusters are empirical response archetypes, "
        "not biological states."
    )

    print(
        "Do not assign pathway names until enrichment is "
        "tested independently."
    )


if __name__ == "__main__":
    main()