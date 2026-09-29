#!/usr/bin/env python3
"""
06_integrated_evidence.py

Phase E: integrate all prior evidence layers without collapsing unlike evidence
into a single opaque score.

The pipeline keeps EXPERIMENTAL SEEDS and EXTERNAL CONNECTORS separate.

Experimental-seed evidence layers:
    1. dataset/mechanism recurrence from Phase 01/02 master table
    2. Phase B module structure and bridge behavior
    3. Phase C/04b directed mechanistic-path robustness

External-connector evidence layers:
    1. Phase C/04b causal robustness
    2. Phase D STRING expansion support
    3. Phase D/05b structural connector audit

Outputs are evidence matrices and evidence classes, not a single "winner".

RUN FROM Project/
-----------------
    py .\\Scripts\\06_integrated_evidence.py

OUTPUT
------
Data/Processed/Network/Integrated/
    experimental_seed_evidence_matrix.csv
    external_connector_evidence_matrix.csv
    integrated_evidence_long.csv
    integrated_diagnostics.json
    seed_evidence_class_summary.csv
    external_evidence_class_summary.csv

Figures/Network/Integrated/
    experimental_seed_evidence_counts.png
    external_connector_evidence_counts.png
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# ============================================================
# Helpers
# ============================================================

def clean_gene(value):
    if pd.isna(value):
        return None
    value = str(value).strip().upper()
    return value or None


def write_json(path, obj):
    path.write_text(
        json.dumps(
            obj,
            indent=2,
            ensure_ascii=False,
            default=str,
        ),
        encoding="utf-8",
    )


def normalize_header(value):
    return (
        str(value)
        .strip()
        .lower()
        .replace("–", "-")
        .replace("—", "-")
        .replace("_", " ")
    )


def find_column(df, candidates):
    lookup = {
        normalize_header(col): col
        for col in df.columns
    }

    for candidate in candidates:
        key = normalize_header(candidate)
        if key in lookup:
            return lookup[key]

    return None


def bool_series(df, candidates):
    col = find_column(df, candidates)

    if col is None:
        return pd.Series(
            False,
            index=df.index,
        )

    values = df[col]

    if pd.api.types.is_bool_dtype(values):
        return values.fillna(False)

    numeric = pd.to_numeric(
        values,
        errors="coerce",
    )

    text = (
        values
        .astype(str)
        .str.strip()
        .str.upper()
    )

    return (
        numeric.eq(1)
        |
        text.isin(
            [
                "YES",
                "TRUE",
                "Y",
            ]
        )
    )


def numeric_series(
    df,
    candidates,
    default=0.0,
):
    col = find_column(
        df,
        candidates,
    )

    if col is None:
        return pd.Series(
            default,
            index=df.index,
            dtype=float,
        )

    return (
        pd.to_numeric(
            df[col],
            errors="coerce",
        )
        .fillna(default)
    )


def text_series(
    df,
    candidates,
    default="",
):
    col = find_column(
        df,
        candidates,
    )

    if col is None:
        return pd.Series(
            default,
            index=df.index,
            dtype=object,
        )

    return (
        df[col]
        .fillna(default)
        .astype(str)
    )


def ensure_gene(df):
    if "Gene" not in df.columns:
        raise ValueError(
            "Input table missing Gene column."
        )

    out = df.copy()

    out["Gene"] = (
        out["Gene"]
        .map(clean_gene)
    )

    out.dropna(
        subset=["Gene"],
        inplace=True,
    )

    out.drop_duplicates(
        "Gene",
        keep="first",
        inplace=True,
    )

    return out


# ============================================================
# Evidence classes: seeds
# ============================================================

def classify_seed(row):
    """
    Descriptive evidence classes.
    They are NOT mutually exclusive biological mechanisms; this is simply
    a compact summary of how many independent analysis layers support a gene.
    """

    expression = bool(
        row["adult_embryo_shared"]
        or row["reversed"]
        or row["andrusiak_overlap"]
        or row["multi_mechanism"]
    )

    structural = bool(
        row["phaseB_nonisolated"]
        and (
            row["phaseB_cross_module_degree"] > 0
            or row["phaseB_participation"] > 0
            or row["phaseB_betweenness"] > 0
        )
    )

    mechanistic = bool(
        row["phaseC_settings_present"] >= 2
        and (
            row["phaseC_max_internal_pair_occurrence"] > 0
            or row["phaseC_required_internal_settings"] > 0
            or row["phaseC_union_nonanchor_bottleneck_settings"] > 0
        )
    )

    strong_mechanistic = bool(
        row["phaseC_settings_present"] == 4
        and (
            row["phaseC_required_internal_settings"] >= 2
            or row["phaseC_union_nonanchor_bottleneck_settings"] >= 2
        )
    )

    if expression and structural and strong_mechanistic:
        return "multi_layer_strong"

    if expression and mechanistic:
        return "experimental_plus_mechanistic"

    if expression and structural:
        return "experimental_plus_structural"

    if structural and mechanistic:
        return "structural_plus_mechanistic"

    if expression:
        return "experimental_only"

    if mechanistic:
        return "mechanistic_only"

    if structural:
        return "structural_only"

    return "limited_network_evidence"


def seed_evidence_count(row):
    flags = [
        row["adult_embryo_shared"],
        row["reversed"],
        row["andrusiak_overlap"],
        row["multi_mechanism"],
        row["phaseB_cross_module_degree"] > 0,
        row["phaseB_participation"] > 0,
        row["phaseC_settings_present"] == 4,
        row["phaseC_required_internal_settings"] > 0,
        row["phaseC_union_nonanchor_bottleneck_settings"] > 0,
    ]

    return int(
        sum(bool(value) for value in flags)
    )


# ============================================================
# Evidence classes: external connectors
# ============================================================

def classify_external(row):
    causal = bool(
        row["causal_settings_present"] == 4
        and (
            row["causal_required_internal_settings"] > 0
            or row["causal_union_nonanchor_bottleneck_settings"] > 0
        )
    )

    structural = bool(
        row["connector_class"]
        in {
            "cross_module_bridge",
            "mixed_module_isolate_bridge",
        }
    )

    unique_structural = bool(
        row["unique_component_pair_count"] > 0
        or row["marginal_seed_component_increase_on_removal"] > 0
    )

    isolate_rescue = bool(
        row["connector_class"]
        == "isolate_rescue"
    )

    if causal and structural and unique_structural:
        return "causal_and_unique_structural"

    if causal and structural:
        return "causal_and_structural"

    if causal:
        return "causal_connector"

    if structural and unique_structural:
        return "unique_structural_bridge"

    if structural:
        return "redundant_structural_bridge"

    if isolate_rescue and unique_structural:
        return "unique_isolate_rescue"

    if isolate_rescue:
        return "redundant_isolate_rescue"

    return "other_external"


def external_evidence_count(row):
    flags = [
        row["causal_settings_present"] == 4,
        row["causal_required_internal_settings"] > 0,
        row["causal_union_nonanchor_bottleneck_settings"] > 0,
        row["connector_class"]
        in {
            "cross_module_bridge",
            "mixed_module_isolate_bridge",
        },
        row["unique_component_pair_count"] > 0,
        row["marginal_seed_component_increase_on_removal"] > 0,
        row["nonisolate_module_count"] >= 2,
    ]

    return int(
        sum(bool(value) for value in flags)
    )


# ============================================================
# Figures
# ============================================================

def plot_evidence_counts(
    df,
    class_col,
    title,
    output_path,
):
    if df.empty:
        return

    counts = (
        df[class_col]
        .value_counts()
        .sort_values(
            ascending=True
        )
    )

    plt.figure(
        figsize=(
            10,
            max(
                5,
                0.55 * len(counts) + 2,
            ),
        )
    )

    plt.barh(
        counts.index,
        counts.values,
    )

    plt.xlabel(
        "Number of genes"
    )

    plt.title(title)

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    project_root = (
        Path(__file__)
        .resolve()
        .parent
        .parent
    )

    processed_network = (
        project_root
        / "Data"
        / "Processed"
        / "Network"
    )

    modules_dir = (
        processed_network
        / "Modules"
    )

    robustness_dir = (
        processed_network
        / "Mechanistic"
        / "Robustness"
    )

    global_dir = (
        processed_network
        / "GlobalExpansion"
    )

    audit_dir = (
        global_dir
        / "ConnectorAudit"
    )

    output_dir = (
        processed_network
        / "Integrated"
    )

    figures_dir = (
        project_root
        / "Figures"
        / "Network"
        / "Integrated"
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    figures_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    assignment_path = (
        modules_dir
        / "module_assignments.csv"
    )

    seed_robustness_path = (
        robustness_dir
        / "robustness_seed_consensus.csv"
    )

    all_robustness_path = (
        robustness_dir
        / "robustness_gene_consensus.csv"
    )

    accepted_path = (
        global_dir
        / "accepted_external_connectors.csv"
    )

    audit_path = (
        audit_dir
        / "connector_audit.csv"
    )

    for path in [
        assignment_path,
        seed_robustness_path,
        all_robustness_path,
        accepted_path,
        audit_path,
    ]:
        if not path.exists():
            raise FileNotFoundError(path)

    assignments = ensure_gene(
        pd.read_csv(
            assignment_path
        )
    )

    seed_robustness = ensure_gene(
        pd.read_csv(
            seed_robustness_path
        )
    )

    all_robustness = ensure_gene(
        pd.read_csv(
            all_robustness_path
        )
    )

    accepted = ensure_gene(
        pd.read_csv(
            accepted_path
        )
    )

    audit = ensure_gene(
        pd.read_csv(
            audit_path
        )
    )

    seed_genes = set(
        assignments[
            "Gene"
        ]
    )

    # ========================================================
    # Experimental seed evidence matrix
    # ========================================================

    seeds = assignments.copy()

    seeds[
        "adult_embryo_shared"
    ] = bool_series(
        seeds,
        [
            "Adult–embryo shared?",
            "Adult-embryo shared?",
        ],
    )

    direction = text_series(
        seeds,
        [
            "Direction across stages"
        ],
    )

    seeds[
        "reversed"
    ] = (
        direction
        .str.upper()
        .eq("REVERSED")
    )

    seeds[
        "andrusiak_overlap"
    ] = bool_series(
        seeds,
        [
            "Andrusiak overlap"
        ],
    )

    mechanism_count = (
        numeric_series(
            seeds,
            [
                "Mechanism count"
            ],
            default=0.0,
        )
    )

    seeds[
        "multi_mechanism"
    ] = (
        mechanism_count > 1
    )

    seeds[
        "mechanism_count"
    ] = mechanism_count

    seeds[
        "phaseB_nonisolated"
    ] = (
        seeds[
            "module_id"
        ]
        .fillna("")
        .ne("ISOLATE")
    )

    seeds[
        "phaseB_cross_module_degree"
    ] = numeric_series(
        seeds,
        [
            "cross_module_degree"
        ],
    )

    seeds[
        "phaseB_participation"
    ] = numeric_series(
        seeds,
        [
            "participation_coefficient"
        ],
    )

    seeds[
        "phaseB_betweenness"
    ] = numeric_series(
        seeds,
        [
            "betweenness_weighted"
        ],
    )

    seed_robust_cols = {
        "settings_present_count":
            "phaseC_settings_present",
        "max_internal_pair_occurrence_count":
            "phaseC_max_internal_pair_occurrence",
        "settings_required_internal_for_at_least_one_pair":
            "phaseC_required_internal_settings",
        "max_required_internal_pair_count":
            "phaseC_max_required_internal_pair_count",
        "settings_union_nonanchor_bottleneck":
            "phaseC_union_nonanchor_bottleneck_settings",
        "max_union_nonanchor_pairs_lost":
            "phaseC_max_union_nonanchor_pairs_lost",
    }

    available = [
        "Gene"
    ] + [
        col
        for col
        in seed_robust_cols
        if col
        in seed_robustness.columns
    ]

    seed_mech = (
        seed_robustness[
            available
        ]
        .rename(
            columns=
                seed_robust_cols
        )
    )

    seeds = seeds.merge(
        seed_mech,
        on="Gene",
        how="left",
    )

    for column in [
        "phaseC_settings_present",
        "phaseC_max_internal_pair_occurrence",
        "phaseC_required_internal_settings",
        "phaseC_max_required_internal_pair_count",
        "phaseC_union_nonanchor_bottleneck_settings",
        "phaseC_max_union_nonanchor_pairs_lost",
    ]:
        if column not in seeds.columns:
            seeds[column] = 0

        seeds[column] = (
            pd.to_numeric(
                seeds[column],
                errors="coerce",
            )
            .fillna(0)
        )

    seeds[
        "integrated_evidence_class"
    ] = seeds.apply(
        classify_seed,
        axis=1,
    )

    seeds[
        "independent_evidence_flag_count"
    ] = seeds.apply(
        seed_evidence_count,
        axis=1,
    )

    # Sort descriptively, not as an overall biological verdict.
    seeds.sort_values(
        [
            "independent_evidence_flag_count",
            "phaseC_union_nonanchor_bottleneck_settings",
            "phaseC_required_internal_settings",
            "adult_embryo_shared",
            "andrusiak_overlap",
            "phaseB_participation",
        ],
        ascending=[
            False,
            False,
            False,
            False,
            False,
            False,
        ],
        inplace=True,
    )

    seed_output_cols = [
        "Gene",
        "module_id",
        "component_id_phaseB",
        "adult_embryo_shared",
        "reversed",
        "andrusiak_overlap",
        "mechanism_count",
        "multi_mechanism",
        "phaseB_nonisolated",
        "phaseB_cross_module_degree",
        "phaseB_participation",
        "phaseB_betweenness",
        "phaseC_settings_present",
        "phaseC_max_internal_pair_occurrence",
        "phaseC_required_internal_settings",
        "phaseC_max_required_internal_pair_count",
        "phaseC_union_nonanchor_bottleneck_settings",
        "phaseC_max_union_nonanchor_pairs_lost",
        "integrated_evidence_class",
        "independent_evidence_flag_count",
    ]

    seed_output_cols = [
        col
        for col
        in seed_output_cols
        if col in seeds.columns
    ]

    seed_matrix = (
        seeds[
            seed_output_cols
        ]
        .copy()
    )

    seed_matrix.to_csv(
        output_dir
        / "experimental_seed_evidence_matrix.csv",
        index=False,
    )

    # ========================================================
    # External connector evidence matrix
    # ========================================================

    external_genes = (
        set(
            accepted["Gene"]
        )
        |
        set(
            audit["Gene"]
        )
    ) - seed_genes

    external = pd.DataFrame({
        "Gene":
            sorted(
                external_genes
            )
    })

    audit_cols = [
        col
        for col in [
            "Gene",
            "connector_class",
            "seed_neighbor_count",
            "isolate_seed_neighbor_count",
            "nonisolate_seed_neighbor_count",
            "nonisolate_module_count",
            "baseline_seed_component_count",
            "unique_component_pair_count",
            "shared_component_pair_count",
            "marginal_seed_component_increase_on_removal",
            "largest_seed_component_drop_on_removal",
            "structurally_distinct",
            "acceptance_reason",
            "max_string_score",
            "max_non_text_evidence",
            "text_mining_dominant_edge_fraction",
        ]
        if col in audit.columns
    ]

    external = external.merge(
        audit[
            audit_cols
        ],
        on="Gene",
        how="left",
    )

    robust_external = (
        all_robustness[
            ~all_robustness[
                "Gene"
            ].isin(
                seed_genes
            )
        ]
        .copy()
    )

    causal_cols = {
        "settings_present_count":
            "causal_settings_present",
        "max_internal_pair_occurrence_count":
            "causal_max_internal_pair_occurrence",
        "settings_required_internal_for_at_least_one_pair":
            "causal_required_internal_settings",
        "max_required_internal_pair_count":
            "causal_max_required_internal_pair_count",
        "settings_union_nonanchor_bottleneck":
            "causal_union_nonanchor_bottleneck_settings",
        "max_union_nonanchor_pairs_lost":
            "causal_max_union_nonanchor_pairs_lost",
    }

    available = [
        "Gene"
    ] + [
        col
        for col
        in causal_cols
        if col in robust_external.columns
    ]

    external = external.merge(
        robust_external[
            available
        ].rename(
            columns=
                causal_cols
        ),
        on="Gene",
        how="left",
    )

    numeric_external_cols = [
        "seed_neighbor_count",
        "isolate_seed_neighbor_count",
        "nonisolate_seed_neighbor_count",
        "nonisolate_module_count",
        "baseline_seed_component_count",
        "unique_component_pair_count",
        "shared_component_pair_count",
        "marginal_seed_component_increase_on_removal",
        "largest_seed_component_drop_on_removal",
        "causal_settings_present",
        "causal_max_internal_pair_occurrence",
        "causal_required_internal_settings",
        "causal_max_required_internal_pair_count",
        "causal_union_nonanchor_bottleneck_settings",
        "causal_max_union_nonanchor_pairs_lost",
    ]

    for column in numeric_external_cols:
        if column not in external.columns:
            external[column] = 0

        external[column] = (
            pd.to_numeric(
                external[column],
                errors="coerce",
            )
            .fillna(0)
        )

    if "connector_class" not in external.columns:
        external[
            "connector_class"
        ] = "other_external"

    external[
        "connector_class"
    ] = (
        external[
            "connector_class"
        ]
        .fillna(
            "causal_only"
        )
    )

    external[
        "integrated_evidence_class"
    ] = external.apply(
        classify_external,
        axis=1,
    )

    external[
        "independent_evidence_flag_count"
    ] = external.apply(
        external_evidence_count,
        axis=1,
    )

    external.sort_values(
        [
            "independent_evidence_flag_count",
            "causal_union_nonanchor_bottleneck_settings",
            "causal_required_internal_settings",
            "unique_component_pair_count",
            "marginal_seed_component_increase_on_removal",
            "nonisolate_module_count",
        ],
        ascending=[
            False,
            False,
            False,
            False,
            False,
            False,
        ],
        inplace=True,
    )

    external_output_cols = [
        "Gene",
        "connector_class",
        "acceptance_reason",
        "seed_neighbor_count",
        "isolate_seed_neighbor_count",
        "nonisolate_seed_neighbor_count",
        "nonisolate_module_count",
        "baseline_seed_component_count",
        "unique_component_pair_count",
        "shared_component_pair_count",
        "marginal_seed_component_increase_on_removal",
        "largest_seed_component_drop_on_removal",
        "max_string_score",
        "max_non_text_evidence",
        "text_mining_dominant_edge_fraction",
        "causal_settings_present",
        "causal_max_internal_pair_occurrence",
        "causal_required_internal_settings",
        "causal_max_required_internal_pair_count",
        "causal_union_nonanchor_bottleneck_settings",
        "causal_max_union_nonanchor_pairs_lost",
        "integrated_evidence_class",
        "independent_evidence_flag_count",
    ]

    external_output_cols = [
        col
        for col
        in external_output_cols
        if col in external.columns
    ]

    external_matrix = (
        external[
            external_output_cols
        ]
        .copy()
    )

    external_matrix.to_csv(
        output_dir
        / "external_connector_evidence_matrix.csv",
        index=False,
    )

    # ========================================================
    # Long-format evidence table
    # ========================================================

    long_rows = []

    for _, row in seed_matrix.iterrows():
        gene = row["Gene"]

        evidence = {
            "adult_embryo_shared":
                bool(
                    row.get(
                        "adult_embryo_shared",
                        False,
                    )
                ),
            "developmental_reversal":
                bool(
                    row.get(
                        "reversed",
                        False,
                    )
                ),
            "andrusiak_overlap":
                bool(
                    row.get(
                        "andrusiak_overlap",
                        False,
                    )
                ),
            "multi_mechanism":
                bool(
                    row.get(
                        "multi_mechanism",
                        False,
                    )
                ),
            "phaseB_cross_module":
                (
                    row.get(
                        "phaseB_cross_module_degree",
                        0,
                    ) > 0
                ),
            "phaseC_required_internal":
                (
                    row.get(
                        "phaseC_required_internal_settings",
                        0,
                    ) > 0
                ),
            "phaseC_union_bottleneck":
                (
                    row.get(
                        "phaseC_union_nonanchor_bottleneck_settings",
                        0,
                    ) > 0
                ),
        }

        for evidence_type, present in evidence.items():
            long_rows.append({
                "Gene": gene,
                "node_origin":
                    "experimental_seed",
                "evidence_type":
                    evidence_type,
                "present":
                    bool(present),
            })

    for _, row in external_matrix.iterrows():
        gene = row["Gene"]

        evidence = {
            "causal_present_all_settings":
                (
                    row.get(
                        "causal_settings_present",
                        0,
                    ) == 4
                ),
            "causal_required_internal":
                (
                    row.get(
                        "causal_required_internal_settings",
                        0,
                    ) > 0
                ),
            "causal_union_bottleneck":
                (
                    row.get(
                        "causal_union_nonanchor_bottleneck_settings",
                        0,
                    ) > 0
                ),
            "cross_module_bridge":
                (
                    row.get(
                        "connector_class",
                        ""
                    )
                    == "cross_module_bridge"
                ),
            "mixed_module_isolate_bridge":
                (
                    row.get(
                        "connector_class",
                        ""
                    )
                    == "mixed_module_isolate_bridge"
                ),
            "unique_component_pair":
                (
                    row.get(
                        "unique_component_pair_count",
                        0,
                    ) > 0
                ),
            "marginal_component_contribution":
                (
                    row.get(
                        "marginal_seed_component_increase_on_removal",
                        0,
                    ) > 0
                ),
        }

        for evidence_type, present in evidence.items():
            long_rows.append({
                "Gene": gene,
                "node_origin":
                    "external_connector",
                "evidence_type":
                    evidence_type,
                "present":
                    bool(present),
            })

    long_df = pd.DataFrame(
        long_rows
    )

    long_df.to_csv(
        output_dir
        / "integrated_evidence_long.csv",
        index=False,
    )

    # ========================================================
    # Class summaries
    # ========================================================

    seed_class_summary = (
        seed_matrix[
            "integrated_evidence_class"
        ]
        .value_counts()
        .rename_axis(
            "integrated_evidence_class"
        )
        .reset_index(
            name="gene_count"
        )
    )

    seed_class_summary.to_csv(
        output_dir
        / "seed_evidence_class_summary.csv",
        index=False,
    )

    external_class_summary = (
        external_matrix[
            "integrated_evidence_class"
        ]
        .value_counts()
        .rename_axis(
            "integrated_evidence_class"
        )
        .reset_index(
            name="gene_count"
        )
    )

    external_class_summary.to_csv(
        output_dir
        / "external_evidence_class_summary.csv",
        index=False,
    )

    # ========================================================
    # Figures
    # ========================================================

    plot_evidence_counts(
        seed_matrix,
        "integrated_evidence_class",
        "Experimental seed genes by integrated evidence class",
        figures_dir
        / "experimental_seed_evidence_counts.png",
    )

    plot_evidence_counts(
        external_matrix,
        "integrated_evidence_class",
        "External connectors by integrated evidence class",
        figures_dir
        / "external_connector_evidence_counts.png",
    )

    # ========================================================
    # Diagnostics
    # ========================================================

    diagnostics = {
        "experimental_seed_count":
            len(
                seed_matrix
            ),
        "external_connector_count":
            len(
                external_matrix
            ),
        "seed_classes":
            seed_class_summary
            .set_index(
                "integrated_evidence_class"
            )[
                "gene_count"
            ]
            .to_dict(),
        "external_classes":
            external_class_summary
            .set_index(
                "integrated_evidence_class"
            )[
                "gene_count"
            ]
            .to_dict(),
        "principle":
            (
                "Experimental seeds and external connectors are "
                "kept separate because their evidence types are "
                "not commensurate."
            ),
    }

    write_json(
        output_dir
        / "integrated_diagnostics.json",
        diagnostics,
    )

    # ========================================================
    # Console
    # ========================================================

    print()
    print(
        "R-Noel integrated evidence - Phase E"
    )
    print(
        "===================================="
    )

    print()
    print(
        "Experimental seed evidence classes"
    )
    print(
        "----------------------------------"
    )

    print(
        seed_class_summary.to_string(
            index=False
        )
    )

    print()
    print(
        "Experimental seeds with the most "
        "independent evidence flags"
    )
    print(
        "-----------------------------------"
    )

    print(
        seed_matrix[
            [
                "Gene",
                "module_id",
                "adult_embryo_shared",
                "reversed",
                "andrusiak_overlap",
                "phaseB_cross_module_degree",
                "phaseC_required_internal_settings",
                "phaseC_union_nonanchor_bottleneck_settings",
                "integrated_evidence_class",
                "independent_evidence_flag_count",
            ]
        ]
        .head(25)
        .to_string(
            index=False
        )
    )

    print()
    print(
        "External connector evidence classes"
    )
    print(
        "-----------------------------------"
    )

    print(
        external_class_summary.to_string(
            index=False
        )
    )

    print()
    print(
        "External connectors with the most "
        "independent evidence flags"
    )
    print(
        "------------------------------------"
    )

    print(
        external_matrix[
            [
                "Gene",
                "connector_class",
                "causal_required_internal_settings",
                "causal_union_nonanchor_bottleneck_settings",
                "nonisolate_module_count",
                "unique_component_pair_count",
                "marginal_seed_component_increase_on_removal",
                "integrated_evidence_class",
                "independent_evidence_flag_count",
            ]
        ]
        .head(25)
        .to_string(
            index=False
        )
    )

    print()
    print(
        "Phase E complete."
    )

    print(
        f"Processed: {output_dir}"
    )

    print(
        f"Figures: {figures_dir}"
    )

    print()
    print(
        "This stage creates evidence classes, not a single "
        "cross-origin gene ranking. Experimental candidates "
        "and external connectors remain separate."
    )


if __name__ == "__main__":
    main()
