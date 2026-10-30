#!/usr/bin/env python3

"""
F6 — Integrated old + new mechanistic synthesis
================================================

Purpose
-------
Integrate the original R-Noel network analysis with the new Fong
downstream pipeline without collapsing heterogeneous evidence into
one numerical "importance score".

Primary evidence layers
-----------------------
OLD PIPELINE
    - curated 126-gene network master
    - adult/embryonic overlap
    - developmental reversal
    - Andrusiak overlap
    - old mechanisms
    - old STRING/module/mechanistic evidence where available

NEW FONG PIPELINE
    - current strict Fong DEG universe
    - F2 cell-fate/stress gene universe
    - F4 CORE and EXTENDED modules
    - F4 cross-module bridge metrics
    - F4 INTERFACE -> CORE bridge structure

DIRECTED MECHANISTIC EVIDENCE
    - corrected F5 support tiers
    - direct source -> target edges
    - shortest experimentally anchored routes
    - exact required-node evidence

PROMOTER EVIDENCE
    - revised Julian E2F4-vs-E2F3 significant promoter-gene set
    - optional Julian binding/expression integration workbook

IMPORTANT
---------
This script DOES NOT calculate a composite biological score.

Each evidence type remains a separate column.

Primary integration classes
---------------------------
    - Old + new + promoter convergence
    - Old + new convergence
    - New downstream CORE
    - New stress-interface
    - Old strict-Fong gene outside F2 cell-fate/stress panel
    - Old candidate absent from current strict Fong/F2

External prior-knowledge connectors are kept in a separate table.


Required inputs
---------------
Project/
    cell death project-analysis-ready.xlsx

Data/Raw/microarray/Fong2022/
    Fong 2022 L2FC1 microarray data.xlsx

Data/Processed/FongDownstream/F2_CellFate/
    cell_fate_master.csv

Data/Processed/FongDownstream/F4_Modules/
    CORE/module_assignments.csv
    CORE/bridge_metrics.csv
    EXTENDED/module_assignments.csv
    EXTENDED/bridge_metrics.csv
    EXTENDED/interface_core_bridges.csv

Data/Processed/FongDownstream/F5_Mechanistic/Corrected/
    anchor_audit.csv
    pair_robustness.csv
    representative_paths.csv
    exact_required_nodes.csv


Optional old-network integration files
--------------------------------------
Data/Processed/Network/Integrated/
    experimental_seed_evidence_matrix.csv
    external_connector_evidence_matrix.csv

Data/Processed/Network/FinalSynthesis/
    final_experimental_candidates.csv
    final_external_connectors.csv


Julian files are discovered recursively, especially:
    E2F4_vs_E2F3_significant_genes.csv
    binding_vs_TKO_expression.xlsx


Outputs
-------
Data/Processed/FongDownstream/F6_Integrated/

    F6_integrated_evidence.xlsx

    experimental_gene_evidence_matrix.csv
    directed_route_summary.csv
    candidate_panels.csv
    mechanistic_hypotheses.csv
    mechanistic_model_edges.csv
    external_prior_knowledge_nodes.csv
    integration_summary.csv
    julian_input_audit.csv
    F6_metadata.json


Run
---
From Project/:

    py .\\Scripts\\F6_integrate_old_new.py
"""

from __future__ import annotations

from pathlib import Path
from collections import Counter, defaultdict
from typing import Any
import json
import re

import pandas as pd

from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]


OLD_MASTER_FILE = (
    PROJECT_ROOT
    / "cell death project-analysis-ready.xlsx"
)


STRICT_FONG_FILE = (
    PROJECT_ROOT
    / "Data"
    / "Raw"
    / "microarray"
    / "Fong2022"
    / "Fong 2022 L2FC1 microarray data.xlsx"
)


OLD_INTEGRATED_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "Network"
    / "Integrated"
)


OLD_FINAL_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "Network"
    / "FinalSynthesis"
)


F2_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F2_CellFate"
)


F4_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F4_Modules"
)


F5_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F5_Mechanistic"
    / "Corrected"
)


OUTPUT_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F6_Integrated"
)


OUTPUT_XLSX = (
    OUTPUT_DIR
    / "F6_integrated_evidence.xlsx"
)


OUTPUT_METADATA = (
    OUTPUT_DIR
    / "F6_metadata.json"
)


# ============================================================
# Required inputs
# ============================================================

F2_MASTER_FILE = (
    F2_DIR
    / "cell_fate_master.csv"
)


CORE_ASSIGN_FILE = (
    F4_DIR
    / "CORE"
    / "module_assignments.csv"
)


CORE_BRIDGE_FILE = (
    F4_DIR
    / "CORE"
    / "bridge_metrics.csv"
)


EXT_ASSIGN_FILE = (
    F4_DIR
    / "EXTENDED"
    / "module_assignments.csv"
)


EXT_BRIDGE_FILE = (
    F4_DIR
    / "EXTENDED"
    / "bridge_metrics.csv"
)


INTERFACE_BRIDGE_FILE = (
    F4_DIR
    / "EXTENDED"
    / "interface_core_bridges.csv"
)


F5_ANCHOR_FILE = (
    F5_DIR
    / "anchor_audit.csv"
)


F5_PAIR_FILE = (
    F5_DIR
    / "pair_robustness.csv"
)


F5_PATH_FILE = (
    F5_DIR
    / "representative_paths.csv"
)


F5_REQUIRED_FILE = (
    F5_DIR
    / "exact_required_nodes.csv"
)


# ============================================================
# Optional old outputs
# ============================================================

OLD_SEED_EVIDENCE_FILE = (
    OLD_INTEGRATED_DIR
    / "experimental_seed_evidence_matrix.csv"
)


OLD_EXTERNAL_EVIDENCE_FILE = (
    OLD_INTEGRATED_DIR
    / "external_connector_evidence_matrix.csv"
)


OLD_FINAL_CANDIDATES_FILE = (
    OLD_FINAL_DIR
    / "final_experimental_candidates.csv"
)


OLD_FINAL_EXTERNAL_FILE = (
    OLD_FINAL_DIR
    / "final_external_connectors.csv"
)


# ============================================================
# F4 module labels
#
# These correspond to the deterministic F4 run using:
#     Louvain
#     seed = 42
#     resolution = 1.0
# ============================================================

CORE_MODULE_LABELS = {

    "M01":
        "Inflammatory / cytokine / immune fate signaling",

    "M02":
        "Trophic / RTK / ECM-support signaling",

    "M03":
        "Mitosis / chromosome segregation",

    "M04":
        "Stress / death-receptor / interferon-associated",

    "M05":
        "DNA repair / homologous recombination",

    "M06":
        "E2F-G1/S / apoptotic-effector module",

    "M07":
        "Regulated necrosis / innate death",

    "M08":
        "NR4A nuclear-receptor response",

    "M09":
        "ABCC9-KCNJ8 potassium-channel pair",

    "M10":
        "EMP3-TGFBR2 pair",

    "M11":
        "MAG-MAL myelin-associated pair",

    "ISOLATE":
        "STRING isolate",
}


EXT_MODULE_LABELS = {

    "M01":
        "Inflammatory / stress / regulated-death program",

    "M02":
        "Homologous recombination / DNA repair",

    "M03":
        "Cell cycle / mitosis / E2F-associated fate machinery",

    "M04":
        "DNA replication / ATR replication-stress program",

    "M05":
        "RTK-MAPK-PI3K trophic / survival signaling",

    "M06":
        "Innate antiviral / interferon response",

    "M07":
        "NR4A nuclear-receptor response",

    "M08":
        "ABCC9-KCNJ8 potassium-channel pair",

    "M09":
        "FOSL2-NFATC2 transcriptional pair",

    "M10":
        "HAP1-SOD3 pair",

    "M11":
        "MAG-MAL myelin-associated pair",

    "ISOLATE":
        "STRING isolate",
}


# ============================================================
# Helpers
# ============================================================

def normalize_gene(
    value: Any,
) -> str | None:

    if pd.isna(value):
        return None

    value = str(value).strip()

    if not value:
        return None

    if value.upper() in {
        "NA",
        "N/A",
        "NONE",
        "NULL",
        "NAN",
    }:
        return None

    return value.upper()


def looks_like_gene_symbol(
    gene: str | None,
) -> bool:

    if gene is None:
        return False

    return (
        re.fullmatch(
            r"[A-Z0-9][A-Z0-9._-]*",
            gene,
        )
        is not None
    )


def truthy(
    value: Any,
) -> bool:

    if isinstance(
        value,
        bool,
    ):
        return value

    if pd.isna(value):
        return False

    if isinstance(
        value,
        (int, float),
    ):
        return value != 0

    return (
        str(value)
        .strip()
        .upper()
        in {
            "TRUE",
            "YES",
            "Y",
            "1",
        }
    )


def safe_number(
    value: Any,
) -> float | None:

    try:

        if pd.isna(value):
            return None

        return float(value)

    except Exception:
        return None


def find_column(
    df: pd.DataFrame,
    candidates: list[str],
) -> str | None:

    lookup = {
        str(column)
        .strip()
        .lower():
            column
        for column
        in df.columns
    }

    for candidate in candidates:

        key = (
            candidate
            .strip()
            .lower()
        )

        if key in lookup:
            return lookup[key]

    return None


def join_unique(
    values,
    separator: str = "; ",
) -> str:

    output = []

    seen = set()

    for value in values:

        if pd.isna(value):
            continue

        value = str(value).strip()

        if not value:
            continue

        if value not in seen:

            seen.add(value)
            output.append(value)

    return separator.join(
        output
    )


def optional_csv(
    path: Path,
) -> pd.DataFrame:

    if not path.exists():
        return pd.DataFrame()

    return pd.read_csv(
        path,
        low_memory=False,
    )


def normalize_gene_column(
    df: pd.DataFrame,
) -> pd.DataFrame:

    result = df.copy()

    gene_col = find_column(
        result,
        [
            "Gene",
            "gene",
            "gene_symbol",
            "Gene Symbol",
            "symbol",
            "SYMBOL",
            "genesymbol",
        ],
    )

    if gene_col is None:

        raise ValueError(
            "Could not identify a gene-symbol column.\n"
            f"Columns: {list(result.columns)}"
        )

    if gene_col != "Gene":

        result = result.rename(
            columns={
                gene_col:
                    "Gene"
            }
        )

    result[
        "Gene"
    ] = (
        result[
            "Gene"
        ]
        .map(
            normalize_gene
        )
    )

    result = result[
        result[
            "Gene"
        ].notna()
    ].copy()

    return result


def merge_if_available(
    left: pd.DataFrame,
    right: pd.DataFrame,
) -> pd.DataFrame:

    if right.empty:
        return left

    if "Gene" not in right.columns:
        return left

    return left.merge(
        right,
        on="Gene",
        how="left",
    )


def recursive_file_search(
    suffix: str,
    required_tokens: list[str],
) -> list[Path]:

    matches = []

    for path in PROJECT_ROOT.rglob(
        f"*{suffix}"
    ):

        name = path.name.lower()

        if all(
            token.lower()
            in name
            for token
            in required_tokens
        ):

            matches.append(path)

    return sorted(
        matches,
        key=lambda path: (
            len(path.parts),
            str(path),
        ),
    )


# ============================================================
# Validate inputs
# ============================================================

print()
print("=" * 82)
print("PHASE F6 — INTEGRATED OLD + NEW MECHANISTIC SYNTHESIS")
print("=" * 82)


required_files = [

    OLD_MASTER_FILE,
    STRICT_FONG_FILE,

    F2_MASTER_FILE,

    CORE_ASSIGN_FILE,
    CORE_BRIDGE_FILE,

    EXT_ASSIGN_FILE,
    EXT_BRIDGE_FILE,

    INTERFACE_BRIDGE_FILE,

    F5_ANCHOR_FILE,
    F5_PAIR_FILE,
    F5_PATH_FILE,
    F5_REQUIRED_FILE,
]


for path in required_files:

    if not path.exists():

        raise FileNotFoundError(
            f"Required input not found:\n"
            f"{path}"
        )


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 1. Old network master
# ============================================================

print()
print(
    "[F6] Loading old curated network master..."
)


old_raw = pd.read_excel(
    OLD_MASTER_FILE,
    sheet_name="network master",
)


old_raw[
    "Gene"
] = (
    old_raw[
        "Gene"
    ]
    .map(
        normalize_gene
    )
)


old = old_raw[
    old_raw[
        "Gene"
    ]
    .map(
        looks_like_gene_symbol
    )
].copy()


old = old.drop_duplicates(
    "Gene"
)


old_annotation = pd.DataFrame(
    {
        "Gene":
            old[
                "Gene"
            ],

        "Old_master":
            True,
    }
)


old_column_map = {

    "fong":
        "Old_Fong_membership",

    "cag":
        "Old_CAG_membership",

    "map2":
        "Old_MAP2_membership",

    "expr_fong":
        "Old_expr_Fong",

    "expr_cag":
        "Old_expr_CAG",

    "expr_map2":
        "Old_expr_MAP2",

    "Mechanism(s)":
        "Old_mechanisms",

    "Adult-embryo shared?":
        "Old_adult_embryo_shared_raw",

    "Embryo overlap":
        "Old_embryo_overlap",

    "Direction across stages":
        "Old_direction_across_stages",

    "Andrusiak overlap":
        "Old_Andrusiak_overlap_raw",

    "Top-distinct":
        "Old_top_distinct_raw",

    "Doctor highlight":
        "Old_doctor_highlight_raw",
}


for old_column, new_column in (
    old_column_map.items()
):

    if old_column in old.columns:

        old_annotation[
            new_column
        ] = (
            old[
                old_column
            ].values
        )


for raw_column, bool_column in [

    (
        "Old_adult_embryo_shared_raw",
        "Old_adult_embryo_shared",
    ),

    (
        "Old_Andrusiak_overlap_raw",
        "Old_Andrusiak_overlap",
    ),

    (
        "Old_top_distinct_raw",
        "Old_top_distinct",
    ),

    (
        "Old_doctor_highlight_raw",
        "Old_doctor_highlight",
    ),
]:

    if raw_column in old_annotation.columns:

        old_annotation[
            bool_column
        ] = (
            old_annotation[
                raw_column
            ]
            .map(
                truthy
            )
        )


if (
    "Old_direction_across_stages"
    in old_annotation.columns
):

    old_annotation[
        "Old_developmental_reversal"
    ] = (
        old_annotation[
            "Old_direction_across_stages"
        ]
        .astype(str)
        .str.upper()
        .eq(
            "REVERSED"
        )
    )


print(
    f"[F6] Valid old-master genes: "
    f"{len(old_annotation):,}"
)


# ============================================================
# 2. Current audited strict-Fong universe
# ============================================================

print()
print(
    "[F6] Loading current strict Fong DEG universe..."
)


strict_fong_raw = pd.read_excel(
    STRICT_FONG_FILE,
    sheet_name="Sig L2FC1",
)


strict_gene_col = find_column(
    strict_fong_raw,
    [
        "gene symbol",
        "Gene",
        "gene",
        "symbol",
        "SYMBOL",
    ],
)


# Known Fong workbook layout:
# third Excel column contains the gene symbol.
if strict_gene_col is None:

    strict_gene_col = (
        strict_fong_raw.columns[
            2
        ]
    )


strict_fong_genes = {
    normalize_gene(value)
    for value
    in strict_fong_raw[
        strict_gene_col
    ]
}


strict_fong_genes = {
    gene
    for gene
    in strict_fong_genes
    if looks_like_gene_symbol(
        gene
    )
}


print(
    f"[F6] Valid strict-Fong genes: "
    f"{len(strict_fong_genes):,}"
)


# ============================================================
# 3. New F2 experimental universe
# ============================================================

print()
print(
    "[F6] Loading F2 Fong cell-fate/stress universe..."
)


f2 = pd.read_csv(
    F2_MASTER_FILE
)


f2 = normalize_gene_column(
    f2
)


f2 = f2.drop_duplicates(
    "Gene"
)


f2_annotation = pd.DataFrame(
    {
        "Gene":
            f2[
                "Gene"
            ],

        "Fong_F2":
            True,
    }
)


f2_map = {

    "Direction":
        "Fong_direction",

    "log2FoldChange":
        "Fong_log2FC",

    "padj":
        "Fong_padj",

    "F2_tier":
        "F2_tier",

    "F2_categories":
        "F2_categories",

    "Old_126_gene_master":
        "F2_old_master_flag",
}


for source_column, target_column in (
    f2_map.items()
):

    if source_column in f2.columns:

        f2_annotation[
            target_column
        ] = (
            f2[
                source_column
            ].values
        )


print(
    f"[F6] F2 genes: "
    f"{len(f2_annotation):,}"
)


# ============================================================
# 4. Experimental integration universe
# ============================================================

old_genes = set(
    old_annotation[
        "Gene"
    ]
)


f2_genes = set(
    f2_annotation[
        "Gene"
    ]
)


experimental_genes = sorted(
    old_genes
    |
    f2_genes
)


master = pd.DataFrame(
    {
        "Gene":
            experimental_genes
    }
)


master = master.merge(
    old_annotation,
    on="Gene",
    how="left",
)


master = master.merge(
    f2_annotation,
    on="Gene",
    how="left",
)


master[
    "Old_master"
] = (
    master[
        "Old_master"
    ]
    .map(
        truthy
    )
)


master[
    "Fong_F2"
] = (
    master[
        "Fong_F2"
    ]
    .map(
        truthy
    )
)


# IMPORTANT:
# current audited strict-Fong membership, not historical workbook flag.
master[
    "Fong_strict_DE"
] = (
    master[
        "Gene"
    ]
    .isin(
        strict_fong_genes
    )
)


# ============================================================
# 5. Old integrated network evidence
# ============================================================

print()
print(
    "[F6] Loading old integrated network evidence where available..."
)


old_seed_evidence = optional_csv(
    OLD_SEED_EVIDENCE_FILE
)


if not old_seed_evidence.empty:

    old_seed_evidence = (
        normalize_gene_column(
            old_seed_evidence
        )
    )


    old_keep_map = {

        "module_id":
            "Old_network_module",

        "component_id_phaseB":
            "Old_component",

        "phaseB_nonisolated":
            "Old_STRING_nonisolated",

        "phaseB_cross_module_degree":
            "Old_cross_module_degree",

        "phaseB_participation":
            "Old_participation",

        "phaseB_betweenness":
            "Old_betweenness",

        "phaseC_settings_present":
            "Old_OmniPath_settings_present",

        "phaseC_required_internal_settings":
            "Old_OmniPath_required_settings",

        "phaseC_union_nonanchor_bottleneck_settings":
            "Old_OmniPath_bottleneck_settings",

        "integrated_evidence_class":
            "Old_integrated_evidence_class",

        "independent_evidence_flag_count":
            "Old_evidence_feature_count_historical",
    }


    keep = [
        "Gene"
    ]


    rename = {}


    for source, target in (
        old_keep_map.items()
    ):

        if source in old_seed_evidence.columns:

            keep.append(
                source
            )

            rename[
                source
            ] = target


    old_seed_evidence = (
        old_seed_evidence[
            keep
        ]
        .rename(
            columns=rename
        )
    )


    master = merge_if_available(
        master,
        old_seed_evidence,
    )


old_final = optional_csv(
    OLD_FINAL_CANDIDATES_FILE
)


if not old_final.empty:

    old_final = (
        normalize_gene_column(
            old_final
        )
    )


    final_columns = [
        column
        for column
        in [
            "Gene",
            "biological_theme",
            "followup_category",
        ]
        if column
        in old_final.columns
    ]


    old_final = (
        old_final[
            final_columns
        ]
        .rename(
            columns={
                "biological_theme":
                    "Old_working_theme",

                "followup_category":
                    "Old_followup_category",
            }
        )
    )


    master = merge_if_available(
        master,
        old_final,
    )


# ============================================================
# 6. F4 module and bridge evidence
# ============================================================

print()
print(
    "[F6] Integrating F4 module and bridge evidence..."
)


def load_f4_network(
    assignment_file: Path,
    bridge_file: Path,
    prefix: str,
    module_labels: dict,
) -> pd.DataFrame:

    assignments = pd.read_csv(
        assignment_file
    )


    assignments = normalize_gene_column(
        assignments
    )


    assignment_keep = [
        column
        for column
        in [
            "Gene",
            "Module",
            "Isolate",
            "Degree",
        ]
        if column
        in assignments.columns
    ]


    assignments = assignments[
        assignment_keep
    ].copy()


    assignments = assignments.rename(
        columns={
            "Module":
                f"{prefix}_module",

            "Isolate":
                f"{prefix}_isolate",

            "Degree":
                f"{prefix}_degree",
        }
    )


    assignments[
        f"{prefix}_module_label"
    ] = (
        assignments[
            f"{prefix}_module"
        ]
        .map(
            module_labels
        )
    )


    bridges = pd.read_csv(
        bridge_file
    )


    bridges = normalize_gene_column(
        bridges
    )


    bridge_map = {

        "External_degree":
            f"{prefix}_external_degree",

        "External_module_count":
            f"{prefix}_external_module_count",

        "Participation_coefficient":
            f"{prefix}_participation",

        "Within_module_degree_z":
            f"{prefix}_within_module_z",

        "Betweenness":
            f"{prefix}_betweenness",

        "PageRank":
            f"{prefix}_pagerank",
    }


    bridge_keep = [
        "Gene"
    ]


    rename = {}


    for source, target in (
        bridge_map.items()
    ):

        if source in bridges.columns:

            bridge_keep.append(
                source
            )

            rename[
                source
            ] = target


    bridges = (
        bridges[
            bridge_keep
        ]
        .rename(
            columns=rename
        )
    )


    return assignments.merge(
        bridges,
        on="Gene",
        how="left",
    )


core_f4 = load_f4_network(
    CORE_ASSIGN_FILE,
    CORE_BRIDGE_FILE,
    "F4_CORE",
    CORE_MODULE_LABELS,
)


ext_f4 = load_f4_network(
    EXT_ASSIGN_FILE,
    EXT_BRIDGE_FILE,
    "F4_EXT",
    EXT_MODULE_LABELS,
)


master = master.merge(
    core_f4,
    on="Gene",
    how="left",
)


master = master.merge(
    ext_f4,
    on="Gene",
    how="left",
)


interface_bridges = pd.read_csv(
    INTERFACE_BRIDGE_FILE
)


interface_bridges = normalize_gene_column(
    interface_bridges
)


interface_map = {

    "CORE_neighbor_count":
        "F4_interface_CORE_neighbor_count",

    "Distinct_CORE_modules_touched":
        "F4_interface_CORE_modules_touched_count",

    "CORE_modules_touched":
        "F4_interface_CORE_modules_touched",

    "Participation_coefficient":
        "F4_interface_participation",

    "Betweenness":
        "F4_interface_betweenness",
}


keep = [
    "Gene"
]


rename = {}


for source, target in (
    interface_map.items()
):

    if source in interface_bridges.columns:

        keep.append(
            source
        )

        rename[
            source
        ] = target


interface_bridges = (
    interface_bridges[
        keep
    ]
    .rename(
        columns=rename
    )
)


master = master.merge(
    interface_bridges,
    on="Gene",
    how="left",
)


# ============================================================
# 7. Corrected F5
# ============================================================

print()
print(
    "[F6] Integrating corrected F5 directed evidence..."
)


anchors = pd.read_csv(
    F5_ANCHOR_FILE
)


pairs = pd.read_csv(
    F5_PAIR_FILE
)


paths = pd.read_csv(
    F5_PATH_FILE
)


required_nodes = pd.read_csv(
    F5_REQUIRED_FILE
)


for frame in [
    anchors,
    pairs,
    paths,
    required_nodes,
]:

    for column in [
        "Gene",
        "Source",
        "Target",
        "Node",
        "Configured_gene",
        "Resolved_graph_gene",
    ]:

        if column in frame.columns:

            frame[
                column
            ] = (
                frame[
                    column
                ]
                .map(
                    normalize_gene
                )
            )


# ============================================================
# 7a. Directed route summary
# ============================================================

route_rows = []


for _, pair in pairs.iterrows():

    source = pair[
        "Source"
    ]


    target = pair[
        "Target"
    ]


    pair_paths = paths[
        (
            paths[
                "Source"
            ]
            ==
            source
        )
        &
        (
            paths[
                "Target"
            ]
            ==
            target
        )
    ].copy()


    direct = pair_paths[
        pair_paths[
            "Edges"
        ]
        ==
        1
    ].copy()


    if pair_paths.empty:

        best_external = None

        best_edges = None

        best_paths = ""

        best_effects = ""

        best_layers = ""

    else:

        pair_paths[
            "External_intermediate_count"
        ] = pd.to_numeric(
            pair_paths[
                "External_intermediate_count"
            ],
            errors="coerce",
        )


        pair_paths[
            "Edges"
        ] = pd.to_numeric(
            pair_paths[
                "Edges"
            ],
            errors="coerce",
        )


        best_external = (
            pair_paths[
                "External_intermediate_count"
            ]
            .min()
        )


        best_subset = pair_paths[
            pair_paths[
                "External_intermediate_count"
            ]
            ==
            best_external
        ].copy()


        best_edges = (
            best_subset[
                "Edges"
            ]
            .min()
        )


        best_subset = best_subset[
            best_subset[
                "Edges"
            ]
            ==
            best_edges
        ].copy()


        best_subset = (
            best_subset
            .drop_duplicates(
                "Path"
            )
            .head(
                3
            )
        )


        best_paths = join_unique(
            best_subset[
                "Path"
            ],
            separator=" || ",
        )


        best_effects = join_unique(
            best_subset[
                "Path_effect"
            ]
        )


        best_layers = join_unique(
            best_subset[
                "Path_layers"
            ]
        )


    required_for_pair = (
        required_nodes[
            (
                required_nodes[
                    "Source"
                ]
                ==
                source
            )
            &
            (
                required_nodes[
                    "Target"
                ]
                ==
                target
            )
            &
            (
                required_nodes[
                    "Exact_required"
                ]
                .map(
                    truthy
                )
            )
        ]
    )


    route_rows.append(
        {
            "Source_family":
                pair.get(
                    "Source_family"
                ),

            "Source":
                source,

            "Target_family":
                pair.get(
                    "Target_family"
                ),

            "Target":
                target,

            "Testable":
                pair.get(
                    "Testable"
                ),

            "Support_tier":
                pair.get(
                    "Support_tier"
                ),

            "Strictest_reachable_setting":
                pair.get(
                    "Strictest_reachable_setting"
                ),

            "Minimum_external_any_setting":
                pair.get(
                    "Minimum_external_any_setting"
                ),

            "Shortest_edges_at_min_external":
                pair.get(
                    "Shortest_edges_at_min_external"
                ),

            "Direct_edge":
                not direct.empty,

            "Direct_effect":
                (
                    join_unique(
                        direct[
                            "Path_effect"
                        ]
                    )
                    if not direct.empty
                    else ""
                ),

            "Direct_layers":
                (
                    join_unique(
                        direct[
                            "Path_layers"
                        ]
                    )
                    if not direct.empty
                    else ""
                ),

            "Best_representative_external":
                best_external,

            "Best_representative_edges":
                best_edges,

            "Best_representative_effect":
                best_effects,

            "Best_representative_layers":
                best_layers,

            "Best_representative_paths":
                best_paths,

            "Exact_required_intermediates":
                (
                    join_unique(
                        required_for_pair[
                            "Node"
                        ]
                    )
                    if not required_for_pair.empty
                    else ""
                ),
        }
    )


route_summary = pd.DataFrame(
    route_rows
)


route_summary.to_csv(
    OUTPUT_DIR
    / "directed_route_summary.csv",
    index=False,
)


# ============================================================
# 7b. Anchor coverage lookup
# ============================================================

anchor_lookup = {}


for _, row in anchors.iterrows():

    gene = row.get(
        "Configured_gene"
    )

    if gene is None:
        continue

    anchor_lookup[
        (
            row.get(
                "Role"
            ),
            gene,
        )
    ] = row


# ============================================================
# 7c. Gene-level F5 summaries
# ============================================================

f5_gene_rows = []


exact_required_only = required_nodes[
    required_nodes[
        "Exact_required"
    ]
    .map(
        truthy
    )
].copy()


for gene in experimental_genes:

    outgoing = route_summary[
        route_summary[
            "Source"
        ]
        ==
        gene
    ]


    incoming = route_summary[
        route_summary[
            "Target"
        ]
        ==
        gene
    ]


    tier_a_out = outgoing[
        outgoing[
            "Support_tier"
        ]
        ==
        "A"
    ]


    tier_b_out = outgoing[
        outgoing[
            "Support_tier"
        ]
        ==
        "B"
    ]


    tier_a_in = incoming[
        incoming[
            "Support_tier"
        ]
        ==
        "A"
    ]


    direct_out = outgoing[
        outgoing[
            "Direct_edge"
        ]
        .map(
            truthy
        )
    ]


    direct_in = incoming[
        incoming[
            "Direct_edge"
        ]
        .map(
            truthy
        )
    ]


    required_here = exact_required_only[
        exact_required_only[
            "Node"
        ]
        ==
        gene
    ]


    source_anchor = anchor_lookup.get(
        (
            "SOURCE",
            gene,
        )
    )


    target_anchor = anchor_lookup.get(
        (
            "TARGET",
            gene,
        )
    )


    direct_out_text = []


    for _, row in direct_out.iterrows():

        direct_out_text.append(
            f"{row['Target']}"
            f"[{row['Direct_effect']}; "
            f"{row['Direct_layers']}]"
        )


    direct_in_text = []


    for _, row in direct_in.iterrows():

        direct_in_text.append(
            f"{row['Source']}"
            f"[{row['Direct_effect']}; "
            f"{row['Direct_layers']}]"
        )


    required_pairs = []


    for _, row in required_here.iterrows():

        required_pairs.append(
            f"{row['Source']}"
            f"->{row['Target']}"
            f"@{row['Setting']}"
        )


    f5_gene_rows.append(
        {
            "Gene":
                gene,

            "F5_configured_source":
                source_anchor
                is not None,

            "F5_source_resolution":
                (
                    source_anchor.get(
                        "Resolution_status"
                    )
                    if source_anchor
                    is not None
                    else None
                ),

            "F5_configured_target":
                target_anchor
                is not None,

            "F5_target_resolution":
                (
                    target_anchor.get(
                        "Resolution_status"
                    )
                    if target_anchor
                    is not None
                    else None
                ),

            "F5_target_in_prior_graph":
                (
                    target_anchor.get(
                        "In_OmniPath_graph"
                    )
                    if target_anchor
                    is not None
                    else None
                ),

            "F5_TierA_outgoing_count":
                len(
                    tier_a_out
                ),

            "F5_TierA_targets":
                join_unique(
                    tier_a_out[
                        "Target"
                    ]
                ),

            "F5_TierB_outgoing_count":
                len(
                    tier_b_out
                ),

            "F5_TierB_targets":
                join_unique(
                    tier_b_out[
                        "Target"
                    ]
                ),

            "F5_TierA_incoming_count":
                len(
                    tier_a_in
                ),

            "F5_TierA_sources":
                join_unique(
                    tier_a_in[
                        "Source"
                    ]
                ),

            "F5_direct_outgoing":
                join_unique(
                    direct_out_text
                ),

            "F5_direct_incoming":
                join_unique(
                    direct_in_text
                ),

            "F5_exact_required_pair_count":
                len(
                    required_here
                ),

            "F5_exact_required_pairs":
                join_unique(
                    required_pairs
                ),
        }
    )


f5_gene = pd.DataFrame(
    f5_gene_rows
)


master = master.merge(
    f5_gene,
    on="Gene",
    how="left",
)


# ============================================================
# 8. Julian / E2F4 promoter evidence
# ============================================================

print()
print(
    "[F6] Searching for revised Julian promoter evidence..."
)


julian_audit_rows = []


julian_sig_candidates = (
    recursive_file_search(
        ".csv",
        [
            "e2f4",
            "e2f3",
            "significant",
            "genes",
        ],
    )
)


julian_sig_candidates = [
    path
    for path
    in julian_sig_candidates
    if "probe"
    not in path.name.lower()
]


julian_sig = pd.DataFrame()


if julian_sig_candidates:

    julian_sig_file = (
        julian_sig_candidates[
            0
        ]
    )


    try:

        julian_sig = pd.read_csv(
            julian_sig_file,
            low_memory=False,
        )


        julian_sig = (
            normalize_gene_column(
                julian_sig
            )
        )


        julian_sig = (
            julian_sig
            .drop_duplicates(
                "Gene"
            )
        )


        julian_audit_rows.append(
            {
                "Evidence":
                    "Revised E2F4-vs-E2F3 significant promoter genes",

                "Status":
                    "loaded",

                "Path":
                    str(
                        julian_sig_file
                    ),

                "Rows":
                    len(
                        julian_sig
                    ),
            }
        )

    except Exception as exc:

        julian_audit_rows.append(
            {
                "Evidence":
                    "Revised E2F4-vs-E2F3 significant promoter genes",

                "Status":
                    f"failed: {exc}",

                "Path":
                    str(
                        julian_sig_file
                    ),

                "Rows":
                    None,
            }
        )

        julian_sig = pd.DataFrame()


else:

    julian_audit_rows.append(
        {
            "Evidence":
                "Revised E2F4-vs-E2F3 significant promoter genes",

            "Status":
                "not found",

            "Path":
                None,

            "Rows":
                None,
        }
    )


if not julian_sig.empty:

    julian_annotation = pd.DataFrame(
        {
            "Gene":
                julian_sig[
                    "Gene"
                ],

            "Julian_E2F4_vs_E2F3_sig_promoter":
                True,
        }
    )


    metric_candidates = {

        "logFC":
            "Julian_binding_logFC",

        "log2FC":
            "Julian_binding_log2FC",

        "AveExpr":
            "Julian_binding_AveExpr",

        "P.Value":
            "Julian_binding_pvalue",

        "adj.P.Val":
            "Julian_binding_adj_pvalue",

        "M_avg":
            "Julian_binding_Mavg",

        "Mavg":
            "Julian_binding_Mavg",
    }


    for source, target in (
        metric_candidates.items()
    ):

        if source in julian_sig.columns:

            julian_annotation[
                target
            ] = (
                julian_sig[
                    source
                ].values
            )


    master = master.merge(
        julian_annotation,
        on="Gene",
        how="left",
    )


else:

    master[
        "Julian_E2F4_vs_E2F3_sig_promoter"
    ] = False


master[
    "Julian_E2F4_vs_E2F3_sig_promoter"
] = (
    master[
        "Julian_E2F4_vs_E2F3_sig_promoter"
    ]
    .map(
        truthy
    )
)


# ============================================================
# 8b. Optional binding-expression workbook
# ============================================================

binding_expression_candidates = (
    recursive_file_search(
        ".xlsx",
        [
            "binding",
            "expression",
        ],
    )
)


julian_binding_overlap = set()

julian_significant_both = set()


if binding_expression_candidates:

    binding_expression_file = (
        binding_expression_candidates[
            0
        ]
    )


    try:

        xl = pd.ExcelFile(
            binding_expression_file
        )


        for sheet in xl.sheet_names:

            sheet_df = pd.read_excel(
                binding_expression_file,
                sheet_name=sheet,
            )


            try:

                sheet_df = (
                    normalize_gene_column(
                        sheet_df
                    )
                )

            except Exception:

                continue


            genes = set(
                sheet_df[
                    "Gene"
                ]
            )


            julian_binding_overlap.update(
                genes
            )


            if (
                "significant"
                in sheet.lower()
                and
                "both"
                in sheet.lower()
            ):

                julian_significant_both.update(
                    genes
                )


        julian_audit_rows.append(
            {
                "Evidence":
                    "Julian binding-expression integration workbook",

                "Status":
                    "loaded",

                "Path":
                    str(
                        binding_expression_file
                    ),

                "Rows":
                    len(
                        julian_binding_overlap
                    ),
            }
        )

    except Exception as exc:

        julian_audit_rows.append(
            {
                "Evidence":
                    "Julian binding-expression integration workbook",

                "Status":
                    f"failed: {exc}",

                "Path":
                    str(
                        binding_expression_file
                    ),

                "Rows":
                    None,
            }
        )


else:

    julian_audit_rows.append(
        {
            "Evidence":
                "Julian binding-expression integration workbook",

            "Status":
                "not found",

            "Path":
                None,

            "Rows":
                None,
        }
    )


master[
    "Julian_binding_expression_overlap"
] = (
    master[
        "Gene"
    ]
    .isin(
        julian_binding_overlap
    )
)


master[
    "Julian_significant_binding_and_expression"
] = (
    master[
        "Gene"
    ]
    .isin(
        julian_significant_both
    )
)


julian_audit = pd.DataFrame(
    julian_audit_rows
)


julian_audit.to_csv(
    OUTPUT_DIR
    / "julian_input_audit.csv",
    index=False,
)


# ============================================================
# 9. Corrected integration classes
# ============================================================

def integration_class(
    row: pd.Series,
) -> str:

    old_present = truthy(
        row.get(
            "Old_master"
        )
    )


    f2_present = truthy(
        row.get(
            "Fong_F2"
        )
    )


    strict_fong = truthy(
        row.get(
            "Fong_strict_DE"
        )
    )


    julian = truthy(
        row.get(
            "Julian_E2F4_vs_E2F3_sig_promoter"
        )
    )


    if (
        old_present
        and
        f2_present
        and
        julian
    ):

        return (
            "Old + new + promoter convergence"
        )


    if (
        old_present
        and
        f2_present
    ):

        return (
            "Old + new convergence"
        )


    if f2_present:

        tier = str(
            row.get(
                "F2_tier",
                ""
            )
        ).upper()


        if tier == "CORE":

            return (
                "New downstream CORE"
            )


        return (
            "New stress-interface"
        )


    if old_present:

        if strict_fong:

            return (
                "Old strict-Fong gene outside F2 cell-fate/stress panel"
            )


        return (
            "Old candidate absent from current strict Fong/F2"
        )


    return (
        "Other"
    )


master[
    "Integration_class"
] = master.apply(
    integration_class,
    axis=1,
)


# ============================================================
# 10. Evidence-layer descriptions
# ============================================================

def evidence_layers(
    row: pd.Series,
) -> str:

    layers = []


    if truthy(
        row.get(
            "Old_master"
        )
    ):

        layers.append(
            "old curated master"
        )


    if truthy(
        row.get(
            "Fong_strict_DE"
        )
    ):

        layers.append(
            "current strict Fong DEG"
        )


    if truthy(
        row.get(
            "Old_adult_embryo_shared"
        )
    ):

        layers.append(
            "adult-embryo shared"
        )


    if truthy(
        row.get(
            "Old_developmental_reversal"
        )
    ):

        layers.append(
            "developmental reversal"
        )


    if truthy(
        row.get(
            "Old_Andrusiak_overlap"
        )
    ):

        layers.append(
            "Andrusiak overlap"
        )


    if truthy(
        row.get(
            "Fong_F2"
        )
    ):

        layers.append(
            "Fong F2"
        )


    if (
        safe_number(
            row.get(
                "F4_CORE_external_degree"
            )
        )
        or 0
    ) > 0:

        layers.append(
            "F4 CORE cross-module"
        )


    if (
        safe_number(
            row.get(
                "F4_EXT_external_degree"
            )
        )
        or 0
    ) > 0:

        layers.append(
            "F4 EXTENDED cross-module"
        )


    if (
        safe_number(
            row.get(
                "F4_interface_CORE_modules_touched_count"
            )
        )
        or 0
    ) >= 2:

        layers.append(
            "F4 multi-CORE interface"
        )


    if truthy(
        row.get(
            "Julian_E2F4_vs_E2F3_sig_promoter"
        )
    ):

        layers.append(
            "Julian E2F4-biased promoter"
        )


    if str(
        row.get(
            "F5_direct_outgoing",
            ""
        )
    ).strip():

        layers.append(
            "F5 direct outgoing"
        )


    if str(
        row.get(
            "F5_direct_incoming",
            ""
        )
    ).strip():

        layers.append(
            "F5 direct incoming"
        )


    if (
        safe_number(
            row.get(
                "F5_exact_required_pair_count"
            )
        )
        or 0
    ) > 0:

        layers.append(
            "F5 exact-required intermediate"
        )


    return "; ".join(
        layers
    )


master[
    "Evidence_layers"
] = master.apply(
    evidence_layers,
    axis=1,
)


# ============================================================
# 11. Explicit mechanistic candidate panels
#
# These are descriptive hypothesis groups.
# They are NOT rankings.
# ============================================================

PANELS = {

    "P1_replication_checkpoint_interface": {

        "title":
            "Replication-stress / checkpoint interface",

        "genes": [
            "PCNA",
            "CLSPN",
            "TOPBP1",
            "CHEK1",
            "CHEK2",
            "BRCA1",
            "H2AFX",
            "CDK1",
            "CDK2",
        ],

        "rationale":
            (
                "F1/F4 identify a strong replication/checkpoint program; "
                "F4 places several members at module interfaces and F5 "
                "connects this state to multiple fate outputs."
            ),
    },


    "P2_E2F_MYC_death_receptor": {

        "title":
            "E2F/MYC to death-receptor branch",

        "genes": [
            "E2F1",
            "E2F2",
            "MYC",
            "DDIT3",
            "TNF",
            "FAS",
            "TNFRSF10B",
        ],

        "rationale":
            (
                "F4 identifies E2F/cell-cycle and stress/death modules; "
                "corrected F5 recovers short Fong-internal routes and "
                "direct MYC-to-death-receptor transcriptional edges."
            ),
    },


    "P3_intrinsic_apoptosis": {

        "title":
            "Intrinsic apoptotic effectors retained outside OmniPath coverage",

        "genes": [
            "PMAIP1",
            "BBC3",
        ],

        "rationale":
            (
                "NOXA/PMAIP1 and PUMA/BBC3 are strongly altered in Fong "
                "and co-localize with the E2F/G1-S fate module, but remain "
                "absent as targets from the directed mouse OmniPath graph."
            ),
    },


    "P4_regulated_necrosis": {

        "title":
            "Regulated-necrosis / inflammatory-death branch",

        "genes": [
            "ZBP1",
            "PYCARD",
            "MLKL",
            "EIF2AK2",
            "TNF",
        ],

        "rationale":
            (
                "F2/F4 recover a regulated-death/innate-response program; "
                "F5 provides directed structural routes to MLKL and PYCARD "
                "although sign resolution is generally weaker."
            ),
    },


    "P5_trophic_survival": {

        "title":
            "Trophic / survival counter-program",

        "genes": [
            "BDNF",
            "KDR",
            "FGF2",
            "IRS2",
            "ANGPT1",
            "MET",
            "CXCL12",
        ],

        "rationale":
            (
                "F4 resolves a distinct RTK/MAPK/PI3K survival module and "
                "F5 identifies direct or short directed routes from the "
                "cell-cycle/checkpoint state into several survival targets."
            ),
    },


    "P6_old_mitotic_recurrence": {

        "title":
            "Old mitotic/checkpoint candidates preserved by the new analysis",

        "genes": [
            "AURKB",
            "BIRC5",
            "BUB1",
            "MELK",
            "KIF14",
            "MYBL2",
        ],

        "rationale":
            (
                "The original network highlighted an aberrant mitotic "
                "program; the independent Fong-wide analysis again resolves "
                "strong mitotic and chromosome-segregation structure."
            ),
    },


    "P7_developmental_context": {

        "title":
            "Developmental-stage context",

        "genes": [
            "FAS",
            "BDNF",
            "CXCL12",
            "NR4A3",
        ],

        "rationale":
            (
                "These genes preserve the old adult-versus-embryonic "
                "comparison and allow the new downstream mechanisms to be "
                "interpreted in developmental context."
            ),
    },
}


panel_rows = []


master_gene_set = set(
    master[
        "Gene"
    ]
)


for panel_id, panel in (
    PANELS.items()
):

    configured = [
        normalize_gene(
            gene
        )
        for gene
        in panel[
            "genes"
        ]
    ]


    present = [
        gene
        for gene
        in configured
        if gene
        in master_gene_set
    ]


    absent = [
        gene
        for gene
        in configured
        if gene
        not in master_gene_set
    ]


    panel_rows.append(
        {
            "Panel":
                panel_id,

            "Title":
                panel[
                    "title"
                ],

            "Genes_present":
                "; ".join(
                    present
                ),

            "Genes_absent":
                "; ".join(
                    absent
                ),

            "Rationale":
                panel[
                    "rationale"
                ],
        }
    )


candidate_panels = pd.DataFrame(
    panel_rows
)


candidate_panels.to_csv(
    OUTPUT_DIR
    / "candidate_panels.csv",
    index=False,
)


gene_panel_membership = defaultdict(
    list
)


for panel_id, panel in (
    PANELS.items()
):

    for gene in panel[
        "genes"
    ]:

        gene_panel_membership[
            normalize_gene(
                gene
            )
        ].append(
            panel_id
        )


master[
    "Candidate_panels"
] = (
    master[
        "Gene"
    ]
    .map(
        lambda gene:
            "; ".join(
                gene_panel_membership.get(
                    gene,
                    [],
                )
            )
    )
)


# ============================================================
# 12. Mechanistic hypotheses
# ============================================================

hypotheses = pd.DataFrame(
    [

        {
            "Hypothesis":
                "H1",

            "Working_model":
                (
                    "Pocket-protein loss produces an aberrant proliferative "
                    "state accompanied by replication stress, checkpoint "
                    "activation, and DNA-repair engagement."
                ),

            "Primary_genes":
                (
                    "PCNA; CLSPN; TOPBP1; CHEK1; CHEK2; BRCA1; "
                    "H2AFX; CDK1; CDK2"
                ),

            "Evidence":
                (
                    "F1 global enrichment + F4 EXT-M02/M03/M04 + "
                    "F4 bridge topology"
                ),

            "Interpretation_limit":
                (
                    "The computational data do not establish whether the "
                    "repair/checkpoint response is protective, damaging, "
                    "or both."
                ),
        },


        {
            "Hypothesis":
                "H2",

            "Working_model":
                (
                    "The aberrant cell-cycle/checkpoint state can feed into "
                    "death-receptor regulation through E2F1/MYC/DDIT3/TNF."
                ),

            "Primary_genes":
                (
                    "E2F1; E2F2; MYC; DDIT3; TNF; FAS; TNFRSF10B"
                ),

            "Evidence":
                (
                    "F4 cross-module topology + corrected F5 X0 routes + "
                    "direct MYC->FAS and MYC->TNFRSF10B prior-knowledge edges"
                ),

            "Interpretation_limit":
                (
                    "OmniPath/CollecTRI edges are prior knowledge and do not "
                    "prove that the same regulation occurs after Rb loss."
                ),
        },


        {
            "Hypothesis":
                "H3",

            "Working_model":
                (
                    "A parallel intrinsic-apoptotic program involving "
                    "NOXA/PMAIP1 and PUMA/BBC3 is engaged in the Fong "
                    "response."
                ),

            "Primary_genes":
                "PMAIP1; BBC3",

            "Evidence":
                (
                    "Strong Fong differential expression + F2 CORE membership "
                    "+ F4 E2F/G1-S fate-module co-localization"
                ),

            "Interpretation_limit":
                (
                    "PMAIP1 and BBC3 are absent from the directed OmniPath "
                    "mouse graph used here, so F5 cannot support or refute "
                    "their directed upstream regulation."
                ),
        },


        {
            "Hypothesis":
                "H4",

            "Working_model":
                (
                    "Regulated-necrosis / inflammatory-death machinery may "
                    "form an additional downstream branch."
                ),

            "Primary_genes":
                "ZBP1; PYCARD; MLKL; EIF2AK2; TNF",

            "Evidence":
                (
                    "F2/F4 regulated-death modules + F5 directed X0 routes "
                    "to MLKL/PYCARD"
                ),

            "Interpretation_limit":
                (
                    "Most directed routes have unresolved or mixed sign, "
                    "and ZBP1 could not be tested as an OmniPath endpoint."
                ),
        },


        {
            "Hypothesis":
                "H5",

            "Working_model":
                (
                    "Rb-family loss simultaneously engages trophic/survival "
                    "signaling rather than producing a purely pro-death state."
                ),

            "Primary_genes":
                "BDNF; KDR; FGF2; IRS2; ANGPT1; MET; CXCL12",

            "Evidence":
                (
                    "F2 anti-death annotations + F4 RTK/MAPK/PI3K module + "
                    "F5 direct/short survival-target routes"
                ),

            "Interpretation_limit":
                (
                    "Genes within this group are not uniformly up- or "
                    "down-regulated, so it should not be treated as one "
                    "unidirectional pathway."
                ),
        },


        {
            "Hypothesis":
                "H6",

            "Working_model":
                (
                    "Developmental stage modifies which downstream fate "
                    "branches dominate after pocket-protein loss."
                ),

            "Primary_genes":
                "FAS; BDNF; CXCL12; NR4A3",

            "Evidence":
                (
                    "Original adult/embryonic overlap and reversal evidence "
                    "+ preservation of several genes in the new Fong network"
                ),

            "Interpretation_limit":
                (
                    "The adult and embryonic datasets differ in biological "
                    "context and experimental design, so developmental "
                    "causality requires direct validation."
                ),
        },
    ]
)


hypotheses.to_csv(
    OUTPUT_DIR
    / "mechanistic_hypotheses.csv",
    index=False,
)


# ============================================================
# 13. Mechanistic model edges
# ============================================================

model_rows = [

    {
        "Source":
            "Pocket-protein loss",

        "Target":
            "E2F / cell-cycle dysregulation",

        "Sign":
            "context-dependent",

        "Evidence_type":
            "project premise / old + new convergence",
    },


    {
        "Source":
            "E2F / cell-cycle dysregulation",

        "Target":
            "Mitotic / replication-stress state",

        "Sign":
            "activation-associated",

        "Evidence_type":
            "F1 enrichment + F4 module structure",
    },


    {
        "Source":
            "Mitotic / replication-stress state",

        "Target":
            "Checkpoint / DNA-repair response",

        "Sign":
            "activation-associated",

        "Evidence_type":
            "F1 + F4 EXT-M02/M04",
    },


    {
        "Source":
            "E2F/MYC/DDIT3/TNF interface",

        "Target":
            "Death-receptor branch",

        "Sign":
            "mixed / route-specific",

        "Evidence_type":
            "F4 + corrected F5",
    },


    {
        "Source":
            "E2F/G1-S fate module",

        "Target":
            "PMAIP1 / BBC3 intrinsic apoptosis",

        "Sign":
            "association only",

        "Evidence_type":
            "Fong expression + F4 co-module evidence",
    },


    {
        "Source":
            "Stress / inflammatory signaling",

        "Target":
            "MLKL / PYCARD regulated-death branch",

        "Sign":
            "mostly unresolved",

        "Evidence_type":
            "F4 + corrected F5",
    },


    {
        "Source":
            "Trophic / RTK signaling",

        "Target":
            "Cell-death susceptibility",

        "Sign":
            "context-dependent",

        "Evidence_type":
            "F2 + F4 + F5",
    },
]


direct_routes = route_summary[
    route_summary[
        "Direct_edge"
    ]
    .map(
        truthy
    )
]


for _, row in direct_routes.iterrows():

    model_rows.append(
        {
            "Source":
                row[
                    "Source"
                ],

            "Target":
                row[
                    "Target"
                ],

            "Sign":
                row[
                    "Direct_effect"
                ],

            "Evidence_type":
                (
                    "Corrected F5 direct prior-knowledge edge: "
                    f"{row['Direct_layers']}"
                ),
        }
    )


model_edges = pd.DataFrame(
    model_rows
)


model_edges.to_csv(
    OUTPUT_DIR
    / "mechanistic_model_edges.csv",
    index=False,
)


# ============================================================
# 14. External prior-knowledge nodes
# ============================================================

external_counter = Counter()

external_pairs = defaultdict(
    set
)

external_settings = defaultdict(
    set
)


if not paths.empty:

    for _, row in paths.iterrows():

        value = row.get(
            "External_intermediates"
        )


        if pd.isna(value):
            continue


        pair_id = (
            f"{row['Source']}"
            f"->{row['Target']}"
        )


        for node in str(
            value
        ).split(
            ";"
        ):

            node = normalize_gene(
                node
            )


            if node is None:
                continue


            external_counter[
                node
            ] += 1


            external_pairs[
                node
            ].add(
                pair_id
            )


            external_settings[
                node
            ].add(
                row[
                    "Setting"
                ]
            )


old_external = optional_csv(
    OLD_FINAL_EXTERNAL_FILE
)


if old_external.empty:

    old_external = optional_csv(
        OLD_EXTERNAL_EVIDENCE_FILE
    )


if not old_external.empty:

    old_external = (
        normalize_gene_column(
            old_external
        )
    )


old_external_genes = (
    set(
        old_external[
            "Gene"
        ]
    )
    if not old_external.empty
    else set()
)


all_external_genes = sorted(
    set(
        external_counter
    )
    |
    old_external_genes
)


external_rows = []


julian_sig_gene_set = (
    set(
        julian_sig[
            "Gene"
        ]
    )
    if not julian_sig.empty
    else set()
)


for gene in all_external_genes:

    external_rows.append(
        {
            "Gene":
                gene,

            "Old_external_connector":
                gene
                in old_external_genes,

            "New_F5_representative_external":
                gene
                in external_counter,

            "New_F5_path_occurrences":
                external_counter.get(
                    gene,
                    0,
                ),

            "New_F5_distinct_pairs":
                len(
                    external_pairs.get(
                        gene,
                        set(),
                    )
                ),

            "New_F5_settings":
                "; ".join(
                    sorted(
                        external_settings.get(
                            gene,
                            set(),
                        )
                    )
                ),

            "Observed_in_current_strict_Fong":
                gene
                in strict_fong_genes,

            "Observed_in_F2":
                gene
                in f2_genes,

            "Observed_in_old_master":
                gene
                in old_genes,

            "Julian_sig_promoter":
                gene
                in julian_sig_gene_set,
        }
    )


external_prior = pd.DataFrame(
    external_rows
)


if not external_prior.empty:

    external_prior = (
        external_prior
        .sort_values(
            by=[
                "Old_external_connector",
                "New_F5_distinct_pairs",
                "New_F5_path_occurrences",
            ],
            ascending=[
                False,
                False,
                False,
            ],
        )
        .reset_index(
            drop=True
        )
    )


external_prior.to_csv(
    OUTPUT_DIR
    / "external_prior_knowledge_nodes.csv",
    index=False,
)


# ============================================================
# 15. Integration summary
# ============================================================

old_new_overlap = (
    old_genes
    &
    f2_genes
)


new_only = (
    f2_genes
    -
    old_genes
)


old_only = (
    old_genes
    -
    f2_genes
)


old_strict_fong = (
    old_genes
    &
    strict_fong_genes
)


old_absent_strict_fong = (
    old_genes
    -
    strict_fong_genes
)


integration_class_counts = (
    master[
        "Integration_class"
    ]
    .value_counts()
)


andrusiak_series = (
    master[
        "Old_Andrusiak_overlap"
    ]
    if "Old_Andrusiak_overlap"
    in master.columns
    else pd.Series(
        False,
        index=master.index,
    )
)


andrusiak_count = int(
    andrusiak_series
    .map(
        truthy
    )
    .sum()
)


interface_touch_series = (
    master[
        "F4_interface_CORE_modules_touched_count"
    ]
    if (
        "F4_interface_CORE_modules_touched_count"
        in master.columns
    )
    else pd.Series(
        0,
        index=master.index,
    )
)


multi_core_interface = int(
    (
        pd.to_numeric(
            interface_touch_series,
            errors="coerce",
        )
        .fillna(
            0
        )
        >=
        2
    )
    .sum()
)


summary_rows = [

    {
        "Metric":
            "Old valid experimental master genes",

        "Value":
            len(
                old_genes
            ),
    },


    {
        "Metric":
            "Current strict Fong genes",

        "Value":
            len(
                strict_fong_genes
            ),
    },


    {
        "Metric":
            "Old master genes in current strict Fong",

        "Value":
            len(
                old_strict_fong
            ),
    },


    {
        "Metric":
            "Old master genes absent from current strict Fong",

        "Value":
            len(
                old_absent_strict_fong
            ),
    },


    {
        "Metric":
            "New F2 genes",

        "Value":
            len(
                f2_genes
            ),
    },


    {
        "Metric":
            "Old-new overlap",

        "Value":
            len(
                old_new_overlap
            ),
    },


    {
        "Metric":
            "New F2 genes absent from old master",

        "Value":
            len(
                new_only
            ),
    },


    {
        "Metric":
            "Old genes absent from F2",

        "Value":
            len(
                old_only
            ),
    },


    {
        "Metric":
            "Integrated experimental universe",

        "Value":
            len(
                experimental_genes
            ),
    },


    {
        "Metric":
            "Julian revised promoter hits in integrated universe",

        "Value":
            int(
                master[
                    "Julian_E2F4_vs_E2F3_sig_promoter"
                ]
                .sum()
            ),
    },


    {
        "Metric":
            "Old Andrusiak overlaps in integrated universe",

        "Value":
            andrusiak_count,
    },


    {
        "Metric":
            "F4 interface genes touching >=2 CORE modules",

        "Value":
            multi_core_interface,
    },


    {
        "Metric":
            "F5 configured source-target pairs",

        "Value":
            len(
                route_summary
            ),
    },


    {
        "Metric":
            "F5 Tier A pairs",

        "Value":
            int(
                (
                    route_summary[
                        "Support_tier"
                    ]
                    ==
                    "A"
                )
                .sum()
            ),
    },


    {
        "Metric":
            "F5 Tier B pairs",

        "Value":
            int(
                (
                    route_summary[
                        "Support_tier"
                    ]
                    ==
                    "B"
                )
                .sum()
            ),
    },


    {
        "Metric":
            "F5 unreachable pairs",

        "Value":
            int(
                (
                    route_summary[
                        "Support_tier"
                    ]
                    ==
                    "UNREACHABLE"
                )
                .sum()
            ),
    },


    {
        "Metric":
            "F5 untestable pairs",

        "Value":
            int(
                (
                    route_summary[
                        "Support_tier"
                    ]
                    ==
                    "UNTESTABLE"
                )
                .sum()
            ),
    },


    {
        "Metric":
            "F5 direct source-target edges",

        "Value":
            int(
                route_summary[
                    "Direct_edge"
                ]
                .map(
                    truthy
                )
                .sum()
            ),
    },


    {
        "Metric":
            "F5 exact-required intermediate events",

        "Value":
            len(
                exact_required_only
            ),
    },
]


for class_name, count in (
    integration_class_counts.items()
):

    summary_rows.append(
        {
            "Metric":
                (
                    "Integration class: "
                    f"{class_name}"
                ),

            "Value":
                int(
                    count
                ),
        }
    )


integration_summary = pd.DataFrame(
    summary_rows
)


integration_summary.to_csv(
    OUTPUT_DIR
    / "integration_summary.csv",
    index=False,
)


# ============================================================
# 16. Final ordering
# ============================================================

class_order = {

    "Old + new + promoter convergence":
        0,

    "Old + new convergence":
        1,

    "New downstream CORE":
        2,

    "New stress-interface":
        3,

    "Old strict-Fong gene outside F2 cell-fate/stress panel":
        4,

    "Old candidate absent from current strict Fong/F2":
        5,

    "Other":
        6,
}


master[
    "_class_order"
] = (
    master[
        "Integration_class"
    ]
    .map(
        class_order
    )
    .fillna(
        99
    )
)


fong_lfc_series = (
    master[
        "Fong_log2FC"
    ]
    if "Fong_log2FC"
    in master.columns
    else pd.Series(
        float("nan"),
        index=master.index,
    )
)


master[
    "_abs_lfc"
] = (
    pd.to_numeric(
        fong_lfc_series,
        errors="coerce",
    )
    .abs()
)


master = (
    master
    .sort_values(
        by=[
            "_class_order",
            "_abs_lfc",
            "Gene",
        ],
        ascending=[
            True,
            False,
            True,
        ],
        na_position="last",
    )
    .drop(
        columns=[
            "_class_order",
            "_abs_lfc",
        ]
    )
    .reset_index(
        drop=True
    )
)


master.to_csv(
    OUTPUT_DIR
    / "experimental_gene_evidence_matrix.csv",
    index=False,
)


# ============================================================
# 17. Excel formatting
# ============================================================

def format_workbook(
    path: Path,
) -> None:

    wb = load_workbook(
        path
    )


    header_fill = PatternFill(
        fill_type="solid",
        fgColor="D9E1F2",
    )


    for ws in wb.worksheets:

        ws.freeze_panes = "A2"


        for cell in ws[1]:

            cell.fill = (
                header_fill
            )


            cell.font = Font(
                bold=True
            )


            cell.alignment = Alignment(
                horizontal="center",
                vertical="center",
            )


        if ws.max_row > 1:

            ws.auto_filter.ref = (
                ws.dimensions
            )


        for column_cells in ws.columns:

            letter = get_column_letter(
                column_cells[
                    0
                ].column
            )


            max_length = 0


            for cell in column_cells[
                :100
            ]:

                if cell.value is None:
                    continue


                max_length = max(
                    max_length,
                    len(
                        str(
                            cell.value
                        )
                    ),
                )


            ws.column_dimensions[
                letter
            ].width = min(
                max(
                    max_length + 2,
                    12,
                ),
                55,
            )


    wb.save(
        path
    )


with pd.ExcelWriter(
    OUTPUT_XLSX,
    engine="openpyxl",
) as writer:

    integration_summary.to_excel(
        writer,
        sheet_name="Summary",
        index=False,
    )


    master.to_excel(
        writer,
        sheet_name="Experimental evidence",
        index=False,
    )


    route_summary.to_excel(
        writer,
        sheet_name="Directed routes",
        index=False,
    )


    candidate_panels.to_excel(
        writer,
        sheet_name="Candidate panels",
        index=False,
    )


    hypotheses.to_excel(
        writer,
        sheet_name="Hypotheses",
        index=False,
    )


    model_edges.to_excel(
        writer,
        sheet_name="Model edges",
        index=False,
    )


    external_prior.to_excel(
        writer,
        sheet_name="External prior nodes",
        index=False,
    )


    julian_audit.to_excel(
        writer,
        sheet_name="Julian audit",
        index=False,
    )


format_workbook(
    OUTPUT_XLSX
)


# ============================================================
# 18. Metadata
# ============================================================

metadata = {

    "phase":
        "F6_integrated_old_new_synthesis_corrected",

    "design":
        (
            "Descriptive multi-layer evidence integration; "
            "no composite biological score."
        ),

    "strict_Fong_definition":
        {
            "source_file":
                str(
                    STRICT_FONG_FILE
                ),

            "sheet":
                "Sig L2FC1",

            "valid_gene_count":
                len(
                    strict_fong_genes
                ),

            "note":
                (
                    "Current audited strict Fong membership is used for "
                    "the old-only integration classification rather than "
                    "the historical Fong flag from the old analysis-ready "
                    "workbook."
                ),
        },

    "experimental_universe":
        {
            "old_master_genes":
                len(
                    old_genes
                ),

            "old_master_in_current_strict_Fong":
                len(
                    old_strict_fong
                ),

            "old_master_absent_current_strict_Fong":
                len(
                    old_absent_strict_fong
                ),

            "new_F2_genes":
                len(
                    f2_genes
                ),

            "old_new_overlap":
                len(
                    old_new_overlap
                ),

            "new_only":
                len(
                    new_only
                ),

            "old_only":
                len(
                    old_only
                ),

            "union":
                len(
                    experimental_genes
                ),
        },

    "primary_integration_classes":
        list(
            class_order.keys()
        ),

    "F4_module_labels":
        {
            "CORE":
                CORE_MODULE_LABELS,

            "EXTENDED":
                EXT_MODULE_LABELS,
        },

    "candidate_panels":
        PANELS,

    "Julian_sources":
        julian_audit.to_dict(
            orient="records"
        ),

    "important_interpretation_notes":
        [

            (
                "Julian E2F4-vs-E2F3 promoter evidence is an annotation "
                "of possible proximal/direct regulation and is not used "
                "as an entry filter for the Fong downstream analysis."
            ),

            (
                "Current strict-Fong membership is derived directly from "
                "the audited Fong L2FC1 workbook and supersedes the "
                "historical Fong-membership flag for final classification."
            ),

            (
                "F5 direct and multi-step routes are prior-knowledge "
                "mechanistic hypotheses, not causal validation in the "
                "Fong experiment."
            ),

            (
                "PMAIP1/BBC3 directed OmniPath absence is treated as "
                "network-coverage limitation rather than negative "
                "biological evidence."
            ),

            (
                "STRING module centrality and OmniPath reachability are "
                "kept as separate evidence dimensions."
            ),

            (
                "No overall gene importance score is calculated."
            ),
        ],
}


OUTPUT_METADATA.write_text(
    json.dumps(
        metadata,
        indent=2,
        ensure_ascii=False,
        default=str,
    ),
    encoding="utf-8",
)


# ============================================================
# 19. Console report
# ============================================================

print()
print("=" * 82)
print("PHASE F6 COMPLETE")
print("=" * 82)


print()
print(
    "Experimental integration"
)


print(
    f"    Old valid master genes              : "
    f"{len(old_genes):,}"
)


print(
    f"    Current strict Fong genes           : "
    f"{len(strict_fong_genes):,}"
)


print(
    f"    Old master genes in strict Fong     : "
    f"{len(old_strict_fong):,}"
)


print(
    f"    Old master absent strict Fong       : "
    f"{len(old_absent_strict_fong):,}"
)


print(
    f"    New F2 genes                        : "
    f"{len(f2_genes):,}"
)


print(
    f"    Old-new overlap                     : "
    f"{len(old_new_overlap):,}"
)


print(
    f"    New F2 genes absent from old master : "
    f"{len(new_only):,}"
)


print(
    f"    Old genes absent from new F2        : "
    f"{len(old_only):,}"
)


print(
    f"    Integrated experimental universe    : "
    f"{len(experimental_genes):,}"
)


print()
print(
    "Integration classes"
)


for class_name in class_order:

    count = int(
        (
            master[
                "Integration_class"
            ]
            ==
            class_name
        )
        .sum()
    )


    if count == 0:
        continue


    print(
        f"    {class_name:<58} "
        f"{count:>4}"
    )


print()
print(
    "Cross-layer annotations"
)


print(
    f"    Julian revised promoter hits        : "
    f"{int(master['Julian_E2F4_vs_E2F3_sig_promoter'].sum()):,}"
)


print(
    f"    Andrusiak overlap                   : "
    f"{andrusiak_count:,}"
)


print(
    f"    F4 INTERFACE genes touching >=2 CORE modules: "
    f"{multi_core_interface:,}"
)


print()
print(
    "Corrected F5"
)


for tier in [
    "A",
    "B",
    "C",
    "D",
    "UNREACHABLE",
    "UNTESTABLE",
]:

    count = int(
        (
            route_summary[
                "Support_tier"
            ]
            ==
            tier
        )
        .sum()
    )


    print(
        f"    {tier:<11}: "
        f"{count:>3}"
    )


print(
    f"    Direct source-target edges          : "
    f"{int(route_summary['Direct_edge'].map(truthy).sum()):,}"
)


print(
    f"    Exact required-node events          : "
    f"{len(exact_required_only):,}"
)


print()
print(
    "Mechanistic panels"
)


for _, row in candidate_panels.iterrows():

    print()

    print(
        f"    {row['Panel']} — "
        f"{row['Title']}"
    )


    print(
        f"        {row['Genes_present']}"
    )


print()
print(
    "Julian input audit"
)


for _, row in julian_audit.iterrows():

    print(
        f"    {row['Evidence']}: "
        f"{row['Status']}"
    )


    if pd.notna(
        row[
            "Path"
        ]
    ):

        print(
            f"        {row['Path']}"
        )


print()
print(
    "Sanity checks"
)


expected_old_strict = (
    len(
        old_strict_fong
    )
)


expected_old_strict_outside_f2 = len(
    (
        old_genes
        &
        strict_fong_genes
    )
    -
    f2_genes
)


expected_old_absent_strict = len(
    old_genes
    -
    strict_fong_genes
)


print(
    f"    Old ∩ strict Fong                  : "
    f"{expected_old_strict:,}"
)


print(
    f"    Old ∩ strict Fong but outside F2   : "
    f"{expected_old_strict_outside_f2:,}"
)


print(
    f"    Old absent from strict Fong        : "
    f"{expected_old_absent_strict:,}"
)


print(
    f"    Partition check                    : "
    f"{len(old_new_overlap)} + "
    f"{expected_old_strict_outside_f2} + "
    f"{expected_old_absent_strict} "
    f"= "
    f"{len(old_new_overlap) + expected_old_strict_outside_f2 + expected_old_absent_strict}"
)


print()
print(
    "Output directory:"
)


print(
    f"    "
    f"{OUTPUT_DIR.relative_to(PROJECT_ROOT)}"
)


print()
print(
    "[NEXT] Interpret the corrected integrated evidence matrix "
    "and prioritize experimental follow-up hypotheses."
)