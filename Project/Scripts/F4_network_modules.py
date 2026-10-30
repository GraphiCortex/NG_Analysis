#!/usr/bin/env python3

"""
Phase F4 — STRING module and bridge analysis
=============================================

Purpose
-------
Analyze the internal organization of the F3 CORE and EXTENDED STRING
networks.

For each network this script:

1. Removes isolates for community detection.
2. Detects weighted Louvain communities using STRING combined scores.
3. Assigns stable module IDs M01, M02, ... by decreasing module size.
4. Calculates module-level structure.
5. Calculates node-level bridge metrics:
       - degree
       - internal degree
       - external/cross-module degree
       - number of neighboring modules
       - participation coefficient
       - within-module degree z-score
       - betweenness centrality
       - PageRank
6. Performs GO:BP + Reactome enrichment for each module.

Important enrichment design
---------------------------
Module enrichment uses the NON-ISOLATED proteins in the same STRING
network as a CUSTOM statistical background.

Therefore the enrichment question is:

    "What biology distinguishes this module from the rest of the
     connected Fong cell-fate network?"

rather than:

    "Is this preselected cell-fate gene set enriched for cell-fate terms?"

For the EXTENDED network, the script also asks:

7. Which INTERFACE proteins directly connect CORE genes?
8. Which INTERFACE proteins touch genes belonging to multiple
   standalone CORE modules?
9. How do standalone CORE modules map into EXTENDED modules?

Input
-----
Data/Processed/FongDownstream/F3_STRING/

    CORE/
        nodes.csv
        edges.csv

    EXTENDED/
        nodes.csv
        edges.csv

Outputs
-------
Data/Processed/FongDownstream/F4_Modules/

    F4_module_analysis.xlsx
    F4_metadata.json

    CORE/
        module_assignments.csv
        module_summary.csv
        bridge_metrics.csv
        cross_module_bridges.csv
        module_enrichment.csv
        module_enrichment_top.csv

    EXTENDED/
        module_assignments.csv
        module_summary.csv
        bridge_metrics.csv
        cross_module_bridges.csv
        interface_core_bridges.csv
        module_enrichment.csv
        module_enrichment_top.csv

    core_extended_correspondence.csv

Run
---
From Project/:

    py .\\Scripts\\F4_network_modules.py
"""

from __future__ import annotations

from pathlib import Path
from collections import Counter, defaultdict
from typing import Any
import json
import math
import time

import pandas as pd
import requests
import networkx as nx

from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

F3_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F3_STRING"
)

CORE_INPUT_DIR = (
    F3_DIR
    / "CORE"
)

EXTENDED_INPUT_DIR = (
    F3_DIR
    / "EXTENDED"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F4_Modules"
)

CORE_OUTPUT_DIR = (
    OUTPUT_DIR
    / "CORE"
)

EXTENDED_OUTPUT_DIR = (
    OUTPUT_DIR
    / "EXTENDED"
)

OUTPUT_XLSX = (
    OUTPUT_DIR
    / "F4_module_analysis.xlsx"
)

OUTPUT_METADATA = (
    OUTPUT_DIR
    / "F4_metadata.json"
)

CORRESPONDENCE_CSV = (
    OUTPUT_DIR
    / "core_extended_correspondence.csv"
)


# ============================================================
# Community-detection configuration
# ============================================================

LOUVAIN_SEED = 42

LOUVAIN_RESOLUTION = 1.0

EDGE_WEIGHT_COLUMN = "combined_score"


# ============================================================
# g:Profiler configuration
# ============================================================

GPROFILER_URL = (
    "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
)

ORGANISM = "mmusculus"

ENRICHMENT_SOURCES = [
    "GO:BP",
    "REAC",
]

FDR_THRESHOLD = 0.05

REQUEST_TIMEOUT = 180

API_SLEEP_SECONDS = 0.5

MIN_MODULE_SIZE_FOR_ENRICHMENT = 3

TOP_TERMS_PER_SOURCE = 5


# ============================================================
# Helpers
# ============================================================

def normalize_gene(
    value: Any,
) -> str | None:

    if pd.isna(value):
        return None

    gene = str(value).strip()

    if not gene:
        return None

    return gene.upper()


def bool_from_value(
    value: Any,
) -> bool:

    if isinstance(value, bool):
        return value

    if pd.isna(value):
        return False

    return (
        str(value)
        .strip()
        .upper()
        in {
            "TRUE",
            "YES",
            "1",
        }
    )


def parse_categories(
    value: Any,
) -> list[str]:

    if pd.isna(value):
        return []

    value = str(value).strip()

    if not value:
        return []

    return [
        x.strip()
        for x in value.split(";")
        if x.strip()
    ]


# ============================================================
# Load F3 network
# ============================================================

def load_network(
    input_dir: Path,
) -> tuple[
    nx.Graph,
    pd.DataFrame,
    pd.DataFrame,
]:

    nodes_file = (
        input_dir
        / "nodes.csv"
    )

    edges_file = (
        input_dir
        / "edges.csv"
    )

    if not nodes_file.exists():
        raise FileNotFoundError(
            f"Missing nodes file:\n{nodes_file}"
        )

    if not edges_file.exists():
        raise FileNotFoundError(
            f"Missing edges file:\n{edges_file}"
        )

    nodes = pd.read_csv(
        nodes_file
    )

    edges = pd.read_csv(
        edges_file
    )

    nodes[
        "Gene"
    ] = (
        nodes[
            "Gene"
        ]
        .map(
            normalize_gene
        )
    )

    G = nx.Graph()

    # Add all mapped proteins, including isolates.
    for _, row in nodes.iterrows():

        gene = row[
            "Gene"
        ]

        if gene is None:
            continue

        attrs = {}

        for column in nodes.columns:

            if column == "Gene":
                continue

            value = row[
                column
            ]

            if pd.isna(value):
                continue

            attrs[
                column
            ] = value

        G.add_node(
            gene,
            **attrs,
        )

    # Add STRING edges.
    for _, row in edges.iterrows():

        gene_a = normalize_gene(
            row[
                "Gene_A"
            ]
        )

        gene_b = normalize_gene(
            row[
                "Gene_B"
            ]
        )

        if (
            gene_a is None
            or
            gene_b is None
            or
            gene_a == gene_b
        ):
            continue

        score = float(
            row.get(
                EDGE_WEIGHT_COLUMN,
                1.0,
            )
        )

        G.add_edge(
            gene_a,
            gene_b,
            score=score,
        )

    return (
        G,
        nodes,
        edges,
    )


# ============================================================
# Community detection
# ============================================================

def detect_modules(
    G: nx.Graph,
) -> tuple[
    nx.Graph,
    list[set[str]],
    dict[str, str],
    float,
]:

    """
    Detect weighted Louvain communities among non-isolated nodes.

    Isolates are intentionally excluded from modules.
    """

    non_isolates = [
        node
        for node in G.nodes()
        if G.degree(
            node
        ) > 0
    ]

    H = G.subgraph(
        non_isolates
    ).copy()

    if H.number_of_nodes() == 0:

        return (
            H,
            [],
            {},
            float("nan"),
        )

    try:

        communities = (
            nx.algorithms.community
            .louvain_communities(
                H,
                weight="score",
                resolution=LOUVAIN_RESOLUTION,
                seed=LOUVAIN_SEED,
            )
        )

    except AttributeError as exc:

        raise RuntimeError(
            "\nYour installed NetworkX version does not provide "
            "louvain_communities().\n"
            "Upgrade it with:\n\n"
            "    py -m pip install --upgrade networkx\n"
        ) from exc

    # Stable ordering:
    # largest community first, then alphabetical tie-break.
    communities = sorted(
        [
            set(
                community
            )
            for community
            in communities
        ],
        key=lambda c: (
            -len(c),
            sorted(c)[0],
        ),
    )

    module_lookup = {}

    for idx, community in enumerate(
        communities,
        start=1,
    ):

        module_id = (
            f"M{idx:02d}"
        )

        for gene in community:

            module_lookup[
                gene
            ] = module_id

    modularity = (
        nx.algorithms.community
        .quality
        .modularity(
            H,
            communities,
            weight="score",
            resolution=LOUVAIN_RESOLUTION,
        )
    )

    return (
        H,
        communities,
        module_lookup,
        modularity,
    )


# ============================================================
# Bridge metrics
# ============================================================

def calculate_bridge_metrics(
    H: nx.Graph,
    module_lookup: dict[str, str],
) -> pd.DataFrame:

    """
    Calculate node-level module participation statistics.

    Participation coefficient:

        P_i = 1 - sum_s (k_is / k_i)^2

    where k_is is the number of edges from node i to module s.

    P ~ 0:
        nearly all neighbors are in one module

    larger P:
        edges are distributed across several modules
    """

    if H.number_of_nodes() == 0:
        return pd.DataFrame()

    # --------------------------------------------------------
    # Betweenness
    #
    # We use unweighted shortest paths here because STRING score
    # represents evidence strength, not physical distance.
    # --------------------------------------------------------

    betweenness = (
        nx.betweenness_centrality(
            H,
            normalized=True,
            weight=None,
        )
    )

    pagerank = nx.pagerank(
        H,
        weight="score",
    )

    node_rows = []

    internal_degree_lookup = {}

    # First pass:
    # calculate internal degree so module distributions exist.
    for node in H.nodes():

        own_module = (
            module_lookup[
                node
            ]
        )

        internal_degree = sum(
            1
            for neighbor
            in H.neighbors(
                node
            )
            if module_lookup[
                neighbor
            ]
            == own_module
        )

        internal_degree_lookup[
            node
        ] = internal_degree

    # Module-wise internal degree statistics.
    module_internal_values = defaultdict(
        list
    )

    for node, value in (
        internal_degree_lookup.items()
    ):

        module_internal_values[
            module_lookup[
                node
            ]
        ].append(
            value
        )

    module_internal_stats = {}

    for module, values in (
        module_internal_values.items()
    ):

        n = len(
            values
        )

        mean_value = (
            sum(values)
            /
            n
        )

        variance = (
            sum(
                (
                    x
                    -
                    mean_value
                )
                ** 2
                for x in values
            )
            /
            n
        )

        sd = math.sqrt(
            variance
        )

        module_internal_stats[
            module
        ] = (
            mean_value,
            sd,
        )

    # --------------------------------------------------------
    # Node metrics
    # --------------------------------------------------------

    for node in H.nodes():

        degree = H.degree(
            node
        )

        own_module = (
            module_lookup[
                node
            ]
        )

        neighbor_module_counts = Counter(
            module_lookup[
                neighbor
            ]
            for neighbor
            in H.neighbors(
                node
            )
        )

        internal_degree = (
            neighbor_module_counts.get(
                own_module,
                0,
            )
        )

        external_degree = (
            degree
            -
            internal_degree
        )

        external_modules = {
            module
            for module
            in neighbor_module_counts
            if module
            != own_module
        }

        if degree > 0:

            participation = (
                1.0
                -
                sum(
                    (
                        count
                        /
                        degree
                    )
                    ** 2
                    for count
                    in neighbor_module_counts.values()
                )
            )

        else:

            participation = 0.0

        mean_internal, sd_internal = (
            module_internal_stats[
                own_module
            ]
        )

        if sd_internal > 0:

            within_module_z = (
                internal_degree
                -
                mean_internal
            ) / sd_internal

        else:

            within_module_z = 0.0

        attrs = H.nodes[
            node
        ]

        node_rows.append(
            {
                "Gene":
                    node,

                "Module":
                    own_module,

                "Degree":
                    degree,

                "Internal_degree":
                    internal_degree,

                "External_degree":
                    external_degree,

                "External_degree_fraction":
                    (
                        external_degree
                        /
                        degree
                        if degree
                        else 0.0
                    ),

                "Neighbor_module_count":
                    len(
                        neighbor_module_counts
                    ),

                "External_module_count":
                    len(
                        external_modules
                    ),

                "Participation_coefficient":
                    participation,

                "Within_module_degree_z":
                    within_module_z,

                "Betweenness":
                    betweenness.get(
                        node,
                        0.0,
                    ),

                "PageRank":
                    pagerank.get(
                        node,
                        0.0,
                    ),

                "Direction":
                    attrs.get(
                        "Direction",
                        "",
                    ),

                "log2FoldChange":
                    attrs.get(
                        "log2FoldChange",
                    ),

                "F2_tier":
                    attrs.get(
                        "F2_tier",
                        "",
                    ),

                "F2_categories":
                    attrs.get(
                        "F2_categories",
                        "",
                    ),

                "Old_126_gene_master":
                    attrs.get(
                        "Old_126_gene_master",
                        False,
                    ),
            }
        )

    metrics = pd.DataFrame(
        node_rows
    )

    metrics = metrics.sort_values(
        by=[
            "Participation_coefficient",
            "External_degree",
            "Betweenness",
            "Degree",
        ],
        ascending=[
            False,
            False,
            False,
            False,
        ],
    ).reset_index(
        drop=True
    )

    return metrics


# ============================================================
# Module summary
# ============================================================

def build_module_summary(
    H: nx.Graph,
    communities: list[set[str]],
    module_lookup: dict[str, str],
) -> pd.DataFrame:

    rows = []

    for idx, genes in enumerate(
        communities,
        start=1,
    ):

        module_id = (
            f"M{idx:02d}"
        )

        sub = H.subgraph(
            genes
        )

        up = 0
        down = 0

        core_count = 0
        interface_count = 0

        old_master_count = 0

        category_counter = Counter()

        lfcs = []

        for gene in genes:

            attrs = H.nodes[
                gene
            ]

            direction = str(
                attrs.get(
                    "Direction",
                    "",
                )
            )

            if direction == "Up":
                up += 1

            elif direction == "Down":
                down += 1

            tier = str(
                attrs.get(
                    "F2_tier",
                    "",
                )
            )

            if tier == "CORE":
                core_count += 1

            elif tier == "INTERFACE":
                interface_count += 1

            if bool_from_value(
                attrs.get(
                    "Old_126_gene_master",
                    False,
                )
            ):
                old_master_count += 1

            categories = (
                parse_categories(
                    attrs.get(
                        "F2_categories",
                        "",
                    )
                )
            )

            category_counter.update(
                categories
            )

            lfc = attrs.get(
                "log2FoldChange"
            )

            try:

                if not pd.isna(
                    lfc
                ):
                    lfcs.append(
                        float(
                            lfc
                        )
                    )

            except Exception:
                pass

        category_text = "; ".join(
            f"{category}:{count}"
            for category, count
            in category_counter.most_common()
        )

        mean_lfc = (
            sum(
                lfcs
            )
            /
            len(
                lfcs
            )
            if lfcs
            else float("nan")
        )

        mean_abs_lfc = (
            sum(
                abs(x)
                for x in lfcs
            )
            /
            len(
                lfcs
            )
            if lfcs
            else float("nan")
        )

        rows.append(
            {
                "Module":
                    module_id,

                "Size":
                    len(
                        genes
                    ),

                "Internal_edges":
                    sub.number_of_edges(),

                "Density":
                    (
                        nx.density(
                            sub
                        )
                        if len(
                            genes
                        ) > 1
                        else 0.0
                    ),

                "Up":
                    up,

                "Down":
                    down,

                "CORE":
                    core_count,

                "INTERFACE":
                    interface_count,

                "Old_master_overlap":
                    old_master_count,

                "Mean_log2FC":
                    mean_lfc,

                "Mean_abs_log2FC":
                    mean_abs_lfc,

                "F2_category_counts":
                    category_text,

                "Genes":
                    ";".join(
                        sorted(
                            genes
                        )
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# g:Profiler module enrichment
# ============================================================

def query_module_enrichment(
    module_genes: list[str],
    background_genes: list[str],
) -> pd.DataFrame:

    payload = {
        "organism":
            ORGANISM,

        "query":
            module_genes,

        "sources":
            ENRICHMENT_SOURCES,

        "user_threshold":
            FDR_THRESHOLD,

        "significance_threshold_method":
            "fdr",

        "domain_scope":
            "custom",

        "background":
            background_genes,

        "ordered":
            False,

        "all_results":
            False,

        "no_evidences":
            True,
    }

    last_error = None

    for attempt in range(
        1,
        4,
    ):

        try:

            response = requests.post(
                GPROFILER_URL,
                json=payload,
                timeout=REQUEST_TIMEOUT,
                headers={
                    "User-Agent":
                        "BIOL_295_FongDownstream_F4/1.0"
                },
            )

            response.raise_for_status()

            data = response.json()

            result = pd.DataFrame(
                data.get(
                    "result",
                    [],
                )
            )

            time.sleep(
                API_SLEEP_SECONDS
            )

            return result

        except Exception as exc:

            last_error = exc

            print(
                f"          [WARNING] enrichment attempt "
                f"{attempt}/3 failed: {exc}"
            )

            time.sleep(
                2
                *
                attempt
            )

    print(
        f"          [WARNING] module enrichment failed completely: "
        f"{last_error}"
    )

    return pd.DataFrame()


def add_enrichment_statistics(
    df: pd.DataFrame,
) -> pd.DataFrame:

    if df.empty:
        return df

    result = df.copy()

    result[
        "neg_log10_FDR"
    ] = (
        result[
            "p_value"
        ]
        .astype(
            float
        )
        .map(
            lambda p:
                -math.log10(
                    max(
                        p,
                        1e-300,
                    )
                )
        )
    )

    required = {
        "intersection_size",
        "query_size",
        "term_size",
        "effective_domain_size",
    }

    if required.issubset(
        result.columns
    ):

        result[
            "fold_enrichment"
        ] = (
            (
                result[
                    "intersection_size"
                ]
                /
                result[
                    "query_size"
                ]
            )
            /
            (
                result[
                    "term_size"
                ]
                /
                result[
                    "effective_domain_size"
                ]
            )
        )

    return result


def enrich_modules(
    network_name: str,
    communities: list[set[str]],
    background_genes: list[str],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
]:

    all_rows = []

    print()
    print(
        f"[{network_name}] Module enrichment "
        f"with {len(background_genes):,}-gene custom background..."
    )

    for idx, genes in enumerate(
        communities,
        start=1,
    ):

        module_id = (
            f"M{idx:02d}"
        )

        if len(
            genes
        ) < MIN_MODULE_SIZE_FOR_ENRICHMENT:

            print(
                f"       {module_id}: "
                f"size={len(genes)} -> skipped"
            )

            continue

        print(
            f"       {module_id}: "
            f"{len(genes)} genes"
        )

        result = query_module_enrichment(
            module_genes=sorted(
                genes
            ),
            background_genes=background_genes,
        )

        if result.empty:
            continue

        result = add_enrichment_statistics(
            result
        )

        result.insert(
            0,
            "Network",
            network_name,
        )

        result.insert(
            1,
            "Module",
            module_id,
        )

        all_rows.append(
            result
        )

    if not all_rows:

        return (
            pd.DataFrame(),
            pd.DataFrame(),
        )

    enrichment = pd.concat(
        all_rows,
        ignore_index=True,
    )

    enrichment = enrichment.sort_values(
        by=[
            "Module",
            "source",
            "p_value",
        ],
        ascending=[
            True,
            True,
            True,
        ],
    ).reset_index(
        drop=True
    )

    # --------------------------------------------------------
    # Compact top-term table
    # --------------------------------------------------------

    top_rows = []

    for (
        module,
        source
    ), group in enrichment.groupby(
        [
            "Module",
            "source",
        ]
    ):

        group = (
            group
            .sort_values(
                by=[
                    "p_value",
                    "fold_enrichment",
                ],
                ascending=[
                    True,
                    False,
                ],
            )
            .head(
                TOP_TERMS_PER_SOURCE
            )
        )

        top_rows.append(
            group
        )

    top = pd.concat(
        top_rows,
        ignore_index=True,
    )

    return (
        enrichment,
        top,
    )


# ============================================================
# Analyze one network
# ============================================================

def analyze_network(
    name: str,
    input_dir: Path,
    output_dir: Path,
) -> dict:

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    print()
    print(
        "=" * 76
    )

    print(
        f"{name} MODULE ANALYSIS"
    )

    print(
        "=" * 76
    )

    G, nodes, edges = load_network(
        input_dir
    )

    (
        H,
        communities,
        module_lookup,
        modularity,
    ) = detect_modules(
        G
    )

    print(
        f"[{name}] Total STRING nodes      : "
        f"{G.number_of_nodes():,}"
    )

    print(
        f"[{name}] Isolates                : "
        f"{len(list(nx.isolates(G))):,}"
    )

    print(
        f"[{name}] Non-isolate nodes       : "
        f"{H.number_of_nodes():,}"
    )

    print(
        f"[{name}] Non-isolate edges       : "
        f"{H.number_of_edges():,}"
    )

    print(
        f"[{name}] Louvain modules         : "
        f"{len(communities):,}"
    )

    print(
        f"[{name}] Weighted modularity Q   : "
        f"{modularity:.4f}"
    )

    # --------------------------------------------------------
    # Assignment table
    # --------------------------------------------------------

    assignment_rows = []

    node_lookup = (
        nodes
        .set_index(
            "Gene"
        )
    )

    for gene in G.nodes():

        is_isolate = (
            G.degree(
                gene
            ) == 0
        )

        module = (
            "ISOLATE"
            if is_isolate
            else module_lookup.get(
                gene
            )
        )

        row = {
            "Gene":
                gene,

            "Module":
                module,

            "Isolate":
                is_isolate,

            "Degree":
                G.degree(
                    gene
                ),
        }

        if gene in node_lookup.index:

            original = (
                node_lookup.loc[
                    gene
                ]
            )

            for column in [
                "Direction",
                "log2FoldChange",
                "padj",
                "F2_tier",
                "F2_categories",
                "Old_126_gene_master",
                "Component",
                "Component_size",
            ]:

                if column in node_lookup.columns:

                    row[
                        column
                    ] = original[
                        column
                    ]

        assignment_rows.append(
            row
        )

    assignments = pd.DataFrame(
        assignment_rows
    )

    assignments = assignments.sort_values(
        by=[
            "Isolate",
            "Module",
            "Degree",
            "Gene",
        ],
        ascending=[
            True,
            True,
            False,
            True,
        ],
    ).reset_index(
        drop=True
    )

    # --------------------------------------------------------
    # Module summary
    # --------------------------------------------------------

    summary = build_module_summary(
        H,
        communities,
        module_lookup,
    )

    # --------------------------------------------------------
    # Bridge metrics
    # --------------------------------------------------------

    bridge_metrics = (
        calculate_bridge_metrics(
            H,
            module_lookup,
        )
    )

    cross_module = (
        bridge_metrics[
            bridge_metrics[
                "External_degree"
            ] > 0
        ]
        .copy()
    )

    cross_module = (
        cross_module
        .sort_values(
            by=[
                "Participation_coefficient",
                "External_degree",
                "Betweenness",
                "Degree",
            ],
            ascending=[
                False,
                False,
                False,
                False,
            ],
        )
        .reset_index(
            drop=True
        )
    )

    # --------------------------------------------------------
    # Module enrichment
    # --------------------------------------------------------

    background_genes = sorted(
        H.nodes()
    )

    enrichment, enrichment_top = (
        enrich_modules(
            network_name=name,
            communities=communities,
            background_genes=background_genes,
        )
    )

    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    assignments.to_csv(
        output_dir
        / "module_assignments.csv",
        index=False,
    )

    summary.to_csv(
        output_dir
        / "module_summary.csv",
        index=False,
    )

    bridge_metrics.to_csv(
        output_dir
        / "bridge_metrics.csv",
        index=False,
    )

    cross_module.to_csv(
        output_dir
        / "cross_module_bridges.csv",
        index=False,
    )

    enrichment.to_csv(
        output_dir
        / "module_enrichment.csv",
        index=False,
    )

    enrichment_top.to_csv(
        output_dir
        / "module_enrichment_top.csv",
        index=False,
    )

    return {
        "name":
            name,

        "graph":
            G,

        "non_isolate_graph":
            H,

        "communities":
            communities,

        "module_lookup":
            module_lookup,

        "modularity":
            modularity,

        "assignments":
            assignments,

        "summary":
            summary,

        "bridge_metrics":
            bridge_metrics,

        "cross_module":
            cross_module,

        "enrichment":
            enrichment,

        "enrichment_top":
            enrichment_top,
    }


# ============================================================
# EXTENDED: interface → CORE bridge analysis
# ============================================================

def build_interface_core_bridge_table(
    extended_result: dict,
    core_result: dict,
) -> pd.DataFrame:

    H = extended_result[
        "non_isolate_graph"
    ]

    ext_metrics = (
        extended_result[
            "bridge_metrics"
        ]
        .set_index(
            "Gene"
        )
    )

    core_module_lookup = (
        core_result[
            "module_lookup"
        ]
    )

    rows = []

    for gene in H.nodes():

        attrs = H.nodes[
            gene
        ]

        if str(
            attrs.get(
                "F2_tier",
                ""
            )
        ) != "INTERFACE":
            continue

        core_neighbors = []

        core_neighbor_modules = set()

        core_neighbor_module_counts = Counter()

        for neighbor in H.neighbors(
            gene
        ):

            if neighbor not in core_module_lookup:
                continue

            core_neighbors.append(
                neighbor
            )

            module = (
                core_module_lookup[
                    neighbor
                ]
            )

            core_neighbor_modules.add(
                module
            )

            core_neighbor_module_counts[
                module
            ] += 1

        metric = (
            ext_metrics.loc[
                gene
            ]
        )

        rows.append(
            {
                "Gene":
                    gene,

                "EXTENDED_module":
                    extended_result[
                        "module_lookup"
                    ].get(
                        gene
                    ),

                "Degree":
                    metric[
                        "Degree"
                    ],

                "External_degree":
                    metric[
                        "External_degree"
                    ],

                "Participation_coefficient":
                    metric[
                        "Participation_coefficient"
                    ],

                "Betweenness":
                    metric[
                        "Betweenness"
                    ],

                "PageRank":
                    metric[
                        "PageRank"
                    ],

                "Direction":
                    metric[
                        "Direction"
                    ],

                "log2FoldChange":
                    metric[
                        "log2FoldChange"
                    ],

                "CORE_neighbor_count":
                    len(
                        core_neighbors
                    ),

                "Distinct_CORE_modules_touched":
                    len(
                        core_neighbor_modules
                    ),

                "CORE_modules_touched":
                    ";".join(
                        sorted(
                            core_neighbor_modules
                        )
                    ),

                "CORE_module_neighbor_counts":
                    "; ".join(
                        f"{module}:{count}"
                        for module, count
                        in sorted(
                            core_neighbor_module_counts.items()
                        )
                    ),

                "CORE_neighbors":
                    ";".join(
                        sorted(
                            core_neighbors
                        )
                    ),
            }
        )

    table = pd.DataFrame(
        rows
    )

    if table.empty:
        return table

    return (
        table
        .sort_values(
            by=[
                "Distinct_CORE_modules_touched",
                "CORE_neighbor_count",
                "Participation_coefficient",
                "Betweenness",
            ],
            ascending=[
                False,
                False,
                False,
                False,
            ],
        )
        .reset_index(
            drop=True
        )
    )


# ============================================================
# CORE ↔ EXTENDED module correspondence
# ============================================================

def build_module_correspondence(
    core_result: dict,
    extended_result: dict,
) -> pd.DataFrame:

    core_lookup = (
        core_result[
            "module_lookup"
        ]
    )

    extended_lookup = (
        extended_result[
            "module_lookup"
        ]
    )

    H_ext = (
        extended_result[
            "non_isolate_graph"
        ]
    )

    rows = []

    for ext_module in sorted(
        set(
            extended_lookup.values()
        )
    ):

        ext_genes = {
            gene
            for gene, module
            in extended_lookup.items()
            if module
            == ext_module
        }

        core_genes = sorted(
            gene
            for gene
            in ext_genes
            if gene
            in core_lookup
        )

        interface_genes = sorted(
            gene
            for gene
            in ext_genes
            if str(
                H_ext.nodes[
                    gene
                ].get(
                    "F2_tier",
                    ""
                )
            )
            == "INTERFACE"
        )

        represented_core_modules = (
            Counter(
                core_lookup[
                    gene
                ]
                for gene
                in core_genes
            )
        )

        rows.append(
            {
                "EXTENDED_module":
                    ext_module,

                "EXTENDED_size":
                    len(
                        ext_genes
                    ),

                "CORE_genes":
                    len(
                        core_genes
                    ),

                "INTERFACE_genes":
                    len(
                        interface_genes
                    ),

                "Distinct_CORE_modules_represented":
                    len(
                        represented_core_modules
                    ),

                "CORE_module_composition":
                    "; ".join(
                        f"{module}:{count}"
                        for module, count
                        in sorted(
                            represented_core_modules.items()
                        )
                    ),

                "CORE_gene_list":
                    ";".join(
                        core_genes
                    ),

                "INTERFACE_gene_list":
                    ";".join(
                        interface_genes
                    ),
            }
        )

    return pd.DataFrame(
        rows
    ).sort_values(
        by=[
            "EXTENDED_size",
            "EXTENDED_module",
        ],
        ascending=[
            False,
            True,
        ],
    )


# ============================================================
# Excel formatting
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

            cell.fill = header_fill

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
                column_cells[0].column
            )

            max_length = 0

            for cell in column_cells[:100]:

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
                50,
            )

    wb.save(
        path
    )


# ============================================================
# Validate inputs
# ============================================================

print()
print("=" * 76)
print("PHASE F4 — STRING MODULE / BRIDGE ANALYSIS")
print("=" * 76)

for required_file in [
    CORE_INPUT_DIR / "nodes.csv",
    CORE_INPUT_DIR / "edges.csv",
    EXTENDED_INPUT_DIR / "nodes.csv",
    EXTENDED_INPUT_DIR / "edges.csv",
]:

    if not required_file.exists():

        raise FileNotFoundError(
            f"Required F3 file not found:\n"
            f"{required_file}"
        )


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

CORE_OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

EXTENDED_OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Analyze CORE
# ============================================================

core_result = analyze_network(
    name="CORE",
    input_dir=CORE_INPUT_DIR,
    output_dir=CORE_OUTPUT_DIR,
)


# ============================================================
# Analyze EXTENDED
# ============================================================

extended_result = analyze_network(
    name="EXTENDED",
    input_dir=EXTENDED_INPUT_DIR,
    output_dir=EXTENDED_OUTPUT_DIR,
)


# ============================================================
# Interface bridge analysis
# ============================================================

print()
print(
    "[EXTENDED] Evaluating INTERFACE proteins against "
    "standalone CORE modules..."
)

interface_core_bridges = (
    build_interface_core_bridge_table(
        extended_result=extended_result,
        core_result=core_result,
    )
)

interface_core_bridges.to_csv(
    EXTENDED_OUTPUT_DIR
    / "interface_core_bridges.csv",
    index=False,
)


# ============================================================
# Module correspondence
# ============================================================

correspondence = (
    build_module_correspondence(
        core_result=core_result,
        extended_result=extended_result,
    )
)

correspondence.to_csv(
    CORRESPONDENCE_CSV,
    index=False,
)


# ============================================================
# Compact network summary
# ============================================================

network_summary = pd.DataFrame(
    [
        {
            "Network":
                "CORE",

            "Total_nodes":
                core_result[
                    "graph"
                ].number_of_nodes(),

            "Non_isolate_nodes":
                core_result[
                    "non_isolate_graph"
                ].number_of_nodes(),

            "Edges":
                core_result[
                    "non_isolate_graph"
                ].number_of_edges(),

            "Modules":
                len(
                    core_result[
                        "communities"
                    ]
                ),

            "Weighted_modularity_Q":
                core_result[
                    "modularity"
                ],

            "Cross_module_nodes":
                len(
                    core_result[
                        "cross_module"
                    ]
                ),
        },

        {
            "Network":
                "EXTENDED",

            "Total_nodes":
                extended_result[
                    "graph"
                ].number_of_nodes(),

            "Non_isolate_nodes":
                extended_result[
                    "non_isolate_graph"
                ].number_of_nodes(),

            "Edges":
                extended_result[
                    "non_isolate_graph"
                ].number_of_edges(),

            "Modules":
                len(
                    extended_result[
                        "communities"
                    ]
                ),

            "Weighted_modularity_Q":
                extended_result[
                    "modularity"
                ],

            "Cross_module_nodes":
                len(
                    extended_result[
                        "cross_module"
                    ]
                ),
        },
    ]
)


# ============================================================
# Export Excel
# ============================================================

with pd.ExcelWriter(
    OUTPUT_XLSX,
    engine="openpyxl",
) as writer:

    network_summary.to_excel(
        writer,
        sheet_name="Network summary",
        index=False,
    )

    core_result[
        "summary"
    ].to_excel(
        writer,
        sheet_name="CORE modules",
        index=False,
    )

    core_result[
        "assignments"
    ].to_excel(
        writer,
        sheet_name="CORE assignments",
        index=False,
    )

    core_result[
        "cross_module"
    ].to_excel(
        writer,
        sheet_name="CORE bridges",
        index=False,
    )

    core_result[
        "enrichment_top"
    ].to_excel(
        writer,
        sheet_name="CORE enrich top",
        index=False,
    )

    extended_result[
        "summary"
    ].to_excel(
        writer,
        sheet_name="EXT modules",
        index=False,
    )

    extended_result[
        "assignments"
    ].to_excel(
        writer,
        sheet_name="EXT assignments",
        index=False,
    )

    extended_result[
        "cross_module"
    ].to_excel(
        writer,
        sheet_name="EXT bridges",
        index=False,
    )

    interface_core_bridges.to_excel(
        writer,
        sheet_name="Interface bridges",
        index=False,
    )

    extended_result[
        "enrichment_top"
    ].to_excel(
        writer,
        sheet_name="EXT enrich top",
        index=False,
    )

    correspondence.to_excel(
        writer,
        sheet_name="Module correspondence",
        index=False,
    )


format_workbook(
    OUTPUT_XLSX
)


# ============================================================
# Metadata
# ============================================================

metadata = {
    "phase":
        "F4_network_modules",

    "community_detection": {
        "algorithm":
            "NetworkX Louvain",

        "weight":
            "STRING combined_score",

        "resolution":
            LOUVAIN_RESOLUTION,

        "seed":
            LOUVAIN_SEED,

        "isolates_excluded":
            True,
    },

    "bridge_metrics": {
        "participation_coefficient":
            (
                "1 - sum_s((k_is / k_i)^2)"
            ),

        "within_module_degree_z":
            (
                "z-score of within-module degree "
                "relative to nodes in the same module"
            ),

        "betweenness":
            (
                "Unweighted shortest-path betweenness"
            ),

        "pagerank":
            (
                "Weighted using STRING combined score"
            ),
    },

    "module_enrichment": {
        "sources":
            ENRICHMENT_SOURCES,

        "correction":
            "FDR",

        "threshold":
            FDR_THRESHOLD,

        "background":
            (
                "All non-isolated proteins in the same "
                "STRING network"
            ),

        "minimum_module_size":
            MIN_MODULE_SIZE_FOR_ENRICHMENT,
    },

    "core": {
        "modules":
            len(
                core_result[
                    "communities"
                ]
            ),

        "weighted_modularity":
            core_result[
                "modularity"
            ],
    },

    "extended": {
        "modules":
            len(
                extended_result[
                    "communities"
                ]
            ),

        "weighted_modularity":
            extended_result[
                "modularity"
            ],
    },
}


with open(
    OUTPUT_METADATA,
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        metadata,
        handle,
        indent=2,
        ensure_ascii=False,
    )


# ============================================================
# Console report
# ============================================================

print()
print("=" * 76)
print("PHASE F4 COMPLETE")
print("=" * 76)

print()

print(
    "CORE"
)

print(
    f"    Non-isolate proteins      : "
    f"{core_result['non_isolate_graph'].number_of_nodes():,}"
)

print(
    f"    Modules                   : "
    f"{len(core_result['communities']):,}"
)

print(
    f"    Weighted modularity Q     : "
    f"{core_result['modularity']:.4f}"
)

print(
    f"    Cross-module proteins     : "
    f"{len(core_result['cross_module']):,}"
)

print()

print(
    "EXTENDED"
)

print(
    f"    Non-isolate proteins      : "
    f"{extended_result['non_isolate_graph'].number_of_nodes():,}"
)

print(
    f"    Modules                   : "
    f"{len(extended_result['communities']):,}"
)

print(
    f"    Weighted modularity Q     : "
    f"{extended_result['modularity']:.4f}"
)

print(
    f"    Cross-module proteins     : "
    f"{len(extended_result['cross_module']):,}"
)

print()

if not interface_core_bridges.empty:

    multi_core = (
        interface_core_bridges[
            interface_core_bridges[
                "Distinct_CORE_modules_touched"
            ] >= 2
        ]
    )

    print(
        "INTERFACE → CORE"
    )

    print(
        f"    Interface proteins with CORE neighbors : "
        f"{(interface_core_bridges['CORE_neighbor_count'] > 0).sum():,}"
    )

    print(
        f"    Touching >=2 CORE modules              : "
        f"{len(multi_core):,}"
    )

    if not multi_core.empty:

        print()
        print(
            "    Top multi-CORE interface candidates:"
        )

        for _, row in (
            multi_core
            .head(
                15
            )
            .iterrows()
        ):

            print(
                f"        "
                f"{row['Gene']:<12} "
                f"CORE neighbors={int(row['CORE_neighbor_count']):>2}  "
                f"CORE modules={int(row['Distinct_CORE_modules_touched']):>2}  "
                f"P={row['Participation_coefficient']:.3f}  "
                f"betweenness={row['Betweenness']:.4f}"
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
    "[NEXT] Inspect the module compositions, enrichment terms, "
    "and cross-module bridge proteins before beginning any "
    "directed mechanistic analysis."
)