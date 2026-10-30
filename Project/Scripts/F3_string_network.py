#!/usr/bin/env python3

"""
Phase F3 — STRING network analysis of the Fong cell-fate universe
=================================================================

Purpose
-------
Construct two high-confidence STRING functional-association networks:

    CORE
        F2 CORE death/survival/autophagy genes

    EXTENDED
        CORE + DNA-damage/checkpoint/repair +
        ER/oxidative-stress interface genes

This allows us to ask two distinct questions:

1. What is the structure of the downstream cell-fate machinery itself?
2. Does the stress/checkpoint interface connect otherwise separated
   parts of that cell-fate network?

Network definition
------------------
STRING version: 12.5
Species: Mus musculus (NCBI taxon 10090)
Network type: functional
Required score: >= 700
Added nodes: 0

IMPORTANT
---------
STRING functional edges represent functional associations, not
necessarily direct physical binding and not causal direction.

Input
-----
Data/Processed/FongDownstream/F2_CellFate/

    core_cell_fate.csv
    cell_fate_master.csv

Outputs
-------
Data/Processed/FongDownstream/F3_STRING/

    F3_STRING_summary.xlsx
    network_comparison.csv
    F3_metadata.json

    CORE/
        mapping.csv
        unmapped.csv
        nodes.csv
        edges.csv
        components.csv
        network.graphml
        largest_component.png
        degree_distribution.png

    EXTENDED/
        mapping.csv
        unmapped.csv
        nodes.csv
        edges.csv
        components.csv
        network.graphml
        largest_component.png
        degree_distribution.png

Run
---
From Project/:

    py .\\Scripts\\F3_string_network.py
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
import json
import math
import time

import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt

from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


# ============================================================
# Project paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

F2_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F2_CellFate"
)

CORE_FILE = (
    F2_DIR
    / "core_cell_fate.csv"
)

EXTENDED_FILE = (
    F2_DIR
    / "cell_fate_master.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F3_STRING"
)

CORE_OUTPUT = OUTPUT_DIR / "CORE"
EXTENDED_OUTPUT = OUTPUT_DIR / "EXTENDED"

SUMMARY_XLSX = (
    OUTPUT_DIR
    / "F3_STRING_summary.xlsx"
)

COMPARISON_CSV = (
    OUTPUT_DIR
    / "network_comparison.csv"
)

METADATA_JSON = (
    OUTPUT_DIR
    / "F3_metadata.json"
)


# ============================================================
# STRING configuration
# ============================================================

STRING_VERSION = "12.5"

STRING_API_BASE = (
    "https://version-12-5.string-db.org/api"
)

SPECIES = 10090  # Mus musculus

REQUIRED_SCORE = 700

NETWORK_TYPE = "functional"

ADD_NODES = 0

CALLER_IDENTITY = (
    "BIOL_295_FongDownstream"
)

REQUEST_TIMEOUT = 180

# STRING asks API users to avoid hammering the public endpoint.
API_SLEEP_SECONDS = 1.0


# ============================================================
# Helpers
# ============================================================

def normalize_gene(value: Any) -> str | None:
    if pd.isna(value):
        return None

    gene = str(value).strip()

    if not gene:
        return None

    return gene.upper()


def bool_from_value(value: Any) -> bool:
    if isinstance(value, bool):
        return value

    if pd.isna(value):
        return False

    return str(value).strip().upper() in {
        "TRUE",
        "YES",
        "1",
    }


def request_json(
    method: str,
    params: dict,
) -> Any:
    """
    POST to one STRING JSON endpoint.
    """

    url = (
        f"{STRING_API_BASE}/json/{method}"
    )

    response = requests.post(
        url,
        data=params,
        timeout=REQUEST_TIMEOUT,
        headers={
            "User-Agent":
                "BIOL_295_FongDownstream/1.0"
        },
    )

    response.raise_for_status()

    time.sleep(
        API_SLEEP_SECONDS
    )

    return response.json()


# ============================================================
# STRING mapping
# ============================================================

def map_genes_to_string(
    genes: list[str],
) -> pd.DataFrame:
    """
    Map gene symbols to STRING identifiers.

    STRING returns one best match per query identifier.
    """

    params = {
        "identifiers":
            "\r".join(
                genes
            ),

        "species":
            SPECIES,

        "echo_query":
            1,

        "caller_identity":
            CALLER_IDENTITY,
    }

    result = request_json(
        "get_string_ids",
        params,
    )

    if not result:
        return pd.DataFrame(
            columns=[
                "queryItem",
                "stringId",
                "preferredName",
                "ncbiTaxonId",
                "taxonName",
                "annotation",
            ]
        )

    mapping = pd.DataFrame(
        result
    )

    if "queryItem" not in mapping.columns:
        raise RuntimeError(
            "STRING mapping response did not contain queryItem."
        )

    mapping[
        "Gene"
    ] = (
        mapping[
            "queryItem"
        ]
        .map(
            normalize_gene
        )
    )

    return mapping


# ============================================================
# STRING network
# ============================================================

def query_network(
    string_ids: list[str],
) -> pd.DataFrame:
    """
    Retrieve high-confidence functional associations among
    the supplied STRING proteins only.
    """

    if len(string_ids) < 2:
        return pd.DataFrame()

    params = {
        "identifiers":
            "\r".join(
                string_ids
            ),

        "species":
            SPECIES,

        "required_score":
            REQUIRED_SCORE,

        "network_type":
            NETWORK_TYPE,

        "add_nodes":
            ADD_NODES,

        "caller_identity":
            CALLER_IDENTITY,
    }

    result = request_json(
        "network",
        params,
    )

    if not result:
        return pd.DataFrame()

    return pd.DataFrame(
        result
    )


# ============================================================
# PPI enrichment
# ============================================================

def query_ppi_enrichment(
    string_ids: list[str],
) -> dict:
    """
    Retrieve STRING PPI-enrichment statistics.
    """

    if len(string_ids) < 2:
        return {}

    params = {
        "identifiers":
            "\r".join(
                string_ids
            ),

        "species":
            SPECIES,

        "required_score":
            REQUIRED_SCORE,

        "caller_identity":
            CALLER_IDENTITY,
    }

    result = request_json(
        "ppi_enrichment",
        params,
    )

    if not result:
        return {}

    return result[0]


# ============================================================
# Deduplicate network edges
# ============================================================

def deduplicate_edges(
    edges: pd.DataFrame,
) -> pd.DataFrame:
    """
    STRING functional networks are undirected for our purposes.

    Store each pair once.
    """

    if edges.empty:
        return edges.copy()

    df = edges.copy()

    df["_pair"] = df.apply(
        lambda row: "||".join(
            sorted(
                [
                    str(row["stringId_A"]),
                    str(row["stringId_B"]),
                ]
            )
        ),
        axis=1,
    )

    if "score" in df.columns:

        df = (
            df.sort_values(
                "score",
                ascending=False,
            )
            .drop_duplicates(
                "_pair",
                keep="first",
            )
        )

    else:

        df = df.drop_duplicates(
            "_pair",
            keep="first",
        )

    return (
        df.drop(
            columns=[
                "_pair",
            ]
        )
        .reset_index(
            drop=True
        )
    )


# ============================================================
# Build NetworkX graph
# ============================================================

def build_graph(
    input_df: pd.DataFrame,
    mapping: pd.DataFrame,
    raw_edges: pd.DataFrame,
) -> tuple[
    nx.Graph,
    pd.DataFrame,
    pd.DataFrame,
]:
    """
    Build an undirected graph using original Fong gene symbols.
    """

    gene_to_string = dict(
        zip(
            mapping[
                "Gene"
            ],
            mapping[
                "stringId"
            ],
        )
    )

    gene_to_preferred = dict(
        zip(
            mapping[
                "Gene"
            ],
            mapping[
                "preferredName"
            ],
        )
    )

    # Detect STRING IDs mapping to more than one query gene.
    collisions = (
        mapping
        .groupby(
            "stringId"
        )[
            "Gene"
        ]
        .agg(
            lambda x:
                sorted(
                    set(x)
                )
        )
    )

    collisions = collisions[
        collisions.map(
            len
        )
        > 1
    ]

    if not collisions.empty:

        print(
            "[WARNING] Several input genes map to the same STRING ID:"
        )

        for string_id, genes in collisions.items():

            print(
                f"          {string_id}: "
                f"{', '.join(genes)}"
            )

    # For ordinary one-to-one mapping.
    string_to_gene = {}

    for _, row in mapping.iterrows():

        string_id = row[
            "stringId"
        ]

        gene = row[
            "Gene"
        ]

        if string_id not in string_to_gene:
            string_to_gene[
                string_id
            ] = gene

    G = nx.Graph()

    # Add every mapped query protein first.
    # This preserves true STRING isolates.
    for _, row in input_df.iterrows():

        gene = normalize_gene(
            row[
                "Gene"
            ]
        )

        if gene not in gene_to_string:
            continue

        attrs = {
            "stringId":
                gene_to_string[
                    gene
                ],

            "preferredName":
                str(
                    gene_to_preferred.get(
                        gene,
                        "",
                    )
                ),

            "Direction":
                str(
                    row.get(
                        "Direction",
                        "",
                    )
                ),

            "F2_tier":
                str(
                    row.get(
                        "F2_tier",
                        "",
                    )
                ),

            "log2FoldChange":
                float(
                    row.get(
                        "log2FoldChange",
                        0.0,
                    )
                ),

            "Old_126_gene_master":
                bool_from_value(
                    row.get(
                        "Old_126_gene_master",
                        False,
                    )
                ),
        }

        G.add_node(
            gene,
            **attrs,
        )

    # Deduplicate STRING edges.
    edges = deduplicate_edges(
        raw_edges
    )

    output_edge_rows = []

    for _, row in edges.iterrows():

        string_a = row[
            "stringId_A"
        ]

        string_b = row[
            "stringId_B"
        ]

        gene_a = string_to_gene.get(
            string_a
        )

        gene_b = string_to_gene.get(
            string_b
        )

        if (
            gene_a is None
            or
            gene_b is None
        ):
            continue

        if gene_a == gene_b:
            continue

        score = float(
            row.get(
                "score",
                0.0,
            )
        )

        G.add_edge(
            gene_a,
            gene_b,
            score=score,
        )

        output_edge_rows.append(
            {
                "Gene_A":
                    gene_a,

                "Gene_B":
                    gene_b,

                "STRING_ID_A":
                    string_a,

                "STRING_ID_B":
                    string_b,

                "PreferredName_A":
                    row.get(
                        "preferredName_A"
                    ),

                "PreferredName_B":
                    row.get(
                        "preferredName_B"
                    ),

                "combined_score":
                    score,

                "score_1000":
                    score * 1000.0,

                "neighborhood_score":
                    row.get(
                        "nscore"
                    ),

                "fusion_score":
                    row.get(
                        "fscore"
                    ),

                "phylogenetic_score":
                    row.get(
                        "pscore"
                    ),

                "coexpression_score":
                    row.get(
                        "ascore"
                    ),

                "experimental_score":
                    row.get(
                        "escore"
                    ),

                "database_score":
                    row.get(
                        "dscore"
                    ),

                "textmining_score":
                    row.get(
                        "tscore"
                    ),
            }
        )

    edge_table = pd.DataFrame(
        output_edge_rows
    )

    # --------------------------------------------------------
    # Component assignments
    # --------------------------------------------------------

    components = list(
        nx.connected_components(
            G
        )
    )

    components = sorted(
        components,
        key=len,
        reverse=True,
    )

    component_id = {}

    for idx, component in enumerate(
        components,
        start=1,
    ):

        for gene in component:

            component_id[
                gene
            ] = idx

    # --------------------------------------------------------
    # Node table
    # --------------------------------------------------------

    rows = []

    input_lookup = (
        input_df
        .copy()
        .assign(
            Gene=lambda x:
                x[
                    "Gene"
                ]
                .map(
                    normalize_gene
                )
        )
        .set_index(
            "Gene"
        )
    )

    for gene in G.nodes():

        original = input_lookup.loc[
            gene
        ]

        degree = G.degree(
            gene
        )

        rows.append(
            {
                "Gene":
                    gene,

                "STRING_ID":
                    gene_to_string.get(
                        gene
                    ),

                "PreferredName":
                    gene_to_preferred.get(
                        gene
                    ),

                "Direction":
                    original.get(
                        "Direction"
                    ),

                "log2FoldChange":
                    original.get(
                        "log2FoldChange"
                    ),

                "padj":
                    original.get(
                        "padj"
                    ),

                "F2_tier":
                    original.get(
                        "F2_tier"
                    ),

                "F2_categories":
                    original.get(
                        "F2_categories"
                    ),

                "Old_126_gene_master":
                    original.get(
                        "Old_126_gene_master"
                    ),

                "Degree":
                    degree,

                "Isolate":
                    degree == 0,

                "Component":
                    component_id.get(
                        gene
                    ),

                "Component_size":
                    len(
                        components[
                            component_id[
                                gene
                            ]
                            - 1
                        ]
                    ),
            }
        )

    node_table = pd.DataFrame(
        rows
    )

    if not node_table.empty:

        node_table = (
            node_table
            .sort_values(
                by=[
                    "Degree",
                    "Gene",
                ],
                ascending=[
                    False,
                    True,
                ],
            )
            .reset_index(
                drop=True
            )
        )

    return (
        G,
        node_table,
        edge_table,
    )


# ============================================================
# Component table
# ============================================================

def component_table(
    G: nx.Graph,
) -> pd.DataFrame:

    components = sorted(
        nx.connected_components(
            G
        ),
        key=len,
        reverse=True,
    )

    rows = []

    for idx, component in enumerate(
        components,
        start=1,
    ):

        sub = G.subgraph(
            component
        )

        genes = sorted(
            component
        )

        rows.append(
            {
                "Component":
                    idx,

                "Size":
                    len(
                        component
                    ),

                "Edges":
                    sub.number_of_edges(),

                "Isolate_component":
                    len(component) == 1,

                "Genes":
                    ";".join(
                        genes
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# Network summary
# ============================================================

def summarize_network(
    name: str,
    input_gene_count: int,
    mapping: pd.DataFrame,
    G: nx.Graph,
    ppi: dict,
) -> dict:

    mapped = G.number_of_nodes()

    edges = G.number_of_edges()

    isolates = len(
        list(
            nx.isolates(
                G
            )
        )
    )

    components = (
        nx.number_connected_components(
            G
        )
        if mapped
        else 0
    )

    non_isolate_components = (
        components
        -
        isolates
    )

    largest_component = 0

    if mapped:

        largest_component = max(
            len(c)
            for c in nx.connected_components(
                G
            )
        )

    average_degree = (
        sum(
            dict(
                G.degree()
            ).values()
        )
        /
        mapped
        if mapped
        else 0.0
    )

    density = (
        nx.density(
            G
        )
        if mapped > 1
        else 0.0
    )

    expected_edges = ppi.get(
        "expected_number_of_edges"
    )

    p_value = ppi.get(
        "p_value"
    )

    enrichment_ratio = None

    if (
        expected_edges is not None
        and
        float(
            expected_edges
        )
        > 0
    ):

        enrichment_ratio = (
            edges
            /
            float(
                expected_edges
            )
        )

    return {
        "Network":
            name,

        "Input_genes":
            input_gene_count,

        "STRING_mapped_genes":
            mapped,

        "Mapping_rate":
            (
                mapped
                /
                input_gene_count
                if input_gene_count
                else 0
            ),

        "Edges":
            edges,

        "Average_degree_local":
            average_degree,

        "Density_local":
            density,

        "Connected_components":
            components,

        "Isolates":
            isolates,

        "Non_isolate_components":
            non_isolate_components,

        "Largest_component":
            largest_component,

        "STRING_expected_edges":
            expected_edges,

        "Observed_expected_ratio":
            enrichment_ratio,

        "STRING_PPI_enrichment_p":
            p_value,

        "STRING_average_node_degree":
            ppi.get(
                "average_node_degree"
            ),

        "STRING_local_clustering_coefficient":
            ppi.get(
                "local_clustering_coefficient"
            ),
    }


# ============================================================
# Figures
# ============================================================

def draw_largest_component(
    G: nx.Graph,
    output_file: Path,
    title: str,
) -> None:

    if (
        G.number_of_nodes()
        == 0
    ):
        return

    largest = max(
        nx.connected_components(
            G
        ),
        key=len,
    )

    sub = G.subgraph(
        largest
    ).copy()

    plt.figure(
        figsize=(
            12,
            10,
        )
    )

    positions = nx.spring_layout(
        sub,
        seed=42,
    )

    degrees = dict(
        sub.degree()
    )

    node_sizes = [
        100
        +
        35
        *
        degrees[
            node
        ]
        for node in sub.nodes()
    ]

    nx.draw_networkx_edges(
        sub,
        positions,
        alpha=0.35,
        width=0.8,
    )

    nx.draw_networkx_nodes(
        sub,
        positions,
        node_size=node_sizes,
        alpha=0.85,
    )

    # Keep labels usable rather than turning a large graph into
    # an unreadable hairball.
    if sub.number_of_nodes() <= 70:

        nx.draw_networkx_labels(
            sub,
            positions,
            font_size=7,
        )

    plt.title(
        f"{title}\n"
        f"Largest connected component: "
        f"{sub.number_of_nodes()} nodes, "
        f"{sub.number_of_edges()} edges"
    )

    plt.axis(
        "off"
    )

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def draw_degree_distribution(
    G: nx.Graph,
    output_file: Path,
    title: str,
) -> None:

    degrees = [
        degree
        for _, degree
        in G.degree()
    ]

    if not degrees:
        return

    max_degree = max(
        degrees
    )

    bins = range(
        0,
        max_degree + 2,
    )

    plt.figure(
        figsize=(
            8,
            5,
        )
    )

    plt.hist(
        degrees,
        bins=bins,
        align="left",
    )

    plt.xlabel(
        "Node degree"
    )

    plt.ylabel(
        "Number of proteins"
    )

    plt.title(
        title
    )

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


# ============================================================
# GraphML export
# ============================================================

def write_graphml(
    G: nx.Graph,
    path: Path,
) -> None:
    """
    Export simple graph attributes compatible with GraphML.
    """

    clean = nx.Graph()

    for node, attrs in G.nodes(
        data=True
    ):

        clean.add_node(
            node,

            stringId=str(
                attrs.get(
                    "stringId",
                    "",
                )
            ),

            preferredName=str(
                attrs.get(
                    "preferredName",
                    "",
                )
            ),

            Direction=str(
                attrs.get(
                    "Direction",
                    "",
                )
            ),

            F2_tier=str(
                attrs.get(
                    "F2_tier",
                    "",
                )
            ),

            log2FoldChange=float(
                attrs.get(
                    "log2FoldChange",
                    0.0,
                )
            ),

            Old_126_gene_master=bool(
                attrs.get(
                    "Old_126_gene_master",
                    False,
                )
            ),
        )

    for a, b, attrs in G.edges(
        data=True
    ):

        clean.add_edge(
            a,
            b,
            score=float(
                attrs.get(
                    "score",
                    0.0,
                )
            ),
        )

    nx.write_graphml(
        clean,
        path,
    )


# ============================================================
# Workbook formatting
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
                45,
            )

    wb.save(
        path
    )


# ============================================================
# Analyze one network
# ============================================================

def analyze_network(
    name: str,
    input_df: pd.DataFrame,
    output_dir: Path,
) -> dict:

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    genes = sorted(
        {
            normalize_gene(
                gene
            )
            for gene
            in input_df[
                "Gene"
            ]
            if normalize_gene(
                gene
            )
        }
    )

    print()
    print(
        "-" * 76
    )

    print(
        f"[{name}] Input genes: "
        f"{len(genes):,}"
    )

    # --------------------------------------------------------
    # Mapping
    # --------------------------------------------------------

    print(
        f"[{name}] Mapping genes to STRING..."
    )

    mapping = map_genes_to_string(
        genes
    )

    mapping.to_csv(
        output_dir
        / "mapping.csv",
        index=False,
    )

    mapped_genes = set(
        mapping[
            "Gene"
        ]
    )

    unmapped_genes = sorted(
        set(
            genes
        )
        -
        mapped_genes
    )

    unmapped = pd.DataFrame(
        {
            "Gene":
                unmapped_genes
        }
    )

    unmapped.to_csv(
        output_dir
        / "unmapped.csv",
        index=False,
    )

    print(
        f"[{name}] STRING mapped: "
        f"{len(mapped_genes):,}/{len(genes):,}"
    )

    if unmapped_genes:

        print(
            f"[{name}] Unmapped genes: "
            f"{len(unmapped_genes):,}"
        )

    string_ids = sorted(
        set(
            mapping[
                "stringId"
            ]
        )
    )

    # --------------------------------------------------------
    # Network query
    # --------------------------------------------------------

    print(
        f"[{name}] Querying STRING network "
        f"(score >= {REQUIRED_SCORE}, add_nodes=0)..."
    )

    raw_edges = query_network(
        string_ids
    )

    # --------------------------------------------------------
    # PPI enrichment
    # --------------------------------------------------------

    print(
        f"[{name}] Querying STRING PPI enrichment..."
    )

    ppi = query_ppi_enrichment(
        string_ids
    )

    # --------------------------------------------------------
    # Local graph
    # --------------------------------------------------------

    G, nodes, edges = build_graph(
        input_df=input_df,
        mapping=mapping,
        raw_edges=raw_edges,
    )

    components = component_table(
        G
    )

    nodes.to_csv(
        output_dir
        / "nodes.csv",
        index=False,
    )

    edges.to_csv(
        output_dir
        / "edges.csv",
        index=False,
    )

    components.to_csv(
        output_dir
        / "components.csv",
        index=False,
    )

    write_graphml(
        G,
        output_dir
        / "network.graphml",
    )

    draw_largest_component(
        G,
        output_dir
        / "largest_component.png",
        title=(
            f"{name} STRING network"
        ),
    )

    draw_degree_distribution(
        G,
        output_dir
        / "degree_distribution.png",
        title=(
            f"{name} STRING degree distribution"
        ),
    )

    summary = summarize_network(
        name=name,
        input_gene_count=len(
            genes
        ),
        mapping=mapping,
        G=G,
        ppi=ppi,
    )

    print(
        f"[{name}] Nodes        : "
        f"{G.number_of_nodes():,}"
    )

    print(
        f"[{name}] Edges        : "
        f"{G.number_of_edges():,}"
    )

    print(
        f"[{name}] Components   : "
        f"{summary['Connected_components']:,}"
    )

    print(
        f"[{name}] Isolates     : "
        f"{summary['Isolates']:,}"
    )

    print(
        f"[{name}] Largest comp : "
        f"{summary['Largest_component']:,}"
    )

    print(
        f"[{name}] Expected edges: "
        f"{summary['STRING_expected_edges']}"
    )

    print(
        f"[{name}] PPI p-value  : "
        f"{summary['STRING_PPI_enrichment_p']}"
    )

    return {
        "name":
            name,

        "input":
            input_df,

        "mapping":
            mapping,

        "unmapped":
            unmapped,

        "graph":
            G,

        "nodes":
            nodes,

        "edges":
            edges,

        "components":
            components,

        "ppi":
            ppi,

        "summary":
            summary,
    }


# ============================================================
# Validate inputs
# ============================================================

print()
print("=" * 76)
print("PHASE F3 — DUAL STRING NETWORK ANALYSIS")
print("=" * 76)

for path in [
    CORE_FILE,
    EXTENDED_FILE,
]:

    if not path.exists():

        raise FileNotFoundError(
            f"Required F2 file not found:\n"
            f"{path}"
        )


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Load F2
# ============================================================

core_df = pd.read_csv(
    CORE_FILE
)

extended_df = pd.read_csv(
    EXTENDED_FILE
)

for name, df in [
    (
        "CORE",
        core_df,
    ),
    (
        "EXTENDED",
        extended_df,
    ),
]:

    if "Gene" not in df.columns:

        raise ValueError(
            f"{name} input does not contain a Gene column."
        )

    df[
        "Gene"
    ] = (
        df[
            "Gene"
        ]
        .map(
            normalize_gene
        )
    )


# ============================================================
# Run CORE and EXTENDED
# ============================================================

core_result = analyze_network(
    name="CORE",
    input_df=core_df,
    output_dir=CORE_OUTPUT,
)

extended_result = analyze_network(
    name="EXTENDED",
    input_df=extended_df,
    output_dir=EXTENDED_OUTPUT,
)


# ============================================================
# Compare networks
# ============================================================

core_summary = (
    core_result[
        "summary"
    ]
)

extended_summary = (
    extended_result[
        "summary"
    ]
)

core_graph = (
    core_result[
        "graph"
    ]
)

extended_graph = (
    extended_result[
        "graph"
    ]
)


# ------------------------------------------------------------
# How many CORE components become connected after adding
# INTERFACE genes?
# ------------------------------------------------------------

mapped_core_genes = set(
    core_graph.nodes()
)

extended_component_lookup = {}

for idx, component in enumerate(
    nx.connected_components(
        extended_graph
    ),
    start=1,
):

    for gene in component:

        extended_component_lookup[
            gene
        ] = idx


extended_components_containing_core = {
    extended_component_lookup[
        gene
    ]
    for gene
    in mapped_core_genes
    if gene
    in extended_component_lookup
}


core_group_count_in_extended = len(
    extended_components_containing_core
)


# Largest number of CORE genes that become mutually connected
# through the extended network.
core_counts_by_extended_component = {}

for gene in mapped_core_genes:

    if gene not in extended_component_lookup:
        continue

    component = (
        extended_component_lookup[
            gene
        ]
    )

    core_counts_by_extended_component[
        component
    ] = (
        core_counts_by_extended_component.get(
            component,
            0,
        )
        + 1
    )


largest_core_group_in_extended = (
    max(
        core_counts_by_extended_component.values()
    )
    if core_counts_by_extended_component
    else 0
)


comparison = pd.DataFrame(
    [
        core_summary,
        extended_summary,
    ]
)


comparison[
    "Core_connectivity_groups_within_extended"
] = None

comparison[
    "Largest_core_group_within_extended"
] = None


comparison.loc[
    comparison[
        "Network"
    ]
    == "EXTENDED",
    "Core_connectivity_groups_within_extended",
] = (
    core_group_count_in_extended
)


comparison.loc[
    comparison[
        "Network"
    ]
    == "EXTENDED",
    "Largest_core_group_within_extended",
] = (
    largest_core_group_in_extended
)


comparison.to_csv(
    COMPARISON_CSV,
    index=False,
)


# ============================================================
# Difference summary
# ============================================================

difference_summary = pd.DataFrame(
    [
        {
            "Metric":
                "Input genes added",

            "CORE":
                core_summary[
                    "Input_genes"
                ],

            "EXTENDED":
                extended_summary[
                    "Input_genes"
                ],

            "Difference":
                (
                    extended_summary[
                        "Input_genes"
                    ]
                    -
                    core_summary[
                        "Input_genes"
                    ]
                ),
        },

        {
            "Metric":
                "Mapped STRING proteins",

            "CORE":
                core_summary[
                    "STRING_mapped_genes"
                ],

            "EXTENDED":
                extended_summary[
                    "STRING_mapped_genes"
                ],

            "Difference":
                (
                    extended_summary[
                        "STRING_mapped_genes"
                    ]
                    -
                    core_summary[
                        "STRING_mapped_genes"
                    ]
                ),
        },

        {
            "Metric":
                "Edges",

            "CORE":
                core_summary[
                    "Edges"
                ],

            "EXTENDED":
                extended_summary[
                    "Edges"
                ],

            "Difference":
                (
                    extended_summary[
                        "Edges"
                    ]
                    -
                    core_summary[
                        "Edges"
                    ]
                ),
        },

        {
            "Metric":
                "Connected components",

            "CORE":
                core_summary[
                    "Connected_components"
                ],

            "EXTENDED":
                extended_summary[
                    "Connected_components"
                ],

            "Difference":
                (
                    extended_summary[
                        "Connected_components"
                    ]
                    -
                    core_summary[
                        "Connected_components"
                    ]
                ),
        },

        {
            "Metric":
                "Isolates",

            "CORE":
                core_summary[
                    "Isolates"
                ],

            "EXTENDED":
                extended_summary[
                    "Isolates"
                ],

            "Difference":
                (
                    extended_summary[
                        "Isolates"
                    ]
                    -
                    core_summary[
                        "Isolates"
                    ]
                ),
        },

        {
            "Metric":
                "Largest component",

            "CORE":
                core_summary[
                    "Largest_component"
                ],

            "EXTENDED":
                extended_summary[
                    "Largest_component"
                ],

            "Difference":
                (
                    extended_summary[
                        "Largest_component"
                    ]
                    -
                    core_summary[
                        "Largest_component"
                    ]
                ),
        },

        {
            "Metric":
                "CORE connectivity groups after allowing INTERFACE paths",

            "CORE":
                core_summary[
                    "Connected_components"
                ],

            "EXTENDED":
                core_group_count_in_extended,

            "Difference":
                (
                    core_group_count_in_extended
                    -
                    core_summary[
                        "Connected_components"
                    ]
                ),
        },

        {
            "Metric":
                "Largest co-connected CORE gene group",

            "CORE":
                core_summary[
                    "Largest_component"
                ],

            "EXTENDED":
                largest_core_group_in_extended,

            "Difference":
                (
                    largest_core_group_in_extended
                    -
                    core_summary[
                        "Largest_component"
                    ]
                ),
        },
    ]
)


# ============================================================
# Excel summary workbook
# ============================================================

with pd.ExcelWriter(
    SUMMARY_XLSX,
    engine="openpyxl",
) as writer:

    comparison.to_excel(
        writer,
        sheet_name="Network summary",
        index=False,
    )

    difference_summary.to_excel(
        writer,
        sheet_name="CORE vs EXTENDED",
        index=False,
    )

    core_result[
        "nodes"
    ].to_excel(
        writer,
        sheet_name="CORE nodes",
        index=False,
    )

    core_result[
        "edges"
    ].to_excel(
        writer,
        sheet_name="CORE edges",
        index=False,
    )

    core_result[
        "components"
    ].to_excel(
        writer,
        sheet_name="CORE components",
        index=False,
    )

    extended_result[
        "nodes"
    ].to_excel(
        writer,
        sheet_name="EXTENDED nodes",
        index=False,
    )

    extended_result[
        "edges"
    ].to_excel(
        writer,
        sheet_name="EXTENDED edges",
        index=False,
    )

    extended_result[
        "components"
    ].to_excel(
        writer,
        sheet_name="EXTENDED components",
        index=False,
    )

    core_result[
        "mapping"
    ].to_excel(
        writer,
        sheet_name="CORE mapping",
        index=False,
    )

    extended_result[
        "mapping"
    ].to_excel(
        writer,
        sheet_name="EXT mapping",
        index=False,
    )


format_workbook(
    SUMMARY_XLSX
)


# ============================================================
# Metadata
# ============================================================

metadata = {
    "phase":
        "F3_STRING_network",

    "string_version":
        STRING_VERSION,

    "string_api_base":
        STRING_API_BASE,

    "species":
        {
            "name":
                "Mus musculus",

            "ncbi_taxon":
                SPECIES,
        },

    "network_type":
        NETWORK_TYPE,

    "required_score":
        REQUIRED_SCORE,

    "add_nodes":
        ADD_NODES,

    "interpretation_note":
        (
            "STRING functional associations are undirected prior-knowledge "
            "associations and do not establish causal direction."
        ),

    "core":
        core_summary,

    "extended":
        extended_summary,

    "core_connectivity_through_extended": {
        "standalone_core_components":
            core_summary[
                "Connected_components"
            ],

        "extended_components_containing_core":
            core_group_count_in_extended,

        "largest_core_group_within_extended":
            largest_core_group_in_extended,
    },
}


with open(
    METADATA_JSON,
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        metadata,
        handle,
        indent=2,
        ensure_ascii=False,
        default=str,
    )


# ============================================================
# Console report
# ============================================================

print()
print("=" * 76)
print("PHASE F3 COMPLETE")
print("=" * 76)

print()

print(
    "CORE"
)

print(
    f"    Input genes          : "
    f"{core_summary['Input_genes']:,}"
)

print(
    f"    STRING mapped        : "
    f"{core_summary['STRING_mapped_genes']:,}"
)

print(
    f"    Edges                : "
    f"{core_summary['Edges']:,}"
)

print(
    f"    Expected edges       : "
    f"{core_summary['STRING_expected_edges']}"
)

print(
    f"    Obs/exp ratio        : "
    f"{core_summary['Observed_expected_ratio']}"
)

print(
    f"    PPI enrichment p     : "
    f"{core_summary['STRING_PPI_enrichment_p']}"
)

print(
    f"    Components           : "
    f"{core_summary['Connected_components']:,}"
)

print(
    f"    Isolates             : "
    f"{core_summary['Isolates']:,}"
)

print(
    f"    Largest component    : "
    f"{core_summary['Largest_component']:,}"
)

print()

print(
    "EXTENDED"
)

print(
    f"    Input genes          : "
    f"{extended_summary['Input_genes']:,}"
)

print(
    f"    STRING mapped        : "
    f"{extended_summary['STRING_mapped_genes']:,}"
)

print(
    f"    Edges                : "
    f"{extended_summary['Edges']:,}"
)

print(
    f"    Expected edges       : "
    f"{extended_summary['STRING_expected_edges']}"
)

print(
    f"    Obs/exp ratio        : "
    f"{extended_summary['Observed_expected_ratio']}"
)

print(
    f"    PPI enrichment p     : "
    f"{extended_summary['STRING_PPI_enrichment_p']}"
)

print(
    f"    Components           : "
    f"{extended_summary['Connected_components']:,}"
)

print(
    f"    Isolates             : "
    f"{extended_summary['Isolates']:,}"
)

print(
    f"    Largest component    : "
    f"{extended_summary['Largest_component']:,}"
)

print()

print(
    "CORE connectivity after adding INTERFACE genes"
)

print(
    f"    CORE standalone components          : "
    f"{core_summary['Connected_components']:,}"
)

print(
    f"    Connectivity groups in EXTENDED     : "
    f"{core_group_count_in_extended:,}"
)

print(
    f"    Largest CORE group in EXTENDED      : "
    f"{largest_core_group_in_extended:,}"
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
    "[NEXT] Phase F4 will perform module/community analysis. "
    "First inspect whether the INTERFACE genes substantially merge "
    "or enlarge the CORE network."
)