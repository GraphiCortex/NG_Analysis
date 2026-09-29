#!/usr/bin/env python3

"""
03_network_modules.py

Phase B of the R-Noel cell-death network analysis.

INPUT
-----
Phase A outputs:

    Data/Processed/Network/string_seed_node_metrics.csv
    Data/Processed/Network/string_seed_edges.csv

QUESTIONS
---------
1. What connected components exist in the 126-gene seed network?
2. Can the larger components be subdivided into graph communities/modules?
3. Which genes connect different modules rather than simply being local hubs?
4. What experimental evidence is concentrated in each module?
5. Are modules functionally enriched relative to the already-selected
   126-gene candidate universe?

IMPORTANT
---------
No outside STRING neighbors are added here.

A graph module is NOT automatically a biological pathway.
Biological interpretation comes only after inspecting:
    - its genes
    - its enrichment
    - adult/embryo recurrence
    - developmental reversal
    - Andrusiak overlap
    - mechanism overlap

OUTPUT
------
Data/Processed/Network/Modules/
    module_assignments.csv
    module_summary.csv
    bridge_metrics.csv
    module_edges.csv
    module_enrichment.csv
    module_diagnostics.json

Data/Raw/Network/Modules/
    gprofiler module-query caches

Figures/Network/Modules/
    module_sizes.png
    bridge_participation_vs_betweenness.png
    Mxx_enrichment.png
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import requests


# ============================================================
# Configuration
# ============================================================

GPROFILER_BASE = "https://biit.cs.ut.ee/gprofiler/api"

GPROFILER_ORGANISM = "mmusculus"

GPROFILER_SOURCES = [
    "GO:BP",
    "KEGG",
    "REAC",
    "WP",
]

DEFAULT_ALPHA = 0.05
DEFAULT_CORRECTION = "fdr"

# Modules smaller than this are not enriched:
DEFAULT_MIN_ENRICHMENT_SIZE = 4

# Components smaller than this remain intact.
# Larger ones are eligible for community subdivision.
DEFAULT_SPLIT_COMPONENT_SIZE = 6

DEFAULT_TOP_TERMS = 10

HTTP_TIMEOUT = 120
HTTP_RETRIES = 3


# ============================================================
# General helpers
# ============================================================

def clean_gene(value):

    if pd.isna(value):
        return None

    value = str(value).strip().upper()

    if not value:
        return None

    return value


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


def read_json(path):

    return json.loads(
        path.read_text(
            encoding="utf-8"
        )
    )


def gene_hash(genes):

    text = "|".join(
        sorted(genes)
    )

    return hashlib.sha256(
        text.encode("utf-8")
    ).hexdigest()[:10]


def normalize_header(value):

    return (
        str(value)
        .strip()
        .lower()
        .replace("–", "-")
        .replace("—", "-")
        .replace("_", " ")
    )


def find_column(df, candidate_names):

    lookup = {

        normalize_header(col):
            col

        for col
        in df.columns
    }

    for name in candidate_names:

        key = normalize_header(
            name
        )

        if key in lookup:
            return lookup[key]

    return None


def is_yes(value):

    if pd.isna(value):
        return False

    if isinstance(
        value,
        bool,
    ):
        return value

    try:

        if float(value) == 1:
            return True

    except (
        TypeError,
        ValueError,
    ):
        pass

    return (
        str(value)
        .strip()
        .upper()
        in {
            "YES",
            "TRUE",
            "Y",
        }
    )


def join_genes(values):

    return "; ".join(
        sorted(
            {
                str(v)
                for v in values
                if not pd.isna(v)
                and str(v).strip()
            }
        )
    )


def request_with_retry(
    session,
    method,
    url,
    **kwargs,
):

    last_error = None

    for attempt in range(
        1,
        HTTP_RETRIES + 1,
    ):

        try:

            response = (
                session.request(
                    method,
                    url,
                    timeout=HTTP_TIMEOUT,
                    **kwargs,
                )
            )

            response.raise_for_status()

            return response

        except requests.RequestException as exc:

            last_error = exc

            if attempt < HTTP_RETRIES:

                wait = (
                    2 ** (attempt - 1)
                )

                print(
                    f"HTTP error: {exc}"
                )

                print(
                    f"Retrying in {wait}s..."
                )

                time.sleep(wait)

    raise RuntimeError(
        "HTTP request failed: "
        f"{last_error}"
    )


# ============================================================
# Input
# ============================================================

def load_phase_a(
    node_path,
    edge_path,
):

    if not node_path.exists():
        raise FileNotFoundError(
            node_path
        )

    if not edge_path.exists():
        raise FileNotFoundError(
            edge_path
        )

    nodes = pd.read_csv(
        node_path
    )

    edges = pd.read_csv(
        edge_path
    )

    required_nodes = {
        "Gene",
        "stringId",
    }

    missing = (
        required_nodes
        - set(nodes.columns)
    )

    if missing:

        raise ValueError(
            "Node table missing: "
            + ", ".join(
                sorted(missing)
            )
        )

    nodes["Gene"] = (
        nodes["Gene"]
        .map(clean_gene)
    )

    nodes.dropna(
        subset=["Gene"],
        inplace=True,
    )

    nodes.drop_duplicates(
        "Gene",
        inplace=True,
    )

    if not edges.empty:

        required_edges = {
            "stringId_A",
            "stringId_B",
            "score",
        }

        missing = (
            required_edges
            - set(edges.columns)
        )

        if missing:

            raise ValueError(
                "Edge table missing: "
                + ", ".join(
                    sorted(missing)
                )
            )

    return nodes, edges


# ============================================================
# Build graph
# ============================================================

def build_graph(
    nodes,
    edges,
):

    G = nx.Graph()

    # Add every mapped candidate,
    # including isolates.
    for _, row in nodes.iterrows():

        sid = row[
            "stringId"
        ]

        if pd.isna(sid):
            continue

        G.add_node(

            str(sid),

            Gene=row["Gene"],
        )

    for _, row in edges.iterrows():

        a = str(
            row["stringId_A"]
        )

        b = str(
            row["stringId_B"]
        )

        if a == b:
            continue

        score = float(
            row["score"]
        )

        G.add_edge(

            a,

            b,

            score=score,

            # STRING score is strength.
            # Path algorithms need distance.
            distance=(
                1.0
                / max(
                    score,
                    1e-12,
                )
            ),
        )

    return G


# ============================================================
# Stable component/module assignments
# ============================================================

def gene_name(
    G,
    node,
):

    return (
        G.nodes[node]
        .get(
            "Gene",
            str(node),
        )
    )


def sorted_components(G):

    return sorted(

        nx.connected_components(G),

        key=lambda comp: (

            -len(comp),

            min(
                gene_name(
                    G,
                    node,
                )
                for node
                in comp
            ),
        ),
    )


def assign_components(G):

    assignment = {}

    for component_id, comp in enumerate(
        sorted_components(G),
        start=1,
    ):

        for node in comp:

            assignment[node] = (
                component_id
            )

    return assignment


def assign_modules(
    G,
    split_component_size,
):

    """

    Components < split threshold:
        one component = one module

    Components >= split threshold:
        weighted greedy-modularity
        subdivision

    Isolates:
        module = ISOLATE

    """

    node_to_module = {}

    module_number = 1

    for comp in sorted_components(G):

        SG = G.subgraph(
            comp
        ).copy()

        n = SG.number_of_nodes()

        # -----------------------
        # Isolated gene
        # -----------------------

        if (
            n == 1
            and SG.number_of_edges()
            == 0
        ):

            node = next(
                iter(comp)
            )

            node_to_module[
                node
            ] = "ISOLATE"

            continue

        # -----------------------
        # Larger component
        # -----------------------

        if (
            n
            >= split_component_size
            and SG.number_of_edges()
            > 0
        ):

            communities = list(

                nx.community
                .greedy_modularity_communities(
                    SG,
                    weight="score",
                )
            )

            communities = sorted(

                communities,

                key=lambda community: (

                    -len(community),

                    min(
                        gene_name(
                            G,
                            node,
                        )
                        for node
                        in community
                    ),
                ),
            )

            if len(
                communities
            ) <= 1:

                communities = [
                    set(comp)
                ]

        else:

            communities = [
                set(comp)
            ]

        # -----------------------
        # Stable module labels
        # -----------------------

        for community in communities:

            module_id = (
                f"M{module_number:02d}"
            )

            for node in community:

                node_to_module[
                    node
                ] = module_id

            module_number += 1

    return node_to_module


# ============================================================
# Bridge metrics
# ============================================================

def participation_coefficient(
    G,
    node,
    module_of,
    weighted=False,
):

    neighbors = list(
        G.neighbors(node)
    )

    if not neighbors:
        return 0.0

    totals = {}

    total_degree = 0.0

    for neighbor in neighbors:

        module = module_of[
            neighbor
        ]

        if weighted:

            value = float(

                G[node][neighbor]
                .get(
                    "score",
                    1.0,
                )
            )

        else:

            value = 1.0

        totals[module] = (
            totals.get(
                module,
                0.0,
            )
            + value
        )

        total_degree += value

    if total_degree == 0:
        return 0.0

    return (

        1.0

        - sum(

            (
                value
                / total_degree
            ) ** 2

            for value
            in totals.values()
        )
    )


def within_module_zscores(
    G,
    module_of,
):

    within_degree = {}

    for node in G.nodes():

        module = module_of[
            node
        ]

        within_degree[
            node
        ] = sum(

            1

            for neighbor
            in G.neighbors(node)

            if module_of[
                neighbor
            ] == module
        )

    module_values = {}

    for node, value in (
        within_degree.items()
    ):

        module = module_of[
            node
        ]

        if module == "ISOLATE":
            continue

        module_values.setdefault(
            module,
            [],
        ).append(value)

    zscores = {}

    for node, value in (
        within_degree.items()
    ):

        module = module_of[
            node
        ]

        if module == "ISOLATE":

            zscores[node] = 0.0

            continue

        values = (
            module_values[module]
        )

        if len(values) < 2:

            zscores[node] = 0.0

            continue

        mean = statistics.mean(
            values
        )

        sd = statistics.pstdev(
            values
        )

        if sd == 0:

            zscores[node] = 0.0

        else:

            zscores[node] = (
                value - mean
            ) / sd

    return zscores


def compute_bridge_metrics(
    G,
    nodes,
    module_of,
    component_of,
):

    zscores = (
        within_module_zscores(
            G,
            module_of,
        )
    )

    # Phase A metrics by STRING ID.
    phase_a = (
        nodes
        .copy()
    )

    phase_a[
        "stringId"
    ] = phase_a[
        "stringId"
    ].astype(str)

    phase_a = (

        phase_a
        .drop_duplicates(
            "stringId"
        )
        .set_index(
            "stringId"
        )
    )

    rows = []

    for node in G.nodes():

        module = module_of[
            node
        ]

        degree = G.degree(
            node
        )

        within_degree = sum(

            1

            for neighbor
            in G.neighbors(node)

            if module_of[
                neighbor
            ] == module
        )

        cross_degree = (
            degree
            - within_degree
        )

        row = {

            "Gene":
                gene_name(
                    G,
                    node,
                ),

            "stringId":
                node,

            "component_id":
                component_of[node],

            "module_id":
                module,

            "degree":
                degree,

            "within_module_degree":
                within_degree,

            "cross_module_degree":
                cross_degree,

            "cross_module_fraction":
                (
                    cross_degree
                    / degree
                    if degree
                    else 0.0
                ),

            "within_module_z":
                zscores[node],

            "participation_coefficient":
                participation_coefficient(
                    G,
                    node,
                    module_of,
                    weighted=False,
                ),

            "weighted_participation_coefficient":
                participation_coefficient(
                    G,
                    node,
                    module_of,
                    weighted=True,
                ),
        }

        if node in phase_a.index:

            source = phase_a.loc[
                node
            ]

            for col in [

                "betweenness_weighted",

                "betweenness_unweighted",

                "pagerank",

                "weighted_degree",

                "degree_centrality",

            ]:

                if col in source.index:

                    row[col] = (
                        source[col]
                    )

        rows.append(row)

    return pd.DataFrame(
        rows
    )


# ============================================================
# g:Profiler module enrichment
# ============================================================

def run_module_enrichment(
    session,
    module_id,
    genes,
    background,
    raw_dir,
    *,
    alpha,
    correction,
    refresh,
):

    tag = (
        f"{module_id}_"
        f"{gene_hash(genes)}"
    )

    request_path = (
        raw_dir
        / (
            f"gprofiler_"
            f"{tag}.request.json"
        )
    )

    response_path = (
        raw_dir
        / (
            f"gprofiler_"
            f"{tag}.response.json"
        )
    )

    payload = {

        "organism":
            GPROFILER_ORGANISM,

        "query":
            genes,

        # Important:
        # compare each module against
        # our 126-gene candidate set.
        "background":
            background,

        "domain_scope":
            "custom_annotated",

        "sources":
            GPROFILER_SOURCES,

        "user_threshold":
            alpha,

        "significance_threshold_method":
            correction,

        "ordered":
            False,

        "all_results":
            False,

        "measure_underrepresentation":
            False,

        "no_iea":
            False,

        "no_evidences":
            False,
    }

    write_json(
        request_path,
        payload,
    )

    if (
        response_path.exists()
        and not refresh
    ):

        raw = read_json(
            response_path
        )

    else:

        response = request_with_retry(

            session,

            "POST",

            (
                f"{GPROFILER_BASE}"
                "/gost/profile/"
            ),

            json=payload,
        )

        raw = response.json()

        write_json(
            response_path,
            raw,
        )

        time.sleep(0.2)

    results = raw.get(
        "result",
        [],
    )

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(
        results
    )

    df.rename(
        columns={
            "p_value":
                "adjusted_p_value"
        },
        inplace=True,
    )

    df.insert(
        0,
        "module_id",
        module_id,
    )

    df.insert(
        1,
        "module_gene_count",
        len(genes),
    )

    # Make intersections easier to read.
    if "intersections" in df.columns:

        def flatten_intersections(x):

            if not isinstance(x, list):
                return x

            flattened = []

            for item in x:

                if isinstance(item, list):
                    flattened.extend(
                        str(gene)
                        for gene in item
                    )

                else:
                    flattened.append(
                        str(item)
                    )

            return ";".join(flattened)

        df["intersections"] = (
            df["intersections"]
            .map(flatten_intersections)
        )
    return df


# ============================================================
# Experimental overlay
# ============================================================

def make_module_summary(
    module_id,
    module_nodes,
    bridge,
    G,
    enrichment,
):

    node_ids = (

        module_nodes[
            "stringId"
        ]
        .astype(str)
        .tolist()
    )

    SG = G.subgraph(
        node_ids
    )

    shared_col = find_column(
        module_nodes,
        [
            "Adult–embryo shared?",
            "Adult-embryo shared?",
        ],
    )

    direction_col = find_column(
        module_nodes,
        [
            "Direction across stages"
        ],
    )

    andrusiak_col = find_column(
        module_nodes,
        [
            "Andrusiak overlap"
        ],
    )

    mechanism_count_col = (
        find_column(
            module_nodes,
            [
                "Mechanism count"
            ],
        )
    )

    top_distinct_col = (
        find_column(
            module_nodes,
            [
                "Top-distinct manual highlight",
                "Top-distinct",
            ],
        )
    )

    yellow_col = (
        find_column(
            module_nodes,
            [
                "Yellow pathway highlight",
                "Doctor highlight",
            ],
        )
    )

    if shared_col:

        shared_mask = (
            module_nodes[
                shared_col
            ].map(is_yes)
        )

    else:

        shared_mask = pd.Series(
            False,
            index=module_nodes.index,
        )

    if direction_col:

        reversed_mask = (

            module_nodes[
                direction_col
            ]
            .astype(str)
            .str.upper()
            .eq("REVERSED")
        )

    else:

        reversed_mask = pd.Series(
            False,
            index=module_nodes.index,
        )

    if andrusiak_col:

        andrusiak_mask = (
            module_nodes[
                andrusiak_col
            ].map(is_yes)
        )

    else:

        andrusiak_mask = pd.Series(
            False,
            index=module_nodes.index,
        )

    if mechanism_count_col:

        multi_mask = (

            pd.to_numeric(
                module_nodes[
                    mechanism_count_col
                ],
                errors="coerce",
            )

            .fillna(0)

            > 1
        )

    else:

        multi_mask = pd.Series(
            False,
            index=module_nodes.index,
        )

    top_mask = (
        module_nodes[
            top_distinct_col
        ].map(is_yes)

        if top_distinct_col

        else pd.Series(
            False,
            index=module_nodes.index,
        )
    )

    yellow_mask = (
        module_nodes[
            yellow_col
        ].map(is_yes)

        if yellow_col

        else pd.Series(
            False,
            index=module_nodes.index,
        )
    )

    # Top graph bridges within this module.
    bridge_sub = (

        bridge[
            bridge["module_id"]
            == module_id
        ]

        .sort_values(

            [
                "participation_coefficient",
                "betweenness_weighted",
                "degree",
            ],

            ascending=[
                False,
                False,
                False,
            ],
        )
    )

    top_bridges = (
        bridge_sub[
            "Gene"
        ]
        .head(5)
        .tolist()
    )

    # Top functional enrichments.
    if enrichment.empty:

        top_terms = []

    else:

        top_terms = (

            enrichment
            .sort_values(
                "adjusted_p_value"
            )
            .head(5)

            .apply(

                lambda row:
                    (
                        f"{row['source']} | "
                        f"{row['name']}"
                    ),

                axis=1,
            )

            .tolist()
        )

    return {

        "module_id":
            module_id,

        "component_id":
            int(
                bridge_sub[
                    "component_id"
                ].iloc[0]
            ),

        "size":
            len(
                module_nodes
            ),

        "internal_edges":
            SG.number_of_edges(),

        "density":
            nx.density(SG)
            if len(module_nodes) > 1
            else 0.0,

        "genes":
            join_genes(
                module_nodes[
                    "Gene"
                ]
            ),

        "adult_embryo_shared_count":
            int(
                shared_mask.sum()
            ),

        "adult_embryo_shared_genes":
            join_genes(
                module_nodes.loc[
                    shared_mask,
                    "Gene",
                ]
            ),

        "reversed_count":
            int(
                reversed_mask.sum()
            ),

        "reversed_genes":
            join_genes(
                module_nodes.loc[
                    reversed_mask,
                    "Gene",
                ]
            ),

        "andrusiak_count":
            int(
                andrusiak_mask.sum()
            ),

        "andrusiak_genes":
            join_genes(
                module_nodes.loc[
                    andrusiak_mask,
                    "Gene",
                ]
            ),

        "multi_mechanism_count":
            int(
                multi_mask.sum()
            ),

        "multi_mechanism_genes":
            join_genes(
                module_nodes.loc[
                    multi_mask,
                    "Gene",
                ]
            ),

        "top_distinct_count":
            int(
                top_mask.sum()
            ),

        "top_distinct_genes":
            join_genes(
                module_nodes.loc[
                    top_mask,
                    "Gene",
                ]
            ),

        "doctor_highlight_count":
            int(
                yellow_mask.sum()
            ),

        "doctor_highlight_genes":
            join_genes(
                module_nodes.loc[
                    yellow_mask,
                    "Gene",
                ]
            ),

        "top_bridge_genes":
            join_genes(
                top_bridges
            ),

        "significant_enrichment_terms":
            len(enrichment),

        "top_enriched_terms":
            "; ".join(
                top_terms
            ),
    }


# ============================================================
# Figures
# ============================================================

def plot_module_sizes(
    summary,
    figures_dir,
):

    if summary.empty:
        return

    df = (
        summary
        .sort_values(
            "size",
            ascending=False,
        )
    )

    plt.figure(
        figsize=(
            max(
                8,
                0.55 * len(df),
            ),
            5,
        )
    )

    plt.bar(
        df["module_id"],
        df["size"],
    )

    plt.xlabel(
        "Module"
    )

    plt.ylabel(
        "Number of genes"
    )

    plt.title(
        "Seed-network module sizes"
    )

    plt.xticks(
        rotation=60
    )

    plt.tight_layout()

    plt.savefig(
        figures_dir
        / "module_sizes.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def plot_bridge_scatter(
    bridge,
    figures_dir,
):

    df = bridge[
        bridge["degree"] > 0
    ].copy()

    if df.empty:
        return

    plt.figure(
        figsize=(8, 6)
    )

    plt.scatter(

        df[
            "participation_coefficient"
        ],

        df[
            "betweenness_weighted"
        ].fillna(0),
    )

    plt.xlabel(
        "Participation coefficient"
    )

    plt.ylabel(
        "Weighted betweenness"
    )

    plt.title(
        "Inter-module participation vs betweenness"
    )

    labels = (

        df.sort_values(

            [
                "participation_coefficient",
                "betweenness_weighted",
            ],

            ascending=[
                False,
                False,
            ],
        )

        .head(12)
    )

    for _, row in labels.iterrows():

        plt.annotate(

            row["Gene"],

            (

                row[
                    "participation_coefficient"
                ],

                row[
                    "betweenness_weighted"
                ],
            ),

            fontsize=8,
        )

    plt.tight_layout()

    plt.savefig(
        figures_dir
        / (
            "bridge_participation_"
            "vs_betweenness.png"
        ),
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def plot_module_enrichment(
    enrichment,
    module_id,
    figures_dir,
    top_n,
):

    if enrichment.empty:
        return

    df = (

        enrichment

        .sort_values(
            "adjusted_p_value"
        )

        .head(top_n)

        .copy()
    )

    df[
        "minus_log10_p"
    ] = (

        df[
            "adjusted_p_value"
        ].map(

            lambda p:
                -math.log10(
                    max(
                        float(p),
                        sys.float_info.min,
                    )
                )
        )
    )

    labels = (

        df.apply(

            lambda row:
                (
                    f"{row['source']} | "
                    f"{row['name']}"
                ),

            axis=1,
        )

        .map(

            lambda x:
                x
                if len(x) <= 80
                else x[:77] + "..."
        )
    )

    height = max(
        4,
        0.45 * len(df)
        + 1.5,
    )

    plt.figure(
        figsize=(
            10,
            height,
        )
    )

    y = range(
        len(df)
    )

    plt.barh(
        y,
        df["minus_log10_p"],
    )

    plt.yticks(
        y,
        labels,
    )

    plt.gca().invert_yaxis()

    plt.xlabel(
        "-log10(adjusted p-value)"
    )

    plt.title(
        (
            f"{module_id}: "
            "enrichment vs master "
            "candidate universe"
        )
    )

    plt.tight_layout()

    plt.savefig(

        figures_dir
        / f"{module_id}_enrichment.png",

        dpi=300,

        bbox_inches="tight",
    )

    plt.close()


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(

        "--min-enrichment-size",

        type=int,

        default=
            DEFAULT_MIN_ENRICHMENT_SIZE,
    )

    parser.add_argument(

        "--split-component-size",

        type=int,

        default=
            DEFAULT_SPLIT_COMPONENT_SIZE,
    )

    parser.add_argument(

        "--alpha",

        type=float,

        default=
            DEFAULT_ALPHA,
    )

    parser.add_argument(

        "--correction",

        choices=[
            "fdr",
            "g_SCS",
            "bonferroni",
        ],

        default=
            DEFAULT_CORRECTION,
    )

    parser.add_argument(

        "--top-terms",

        type=int,

        default=
            DEFAULT_TOP_TERMS,
    )

    parser.add_argument(

        "--refresh",

        action="store_true",
    )

    args = parser.parse_args()

    # Project/
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

    node_path = (

        processed_network
        / "string_seed_node_metrics.csv"
    )

    edge_path = (

        processed_network
        / "string_seed_edges.csv"
    )

    processed_modules = (

        processed_network
        / "Modules"
    )

    raw_modules = (

        project_root
        / "Data"
        / "Raw"
        / "Network"
        / "Modules"
    )

    figures_modules = (

        project_root
        / "Figures"
        / "Network"
        / "Modules"
    )

    for directory in [

        processed_modules,
        raw_modules,
        figures_modules,

    ]:

        directory.mkdir(
            parents=True,
            exist_ok=True,
        )

    print()
    print(
        "R-Noel network analysis - Phase B"
    )
    print(
        "================================="
    )

    nodes, edges = (
        load_phase_a(
            node_path,
            edge_path,
        )
    )

    G = build_graph(
        nodes,
        edges,
    )

    print(
        f"Seed genes: "
        f"{G.number_of_nodes()}"
    )

    print(
        f"Seed edges: "
        f"{G.number_of_edges()}"
    )

    # ========================================================
    # Components + modules
    # ========================================================

    component_of = (
        assign_components(G)
    )

    module_of = (
        assign_modules(

            G,

            args.split_component_size,
        )
    )

    modules = sorted(

        {
            module
            for module
            in module_of.values()
            if module != "ISOLATE"
        }
    )

    isolates = sum(

        module == "ISOLATE"

        for module
        in module_of.values()
    )

    print(
        f"Non-isolate modules: "
        f"{len(modules)}"
    )

    print(
        f"Isolated genes: "
        f"{isolates}"
    )

    # ========================================================
    # Bridge metrics
    # ========================================================

    bridge = (
        compute_bridge_metrics(

            G,
            nodes,
            module_of,
            component_of,
        )
    )

    bridge.sort_values(

        [
            "participation_coefficient",
            "betweenness_weighted",
            "degree",
        ],

        ascending=[
            False,
            False,
            False,
        ],

        inplace=True,
    )

    bridge.to_csv(

        processed_modules
        / "bridge_metrics.csv",

        index=False,
    )

    # ========================================================
    # Node/module assignments
    # ========================================================

    assignments = (
        nodes.copy()
    )

    assignments[
        "stringId"
    ] = assignments[
        "stringId"
    ].astype(str)

    assignments[
        "component_id_phaseB"
    ] = (

        assignments[
            "stringId"
        ].map(
            component_of
        )
    )

    assignments[
        "module_id"
    ] = (

        assignments[
            "stringId"
        ].map(
            module_of
        )
    )

    bridge_columns = [

        "Gene",

        "within_module_degree",

        "cross_module_degree",

        "cross_module_fraction",

        "within_module_z",

        "participation_coefficient",

        "weighted_participation_coefficient",
    ]

    assignments = (

        assignments.merge(

            bridge[
                bridge_columns
            ],

            on="Gene",

            how="left",
        )
    )

    assignments.to_csv(

        processed_modules
        / "module_assignments.csv",

        index=False,
    )

    # ========================================================
    # Label edges by module
    # ========================================================

    module_edges = (
        edges.copy()
    )

    module_edges[
        "stringId_A"
    ] = module_edges[
        "stringId_A"
    ].astype(str)

    module_edges[
        "stringId_B"
    ] = module_edges[
        "stringId_B"
    ].astype(str)

    module_edges[
        "module_A"
    ] = (

        module_edges[
            "stringId_A"
        ].map(
            module_of
        )
    )

    module_edges[
        "module_B"
    ] = (

        module_edges[
            "stringId_B"
        ].map(
            module_of
        )
    )

    module_edges[
        "cross_module"
    ] = (

        module_edges[
            "module_A"
        ]

        !=

        module_edges[
            "module_B"
        ]
    )

    module_edges.to_csv(

        processed_modules
        / "module_edges.csv",

        index=False,
    )

    # ========================================================
    # Gene sets by module
    # ========================================================

    gene_sets = {}

    for module_id in modules:

        genes = (

            assignments[

                assignments[
                    "module_id"
                ]
                == module_id

            ]["Gene"]

            .dropna()

            .astype(str)

            .tolist()
        )

        gene_sets[
            module_id
        ] = sorted(genes)

    write_json(

        processed_modules
        / "module_gene_sets.json",

        gene_sets,
    )

    # ========================================================
    # Module enrichment
    # ========================================================

    session = (
        requests.Session()
    )

    session.headers.update(
        {
            "User-Agent":
                (
                    "R-Noel_BIOL295_"
                    "module_pipeline/1.0"
                )
        }
    )

    background = (

        nodes[
            "Gene"
        ]

        .dropna()

        .astype(str)

        .tolist()
    )

    enrichment_by_module = {}

    enrichment_frames = []

    print()
    print(
        "Module enrichment"
    )
    print(
        "Background = 126-gene "
        "master candidate universe"
    )
    print(
        "--------------------------------"
    )

    for module_id in modules:

        genes = gene_sets[
            module_id
        ]

        if (
            len(genes)
            < args.min_enrichment_size
        ):

            print(

                f"{module_id}: "
                f"size={len(genes)} "
                "-> skipped"
            )

            enrichment_by_module[
                module_id
            ] = pd.DataFrame()

            continue

        print(
            f"{module_id}: "
            f"size={len(genes)}"
        )

        enrichment = (
            run_module_enrichment(

                session,

                module_id,

                genes,

                background,

                raw_modules,

                alpha=args.alpha,

                correction=
                    args.correction,

                refresh=
                    args.refresh,
            )
        )

        enrichment_by_module[
            module_id
        ] = enrichment

        if not enrichment.empty:

            enrichment_frames.append(
                enrichment
            )

            plot_module_enrichment(

                enrichment,

                module_id,

                figures_modules,

                args.top_terms,
            )

    if enrichment_frames:

        enrichment_all = pd.concat(

            enrichment_frames,

            ignore_index=True,
        )

        enrichment_all.sort_values(

            [
                "module_id",
                "adjusted_p_value",
            ],

            inplace=True,
        )

    else:

        enrichment_all = (
            pd.DataFrame()
        )

    enrichment_all.to_csv(

        processed_modules
        / "module_enrichment.csv",

        index=False,
    )

    # ========================================================
    # Module summary
    # ========================================================

    summary_rows = []

    for module_id in modules:

        subset = (

            assignments[

                assignments[
                    "module_id"
                ]
                == module_id

            ].copy()
        )

        summary_rows.append(

            make_module_summary(

                module_id,

                subset,

                bridge,

                G,

                enrichment_by_module[
                    module_id
                ],
            )
        )

    summary = pd.DataFrame(
        summary_rows
    )

    summary.sort_values(

        [
            "size",
            "adult_embryo_shared_count",
        ],

        ascending=[
            False,
            False,
        ],

        inplace=True,
    )

    summary.to_csv(

        processed_modules
        / "module_summary.csv",

        index=False,
    )

    diagnostics = {

        "seed_nodes":
            G.number_of_nodes(),

        "seed_edges":
            G.number_of_edges(),

        "connected_components":
            nx.number_connected_components(
                G
            ),

        "isolated_genes":
            isolates,

        "non_isolate_modules":
            len(modules),

        "split_component_size":
            args.split_component_size,

        "minimum_enrichment_size":
            args.min_enrichment_size,

        "enrichment_background":
            (
                "126-gene master "
                "candidate universe"
            ),

        "module_sizes":
            {
                module:
                    len(
                        gene_sets[module]
                    )
                for module
                in modules
            },
    }

    write_json(

        processed_modules
        / "module_diagnostics.json",

        diagnostics,
    )

    # ========================================================
    # Figures
    # ========================================================

    plot_module_sizes(
        summary,
        figures_modules,
    )

    plot_bridge_scatter(
        bridge,
        figures_modules,
    )

    # ========================================================
    # Console summary
    # ========================================================

    print()
    print(
        "Largest modules"
    )
    print(
        "---------------"
    )

    print(

        summary[

            [
                "module_id",
                "component_id",
                "size",
                "adult_embryo_shared_count",
                "reversed_count",
                "top_bridge_genes",
                "top_enriched_terms",
            ]

        ]

        .head(12)

        .to_string(
            index=False
        )
    )

    print()
    print(
        "Highest inter-module participation"
    )
    print(
        "----------------------------------"
    )

    display = (

        bridge[
            bridge["degree"] > 0
        ]

        .head(15)
    )

    print(

        display[

            [
                "Gene",
                "module_id",
                "degree",
                "cross_module_degree",
                "participation_coefficient",
                "betweenness_weighted",
            ]

        ]

        .to_string(
            index=False
        )
    )

    print()
    print(
        "Phase B complete."
    )

    print(
        "Processed outputs:",
        processed_modules,
    )

    print(
        "Figures:",
        figures_modules,
    )

    print()
    print(
        "Do not interpret a graph module "
        "as a biological pathway until "
        "we inspect its genes and enrichment."
    )


if __name__ == "__main__":

    main()