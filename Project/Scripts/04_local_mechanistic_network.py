#!/usr/bin/env python3

"""
04_local_mechanistic_network.py

Phase C: directed local mechanistic network.

STRING (Phases A/B):
    undirected functional-association context

OmniPath (Phase C):
    directed / signed activity-flow interactions

Reactome:
    independent pathway-context validation

We search for short directed paths from an Rb/E2F/cell-cycle source set
to death, autophagy, and trophic/survival endpoints.

External intermediate genes are allowed, but tightly constrained.

IMPORTANT
---------
OmniPath mouse interactions are largely homology-translated from human.
They are mechanistic hypotheses, not proof that the interaction occurs in
the exact neural system studied here.

Bottleneck analysis is performed only on the retained short-path network.

RUN
---
py .\Scripts\04_local_mechanistic_network.py
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import requests


# ============================================================
# Configuration
# ============================================================

OMNIPATH_URL = "https://omnipathdb.org/interactions/"

REACTOME_ANALYSIS_URL = (
    "https://reactome.org/"
    "AnalysisService/identifiers/projection/"
)

REACTOME_VERSION_URL = (
    "https://reactome.org/"
    "ContentService/data/database/version"
)

MOUSE_TAXON = 10090


# Rb/E2F / cell-cycle entry points.
DEFAULT_SOURCES = [
    "RB1",
    "E2F3",
    "CCND1",
    "CDKN2A",
    "CDKN2D",
]


DEFAULT_TARGET_GROUPS = {

    "death": [
        "FAS",
        "FADD",
        "MLKL",
        "BCL2L11",
        "BOK",
    ],

    "autophagy": [
        "ATG7",
        "GABARAPL1",
        "BNIP3L",
    ],

    "trophic_survival": [
        "BDNF",
        "CXCL12",
        "KDR",
    ],
}


# Start conservatively.
DEFAULT_MAX_EDGES = 4
DEFAULT_MAX_EXTERNAL = 2

DEFAULT_PATHS_PER_PAIR = 25

# Safety cap before filtering paths.
DEFAULT_ENUMERATION_CAP = 5000

HTTP_TIMEOUT = 180
HTTP_RETRIES = 3


# ============================================================
# General helpers
# ============================================================

def clean_gene(value):

    if pd.isna(value):
        return None

    value = (
        str(value)
        .strip()
        .upper()
    )

    if not value:
        return None

    value = re.sub(
        r"\s*\([^)]*\)\s*$",
        "",
        value,
    )

    if any(
        ch.isspace()
        for ch in value
    ):
        return None

    return value


def unique_ordered(values):

    seen = set()
    output = []

    for value in values:

        if value not in seen:

            seen.add(value)
            output.append(value)

    return output


def parse_gene_list(
    text,
    default,
):

    if text is None:
        return list(default)

    output = []

    for item in text.split(","):

        gene = clean_gene(item)

        if gene:
            output.append(gene)

    return unique_ordered(output)


def write_json(
    path,
    obj,
):

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


def split_semicolon(value):

    if pd.isna(value):
        return set()

    value = str(value).strip()

    if not value:
        return set()

    return {

        item.strip()

        for item
        in value.split(";")

        if item.strip()
    }


def flag(value):

    if pd.isna(value):
        return 0

    if isinstance(
        value,
        bool,
    ):
        return int(value)

    try:

        return (
            1
            if float(value) == 1
            else 0
        )

    except (
        TypeError,
        ValueError,
    ):

        return (
            1
            if str(value)
            .strip()
            .lower()
            in {
                "true",
                "yes",
                "y",
            }
            else 0
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

            response = session.request(

                method,

                url,

                timeout=HTTP_TIMEOUT,

                **kwargs,
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
                    f"Retrying in "
                    f"{wait}s..."
                )

                time.sleep(wait)

    raise RuntimeError(
        "HTTP request failed: "
        f"{last_error}"
    )


# ============================================================
# Load Phase B candidates
# ============================================================

def load_candidates(path):

    if not path.exists():

        raise FileNotFoundError(
            path
        )

    df = pd.read_csv(path)

    if "Gene" not in df.columns:

        raise ValueError(
            "module_assignments.csv "
            "has no Gene column."
        )

    df["Gene"] = (
        df["Gene"]
        .map(clean_gene)
    )

    df.dropna(
        subset=["Gene"],
        inplace=True,
    )

    df.drop_duplicates(
        "Gene",
        inplace=True,
    )

    return df


# ============================================================
# OmniPath
# ============================================================

def download_omnipath(
    session,
    raw_dir,
    refresh,
):

    raw_path = (

        raw_dir
        / "omnipath_mouse_directed.tsv"
    )

    request_path = (

        raw_dir
        / "omnipath_request.json"
    )

    params = {

        "organisms":
            MOUSE_TAXON,

        "datasets":
            "omnipath",

        "genesymbols":
            1,

        "fields":
            (
                "sources,"
                "references,"
                "curation_effort,"
                "type"
            ),
    }

    write_json(
        request_path,
        params,
    )

    if (
        raw_path.exists()
        and not refresh
    ):

        print(
            "Using cached "
            "OmniPath network."
        )

    else:

        print(
            "Downloading OmniPath "
            "mouse directed network..."
        )

        response = request_with_retry(

            session,

            "GET",

            OMNIPATH_URL,

            params=params,
        )

        raw_path.write_text(

            response.text,

            encoding="utf-8",
        )

    return pd.read_csv(

        raw_path,

        sep="\t",

        low_memory=False,
    )


def symbol_column(
    df,
    side,
):

    candidates = [

        f"{side}_genesymbol",

        f"{side}_genesymb",

        f"{side}_gene_symbol",
    ]

    for col in candidates:

        if col in df.columns:
            return col

    raise ValueError(

        f"Could not identify "
        f"{side} gene-symbol "
        f"column.\n\n"
        f"Columns returned:\n"
        f"{list(df.columns)}"
    )


def row_sign(row):

    if (
        "consensus_stimulation"
        in row.index
    ):

        stimulation = flag(
            row[
                "consensus_stimulation"
            ]
        )

    else:

        stimulation = flag(
            row.get(
                "is_stimulation"
            )
        )

    if (
        "consensus_inhibition"
        in row.index
    ):

        inhibition = flag(
            row[
                "consensus_inhibition"
            ]
        )

    else:

        inhibition = flag(
            row.get(
                "is_inhibition"
            )
        )

    if (
        stimulation
        and not inhibition
    ):

        return (
            1,
            "activation",
        )

    if (
        inhibition
        and not stimulation
    ):

        return (
            -1,
            "inhibition",
        )

    if (
        stimulation
        and inhibition
    ):

        return (
            0,
            "conflicting",
        )

    return (
        0,
        "unsigned",
    )


def normalize_omnipath(raw):

    source_col = symbol_column(
        raw,
        "source",
    )

    target_col = symbol_column(
        raw,
        "target",
    )

    df = raw.copy()

    # Explicitly keep directed edges.
    if "is_directed" in df.columns:

        df = df[
            df[
                "is_directed"
            ].map(flag) == 1
        ].copy()

    # This stage is gene-level,
    # so exclude complex entities.
    if "source" in df.columns:

        df = df[

            ~df["source"]
            .astype(str)
            .str.startswith(
                "COMPLEX:"
            )

        ].copy()

    if "target" in df.columns:

        df = df[

            ~df["target"]
            .astype(str)
            .str.startswith(
                "COMPLEX:"
            )

        ].copy()

    df[
        "source_gene"
    ] = df[
        source_col
    ].map(clean_gene)

    df[
        "target_gene"
    ] = df[
        target_col
    ].map(clean_gene)

    df.dropna(

        subset=[
            "source_gene",
            "target_gene",
        ],

        inplace=True,
    )

    df = df[

        df["source_gene"]
        !=
        df["target_gene"]

    ].copy()

    signs = df.apply(
        row_sign,
        axis=1,
    )

    df["row_sign"] = [
        x[0]
        for x in signs
    ]

    df["row_sign_label"] = [
        x[1]
        for x in signs
    ]

    # Collapse multiple database records
    # for the same directed gene pair.
    rows = []

    for (
        source,
        target,
    ), group in df.groupby(

        [
            "source_gene",
            "target_gene",
        ]
    ):

        resources = set()
        references = set()
        interaction_types = set()

        if "sources" in group.columns:

            for value in group[
                "sources"
            ]:

                resources |= (
                    split_semicolon(
                        value
                    )
                )

        if "references" in group.columns:

            for value in group[
                "references"
            ]:

                references |= (
                    split_semicolon(
                        value
                    )
                )

        if "type" in group.columns:

            for value in group[
                "type"
            ]:

                interaction_types |= (
                    split_semicolon(
                        value
                    )
                )

        signs_present = (

            set(
                group[
                    "row_sign"
                ]
            )
            -
            {0}
        )

        labels_present = set(
            group[
                "row_sign_label"
            ]
        )

        if signs_present == {1}:

            sign = 1
            sign_label = (
                "activation"
            )

        elif signs_present == {-1}:

            sign = -1
            sign_label = (
                "inhibition"
            )

        elif (
            signs_present
            == {1, -1}
            or "conflicting"
            in labels_present
        ):

            sign = 0
            sign_label = (
                "conflicting"
            )

        else:

            sign = 0
            sign_label = (
                "unsigned"
            )

        if (
            "curation_effort"
            in group.columns
        ):

            curation = (
                pd.to_numeric(
                    group[
                        "curation_effort"
                    ],
                    errors="coerce",
                )
            )

            curation_effort = (

                float(
                    curation.max()
                )

                if curation
                .notna()
                .any()

                else 0.0
            )

        else:

            curation_effort = 0.0

        n_resources = len(
            resources
        )

        n_references = len(
            references
        )

        # Descriptive evidence measure.
        evidence_score = (

            1

            + math.log1p(
                n_resources
            )

            + math.log1p(
                n_references
            )

            + 0.25
            * math.log1p(
                curation_effort
            )
        )

        # Edge count remains dominant.
        evidence_cost = (

            1

            + 1
            / evidence_score
        )

        rows.append({

            "source_gene":
                source,

            "target_gene":
                target,

            "sign":
                sign,

            "sign_label":
                sign_label,

            "n_resources":
                n_resources,

            "n_references":
                n_references,

            "curation_effort":
                curation_effort,

            "resources":
                ";".join(
                    sorted(resources)
                ),

            "references":
                ";".join(
                    sorted(references)
                ),

            "interaction_types":
                ";".join(
                    sorted(
                        interaction_types
                    )
                ),

            "evidence_score":
                evidence_score,

            "evidence_cost":
                evidence_cost,

            "raw_rows_collapsed":
                len(group),
        })

    return pd.DataFrame(
        rows
    )


def build_causal_graph(edges):

    G = nx.DiGraph()

    for _, row in edges.iterrows():

        G.add_edge(

            row["source_gene"],

            row["target_gene"],

            sign=
                int(
                    row["sign"]
                ),

            sign_label=
                row[
                    "sign_label"
                ],

            n_resources=
                int(
                    row[
                        "n_resources"
                    ]
                ),

            n_references=
                int(
                    row[
                        "n_references"
                    ]
                ),

            curation_effort=
                float(
                    row[
                        "curation_effort"
                    ]
                ),

            evidence_score=
                float(
                    row[
                        "evidence_score"
                    ]
                ),

            evidence_cost=
                float(
                    row[
                        "evidence_cost"
                    ]
                ),

            resources=
                row["resources"],

            references=
                row["references"],

            interaction_types=
                row[
                    "interaction_types"
                ],
        )

    return G


# ============================================================
# Path search
# ============================================================

def target_group_lookup(
    groups,
):

    output = {}

    for group, genes in (
        groups.items()
    ):

        for gene in genes:

            output[gene] = group

    return output


def corridor(
    G,
    source,
    target,
    max_edges,
):

    forward = (
        nx.single_source_shortest_path_length(

            G,

            source,

            cutoff=max_edges,
        )
    )

    reverse = (
        nx.single_source_shortest_path_length(

            G.reverse(
                copy=False
            ),

            target,

            cutoff=max_edges,
        )
    )

    nodes = {

        node

        for node
        in set(forward)
        .intersection(reverse)

        if (
            forward[node]
            + reverse[node]
            <= max_edges
        )
    }

    return G.subgraph(
        nodes
    ).copy()


def path_sign(
    G,
    path,
):

    product = 1

    unsigned = False
    conflicting = False

    for source, target in zip(

        path[:-1],
        path[1:],

    ):

        attrs = G[
            source
        ][
            target
        ]

        sign = int(
            attrs.get(
                "sign",
                0,
            )
        )

        label = attrs.get(
            "sign_label",
            "unsigned",
        )

        if label == "conflicting":

            conflicting = True

        elif sign == 0:

            unsigned = True

        else:

            product *= sign

    if conflicting:

        return (
            0,
            "conflicting",
        )

    if unsigned:

        return (
            0,
            "unknown",
        )

    if product > 0:

        return (
            1,
            "net_activation",
        )

    return (
        -1,
        "net_inhibition",
    )


def summarize_path(
    G,
    path,
    seed_genes,
    source,
    target,
    group,
):

    intermediates = (
        path[1:-1]
    )

    external = [

        gene

        for gene
        in intermediates

        if gene
        not in seed_genes
    ]

    attrs = [

        G[a][b]

        for a, b
        in zip(
            path[:-1],
            path[1:],
        )
    ]

    sign, sign_label = (
        path_sign(
            G,
            path,
        )
    )

    evidence = [

        float(
            x[
                "evidence_score"
            ]
        )

        for x in attrs
    ]

    resources = [

        int(
            x[
                "n_resources"
            ]
        )

        for x in attrs
    ]

    references = [

        int(
            x[
                "n_references"
            ]
        )

        for x in attrs
    ]

    return {

        "source_anchor":
            source,

        "target_anchor":
            target,

        "target_group":
            group,

        "num_edges":
            len(path) - 1,

        "num_intermediate_nodes":
            len(intermediates),

        "num_external_intermediates":
            len(external),

        "external_intermediates":
            ";".join(external),

        "path_nodes":
            ";".join(path),

        "path_sign":
            sign,

        "path_sign_label":
            sign_label,

        "path_cost":
            sum(
                float(
                    x[
                        "evidence_cost"
                    ]
                )
                for x in attrs
            ),

        "min_edge_evidence_score":
            min(evidence),

        "mean_edge_evidence_score":
            (
                sum(evidence)
                / len(evidence)
            ),

        "min_edge_resources":
            min(resources),

        "min_edge_references":
            min(references),
    }


def find_pair_paths(

    G,

    source,

    target,

    group,

    seed_genes,

    max_edges,

    max_external,

    paths_per_pair,

    enumeration_cap,
):

    if (
        source not in G
        or target not in G
    ):

        return []

    try:

        shortest = (
            nx.shortest_path_length(

                G,

                source,

                target,
            )
        )

    except nx.NetworkXNoPath:

        return []

    if shortest > max_edges:

        return []

    local = corridor(

        G,

        source,

        target,

        max_edges,
    )

    generator = (

        nx.all_simple_paths(

            local,

            source=source,

            target=target,

            cutoff=max_edges,
        )
    )

    rows = []

    for path in itertools.islice(

        generator,

        enumeration_cap,
    ):

        external_count = sum(

            gene
            not in seed_genes

            for gene
            in path[1:-1]
        )

        if (
            external_count
            > max_external
        ):

            continue

        rows.append(

            summarize_path(

                G,

                path,

                seed_genes,

                source,

                target,

                group,
            )
        )

    rows.sort(

        key=lambda row: (

            row["num_edges"],

            row[
                "num_external_intermediates"
            ],

            row["path_cost"],

            -row[
                "min_edge_evidence_score"
            ],
        )
    )

    return rows[
        :paths_per_pair
    ]


def search_paths(

    G,

    sources,

    target_groups,

    seed_genes,

    max_edges,

    max_external,

    paths_per_pair,

    enumeration_cap,
):

    group_lookup = (
        target_group_lookup(
            target_groups
        )
    )

    targets = (
        unique_ordered(

            gene

            for genes
            in target_groups.values()

            for gene
            in genes
        )
    )

    rows = []

    for source in sources:

        for target in targets:

            if source == target:
                continue

            found = find_pair_paths(

                G,

                source,

                target,

                group_lookup[
                    target
                ],

                seed_genes,

                max_edges,

                max_external,

                paths_per_pair,

                enumeration_cap,
            )

            print(

                f"{source:7s}"
                f" -> "
                f"{target:10s}"
                f" : "
                f"{len(found)} path(s)"
            )

            rows.extend(found)

    if not rows:

        return pd.DataFrame()

    df = pd.DataFrame(rows)

    df.insert(

        0,

        "path_id",

        [
            f"P{i:05d}"

            for i in range(
                1,
                len(df) + 1,
            )
        ],
    )

    return df


# ============================================================
# Build retained causal network
# ============================================================

def make_path_edges(
    full_graph,
    paths,
    seed_genes,
):

    edge_paths = defaultdict(set)
    edge_pairs = defaultdict(set)

    for _, row in paths.iterrows():

        nodes = (
            row["path_nodes"]
            .split(";")
        )

        pair = (

            row[
                "source_anchor"
            ],

            row[
                "target_anchor"
            ],
        )

        for source, target in zip(

            nodes[:-1],
            nodes[1:],

        ):

            edge_paths[
                (source, target)
            ].add(
                row["path_id"]
            )

            edge_pairs[
                (source, target)
            ].add(pair)

    rows = []

    for (
        source,
        target,
    ), path_ids in (
        edge_paths.items()
    ):

        attrs = (
            full_graph[
                source
            ][
                target
            ]
        )

        rows.append({

            "source_gene":
                source,

            "target_gene":
                target,

            "source_is_seed":
                source
                in seed_genes,

            "target_is_seed":
                target
                in seed_genes,

            "sign":
                attrs[
                    "sign"
                ],

            "sign_label":
                attrs[
                    "sign_label"
                ],

            "n_resources":
                attrs[
                    "n_resources"
                ],

            "n_references":
                attrs[
                    "n_references"
                ],

            "curation_effort":
                attrs[
                    "curation_effort"
                ],

            "evidence_score":
                attrs[
                    "evidence_score"
                ],

            "resources":
                attrs[
                    "resources"
                ],

            "references":
                attrs[
                    "references"
                ],

            "interaction_types":
                attrs[
                    "interaction_types"
                ],

            "path_count":
                len(path_ids),

            "source_target_pair_count":
                len(
                    edge_pairs[
                        (
                            source,
                            target,
                        )
                    ]
                ),
        })

    return pd.DataFrame(
        rows
    )


def retained_graph(
    edge_table,
):

    G = nx.DiGraph()

    for _, row in (
        edge_table.iterrows()
    ):

        G.add_edge(

            row[
                "source_gene"
            ],

            row[
                "target_gene"
            ],

            sign=
                int(
                    row["sign"]
                ),

            sign_label=
                row[
                    "sign_label"
                ],

            evidence_score=
                float(
                    row[
                        "evidence_score"
                    ]
                ),

            resources=
                row["resources"],

            references=
                row["references"],
        )

    return G


# ============================================================
# Bottleneck analysis
# ============================================================

def reachable_pairs(
    G,
    sources,
    targets,
):

    pairs = set()

    for source in sources:

        if source not in G:
            continue

        reachable = (
            nx.descendants(
                G,
                source,
            )
        )

        for target in targets:

            if target in reachable:

                pairs.add(
                    (
                        source,
                        target,
                    )
                )

    return pairs


def bottleneck_analysis(
    G,
    sources,
    targets,
):

    baseline = reachable_pairs(

        G,

        sources,

        targets,
    )

    rows = []

    for node in G.nodes():

        H = G.copy()

        H.remove_node(
            node
        )

        remaining = (
            reachable_pairs(

                H,

                sources,

                targets,
            )
        )

        lost = (
            baseline
            - remaining
        )

        rows.append({

            "Gene":
                node,

            "baseline_reachable_pairs":
                len(baseline),

            "pairs_lost_on_removal":
                len(lost),

            "pair_loss_fraction":
                (
                    len(lost)
                    / len(baseline)

                    if baseline
                    else 0
                ),

            "lost_pairs":
                ";".join(

                    f"{a}->{b}"

                    for a, b
                    in sorted(lost)
                ),

            "is_source_anchor":
                node in sources,

            "is_target_anchor":
                node in targets,
        })

    return pd.DataFrame(
        rows
    )


# ============================================================
# Node metrics
# ============================================================

def node_metrics(
    G,
    paths,
    assignments,
    sources,
    targets,
    bottlenecks,
):

    seed_genes = set(
        assignments[
            "Gene"
        ]
    )

    occurrence = defaultdict(int)

    pair_occurrence = (
        defaultdict(set)
    )

    for _, row in paths.iterrows():

        nodes = set(

            row[
                "path_nodes"
            ].split(";")
        )

        pair = (

            row[
                "source_anchor"
            ],

            row[
                "target_anchor"
            ],
        )

        for node in nodes:

            occurrence[
                node
            ] += 1

            pair_occurrence[
                node
            ].add(pair)

    bottleneck_lookup = (

        bottlenecks

        .set_index(
            "Gene"
        )
    )

    rows = []

    for node in G.nodes():

        row = {

            "Gene":
                node,

            "is_seed_candidate":
                node
                in seed_genes,

            "is_external_intermediate":
                node
                not in seed_genes,

            "is_source_anchor":
                node
                in sources,

            "is_target_anchor":
                node
                in targets,

            "path_occurrence_count":
                occurrence[node],

            "source_target_pair_occurrence_count":
                len(
                    pair_occurrence[
                        node
                    ]
                ),

            "causal_in_degree":
                G.in_degree(node),

            "causal_out_degree":
                G.out_degree(node),
        }

        if node in (
            bottleneck_lookup.index
        ):

            b = (
                bottleneck_lookup
                .loc[node]
            )

            row[
                "pairs_lost_on_removal"
            ] = int(

                b[
                    "pairs_lost_on_removal"
                ]
            )

            row[
                "pair_loss_fraction"
            ] = float(

                b[
                    "pair_loss_fraction"
                ]
            )

        rows.append(row)

    return pd.DataFrame(
        rows
    )


# ============================================================
# Reactome context validation
# ============================================================

def reactome_version(
    session,
    raw_dir,
    refresh,
):

    path = (

        raw_dir
        / "reactome_database_version.txt"
    )

    if (
        path.exists()
        and not refresh
    ):

        return (

            path
            .read_text(
                encoding="utf-8"
            )
            .strip()
        )

    response = request_with_retry(

        session,

        "GET",

        REACTOME_VERSION_URL,

        headers={
            "Accept":
                "text/plain"
        },
    )

    version = (
        response.text.strip()
    )

    path.write_text(

        version + "\n",

        encoding="utf-8",
    )

    return version


def reactome_analysis(

    session,

    genes,

    label,

    raw_dir,

    processed_dir,

    refresh,
):

    raw_path = (

        raw_dir
        / (
            f"reactome_"
            f"{label}_projection.json"
        )
    )

    if (
        raw_path.exists()
        and not refresh
    ):

        raw = read_json(
            raw_path
        )

    else:

        body = (

            "#Genes\n"

            + "\n".join(
                genes
            )

            + "\n"
        )

        response = request_with_retry(

            session,

            "POST",

            REACTOME_ANALYSIS_URL,

            params={
                "pageSize":
                    100,

                "page":
                    1,
            },

            data=body,

            headers={

                "Content-Type":
                    "text/plain",

                "Accept":
                    "application/json",
            },
        )

        raw = (
            response.json()
        )

        write_json(
            raw_path,
            raw,
        )

    rows = []

    for pathway in raw.get(
        "pathways",
        [],
    ):

        entities = (
            pathway.get(
                "entities",
                {},
            )
            or {}
        )

        reactions = (
            pathway.get(
                "reactions",
                {},
            )
            or {}
        )

        species = pathway.get(
            "species"
        )

        if isinstance(
            species,
            dict,
        ):

            species = species.get(
                "name"
            )

        rows.append({

            "stId":
                pathway.get(
                    "stId"
                ),

            "name":
                pathway.get(
                    "name"
                ),

            "species":
                species,

            "entities_found":
                entities.get(
                    "found"
                ),

            "entities_total":
                entities.get(
                    "total"
                ),

            "p_value":
                entities.get(
                    "pValue"
                ),

            "fdr":
                entities.get(
                    "fdr"
                ),

            "reactions_found":
                reactions.get(
                    "found"
                ),

            "reactions_total":
                reactions.get(
                    "total"
                ),
        })

    df = pd.DataFrame(
        rows
    )

    if (
        not df.empty
        and "fdr"
        in df.columns
    ):

        df.sort_values(

            [
                "fdr",
                "p_value",
            ],

            inplace=True,
        )

    df.to_csv(

        processed_dir
        / (
            f"reactome_"
            f"{label}_pathways.csv"
        ),

        index=False,
    )

    return df


# ============================================================
# Main
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sources",
        default=None,
    )

    parser.add_argument(
        "--death-targets",
        default=None,
    )

    parser.add_argument(
        "--autophagy-targets",
        default=None,
    )

    parser.add_argument(
        "--survival-targets",
        default=None,
    )

    parser.add_argument(
        "--max-edges",
        type=int,
        default=
            DEFAULT_MAX_EDGES,
    )

    parser.add_argument(
        "--max-external",
        type=int,
        default=
            DEFAULT_MAX_EXTERNAL,
    )

    parser.add_argument(
        "--paths-per-pair",
        type=int,
        default=
            DEFAULT_PATHS_PER_PAIR,
    )

    parser.add_argument(
        "--enumeration-cap",
        type=int,
        default=
            DEFAULT_ENUMERATION_CAP,
    )

    parser.add_argument(
        "--refresh",
        action="store_true",
    )

    args = parser.parse_args()

    project_root = (

        Path(__file__)
        .resolve()
        .parent
        .parent
    )

    assignment_path = (

        project_root

        / "Data"
        / "Processed"
        / "Network"
        / "Modules"
        / "module_assignments.csv"
    )

    raw_dir = (

        project_root

        / "Data"
        / "Raw"
        / "Network"
        / "Mechanistic"
    )

    processed_dir = (

        project_root

        / "Data"
        / "Processed"
        / "Network"
        / "Mechanistic"
    )

    figures_dir = (

        project_root

        / "Figures"
        / "Network"
        / "Mechanistic"
    )

    for directory in [

        raw_dir,

        processed_dir,

        figures_dir,

    ]:

        directory.mkdir(

            parents=True,

            exist_ok=True,
        )

    sources = parse_gene_list(

        args.sources,

        DEFAULT_SOURCES,
    )

    target_groups = {

        "death":
            parse_gene_list(

                args.death_targets,

                DEFAULT_TARGET_GROUPS[
                    "death"
                ],
            ),

        "autophagy":
            parse_gene_list(

                args.autophagy_targets,

                DEFAULT_TARGET_GROUPS[
                    "autophagy"
                ],
            ),

        "trophic_survival":
            parse_gene_list(

                args.survival_targets,

                DEFAULT_TARGET_GROUPS[
                    "trophic_survival"
                ],
            ),
    }

    targets = unique_ordered(

        gene

        for genes
        in target_groups.values()

        for gene
        in genes
    )

    print()
    print(
        "R-Noel network analysis "
        "- Phase C"
    )
    print(
        "========================="
    )

    print(
        "Sources:",
        ", ".join(sources),
    )

    print(
        "Targets:",
        ", ".join(targets),
    )

    print(
        "Max edges:",
        args.max_edges,
    )

    print(
        "Max external intermediates:",
        args.max_external,
    )

    assignments = (
        load_candidates(
            assignment_path
        )
    )

    seed_genes = set(
        assignments[
            "Gene"
        ]
    )

    session = (
        requests.Session()
    )

    session.headers.update({

        "User-Agent":
            (
                "R-Noel_BIOL295_"
                "mechanistic_pipeline/1.0"
            )
    })

    # ========================================================
    # OmniPath
    # ========================================================

    raw = download_omnipath(

        session,

        raw_dir,

        args.refresh,
    )

    print(
        "Raw OmniPath rows:",
        len(raw),
    )

    edges = normalize_omnipath(
        raw
    )

    edges.to_csv(

        processed_dir
        / "omnipath_directed_edges.csv",

        index=False,
    )

    print(
        "Collapsed directed edges:",
        len(edges),
    )

    G = build_causal_graph(
        edges
    )

    # Anchor diagnostics.
    anchor_rows = []

    group_lookup = (
        target_group_lookup(
            target_groups
        )
    )

    for gene in sources:

        anchor_rows.append({

            "Gene":
                gene,

            "anchor_role":
                "source",

            "is_seed_candidate":
                gene
                in seed_genes,

            "present_in_omnipath":
                gene in G,

            "in_degree":
                G.in_degree(gene)
                if gene in G
                else 0,

            "out_degree":
                G.out_degree(gene)
                if gene in G
                else 0,
        })

    for gene in targets:

        anchor_rows.append({

            "Gene":
                gene,

            "anchor_role":
                group_lookup[
                    gene
                ],

            "is_seed_candidate":
                gene
                in seed_genes,

            "present_in_omnipath":
                gene in G,

            "in_degree":
                G.in_degree(gene)
                if gene in G
                else 0,

            "out_degree":
                G.out_degree(gene)
                if gene in G
                else 0,
        })

    pd.DataFrame(
        anchor_rows
    ).to_csv(

        processed_dir
        / "anchor_status.csv",

        index=False,
    )

    # ========================================================
    # Directed paths
    # ========================================================

    print()
    print(
        "Searching short "
        "directed paths"
    )
    print(
        "-------------------------"
    )

    paths = search_paths(

        G,

        sources,

        target_groups,

        seed_genes,

        args.max_edges,

        args.max_external,

        args.paths_per_pair,

        args.enumeration_cap,
    )

    paths.to_csv(

        processed_dir
        / "causal_paths.csv",

        index=False,
    )

    if paths.empty:

        print()
        print(
            "No paths passed "
            "the current constraints."
        )

        print(
            "Next try:"
        )

        print(
            "py .\\Scripts\\"
            "04_local_mechanistic_network.py "
            "--max-edges 5"
        )

        return

    # ========================================================
    # Retained path network
    # ========================================================

    path_edges = make_path_edges(

        G,

        paths,

        seed_genes,
    )

    path_edges.to_csv(

        processed_dir
        / "causal_path_edges.csv",

        index=False,
    )

    local_G = retained_graph(
        path_edges
    )

    nx.write_graphml(

        local_G,

        processed_dir
        / (
            "local_mechanistic_"
            "network.graphml"
        ),
    )

    bottlenecks = (
        bottleneck_analysis(

            local_G,

            sources,

            targets,
        )
    )

    bottlenecks.to_csv(

        processed_dir
        / "bottleneck_analysis.csv",

        index=False,
    )

    metrics = node_metrics(

        local_G,

        paths,

        assignments,

        sources,

        targets,

        bottlenecks,
    )

    metrics.to_csv(

        processed_dir
        / "causal_node_metrics.csv",

        index=False,
    )

    # All 126 seeds, including those not
    # represented on a retained path.
    seed_metrics = (
        assignments.merge(

            metrics[
                metrics[
                    "is_seed_candidate"
                ]
            ],

            on="Gene",

            how="left",
        )
    )

    seed_metrics.to_csv(

        processed_dir
        / (
            "seed_candidate_"
            "mechanistic_metrics.csv"
        ),

        index=False,
    )

    # ========================================================
    # Reactome validation
    # ========================================================

    version = reactome_version(

        session,

        raw_dir,

        args.refresh,
    )

    local_seed_nodes = sorted(

        set(
            local_G.nodes()
        )

        .intersection(
            seed_genes
        )
    )

    local_all_nodes = sorted(
        local_G.nodes()
    )

    reactome_seed = (
        reactome_analysis(

            session,

            local_seed_nodes,

            "local_seed",

            raw_dir,

            processed_dir,

            args.refresh,
        )
    )

    reactome_analysis(

        session,

        local_all_nodes,

        "local_all",

        raw_dir,

        processed_dir,

        args.refresh,
    )

    # ========================================================
    # Diagnostics
    # ========================================================

    pair_count = (

        paths[

            [
                "source_anchor",
                "target_anchor",
            ]

        ]

        .drop_duplicates()

        .shape[0]
    )

    diagnostics = {

        "seed_gene_count":
            len(seed_genes),

        "omnipath_taxon":
            MOUSE_TAXON,

        "source_anchors":
            sources,

        "target_groups":
            target_groups,

        "max_edges":
            args.max_edges,

        "max_external_intermediates":
            args.max_external,

        "retained_paths":
            len(paths),

        "reachable_source_target_pairs":
            pair_count,

        "local_network_nodes":
            local_G.number_of_nodes(),

        "local_network_edges":
            local_G.number_of_edges(),

        "local_seed_nodes":
            len(
                local_seed_nodes
            ),

        "external_intermediate_nodes":
            (
                len(local_all_nodes)
                - len(
                    local_seed_nodes
                )
            ),

        "reactome_version":
            version,
    }

    write_json(

        processed_dir
        / (
            "mechanistic_"
            "diagnostics.json"
        ),

        diagnostics,
    )

    # ========================================================
    # Console report
    # ========================================================

    print()
    print(
        "Phase C summary"
    )
    print(
        "---------------"
    )

    print(
        "Retained paths:",
        len(paths),
    )

    print(
        "Reachable source-target pairs:",
        pair_count,
    )

    print(
        "Local network nodes:",
        local_G.number_of_nodes(),
    )

    print(
        "Local network edges:",
        local_G.number_of_edges(),
    )

    print(
        "Seed nodes on paths:",
        len(local_seed_nodes),
    )

    print(
        "External intermediates:",
        (
            len(local_all_nodes)
            - len(local_seed_nodes)
        ),
    )

    print(
        "Reactome version:",
        version,
    )

    print()
    print(
        "Most recurrent seed candidates"
    )
    print(
        "------------------------------"
    )

    display = (

        metrics[

            metrics[
                "is_seed_candidate"
            ]

        ]

        .sort_values(

            [
                "source_target_pair_occurrence_count",
                "path_occurrence_count",
                "pairs_lost_on_removal",
            ],

            ascending=False,
        )

        .head(20)
    )

    print(

        display[

            [
                "Gene",
                "source_target_pair_occurrence_count",
                "path_occurrence_count",
                "pairs_lost_on_removal",
            ]

        ]

        .to_string(
            index=False
        )
    )

    print()
    print(
        "Highest non-anchor bottlenecks"
    )
    print(
        "------------------------------"
    )

    b_display = (

        bottlenecks[

            ~bottlenecks[
                "is_source_anchor"
            ]

            &

            ~bottlenecks[
                "is_target_anchor"
            ]

        ]

        .sort_values(

            [
                "pairs_lost_on_removal",
                "pair_loss_fraction",
            ],

            ascending=False,
        )

        .head(20)
    )

    print(

        b_display[

            [
                "Gene",
                "pairs_lost_on_removal",
                "pair_loss_fraction",
            ]

        ]

        .to_string(
            index=False
        )
    )

    print()
    print(
        "Top Reactome context pathways"
    )
    print(
        "-----------------------------"
    )

    if reactome_seed.empty:

        print(
            "No pathways returned."
        )

    else:

        print(

            reactome_seed[

                [
                    "stId",
                    "name",
                    "entities_found",
                    "fdr",
                ]

            ]

            .head(15)

            .to_string(
                index=False
            )
        )

    print()
    print(
        "Phase C complete."
    )

    print(
        "Processed:",
        processed_dir,
    )

    print()
    print(
        "External intermediates are "
        "causal connectors, not automatically "
        "new experimental candidate genes."
    )


if __name__ == "__main__":

    main()