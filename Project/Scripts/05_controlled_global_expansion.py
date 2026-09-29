#!/usr/bin/env python3
"""
05_controlled_global_expansion.py

Phase D: controlled STRING expansion for the R-Noel cell-death network.

WHY THIS EXISTS
---------------
Phases 02-03 established the high-confidence seed-only STRING structure.
Phases 04-04b established a directed causal layer from OmniPath.

This phase asks a different question:

    Are there proteins OUTSIDE the 126 experimental candidates that
    reproducibly connect otherwise separate experimental modules/components?

The script deliberately avoids a single unrestricted "add lots of nodes"
network expansion.

Instead it:

1. Expands each non-isolate Phase B module independently by a small number
   of STRING neighbors.
2. Expands the isolate pool separately, so the 59 high-confidence STRING
   isolates are not ignored.
3. Scores every external candidate by:
       - number of independent expansion queries supporting it
       - number of experimental seed neighbors
       - number of seed modules touched
       - number of ORIGINAL seed-only connected components touched
       - STRING evidence composition
4. Imports the strongest robust external causal connectors from 04b.
5. Accepts an external node only if it satisfies a transparent rule:
       A) bridges >=2 original seed components with >=2 seed neighbors, OR
       B) recurs in >=2 independent expansion queries with >=2 seed neighbors,
       OR
       C) is a robust 04b causal bottleneck connector.
6. Re-queries STRING with ONLY:
       126 experimental seeds + accepted connectors
   and add_nodes=0.
   This creates a fixed expanded network and prevents a second uncontrolled
   expansion wave.
7. Measures how much accepted external nodes actually compress/bridge the
   original seed-only component structure.

NODE ORIGINS ARE NEVER MIXED:
    experimental_seed
    robust_causal_connector
    new_string_neighbor

IMPORTANT
---------
STRING is a functional-association network, not a causal graph.
A newly accepted STRING connector is a hypothesis-generating network node,
not a newly validated experimental gene.

RUN FROM Project/
-----------------
    py .\Scripts\05_controlled_global_expansion.py

Useful options
--------------
    py .\Scripts\05_controlled_global_expansion.py --refresh
    py .\Scripts\05_controlled_global_expansion.py --module-add-nodes 5 --isolate-add-nodes 20
    py .\Scripts\05_controlled_global_expansion.py --required-score 700

INPUTS
------
Data/Processed/Network/Modules/module_assignments.csv
Data/Processed/Network/string_seed_edges.csv
Data/Processed/Network/Mechanistic/Robustness/robustness_gene_consensus.csv

OUTPUTS
-------
Data/Raw/Network/GlobalExpansion/
    STRING query metadata + raw TSV responses

Data/Processed/Network/GlobalExpansion/
    expansion_query_summary.csv
    candidate_external_neighbors.csv
    accepted_external_connectors.csv
    expanded_edges.csv
    expanded_node_metrics.csv
    component_bridge_summary.csv
    global_expansion_diagnostics.json
    controlled_expanded_network.graphml

Figures/Network/GlobalExpansion/
    accepted_connector_support.png
    component_compression.png
"""

from __future__ import annotations

import argparse
import io
import json
import math
import re
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import requests


# ============================================================
# Configuration
# ============================================================

# Keep the same STRING release used in the earlier pipeline.
STRING_API = "https://version-12-0.string-db.org/api"
STRING_SPECIES = 10090
STRING_NETWORK_TYPE = "functional"

DEFAULT_REQUIRED_SCORE = 700
DEFAULT_MODULE_ADD_NODES = 5
DEFAULT_ISOLATE_ADD_NODES = 20

HTTP_TIMEOUT = 180
HTTP_RETRIES = 3

# Transparent acceptance criteria.
DEFAULT_MIN_SEED_NEIGHBORS = 2
DEFAULT_MIN_QUERY_SUPPORT = 2
DEFAULT_MIN_COMPONENTS_BRIDGED = 2

# 04b robust causal-core criterion:
# present in every tested setting AND a non-anchor union bottleneck
# in at least this many settings.
DEFAULT_MIN_CAUSAL_BOTTLENECK_SETTINGS = 2


# ============================================================
# General helpers
# ============================================================

def clean_gene(value: Any) -> str | None:
    if pd.isna(value):
        return None

    value = str(value).strip().upper()

    if not value:
        return None

    value = re.sub(
        r"\s*\([^)]*\)\s*$",
        "",
        value,
    )

    return value or None


def write_json(path: Path, obj: Any) -> None:
    path.write_text(
        json.dumps(
            obj,
            indent=2,
            ensure_ascii=False,
            default=str,
        ),
        encoding="utf-8",
    )


def sanitize_label(value: str) -> str:
    return re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        str(value),
    )


def split_semicolon(value: Any) -> set[str]:
    if pd.isna(value):
        return set()

    value = str(value).strip()

    if not value:
        return set()

    return {
        item.strip()
        for item in value.split(";")
        if item.strip()
    }


def request_with_retry(
    session: requests.Session,
    method: str,
    url: str,
    **kwargs,
) -> requests.Response:
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
                wait = 2 ** (attempt - 1)

                print(
                    f"  HTTP error: {exc}"
                )
                print(
                    f"  retrying in {wait}s..."
                )

                time.sleep(wait)

    raise RuntimeError(
        "HTTP request failed after "
        f"{HTTP_RETRIES} attempts: "
        f"{last_error}"
    )


# ============================================================
# Input loading
# ============================================================

def load_assignments(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path)

    required = {
        "Gene",
        "module_id",
        "stringId",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            "module_assignments.csv missing: "
            + ", ".join(
                sorted(missing)
            )
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
        keep="first",
        inplace=True,
    )

    return df


def load_seed_edges(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)

    df = pd.read_csv(path)

    required = {
        "stringId_A",
        "stringId_B",
        "score",
    }

    missing = required - set(df.columns)

    if missing:
        raise ValueError(
            "Seed edge table missing: "
            + ", ".join(
                sorted(missing)
            )
        )

    return df


def load_robustness_consensus(
    path: Path,
) -> pd.DataFrame:
    if not path.exists():
        print(
            "WARNING: robustness_gene_consensus.csv "
            "not found. Continuing without 04b "
            "causal-core connectors."
        )
        return pd.DataFrame()

    df = pd.read_csv(path)

    if "Gene" not in df.columns:
        raise ValueError(
            "robustness_gene_consensus.csv "
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

    return df


# ============================================================
# Baseline seed graph/components
# ============================================================

def build_baseline_seed_graph(
    assignments: pd.DataFrame,
    seed_edges: pd.DataFrame,
) -> tuple[nx.Graph, dict[str, int]]:
    string_to_gene = (
        assignments
        .assign(
            stringId=lambda x:
                x["stringId"]
                .astype(str)
        )
        .drop_duplicates(
            "stringId"
        )
        .set_index(
            "stringId"
        )["Gene"]
        .to_dict()
    )

    G = nx.Graph()

    for gene in assignments["Gene"]:
        G.add_node(gene)

    for _, row in seed_edges.iterrows():
        a_id = str(
            row["stringId_A"]
        )
        b_id = str(
            row["stringId_B"]
        )

        a = string_to_gene.get(
            a_id
        )
        b = string_to_gene.get(
            b_id
        )

        if (
            a is None
            or b is None
            or a == b
        ):
            continue

        G.add_edge(
            a,
            b,
            score=float(
                row["score"]
            ),
        )

    components = sorted(
        nx.connected_components(G),
        key=lambda comp: (
            -len(comp),
            min(comp),
        ),
    )

    component_of = {}

    for component_id, comp in enumerate(
        components,
        start=1,
    ):
        for gene in comp:
            component_of[
                gene
            ] = component_id

    return G, component_of


# ============================================================
# STRING API
# ============================================================

STRING_EXPECTED_COLUMNS = [
    "stringId_A",
    "stringId_B",
    "preferredName_A",
    "preferredName_B",
    "ncbiTaxonId",
    "score",
    "nscore",
    "fscore",
    "pscore",
    "ascore",
    "escore",
    "dscore",
    "tscore",
]


def query_string_network(
    session: requests.Session,
    genes: list[str],
    label: str,
    raw_dir: Path,
    *,
    required_score: int,
    add_nodes: int,
    refresh: bool,
) -> pd.DataFrame:
    genes = sorted(
        {
            clean_gene(gene)
            for gene in genes
            if clean_gene(gene)
        }
    )

    if not genes:
        return pd.DataFrame(
            columns=STRING_EXPECTED_COLUMNS
        )

    safe_label = sanitize_label(
        label
    )

    response_path = (
        raw_dir
        / f"{safe_label}.tsv"
    )

    request_path = (
        raw_dir
        / f"{safe_label}.request.json"
    )

    payload = {
        "identifiers": "\r".join(
            genes
        ),
        "species": STRING_SPECIES,
        "required_score": required_score,
        "network_type": STRING_NETWORK_TYPE,
        "add_nodes": add_nodes,
        "caller_identity":
            "R-Noel_BIOL295_global_expansion",
    }

    write_json(
        request_path,
        {
            **payload,
            "identifiers": genes,
            "string_api": STRING_API,
        },
    )

    if (
        response_path.exists()
        and not refresh
    ):
        text = (
            response_path
            .read_text(
                encoding="utf-8"
            )
        )

    else:
        response = request_with_retry(
            session,
            "POST",
            f"{STRING_API}/tsv/network",
            data=payload,
        )

        text = response.text

        response_path.write_text(
            text,
            encoding="utf-8",
        )

        time.sleep(0.25)

    if not text.strip():
        return pd.DataFrame(
            columns=STRING_EXPECTED_COLUMNS
        )

    try:
        df = pd.read_csv(
            io.StringIO(text),
            sep="\t",
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame(
            columns=STRING_EXPECTED_COLUMNS
        )

    for column in [
        "preferredName_A",
        "preferredName_B",
    ]:
        if column in df.columns:
            df[column] = (
                df[column]
                .map(clean_gene)
            )

    return df


# ============================================================
# Expansion query definitions
# ============================================================

def build_expansion_queries(
    assignments: pd.DataFrame,
    *,
    module_add_nodes: int,
    isolate_add_nodes: int,
) -> list[dict[str, Any]]:
    queries = []

    non_isolate = (
        assignments[
            assignments[
                "module_id"
            ] != "ISOLATE"
        ]
    )

    for module_id, group in (
        non_isolate.groupby(
            "module_id"
        )
    ):
        genes = sorted(
            group[
                "Gene"
            ].tolist()
        )

        queries.append({
            "query_id":
                f"MODULE_{module_id}",
            "query_type":
                "module",
            "module_id":
                module_id,
            "genes":
                genes,
            "add_nodes":
                module_add_nodes,
        })

    isolate_genes = sorted(
        assignments.loc[
            assignments[
                "module_id"
            ] == "ISOLATE",
            "Gene",
        ].tolist()
    )

    if isolate_genes:
        queries.append({
            "query_id":
                "ISOLATE_POOL",
            "query_type":
                "isolate_pool",
            "module_id":
                "ISOLATE",
            "genes":
                isolate_genes,
            "add_nodes":
                isolate_add_nodes,
        })

    return queries


# ============================================================
# Candidate extraction
# ============================================================

def edge_non_text_max(
    row: pd.Series,
) -> float:
    channels = [
        "nscore",
        "fscore",
        "pscore",
        "ascore",
        "escore",
        "dscore",
    ]

    values = []

    for channel in channels:
        if channel in row.index:
            try:
                values.append(
                    float(row[channel])
                )
            except (
                TypeError,
                ValueError,
            ):
                pass

    return max(
        values,
        default=0.0,
    )


def extract_candidates_from_query(
    network: pd.DataFrame,
    query_info: dict[str, Any],
    seed_genes: set[str],
    seed_module_of: dict[str, str],
    baseline_component_of: dict[str, int],
) -> tuple[
    list[dict[str, Any]],
    dict[str, Any],
]:
    query_id = query_info[
        "query_id"
    ]

    query_seed_set = set(
        query_info["genes"]
    )

    nodes_returned = set()

    candidate_records = []

    for _, row in network.iterrows():
        a = clean_gene(
            row.get(
                "preferredName_A"
            )
        )
        b = clean_gene(
            row.get(
                "preferredName_B"
            )
        )

        if a:
            nodes_returned.add(a)

        if b:
            nodes_returned.add(b)

        if (
            a is None
            or b is None
        ):
            continue

        score = float(
            row.get(
                "score",
                0.0,
            )
        )

        tscore = float(
            row.get(
                "tscore",
                0.0,
            )
            or 0.0
        )

        non_text_max = (
            edge_non_text_max(
                row
            )
        )

        # External A connected to seed B.
        if (
            a not in seed_genes
            and b in seed_genes
        ):
            candidate_records.append({
                "candidate_gene":
                    a,
                "seed_neighbor":
                    b,
                "seed_module":
                    seed_module_of.get(
                        b
                    ),
                "seed_component":
                    baseline_component_of.get(
                        b
                    ),
                "query_id":
                    query_id,
                "query_type":
                    query_info[
                        "query_type"
                    ],
                "query_module":
                    query_info[
                        "module_id"
                    ],
                "score":
                    score,
                "tscore":
                    tscore,
                "non_text_max":
                    non_text_max,
                "text_mining_dominant":
                    (
                        tscore
                        > non_text_max
                    ),
            })

        # External B connected to seed A.
        if (
            b not in seed_genes
            and a in seed_genes
        ):
            candidate_records.append({
                "candidate_gene":
                    b,
                "seed_neighbor":
                    a,
                "seed_module":
                    seed_module_of.get(
                        a
                    ),
                "seed_component":
                    baseline_component_of.get(
                        a
                    ),
                "query_id":
                    query_id,
                "query_type":
                    query_info[
                        "query_type"
                    ],
                "query_module":
                    query_info[
                        "module_id"
                    ],
                "score":
                    score,
                "tscore":
                    tscore,
                "non_text_max":
                    non_text_max,
                "text_mining_dominant":
                    (
                        tscore
                        > non_text_max
                    ),
            })

    external_nodes = (
        nodes_returned
        - seed_genes
    )

    summary = {
        "query_id":
            query_id,
        "query_type":
            query_info[
                "query_type"
            ],
        "module_id":
            query_info[
                "module_id"
            ],
        "input_seed_count":
            len(query_seed_set),
        "add_nodes":
            query_info[
                "add_nodes"
            ],
        "returned_edge_count":
            len(network),
        "returned_node_count":
            len(nodes_returned),
        "returned_external_node_count":
            len(external_nodes),
        "external_nodes":
            ";".join(
                sorted(
                    external_nodes
                )
            ),
    }

    return (
        candidate_records,
        summary,
    )


def aggregate_candidates(
    records: list[dict[str, Any]],
) -> pd.DataFrame:
    if not records:
        return pd.DataFrame()

    long_df = pd.DataFrame(
        records
    )

    rows = []

    for candidate, group in (
        long_df.groupby(
            "candidate_gene"
        )
    ):
        seed_neighbors = sorted(
            set(
                group[
                    "seed_neighbor"
                ]
            )
        )

        seed_modules = sorted(
            {
                value
                for value
                in group[
                    "seed_module"
                ]
                if pd.notna(value)
            }
        )

        seed_components = sorted(
            {
                int(value)
                for value
                in group[
                    "seed_component"
                ]
                if pd.notna(value)
            }
        )

        query_ids = sorted(
            set(
                group[
                    "query_id"
                ]
            )
        )

        query_types = sorted(
            set(
                group[
                    "query_type"
                ]
            )
        )

        rows.append({
            "Gene":
                candidate,
            "query_support_count":
                len(query_ids),
            "supporting_queries":
                ";".join(
                    query_ids
                ),
            "supporting_query_types":
                ";".join(
                    query_types
                ),
            "seed_neighbor_count":
                len(
                    seed_neighbors
                ),
            "seed_neighbors":
                ";".join(
                    seed_neighbors
                ),
            "seed_module_count":
                len(
                    seed_modules
                ),
            "seed_modules":
                ";".join(
                    seed_modules
                ),
            "baseline_seed_component_count":
                len(
                    seed_components
                ),
            "baseline_seed_components":
                ";".join(
                    str(value)
                    for value
                    in seed_components
                ),
            "max_string_score":
                float(
                    group[
                        "score"
                    ].max()
                ),
            "mean_string_score":
                float(
                    group[
                        "score"
                    ].mean()
                ),
            "max_non_text_evidence":
                float(
                    group[
                        "non_text_max"
                    ].max()
                ),
            "mean_non_text_evidence":
                float(
                    group[
                        "non_text_max"
                    ].mean()
                ),
            "text_mining_dominant_edge_fraction":
                float(
                    group[
                        "text_mining_dominant"
                    ].mean()
                ),
        })

    df = pd.DataFrame(
        rows
    )

    df.sort_values(
        [
            "baseline_seed_component_count",
            "seed_module_count",
            "seed_neighbor_count",
            "query_support_count",
            "max_string_score",
        ],
        ascending=[
            False,
            False,
            False,
            False,
            False,
        ],
        inplace=True,
    )

    return df


# ============================================================
# Causal-core connectors from 04b
# ============================================================

def identify_robust_causal_core(
    consensus: pd.DataFrame,
    seed_genes: set[str],
    min_bottleneck_settings: int,
) -> pd.DataFrame:
    if consensus.empty:
        return pd.DataFrame(
            columns=[
                "Gene",
                "causal_core",
            ]
        )

    required = {
        "Gene",
        "settings_present_count",
        "settings_union_nonanchor_bottleneck",
    }

    missing = (
        required
        - set(
            consensus.columns
        )
    )

    if missing:
        print(
            "WARNING: robustness consensus "
            "does not contain corrected v2 "
            "anchor-excluded columns. "
            "No causal-core connectors imported."
        )

        return pd.DataFrame(
            columns=[
                "Gene",
                "causal_core",
            ]
        )

    max_settings = int(
        consensus[
            "settings_present_count"
        ].max()
    )

    core = consensus[
        (
            consensus[
                "settings_present_count"
            ]
            == max_settings
        )
        &
        (
            consensus[
                "settings_union_nonanchor_bottleneck"
            ]
            >= min_bottleneck_settings
        )
        &
        (
            ~consensus[
                "Gene"
            ].isin(
                seed_genes
            )
        )
    ].copy()

    core["causal_core"] = True

    return core


# ============================================================
# Candidate acceptance
# ============================================================

def apply_acceptance_rules(
    candidates: pd.DataFrame,
    causal_core: pd.DataFrame,
    *,
    min_seed_neighbors: int,
    min_query_support: int,
    min_components_bridged: int,
) -> pd.DataFrame:
    if candidates.empty:
        candidates = pd.DataFrame(
            columns=[
                "Gene",
                "query_support_count",
                "seed_neighbor_count",
                "baseline_seed_component_count",
            ]
        )

    causal_genes = set(
        causal_core["Gene"]
        if not causal_core.empty
        else []
    )

    output = (
        candidates.copy()
    )

    if "Gene" not in output.columns:
        output["Gene"] = []

    output[
        "robust_causal_core"
    ] = output[
        "Gene"
    ].isin(
        causal_genes
    )

    output[
        "rule_cross_component"
    ] = (
        (
            output[
                "seed_neighbor_count"
            ]
            >= min_seed_neighbors
        )
        &
        (
            output[
                "baseline_seed_component_count"
            ]
            >= min_components_bridged
        )
    )

    output[
        "rule_recurrent_neighbor"
    ] = (
        (
            output[
                "query_support_count"
            ]
            >= min_query_support
        )
        &
        (
            output[
                "seed_neighbor_count"
            ]
            >= min_seed_neighbors
        )
    )

    output[
        "rule_causal_core"
    ] = output[
        "robust_causal_core"
    ]

    output[
        "accepted"
    ] = (
        output[
            "rule_cross_component"
        ]
        |
        output[
            "rule_recurrent_neighbor"
        ]
        |
        output[
            "rule_causal_core"
        ]
    )

    def reasons(row):
        values = []

        if row[
            "rule_cross_component"
        ]:
            values.append(
                "bridges_seed_components"
            )

        if row[
            "rule_recurrent_neighbor"
        ]:
            values.append(
                "recurrent_module_support"
            )

        if row[
            "rule_causal_core"
        ]:
            values.append(
                "robust_causal_bottleneck"
            )

        return ";".join(
            values
        )

    output[
        "acceptance_reason"
    ] = output.apply(
        reasons,
        axis=1,
    )

    # Robust causal connectors may not have appeared in
    # the module-wise STRING expansion. Add them explicitly.
    missing_causal = sorted(
        causal_genes
        - set(
            output[
                "Gene"
            ]
        )
    )

    if missing_causal:
        extra = []

        causal_lookup = (
            causal_core
            .set_index(
                "Gene"
            )
        )

        for gene in missing_causal:
            row = {
                column:
                    None
                for column
                in output.columns
            }

            row["Gene"] = gene
            row["robust_causal_core"] = True
            row["rule_cross_component"] = False
            row["rule_recurrent_neighbor"] = False
            row["rule_causal_core"] = True
            row["accepted"] = True
            row[
                "acceptance_reason"
            ] = (
                "robust_causal_bottleneck"
            )

            extra.append(row)

        output = pd.concat(
            [
                output,
                pd.DataFrame(extra),
            ],
            ignore_index=True,
        )

    # Attach selected 04b metrics for provenance.
    if not causal_core.empty:
        causal_cols = [
            column
            for column in [
                "Gene",
                "settings_present_count",
                "max_internal_pair_occurrence_count",
                "settings_required_internal_for_at_least_one_pair",
                "max_required_internal_pair_count",
                "settings_union_nonanchor_bottleneck",
                "max_union_nonanchor_pairs_lost",
            ]
            if column
            in causal_core.columns
        ]

        output = output.merge(
            causal_core[
                causal_cols
            ],
            on="Gene",
            how="left",
            suffixes=(
                "",
                "_04b",
            ),
        )

    return output


# ============================================================
# Final fixed expanded network
# ============================================================

def build_expanded_graph(
    fixed_network: pd.DataFrame,
    all_nodes: set[str],
) -> nx.Graph:
    G = nx.Graph()

    for gene in all_nodes:
        G.add_node(
            gene
        )

    for _, row in (
        fixed_network.iterrows()
    ):
        a = clean_gene(
            row.get(
                "preferredName_A"
            )
        )
        b = clean_gene(
            row.get(
                "preferredName_B"
            )
        )

        if (
            a is None
            or b is None
            or a == b
        ):
            continue

        G.add_edge(
            a,
            b,
            score=float(
                row.get(
                    "score",
                    0.0,
                )
            ),
            nscore=float(
                row.get(
                    "nscore",
                    0.0,
                )
                or 0.0
            ),
            fscore=float(
                row.get(
                    "fscore",
                    0.0,
                )
                or 0.0
            ),
            pscore=float(
                row.get(
                    "pscore",
                    0.0,
                )
                or 0.0
            ),
            ascore=float(
                row.get(
                    "ascore",
                    0.0,
                )
                or 0.0
            ),
            escore=float(
                row.get(
                    "escore",
                    0.0,
                )
                or 0.0
            ),
            dscore=float(
                row.get(
                    "dscore",
                    0.0,
                )
                or 0.0
            ),
            tscore=float(
                row.get(
                    "tscore",
                    0.0,
                )
                or 0.0
            ),
        )

    return G


def seed_component_compression(
    G: nx.Graph,
    seed_genes: set[str],
) -> dict[str, Any]:
    seed_component_sizes = []

    for comp in nx.connected_components(
        G
    ):
        seed_count = len(
            set(comp)
            .intersection(
                seed_genes
            )
        )

        if seed_count > 0:
            seed_component_sizes.append(
                seed_count
            )

    return {
        "seed_connected_component_count":
            len(
                seed_component_sizes
            ),
        "largest_seed_connected_component":
            max(
                seed_component_sizes,
                default=0,
            ),
        "seed_component_sizes":
            sorted(
                seed_component_sizes,
                reverse=True,
            ),
    }


def expanded_node_metrics(
    G: nx.Graph,
    assignments: pd.DataFrame,
    accepted: pd.DataFrame,
    baseline_component_of: dict[str, int],
    causal_core_genes: set[str],
) -> pd.DataFrame:
    seed_genes = set(
        assignments[
            "Gene"
        ]
    )

    seed_module_of = (
        assignments
        .set_index(
            "Gene"
        )["module_id"]
        .to_dict()
    )

    accepted_genes = set(
        accepted.loc[
            accepted[
                "accepted"
            ],
            "Gene",
        ].dropna()
    )

    # Use inverse score as a distance for betweenness.
    H = nx.Graph()

    for node in G.nodes():
        H.add_node(node)

    for a, b, attrs in (
        G.edges(
            data=True
        )
    ):
        score = float(
            attrs.get(
                "score",
                0.0,
            )
        )

        H.add_edge(
            a,
            b,
            distance=(
                1.0
                / max(
                    score,
                    1e-12,
                )
            ),
        )

    betweenness = (
        nx.betweenness_centrality(
            H,
            weight="distance",
            normalized=True,
        )
        if G.number_of_nodes()
        else {}
    )

    pagerank = (
        nx.pagerank(
            G,
            weight="score",
        )
        if G.number_of_nodes()
        else {}
    )

    rows = []

    for node in G.nodes():
        if node in seed_genes:
            origin = (
                "experimental_seed"
            )
        elif node in causal_core_genes:
            origin = (
                "robust_causal_connector"
            )
        elif node in accepted_genes:
            origin = (
                "new_string_neighbor"
            )
        else:
            origin = (
                "external_other"
            )

        seed_neighbors = sorted(
            {
                neighbor
                for neighbor
                in G.neighbors(
                    node
                )
                if neighbor
                in seed_genes
            }
        )

        seed_modules = sorted(
            {
                seed_module_of.get(
                    gene
                )
                for gene
                in seed_neighbors
                if seed_module_of.get(
                    gene
                )
                is not None
            }
        )

        seed_components = sorted(
            {
                baseline_component_of.get(
                    gene
                )
                for gene
                in seed_neighbors
                if baseline_component_of.get(
                    gene
                )
                is not None
            }
        )

        rows.append({
            "Gene":
                node,
            "node_origin":
                origin,
            "seed_module":
                seed_module_of.get(
                    node
                ),
            "baseline_seed_component":
                baseline_component_of.get(
                    node
                ),
            "degree":
                G.degree(
                    node
                ),
            "weighted_degree":
                sum(
                    float(
                        attrs.get(
                            "score",
                            0.0,
                        )
                    )
                    for _, _, attrs
                    in G.edges(
                        node,
                        data=True,
                    )
                ),
            "betweenness_weighted":
                betweenness.get(
                    node,
                    0.0,
                ),
            "pagerank":
                pagerank.get(
                    node,
                    0.0,
                ),
            "seed_neighbor_count":
                len(
                    seed_neighbors
                ),
            "seed_neighbors":
                ";".join(
                    seed_neighbors
                ),
            "seed_module_neighbor_count":
                len(
                    seed_modules
                ),
            "seed_modules_neighbors":
                ";".join(
                    seed_modules
                ),
            "baseline_seed_component_neighbor_count":
                len(
                    seed_components
                ),
            "baseline_seed_components_neighbors":
                ";".join(
                    str(value)
                    for value
                    in seed_components
                ),
        })

    return pd.DataFrame(
        rows
    )


def component_bridge_summary(
    G: nx.Graph,
    external_genes: set[str],
    seed_genes: set[str],
    baseline_component_of: dict[str, int],
) -> pd.DataFrame:
    rows = []

    for gene in sorted(
        external_genes
    ):
        if gene not in G:
            continue

        seed_neighbors = sorted(
            {
                neighbor
                for neighbor
                in G.neighbors(
                    gene
                )
                if neighbor
                in seed_genes
            }
        )

        components = sorted(
            {
                baseline_component_of[
                    seed
                ]
                for seed
                in seed_neighbors
                if seed
                in baseline_component_of
            }
        )

        rows.append({
            "Gene":
                gene,
            "seed_neighbor_count":
                len(
                    seed_neighbors
                ),
            "seed_neighbors":
                ";".join(
                    seed_neighbors
                ),
            "baseline_seed_components_bridged":
                len(
                    components
                ),
            "baseline_component_ids":
                ";".join(
                    str(value)
                    for value
                    in components
                ),
            "is_actual_cross_component_bridge":
                len(
                    components
                ) >= 2,
        })

    df = pd.DataFrame(
        rows
    )

    if not df.empty:
        df.sort_values(
            [
                "baseline_seed_components_bridged",
                "seed_neighbor_count",
            ],
            ascending=[
                False,
                False,
            ],
            inplace=True,
        )

    return df


# ============================================================
# Figures
# ============================================================

def plot_accepted_connector_support(
    accepted: pd.DataFrame,
    figures_dir: Path,
) -> None:
    df = accepted[
        accepted[
            "accepted"
        ]
    ].copy()

    if df.empty:
        return

    for column in [
        "baseline_seed_component_count",
        "seed_neighbor_count",
    ]:
        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        ).fillna(0)

    df.sort_values(
        [
            "baseline_seed_component_count",
            "seed_neighbor_count",
        ],
        ascending=[
            False,
            False,
        ],
        inplace=True,
    )

    df = df.head(
        30
    )

    df = df.iloc[::-1]

    plt.figure(
        figsize=(
            10,
            max(
                6,
                0.32 * len(df)
                + 2,
            ),
        )
    )

    y = range(
        len(df)
    )

    plt.barh(
        y,
        df[
            "baseline_seed_component_count"
        ],
    )

    plt.yticks(
        y,
        df["Gene"],
    )

    plt.xlabel(
        "Original seed components touched"
    )

    plt.title(
        "Accepted external STRING connector support"
    )

    plt.tight_layout()

    plt.savefig(
        figures_dir
        / "accepted_connector_support.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def plot_component_compression(
    baseline_count: int,
    expanded_count: int,
    figures_dir: Path,
) -> None:
    plt.figure(
        figsize=(7, 5)
    )

    plt.bar(
        [
            "Seed-only",
            "With accepted\nconnectors",
        ],
        [
            baseline_count,
            expanded_count,
        ],
    )

    plt.ylabel(
        "Connected components containing seed genes"
    )

    plt.title(
        "Compression of seed-network components"
    )

    plt.tight_layout()

    plt.savefig(
        figures_dir
        / "component_compression.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


# ============================================================
# Main
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Controlled high-confidence STRING "
            "expansion of the R-Noel seed network."
        )
    )

    parser.add_argument(
        "--required-score",
        type=int,
        default=
            DEFAULT_REQUIRED_SCORE,
    )

    parser.add_argument(
        "--module-add-nodes",
        type=int,
        default=
            DEFAULT_MODULE_ADD_NODES,
    )

    parser.add_argument(
        "--isolate-add-nodes",
        type=int,
        default=
            DEFAULT_ISOLATE_ADD_NODES,
    )

    parser.add_argument(
        "--min-seed-neighbors",
        type=int,
        default=
            DEFAULT_MIN_SEED_NEIGHBORS,
    )

    parser.add_argument(
        "--min-query-support",
        type=int,
        default=
            DEFAULT_MIN_QUERY_SUPPORT,
    )

    parser.add_argument(
        "--min-components-bridged",
        type=int,
        default=
            DEFAULT_MIN_COMPONENTS_BRIDGED,
    )

    parser.add_argument(
        "--min-causal-bottleneck-settings",
        type=int,
        default=
            DEFAULT_MIN_CAUSAL_BOTTLENECK_SETTINGS,
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

    processed_network = (
        project_root
        / "Data"
        / "Processed"
        / "Network"
    )

    assignments_path = (
        processed_network
        / "Modules"
        / "module_assignments.csv"
    )

    seed_edges_path = (
        processed_network
        / "string_seed_edges.csv"
    )

    robustness_path = (
        processed_network
        / "Mechanistic"
        / "Robustness"
        / "robustness_gene_consensus.csv"
    )

    raw_dir = (
        project_root
        / "Data"
        / "Raw"
        / "Network"
        / "GlobalExpansion"
    )

    processed_dir = (
        processed_network
        / "GlobalExpansion"
    )

    figures_dir = (
        project_root
        / "Figures"
        / "Network"
        / "GlobalExpansion"
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

    print()
    print(
        "R-Noel network analysis - Phase D"
    )
    print(
        "================================="
    )

    print(
        f"STRING endpoint: {STRING_API}"
    )

    print(
        f"Required score: "
        f"{args.required_score}"
    )

    # ========================================================
    # Inputs
    # ========================================================

    assignments = (
        load_assignments(
            assignments_path
        )
    )

    seed_edges = (
        load_seed_edges(
            seed_edges_path
        )
    )

    robustness = (
        load_robustness_consensus(
            robustness_path
        )
    )

    seed_genes = set(
        assignments[
            "Gene"
        ]
    )

    seed_module_of = (
        assignments
        .set_index(
            "Gene"
        )["module_id"]
        .to_dict()
    )

    baseline_G, baseline_component_of = (
        build_baseline_seed_graph(
            assignments,
            seed_edges,
        )
    )

    baseline_stats = (
        seed_component_compression(
            baseline_G,
            seed_genes,
        )
    )

    print(
        f"Experimental seeds: "
        f"{len(seed_genes)}"
    )

    print(
        f"Baseline seed components: "
        f"{baseline_stats['seed_connected_component_count']}"
    )

    print(
        f"Baseline largest seed component: "
        f"{baseline_stats['largest_seed_connected_component']}"
    )

    # ========================================================
    # Robust causal external core
    # ========================================================

    causal_core = (
        identify_robust_causal_core(
            robustness,
            seed_genes,
            args.min_causal_bottleneck_settings,
        )
    )

    causal_core_genes = set(
        causal_core[
            "Gene"
        ]
        if not causal_core.empty
        else []
    )

    print(
        "Robust external causal core:",
        (
            ", ".join(
                sorted(
                    causal_core_genes
                )
            )
            if causal_core_genes
            else "none"
        ),
    )

    # ========================================================
    # Independent controlled expansion queries
    # ========================================================

    queries = (
        build_expansion_queries(
            assignments,
            module_add_nodes=
                args.module_add_nodes,
            isolate_add_nodes=
                args.isolate_add_nodes,
        )
    )

    print()
    print(
        "Independent STRING expansion queries"
    )
    print(
        "------------------------------------"
    )

    print(
        f"Queries: {len(queries)}"
    )

    session = requests.Session()

    session.headers.update({
        "User-Agent":
            "R-Noel_BIOL295_global_expansion/1.0"
    })

    candidate_records = []
    query_summaries = []

    for query in queries:
        print(
            f"{query['query_id']}: "
            f"{len(query['genes'])} seed(s), "
            f"+{query['add_nodes']} node allowance"
        )

        network = query_string_network(
            session,
            query["genes"],
            query["query_id"],
            raw_dir,
            required_score=
                args.required_score,
            add_nodes=
                query["add_nodes"],
            refresh=
                args.refresh,
        )

        records, summary = (
            extract_candidates_from_query(
                network,
                query,
                seed_genes,
                seed_module_of,
                baseline_component_of,
            )
        )

        candidate_records.extend(
            records
        )

        query_summaries.append(
            summary
        )

    query_summary_df = (
        pd.DataFrame(
            query_summaries
        )
    )

    query_summary_df.to_csv(
        processed_dir
        / "expansion_query_summary.csv",
        index=False,
    )

    candidates = (
        aggregate_candidates(
            candidate_records
        )
    )

    # ========================================================
    # Acceptance rules
    # ========================================================

    accepted_all = (
        apply_acceptance_rules(
            candidates,
            causal_core,
            min_seed_neighbors=
                args.min_seed_neighbors,
            min_query_support=
                args.min_query_support,
            min_components_bridged=
                args.min_components_bridged,
        )
    )

    candidates_out = (
        accepted_all.copy()
    )

    candidates_out.to_csv(
        processed_dir
        / "candidate_external_neighbors.csv",
        index=False,
    )

    accepted = (
        accepted_all[
            accepted_all[
                "accepted"
            ]
        ].copy()
    )

    accepted.sort_values(
        [
            "rule_causal_core",
            "baseline_seed_component_count",
            "seed_neighbor_count",
            "query_support_count",
        ],
        ascending=[
            False,
            False,
            False,
            False,
        ],
        inplace=True,
        na_position="last",
    )

    accepted.to_csv(
        processed_dir
        / "accepted_external_connectors.csv",
        index=False,
    )

    accepted_genes = set(
        accepted[
            "Gene"
        ].dropna()
    )

    print()
    print(
        f"External candidates observed: "
        f"{len(candidates)}"
    )

    print(
        f"Accepted external connectors: "
        f"{len(accepted_genes)}"
    )

    # ========================================================
    # Fixed expanded graph: seeds + accepted only
    # ========================================================

    fixed_nodes = (
        seed_genes
        | accepted_genes
    )

    print()
    print(
        "Building fixed expanded network "
        "(add_nodes=0)..."
    )

    fixed_network = (
        query_string_network(
            session,
            sorted(
                fixed_nodes
            ),
            "FINAL_FIXED_EXPANDED_NETWORK",
            raw_dir,
            required_score=
                args.required_score,
            add_nodes=0,
            refresh=
                args.refresh,
        )
    )

    fixed_network.to_csv(
        processed_dir
        / "expanded_edges.csv",
        index=False,
    )

    expanded_G = (
        build_expanded_graph(
            fixed_network,
            fixed_nodes,
        )
    )

    nx.write_graphml(
        expanded_G,
        processed_dir
        / "controlled_expanded_network.graphml",
    )

    expanded_stats = (
        seed_component_compression(
            expanded_G,
            seed_genes,
        )
    )

    node_metrics = (
        expanded_node_metrics(
            expanded_G,
            assignments,
            accepted,
            baseline_component_of,
            causal_core_genes,
        )
    )

    node_metrics.sort_values(
        [
            "node_origin",
            "baseline_seed_component_neighbor_count",
            "seed_neighbor_count",
            "betweenness_weighted",
        ],
        ascending=[
            True,
            False,
            False,
            False,
        ],
        inplace=True,
    )

    node_metrics.to_csv(
        processed_dir
        / "expanded_node_metrics.csv",
        index=False,
    )

    bridges = (
        component_bridge_summary(
            expanded_G,
            accepted_genes,
            seed_genes,
            baseline_component_of,
        )
    )

    bridges.to_csv(
        processed_dir
        / "component_bridge_summary.csv",
        index=False,
    )

    # ========================================================
    # Diagnostics + figures
    # ========================================================

    actual_cross_component = (
        int(
            bridges[
                "is_actual_cross_component_bridge"
            ].sum()
        )
        if not bridges.empty
        else 0
    )

    diagnostics = {
        "string_api":
            STRING_API,
        "string_species":
            STRING_SPECIES,
        "network_type":
            STRING_NETWORK_TYPE,
        "required_score":
            args.required_score,
        "experimental_seed_count":
            len(
                seed_genes
            ),
        "non_isolate_module_add_nodes":
            args.module_add_nodes,
        "isolate_pool_add_nodes":
            args.isolate_add_nodes,
        "expansion_query_count":
            len(
                queries
            ),
        "external_candidate_count":
            len(
                candidates
            ),
        "accepted_external_connector_count":
            len(
                accepted_genes
            ),
        "accepted_external_connectors":
            sorted(
                accepted_genes
            ),
        "robust_causal_core":
            sorted(
                causal_core_genes
            ),
        "acceptance_rules": {
            "minimum_seed_neighbors":
                args.min_seed_neighbors,
            "minimum_query_support":
                args.min_query_support,
            "minimum_original_components_bridged":
                args.min_components_bridged,
            "minimum_04b_union_bottleneck_settings":
                args.min_causal_bottleneck_settings,
        },
        "baseline_seed_component_count":
            baseline_stats[
                "seed_connected_component_count"
            ],
        "baseline_largest_seed_component":
            baseline_stats[
                "largest_seed_connected_component"
            ],
        "expanded_seed_component_count":
            expanded_stats[
                "seed_connected_component_count"
            ],
        "expanded_largest_seed_component":
            expanded_stats[
                "largest_seed_connected_component"
            ],
        "seed_component_count_reduction":
            (
                baseline_stats[
                    "seed_connected_component_count"
                ]
                -
                expanded_stats[
                    "seed_connected_component_count"
                ]
            ),
        "actual_external_cross_component_bridges":
            actual_cross_component,
        "expanded_graph_nodes":
            expanded_G.number_of_nodes(),
        "expanded_graph_edges":
            expanded_G.number_of_edges(),
    }

    write_json(
        processed_dir
        / "global_expansion_diagnostics.json",
        diagnostics,
    )

    plot_accepted_connector_support(
        accepted,
        figures_dir,
    )

    plot_component_compression(
        baseline_stats[
            "seed_connected_component_count"
        ],
        expanded_stats[
            "seed_connected_component_count"
        ],
        figures_dir,
    )

    # ========================================================
    # Console report
    # ========================================================

    print()
    print(
        "Phase D summary"
    )
    print(
        "---------------"
    )

    print(
        "Baseline seed components:",
        baseline_stats[
            "seed_connected_component_count"
        ],
    )

    print(
        "Expanded seed components:",
        expanded_stats[
            "seed_connected_component_count"
        ],
    )

    print(
        "Largest seed component:",
        (
            f"{baseline_stats['largest_seed_connected_component']}"
            " -> "
            f"{expanded_stats['largest_seed_connected_component']}"
        ),
    )

    print(
        "Actual accepted connectors touching "
        ">=2 original seed components:",
        actual_cross_component,
    )

    print()
    print(
        "Accepted external connectors"
    )
    print(
        "----------------------------"
    )

    if accepted.empty:
        print(
            "None passed the acceptance rules."
        )
    else:
        display_cols = [
            column
            for column in [
                "Gene",
                "acceptance_reason",
                "query_support_count",
                "seed_neighbor_count",
                "seed_module_count",
                "baseline_seed_component_count",
                "max_string_score",
                "max_non_text_evidence",
                "text_mining_dominant_edge_fraction",
                "settings_union_nonanchor_bottleneck",
                "max_union_nonanchor_pairs_lost",
            ]
            if column
            in accepted.columns
        ]

        print(
            accepted[
                display_cols
            ]
            .head(30)
            .to_string(
                index=False
            )
        )

    print()
    print(
        "Strongest ACTUAL cross-component bridges "
        "in the final fixed network"
    )
    print(
        "------------------------------------------"
    )

    if bridges.empty:
        print("None.")
    else:
        print(
            bridges[
                [
                    "Gene",
                    "seed_neighbor_count",
                    "baseline_seed_components_bridged",
                    "seed_neighbors",
                ]
            ]
            .head(25)
            .to_string(
                index=False
            )
        )

    print()
    print(
        "Phase D complete."
    )

    print(
        f"Processed: {processed_dir}"
    )

    print(
        f"Figures: {figures_dir}"
    )

    print()
    print(
        "Interpretation rule: accepted STRING neighbors are "
        "network hypotheses. They remain distinct from the "
        "126 experimental seed genes and from 04b robust "
        "causal connectors."
    )


if __name__ == "__main__":
    main()
