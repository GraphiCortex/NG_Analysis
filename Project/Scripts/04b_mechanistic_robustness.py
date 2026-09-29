#!/usr/bin/env python3
"""
04b_mechanistic_robustness.py

Sensitivity/robustness analysis for Phase C of the R-Noel cell-death project.

This script DOES NOT re-download OmniPath. It begins from the normalized directed
edge table written by 04_local_mechanistic_network.py:

    Data/Processed/Network/Mechanistic/omnipath_directed_edges.csv

and from the Phase B candidate table:

    Data/Processed/Network/Modules/module_assignments.csv

It repeats the short-path search under four nearby parameter settings:

    max_edges  max_external_intermediates
        4                 2
        5                 2
        4                 3
        5                 3

For every setting it distinguishes two different notions of bottleneck, and
reports anchor-inclusive and anchor-excluded versions separately:

1. PATH-SET REQUIREMENT
   A node is "required" for a source->target pair if it occurs in EVERY retained
   path for that explicit pair.

2. UNION-GRAPH BOTTLENECK
   All retained paths are merged into one directed graph. Edges from different
   paths can recombine, so this graph may connect more source->target pairs than
   were explicitly retained. Node-removal connectivity is measured separately
   on this union graph.

The script also flags UniProt-accession-like external nodes (e.g. A0A..., D3...)
and, by default, excludes them from the robustness path search. They are written
to a separate audit file rather than silently discarded.

RUN FROM Project/
-----------------
    py .\Scripts\04b_mechanistic_robustness.py

Useful:
    py .\Scripts\04b_mechanistic_robustness.py --keep-accession-like

OUTPUT
------
Data/Processed/Network/Mechanistic/Robustness/
    robustness_setting_summary.csv
    robustness_paths.csv
    robustness_pair_metrics.csv
    robustness_gene_by_setting.csv
    robustness_gene_consensus.csv
    robustness_seed_consensus.csv
    accession_like_external_nodes.csv
    robustness_diagnostics.json

Figures/Network/Mechanistic/Robustness/
    gene_setting_stability.png
    source_target_reachability.png
    paths_per_setting.png
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


# ============================================================
# Configuration
# ============================================================

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

SETTINGS = [
    ("E4_X2", 4, 2),
    ("E5_X2", 5, 2),
    ("E4_X3", 4, 3),
    ("E5_X3", 5, 3),
]

DEFAULT_PATHS_PER_PAIR = 25
DEFAULT_ENUMERATION_CAP = 10000


# ============================================================
# Helpers
# ============================================================

def clean_gene(value):
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


def unique_ordered(values):
    seen = set()
    out = []

    for value in values:
        if value not in seen:
            seen.add(value)
            out.append(value)

    return out


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


def parse_gene_list(text, default):
    if text is None:
        return list(default)

    genes = []

    for item in text.split(","):
        gene = clean_gene(item)
        if gene:
            genes.append(gene)

    return unique_ordered(genes)


def target_group_lookup(groups):
    lookup = {}

    for group, genes in groups.items():
        for gene in genes:
            lookup[gene] = group

    return lookup


# UniProt accession patterns.
# This catches examples such as D3YVY2 and A0A087WSP5 without treating ordinary
# gene symbols such as BCL2L11 or CDKN2A as accessions.
UNIPROT_RE = re.compile(
    r"^(?:"
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]"
    r"|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
    r")$"
)


def looks_uniprot_accession(value):
    if value is None:
        return False
    return bool(
        UNIPROT_RE.fullmatch(
            str(value).strip().upper()
        )
    )


# ============================================================
# Input
# ============================================================

def load_inputs(edge_path, assignment_path):
    if not edge_path.exists():
        raise FileNotFoundError(edge_path)

    if not assignment_path.exists():
        raise FileNotFoundError(assignment_path)

    edges = pd.read_csv(edge_path)
    assignments = pd.read_csv(assignment_path)

    required_edges = {
        "source_gene",
        "target_gene",
        "sign",
        "sign_label",
        "evidence_score",
        "evidence_cost",
    }

    missing = required_edges - set(edges.columns)

    if missing:
        raise ValueError(
            "OmniPath edge table missing: "
            + ", ".join(sorted(missing))
        )

    if "Gene" not in assignments.columns:
        raise ValueError(
            "module_assignments.csv is missing the Gene column."
        )

    edges["source_gene"] = edges["source_gene"].map(clean_gene)
    edges["target_gene"] = edges["target_gene"].map(clean_gene)

    edges.dropna(
        subset=["source_gene", "target_gene"],
        inplace=True,
    )

    assignments["Gene"] = assignments["Gene"].map(clean_gene)

    assignments.dropna(
        subset=["Gene"],
        inplace=True,
    )

    assignments.drop_duplicates(
        "Gene",
        keep="first",
        inplace=True,
    )

    return edges, assignments


def audit_accession_like_nodes(edges, seed_genes):
    all_nodes = set(edges["source_gene"]) | set(edges["target_gene"])

    accession_nodes = sorted(
        node
        for node in all_nodes
        if node not in seed_genes
        and looks_uniprot_accession(node)
    )

    rows = []

    for node in accession_nodes:
        outgoing = int(
            (edges["source_gene"] == node).sum()
        )
        incoming = int(
            (edges["target_gene"] == node).sum()
        )

        rows.append({
            "identifier": node,
            "is_seed_candidate": False,
            "incoming_edges": incoming,
            "outgoing_edges": outgoing,
            "total_incident_edges": incoming + outgoing,
        })

    return pd.DataFrame(rows)


# ============================================================
# Graph + path search
# ============================================================

def build_graph(edges, excluded_nodes=None):
    excluded_nodes = excluded_nodes or set()

    G = nx.DiGraph()

    for _, row in edges.iterrows():
        source = row["source_gene"]
        target = row["target_gene"]

        if source in excluded_nodes or target in excluded_nodes:
            continue

        G.add_edge(
            source,
            target,
            sign=int(row["sign"]),
            sign_label=row["sign_label"],
            evidence_score=float(row["evidence_score"]),
            evidence_cost=float(row["evidence_cost"]),
            n_resources=int(row.get("n_resources", 0)),
            n_references=int(row.get("n_references", 0)),
            resources=row.get("resources", ""),
            references=row.get("references", ""),
        )

    return G


def corridor(G, source, target, max_edges):
    forward = nx.single_source_shortest_path_length(
        G,
        source,
        cutoff=max_edges,
    )

    reverse = nx.single_source_shortest_path_length(
        G.reverse(copy=False),
        target,
        cutoff=max_edges,
    )

    nodes = {
        node
        for node in set(forward).intersection(reverse)
        if forward[node] + reverse[node] <= max_edges
    }

    return G.subgraph(nodes).copy()


def path_sign(G, path):
    product = 1
    unsigned = False
    conflicting = False

    for source, target in zip(
        path[:-1],
        path[1:],
    ):
        attrs = G[source][target]
        sign = int(attrs.get("sign", 0))
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
        return 0, "conflicting"

    if unsigned:
        return 0, "unknown"

    if product > 0:
        return 1, "net_activation"

    return -1, "net_inhibition"


def summarize_path(
    G,
    path,
    seed_genes,
    source,
    target,
    target_group,
):
    intermediates = path[1:-1]

    external = [
        gene
        for gene in intermediates
        if gene not in seed_genes
    ]

    edge_attrs = [
        G[a][b]
        for a, b in zip(
            path[:-1],
            path[1:],
        )
    ]

    sign, sign_label = path_sign(
        G,
        path,
    )

    return {
        "source_anchor": source,
        "target_anchor": target,
        "target_group": target_group,
        "num_edges": len(path) - 1,
        "num_intermediate_nodes": len(intermediates),
        "num_external_intermediates": len(external),
        "external_intermediates": ";".join(external),
        "path_nodes": ";".join(path),
        "path_sign": sign,
        "path_sign_label": sign_label,
        "path_cost": sum(
            float(attrs["evidence_cost"])
            for attrs in edge_attrs
        ),
        "min_edge_evidence_score": min(
            float(attrs["evidence_score"])
            for attrs in edge_attrs
        ),
        "mean_edge_evidence_score": sum(
            float(attrs["evidence_score"])
            for attrs in edge_attrs
        ) / len(edge_attrs),
    }


def find_pair_paths(
    G,
    source,
    target,
    target_group,
    seed_genes,
    max_edges,
    max_external,
    paths_per_pair,
    enumeration_cap,
):
    if source not in G or target not in G:
        return []

    try:
        shortest = nx.shortest_path_length(
            G,
            source,
            target,
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

    generator = nx.all_simple_paths(
        local,
        source=source,
        target=target,
        cutoff=max_edges,
    )

    rows = []

    for path in itertools.islice(
        generator,
        enumeration_cap,
    ):
        external_count = sum(
            node not in seed_genes
            for node in path[1:-1]
        )

        if external_count > max_external:
            continue

        rows.append(
            summarize_path(
                G,
                path,
                seed_genes,
                source,
                target,
                target_group,
            )
        )

    rows.sort(
        key=lambda row: (
            row["num_edges"],
            row["num_external_intermediates"],
            row["path_cost"],
            -row["min_edge_evidence_score"],
            row["path_nodes"],
        )
    )

    return rows[:paths_per_pair]


def search_setting(
    G,
    setting_id,
    max_edges,
    max_external,
    sources,
    target_groups,
    seed_genes,
    paths_per_pair,
    enumeration_cap,
):
    group_lookup = target_group_lookup(
        target_groups
    )

    targets = unique_ordered(
        gene
        for genes in target_groups.values()
        for gene in genes
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
                group_lookup[target],
                seed_genes,
                max_edges,
                max_external,
                paths_per_pair,
                enumeration_cap,
            )

            for row in found:
                row["setting_id"] = setting_id
                row["max_edges"] = max_edges
                row["max_external"] = max_external

            rows.extend(found)

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)

    df.insert(
        0,
        "path_id",
        [
            f"{setting_id}_P{i:05d}"
            for i in range(1, len(df) + 1)
        ],
    )

    return df


# ============================================================
# Explicit path-set metrics
# ============================================================

def pair_metrics_for_setting(paths):
    if paths.empty:
        return pd.DataFrame()

    rows = []

    for (
        source,
        target,
        group,
    ), pair_df in paths.groupby(
        [
            "source_anchor",
            "target_anchor",
            "target_group",
        ]
    ):
        path_sets = [
            set(
                str(value).split(";")
            )
            for value in pair_df["path_nodes"]
        ]

        required = set.intersection(
            *path_sets
        )

        # Anchors are always "required" by definition; keep both fields,
        # but distinguish internal required nodes.
        required_internal = (
            required
            - {source, target}
        )

        all_nodes = set.union(
            *path_sets
        )

        rows.append({
            "setting_id": pair_df["setting_id"].iloc[0],
            "source_anchor": source,
            "target_anchor": target,
            "target_group": group,
            "retained_path_count": len(pair_df),
            "required_nodes_including_anchors": ";".join(
                sorted(required)
            ),
            "required_internal_nodes": ";".join(
                sorted(required_internal)
            ),
            "required_internal_count": len(required_internal),
            "all_nodes_on_retained_paths": ";".join(
                sorted(all_nodes)
            ),
        })

    return pd.DataFrame(rows)


def explicit_gene_metrics(paths, seed_genes):
    """
    Per-gene metrics on the EXPLICIT retained path set.

    Two kinds of counts are kept separate:
      * anchor-inclusive counts: a source/target is trivially present on paths
        that start/end at that anchor;
      * internal counts: the gene is neither the source nor the target for the
        source->target pair. These are the relevant connector/bottleneck metrics.
    """
    if paths.empty:
        return pd.DataFrame()

    setting_id = paths["setting_id"].iloc[0]

    pair_groups = list(
        paths.groupby(
            ["source_anchor", "target_anchor"]
        )
    )

    all_genes = sorted(
        {
            node
            for value in paths["path_nodes"]
            for node in str(value).split(";")
        }
    )

    rows = []

    for gene in all_genes:
        path_occurrence = 0
        pair_occurrence = 0
        required_pair_count = 0

        internal_path_occurrence = 0
        internal_pair_occurrence = 0
        required_internal_pair_count = 0
        eligible_internal_pair_count = 0

        for (source, target), pair_df in pair_groups:
            path_sets = [
                set(str(value).split(";"))
                for value in pair_df["path_nodes"]
            ]

            in_flags = [
                gene in path_set
                for path_set in path_sets
            ]

            path_occurrence += sum(in_flags)

            if any(in_flags):
                pair_occurrence += 1

            if all(in_flags):
                required_pair_count += 1

            # Anchor-excluded connector metrics.
            if gene not in {source, target}:
                eligible_internal_pair_count += 1
                internal_path_occurrence += sum(in_flags)

                if any(in_flags):
                    internal_pair_occurrence += 1

                if all(in_flags):
                    required_internal_pair_count += 1

        rows.append({
            "setting_id": setting_id,
            "Gene": gene,
            "is_seed_candidate": gene in seed_genes,
            "path_occurrence_count": path_occurrence,
            "pair_occurrence_count": pair_occurrence,
            "required_pair_count": required_pair_count,
            "internal_path_occurrence_count": internal_path_occurrence,
            "internal_pair_occurrence_count": internal_pair_occurrence,
            "required_internal_pair_count": required_internal_pair_count,
            "eligible_internal_pair_count": eligible_internal_pair_count,
            "path_occurrence_fraction": (
                path_occurrence / len(paths)
                if len(paths)
                else 0.0
            ),
            "pair_occurrence_fraction": (
                pair_occurrence / len(pair_groups)
                if pair_groups
                else 0.0
            ),
            "required_pair_fraction": (
                required_pair_count / len(pair_groups)
                if pair_groups
                else 0.0
            ),
            "internal_pair_occurrence_fraction": (
                internal_pair_occurrence / eligible_internal_pair_count
                if eligible_internal_pair_count
                else 0.0
            ),
            "required_internal_pair_fraction": (
                required_internal_pair_count / eligible_internal_pair_count
                if eligible_internal_pair_count
                else 0.0
            ),
        })

    return pd.DataFrame(rows)


# ============================================================
# Union-graph bottleneck metrics
# ============================================================

def union_graph_from_paths(G, paths):
    H = nx.DiGraph()

    if paths.empty:
        return H

    edge_set = set()

    for value in paths["path_nodes"]:
        nodes = str(value).split(";")

        for source, target in zip(
            nodes[:-1],
            nodes[1:],
        ):
            edge_set.add(
                (source, target)
            )

    for source, target in sorted(edge_set):
        H.add_edge(
            source,
            target,
            **G[source][target],
        )

    return H


def reachable_pairs(G, sources, targets):
    pairs = set()

    for source in sources:
        if source not in G:
            continue

        descendants = nx.descendants(
            G,
            source,
        )

        for target in targets:
            if target in descendants:
                pairs.add(
                    (source, target)
                )

    return pairs


def union_bottlenecks(
    G,
    sources,
    targets,
    seed_genes,
    setting_id,
):
    """
    Node-removal bottlenecks on the UNION graph.

    Anchor-inclusive and anchor-excluded losses are reported separately. The
    latter prevents a source or target from looking like a connector simply
    because deleting an endpoint necessarily destroys its own pairs.
    """
    baseline = reachable_pairs(
        G,
        sources,
        targets,
    )

    rows = []

    for node in G.nodes():
        H = G.copy()
        H.remove_node(node)

        remaining = reachable_pairs(
            H,
            sources,
            targets,
        )

        lost = baseline - remaining

        eligible_nonanchor = {
            pair
            for pair in baseline
            if node not in pair
        }

        lost_nonanchor = {
            pair
            for pair in lost
            if node not in pair
        }

        rows.append({
            "setting_id": setting_id,
            "Gene": node,
            "is_seed_candidate": node in seed_genes,
            "union_graph_reachable_pair_count": len(baseline),
            "union_pairs_lost_on_removal": len(lost),
            "union_pair_loss_fraction": (
                len(lost) / len(baseline)
                if baseline
                else 0.0
            ),
            "union_lost_pairs": ";".join(
                f"{source}->{target}"
                for source, target in sorted(lost)
            ),
            "union_nonanchor_eligible_pair_count": len(eligible_nonanchor),
            "union_nonanchor_pairs_lost_on_removal": len(lost_nonanchor),
            "union_nonanchor_pair_loss_fraction": (
                len(lost_nonanchor) / len(eligible_nonanchor)
                if eligible_nonanchor
                else 0.0
            ),
            "union_nonanchor_lost_pairs": ";".join(
                f"{source}->{target}"
                for source, target in sorted(lost_nonanchor)
            ),
        })

    return pd.DataFrame(rows), baseline


# ============================================================
# Consensus across settings
# ============================================================

def build_consensus(
    gene_by_setting,
    settings,
):
    if gene_by_setting.empty:
        return pd.DataFrame()

    n_settings = len(settings)
    rows = []

    for gene, group in gene_by_setting.groupby("Gene"):
        present_settings = sorted(set(group["setting_id"]))
        is_seed = bool(group["is_seed_candidate"].max())

        rows.append({
            "Gene": gene,
            "is_seed_candidate": is_seed,
            "is_external_intermediate": not is_seed,
            "settings_present_count": len(present_settings),
            "settings_present_fraction": len(present_settings) / n_settings,
            "settings_present": ";".join(present_settings),
            "mean_pair_occurrence_count": group["pair_occurrence_count"].mean(),
            "max_pair_occurrence_count": group["pair_occurrence_count"].max(),
            "mean_internal_pair_occurrence_count": group[
                "internal_pair_occurrence_count"
            ].mean(),
            "max_internal_pair_occurrence_count": group[
                "internal_pair_occurrence_count"
            ].max(),
            "mean_required_pair_count": group["required_pair_count"].mean(),
            "max_required_pair_count": group["required_pair_count"].max(),
            "settings_required_for_at_least_one_pair": int(
                (group["required_pair_count"] > 0).sum()
            ),
            "mean_required_internal_pair_count": group[
                "required_internal_pair_count"
            ].mean(),
            "max_required_internal_pair_count": group[
                "required_internal_pair_count"
            ].max(),
            "settings_required_internal_for_at_least_one_pair": int(
                (group["required_internal_pair_count"] > 0).sum()
            ),
            "mean_union_pairs_lost": group[
                "union_pairs_lost_on_removal"
            ].mean(),
            "max_union_pairs_lost": group[
                "union_pairs_lost_on_removal"
            ].max(),
            "settings_union_bottleneck": int(
                (group["union_pairs_lost_on_removal"] > 0).sum()
            ),
            "mean_union_nonanchor_pairs_lost": group[
                "union_nonanchor_pairs_lost_on_removal"
            ].mean(),
            "max_union_nonanchor_pairs_lost": group[
                "union_nonanchor_pairs_lost_on_removal"
            ].max(),
            "settings_union_nonanchor_bottleneck": int(
                (group["union_nonanchor_pairs_lost_on_removal"] > 0).sum()
            ),
            "mean_path_occurrence_fraction": group[
                "path_occurrence_fraction"
            ].mean(),
        })

    consensus = pd.DataFrame(rows)

    # Connector-oriented ordering: robustness first, then anchor-excluded
    # requiredness and union-graph bottleneck evidence.
    consensus.sort_values(
        [
            "settings_present_count",
            "settings_required_internal_for_at_least_one_pair",
            "max_required_internal_pair_count",
            "settings_union_nonanchor_bottleneck",
            "max_union_nonanchor_pairs_lost",
            "mean_internal_pair_occurrence_count",
        ],
        ascending=[False, False, False, False, False, False],
        inplace=True,
    )

    return consensus


# ============================================================
# Figures
# ============================================================

def plot_setting_stability(
    consensus,
    figures_dir,
):
    if consensus.empty:
        return

    df = (
        consensus
        .sort_values(
            [
                "settings_present_count",
                "settings_required_for_at_least_one_pair",
                "mean_pair_occurrence_count",
            ],
            ascending=[
                False,
                False,
                False,
            ],
        )
        .head(30)
        .copy()
    )

    df = df.iloc[::-1]

    plt.figure(
        figsize=(10, max(6, 0.30 * len(df) + 2))
    )

    y = range(len(df))

    plt.barh(
        y,
        df["settings_present_count"],
    )

    plt.yticks(
        y,
        df["Gene"],
    )

    plt.xlim(
        0,
        len(SETTINGS) + 0.2,
    )

    plt.xlabel(
        "Number of robustness settings containing gene"
    )

    plt.title(
        "Mechanistic-path stability across parameter settings"
    )

    plt.tight_layout()

    plt.savefig(
        figures_dir / "gene_setting_stability.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def plot_paths_per_setting(
    setting_summary,
    figures_dir,
):
    if setting_summary.empty:
        return

    plt.figure(
        figsize=(8, 5)
    )

    plt.bar(
        setting_summary["setting_id"],
        setting_summary["retained_path_count"],
    )

    plt.xlabel("Setting")
    plt.ylabel("Retained paths")
    plt.title(
        "Retained causal paths across robustness settings"
    )
    plt.tight_layout()

    plt.savefig(
        figures_dir / "paths_per_setting.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def plot_reachability(
    pair_all,
    sources,
    targets,
    figures_dir,
):
    matrix = pd.DataFrame(
        0,
        index=sources,
        columns=targets,
        dtype=int,
    )

    if not pair_all.empty:
        for (
            source,
            target,
        ), group in pair_all.groupby(
            [
                "source_anchor",
                "target_anchor",
            ]
        ):
            matrix.loc[
                source,
                target,
            ] = group[
                "setting_id"
            ].nunique()

    plt.figure(
        figsize=(
            max(9, 0.75 * len(targets) + 3),
            max(4, 0.75 * len(sources) + 2),
        )
    )

    image = plt.imshow(
        matrix.values,
        aspect="auto",
        vmin=0,
        vmax=len(SETTINGS),
    )

    plt.colorbar(
        image,
        label="Number of settings with ≥1 retained path",
    )

    plt.xticks(
        range(len(targets)),
        targets,
        rotation=60,
        ha="right",
    )

    plt.yticks(
        range(len(sources)),
        sources,
    )

    plt.xlabel("Endpoint")
    plt.ylabel("Source anchor")
    plt.title(
        "Source→endpoint robustness across path constraints"
    )

    for i in range(len(sources)):
        for j in range(len(targets)):
            plt.text(
                j,
                i,
                str(matrix.iloc[i, j]),
                ha="center",
                va="center",
            )

    plt.tight_layout()

    plt.savefig(
        figures_dir / "source_target_reachability.png",
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Robustness analysis for the Phase C directed mechanistic network."
        )
    )

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
        "--paths-per-pair",
        type=int,
        default=DEFAULT_PATHS_PER_PAIR,
    )

    parser.add_argument(
        "--enumeration-cap",
        type=int,
        default=DEFAULT_ENUMERATION_CAP,
    )

    parser.add_argument(
        "--keep-accession-like",
        action="store_true",
        help=(
            "Keep UniProt-accession-like external nodes in the path search. "
            "By default they are excluded after being audited."
        ),
    )

    args = parser.parse_args()

    project_root = (
        Path(__file__)
        .resolve()
        .parent
        .parent
    )

    mechanistic_dir = (
        project_root
        / "Data"
        / "Processed"
        / "Network"
        / "Mechanistic"
    )

    edge_path = (
        mechanistic_dir
        / "omnipath_directed_edges.csv"
    )

    assignment_path = (
        project_root
        / "Data"
        / "Processed"
        / "Network"
        / "Modules"
        / "module_assignments.csv"
    )

    processed_dir = (
        mechanistic_dir
        / "Robustness"
    )

    figures_dir = (
        project_root
        / "Figures"
        / "Network"
        / "Mechanistic"
        / "Robustness"
    )

    processed_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    figures_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    sources = parse_gene_list(
        args.sources,
        DEFAULT_SOURCES,
    )

    target_groups = {
        "death": parse_gene_list(
            args.death_targets,
            DEFAULT_TARGET_GROUPS["death"],
        ),
        "autophagy": parse_gene_list(
            args.autophagy_targets,
            DEFAULT_TARGET_GROUPS["autophagy"],
        ),
        "trophic_survival": parse_gene_list(
            args.survival_targets,
            DEFAULT_TARGET_GROUPS["trophic_survival"],
        ),
    }

    targets = unique_ordered(
        gene
        for genes in target_groups.values()
        for gene in genes
    )

    edges, assignments = load_inputs(
        edge_path,
        assignment_path,
    )

    seed_genes = set(
        assignments["Gene"]
    )

    accession_audit = audit_accession_like_nodes(
        edges,
        seed_genes,
    )

    accession_audit.to_csv(
        processed_dir / "accession_like_external_nodes.csv",
        index=False,
    )

    excluded_nodes = set()

    if not args.keep_accession_like:
        excluded_nodes = set(
            accession_audit["identifier"]
        )

    G = build_graph(
        edges,
        excluded_nodes=excluded_nodes,
    )

    print()
    print(
        "R-Noel mechanistic robustness"
    )
    print(
        "============================="
    )

    print(
        f"Seed genes: {len(seed_genes)}"
    )

    print(
        f"Directed OmniPath nodes after cleanup: {G.number_of_nodes()}"
    )

    print(
        f"Directed OmniPath edges after cleanup: {G.number_of_edges()}"
    )

    print(
        f"Accession-like external nodes audited: {len(accession_audit)}"
    )

    print(
        "Accession-like nodes excluded from robustness search:"
        f" {len(excluded_nodes)}"
    )

    all_paths = []
    all_pairs = []
    all_gene_setting = []
    setting_rows = []

    for setting_id, max_edges, max_external in SETTINGS:
        print()
        print(
            f"{setting_id}: max_edges={max_edges}, "
            f"max_external={max_external}"
        )

        paths = search_setting(
            G,
            setting_id,
            max_edges,
            max_external,
            sources,
            target_groups,
            seed_genes,
            args.paths_per_pair,
            args.enumeration_cap,
        )

        pair_df = pair_metrics_for_setting(
            paths
        )

        explicit_gene = explicit_gene_metrics(
            paths,
            seed_genes,
        )

        union_G = union_graph_from_paths(
            G,
            paths,
        )

        union_gene, union_pairs = union_bottlenecks(
            union_G,
            sources,
            targets,
            seed_genes,
            setting_id,
        )

        explicit_pairs = set()

        if not pair_df.empty:
            explicit_pairs = set(
                zip(
                    pair_df["source_anchor"],
                    pair_df["target_anchor"],
                )
            )

        if not explicit_gene.empty:
            gene_setting = explicit_gene.merge(
                union_gene[
                    [
                        "Gene",
                        "union_graph_reachable_pair_count",
                        "union_pairs_lost_on_removal",
                        "union_pair_loss_fraction",
                        "union_nonanchor_eligible_pair_count",
                        "union_nonanchor_pairs_lost_on_removal",
                        "union_nonanchor_pair_loss_fraction",
                    ]
                ],
                on="Gene",
                how="left",
            )

            gene_setting[
                "union_graph_reachable_pair_count"
            ] = gene_setting[
                "union_graph_reachable_pair_count"
            ].fillna(
                len(union_pairs)
            )

            gene_setting[
                "union_pairs_lost_on_removal"
            ] = gene_setting[
                "union_pairs_lost_on_removal"
            ].fillna(0)

            gene_setting[
                "union_pair_loss_fraction"
            ] = gene_setting[
                "union_pair_loss_fraction"
            ].fillna(0.0)

            gene_setting[
                "union_nonanchor_eligible_pair_count"
            ] = gene_setting[
                "union_nonanchor_eligible_pair_count"
            ].fillna(0)

            gene_setting[
                "union_nonanchor_pairs_lost_on_removal"
            ] = gene_setting[
                "union_nonanchor_pairs_lost_on_removal"
            ].fillna(0)

            gene_setting[
                "union_nonanchor_pair_loss_fraction"
            ] = gene_setting[
                "union_nonanchor_pair_loss_fraction"
            ].fillna(0.0)

            all_gene_setting.append(
                gene_setting
            )

        if not paths.empty:
            all_paths.append(
                paths
            )

        if not pair_df.empty:
            all_pairs.append(
                pair_df
            )

        setting_rows.append({
            "setting_id": setting_id,
            "max_edges": max_edges,
            "max_external": max_external,
            "retained_path_count": len(paths),
            "explicit_reachable_pair_count": len(explicit_pairs),
            "union_graph_reachable_pair_count": len(union_pairs),
            "union_graph_nodes": union_G.number_of_nodes(),
            "union_graph_edges": union_G.number_of_edges(),
            "accession_like_nodes_excluded": len(excluded_nodes),
        })

        print(
            f"  retained paths: {len(paths)}"
        )

        print(
            f"  explicit source-target pairs: {len(explicit_pairs)}"
        )

        print(
            f"  union-graph reachable pairs: {len(union_pairs)}"
        )

        print(
            f"  union graph: "
            f"{union_G.number_of_nodes()} nodes, "
            f"{union_G.number_of_edges()} edges"
        )

    setting_summary = pd.DataFrame(
        setting_rows
    )

    paths_all = (
        pd.concat(
            all_paths,
            ignore_index=True,
        )
        if all_paths
        else pd.DataFrame()
    )

    pairs_all = (
        pd.concat(
            all_pairs,
            ignore_index=True,
        )
        if all_pairs
        else pd.DataFrame()
    )

    gene_by_setting = (
        pd.concat(
            all_gene_setting,
            ignore_index=True,
        )
        if all_gene_setting
        else pd.DataFrame()
    )

    consensus = build_consensus(
        gene_by_setting,
        SETTINGS,
    )

    seed_consensus = (
        consensus[
            consensus["is_seed_candidate"]
        ].copy()
        if not consensus.empty
        else pd.DataFrame()
    )

    setting_summary.to_csv(
        processed_dir / "robustness_setting_summary.csv",
        index=False,
    )

    paths_all.to_csv(
        processed_dir / "robustness_paths.csv",
        index=False,
    )

    pairs_all.to_csv(
        processed_dir / "robustness_pair_metrics.csv",
        index=False,
    )

    gene_by_setting.to_csv(
        processed_dir / "robustness_gene_by_setting.csv",
        index=False,
    )

    consensus.to_csv(
        processed_dir / "robustness_gene_consensus.csv",
        index=False,
    )

    seed_consensus.to_csv(
        processed_dir / "robustness_seed_consensus.csv",
        index=False,
    )

    diagnostics = {
        "settings": [
            {
                "setting_id": setting_id,
                "max_edges": max_edges,
                "max_external": max_external,
            }
            for setting_id, max_edges, max_external in SETTINGS
        ],
        "sources": sources,
        "target_groups": target_groups,
        "seed_gene_count": len(seed_genes),
        "accession_like_external_nodes_audited": (
            accession_audit["identifier"].tolist()
            if not accession_audit.empty
            else []
        ),
        "exclude_accession_like_default": (
            not args.keep_accession_like
        ),
        "paths_per_pair": args.paths_per_pair,
        "enumeration_cap": args.enumeration_cap,
    }

    write_json(
        processed_dir / "robustness_diagnostics.json",
        diagnostics,
    )

    plot_setting_stability(
        consensus,
        figures_dir,
    )

    plot_paths_per_setting(
        setting_summary,
        figures_dir,
    )

    plot_reachability(
        pairs_all,
        sources,
        targets,
        figures_dir,
    )

    print()
    print(
        "Robust seed candidates"
    )
    print(
        "----------------------"
    )

    if seed_consensus.empty:
        print(
            "No seed candidates occurred on retained paths."
        )
    else:
        display_cols = [
            "Gene",
            "settings_present_count",
            "max_internal_pair_occurrence_count",
            "settings_required_internal_for_at_least_one_pair",
            "max_required_internal_pair_count",
            "settings_union_nonanchor_bottleneck",
            "max_union_nonanchor_pairs_lost",
        ]

        print(
            seed_consensus[
                display_cols
            ]
            .head(25)
            .to_string(
                index=False
            )
        )

    print()
    print(
        "Robust external connectors"
    )
    print(
        "--------------------------"
    )

    if consensus.empty:
        print("None.")
    else:
        external = consensus[
            ~consensus["is_seed_candidate"]
        ]

        print(
            external[
                [
                    "Gene",
                    "settings_present_count",
                    "max_internal_pair_occurrence_count",
                    "settings_required_internal_for_at_least_one_pair",
                    "max_required_internal_pair_count",
                    "settings_union_nonanchor_bottleneck",
                    "max_union_nonanchor_pairs_lost",
                ]
            ]
            .head(25)
            .to_string(
                index=False
            )
        )

    print()
    print(
        "Robustness analysis complete."
    )

    print(
        f"Processed: {processed_dir}"
    )

    print(
        f"Figures: {figures_dir}"
    )

    print()
    print(
        "Primary connector interpretation should use the anchor-excluded fields: "
        "'required_internal_pair_count' and "
        "'union_nonanchor_pairs_lost_on_removal'. Anchor-inclusive fields are "
        "retained for provenance only."
    )


if __name__ == "__main__":
    main()
