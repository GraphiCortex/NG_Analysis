#!/usr/bin/env python3
"""
05b_connector_audit.py

Audit Phase D accepted external connectors before final prioritization.

This script does NOT call STRING or any other web service.
It works entirely from Phase D outputs.

Main goals
----------
1. Distinguish:
      - cross-module bridges
      - mixed module/isolate bridges
      - isolate-rescue connectors
      - causal-only connectors
2. Measure each connector's marginal contribution to seed-network connectivity.
3. Detect redundant connectors that bridge the same original seed components.
4. Separate "many alternative proteins for the same gap" from unique bridges.

Run from Project/
-----------------
    py .\Scripts\05b_connector_audit.py
"""

from __future__ import annotations

import itertools
import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


def write_json(path, obj):
    path.write_text(
        json.dumps(obj, indent=2, ensure_ascii=False, default=str),
        encoding="utf-8",
    )


def clean_gene(value):
    if pd.isna(value):
        return None
    value = str(value).strip().upper()
    return value or None


def build_graph_from_expanded_edges(edges, all_nodes):
    G = nx.Graph()

    for gene in all_nodes:
        G.add_node(gene)

    for _, row in edges.iterrows():
        a = clean_gene(row.get("preferredName_A"))
        b = clean_gene(row.get("preferredName_B"))

        if a is None or b is None or a == b:
            continue

        G.add_edge(
            a,
            b,
            score=float(row.get("score", 0.0)),
        )

    return G


def seed_component_stats(G, seed_genes):
    sizes = []

    for comp in nx.connected_components(G):
        n_seed = len(set(comp).intersection(seed_genes))
        if n_seed:
            sizes.append(n_seed)

    sizes.sort(reverse=True)

    return {
        "component_count": len(sizes),
        "largest_seed_component": max(sizes, default=0),
        "sizes": sizes,
    }


def baseline_components_from_assignments(assignments):
    """
    Phase B already stored the original seed-only connected-component label.
    Prefer component_id_phaseB if available.
    """
    if "component_id_phaseB" in assignments.columns:
        values = assignments[
            ["Gene", "component_id_phaseB"]
        ].copy()

        values["component_id_phaseB"] = pd.to_numeric(
            values["component_id_phaseB"],
            errors="coerce",
        )

        return (
            values.dropna(subset=["component_id_phaseB"])
            .assign(
                component_id_phaseB=lambda x:
                    x["component_id_phaseB"].astype(int)
            )
            .set_index("Gene")["component_id_phaseB"]
            .to_dict()
        )

    raise ValueError(
        "module_assignments.csv has no component_id_phaseB column."
    )


def classify_connector(
    robust_causal,
    seed_neighbor_count,
    isolate_seed_neighbors,
    nonisolate_seed_neighbors,
    nonisolate_module_count,
    component_count,
):
    if robust_causal and seed_neighbor_count == 0:
        return "causal_only"

    if nonisolate_module_count >= 2:
        return "cross_module_bridge"

    if (
        nonisolate_seed_neighbors >= 1
        and isolate_seed_neighbors >= 1
        and component_count >= 2
    ):
        return "mixed_module_isolate_bridge"

    if (
        nonisolate_seed_neighbors == 0
        and isolate_seed_neighbors >= 2
        and component_count >= 2
    ):
        return "isolate_rescue"

    if robust_causal:
        return "causal_plus_string"

    if nonisolate_module_count == 1:
        return "local_module_neighbor"

    return "other"


def component_pairs(component_ids):
    ids = sorted(set(component_ids))

    return {
        tuple(pair)
        for pair in itertools.combinations(ids, 2)
    }


def main():
    project_root = Path(__file__).resolve().parent.parent

    network_dir = (
        project_root / "Data" / "Processed" / "Network"
    )

    module_path = (
        network_dir / "Modules" / "module_assignments.csv"
    )

    expansion_dir = (
        network_dir / "GlobalExpansion"
    )

    accepted_path = (
        expansion_dir / "accepted_external_connectors.csv"
    )

    edges_path = (
        expansion_dir / "expanded_edges.csv"
    )

    output_dir = (
        expansion_dir / "ConnectorAudit"
    )

    figures_dir = (
        project_root
        / "Figures"
        / "Network"
        / "GlobalExpansion"
        / "ConnectorAudit"
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    assignments = pd.read_csv(module_path)
    accepted = pd.read_csv(accepted_path)
    expanded_edges = pd.read_csv(edges_path)

    assignments["Gene"] = assignments["Gene"].map(clean_gene)
    accepted["Gene"] = accepted["Gene"].map(clean_gene)

    assignments.dropna(subset=["Gene"], inplace=True)
    accepted.dropna(subset=["Gene"], inplace=True)

    seed_genes = set(assignments["Gene"])
    accepted_genes = set(accepted["Gene"])

    module_of = (
        assignments
        .set_index("Gene")["module_id"]
        .to_dict()
    )

    baseline_component_of = baseline_components_from_assignments(
        assignments
    )

    all_nodes = seed_genes | accepted_genes

    G = build_graph_from_expanded_edges(
        expanded_edges,
        all_nodes,
    )

    full_stats = seed_component_stats(
        G,
        seed_genes,
    )

    accepted_lookup = (
        accepted
        .drop_duplicates("Gene")
        .set_index("Gene")
    )

    # --------------------------------------------------------
    # First pass: connector neighborhoods and component pairs
    # --------------------------------------------------------

    connector_info = {}
    pair_to_connectors = defaultdict(set)

    for gene in sorted(accepted_genes):
        if gene not in G:
            seed_neighbors = []
        else:
            seed_neighbors = sorted(
                nbr
                for nbr in G.neighbors(gene)
                if nbr in seed_genes
            )

        isolate_neighbors = [
            nbr
            for nbr in seed_neighbors
            if module_of.get(nbr) == "ISOLATE"
        ]

        nonisolate_neighbors = [
            nbr
            for nbr in seed_neighbors
            if module_of.get(nbr) != "ISOLATE"
        ]

        nonisolate_modules = sorted(
            {
                module_of[nbr]
                for nbr in nonisolate_neighbors
                if module_of.get(nbr) is not None
            }
        )

        baseline_components = sorted(
            {
                baseline_component_of[nbr]
                for nbr in seed_neighbors
                if nbr in baseline_component_of
            }
        )

        pairs = component_pairs(
            baseline_components
        )

        for pair in pairs:
            pair_to_connectors[pair].add(gene)

        robust_causal = False

        if gene in accepted_lookup.index:
            value = accepted_lookup.loc[
                gene
            ].get(
                "rule_causal_core",
                False,
            )

            if pd.notna(value):
                if isinstance(value, str):
                    robust_causal = (
                        value.strip().lower()
                        in {"true", "1", "yes"}
                    )
                else:
                    robust_causal = bool(value)

        connector_info[gene] = {
            "Gene": gene,
            "robust_causal_core": robust_causal,
            "seed_neighbor_count": len(seed_neighbors),
            "seed_neighbors": ";".join(seed_neighbors),
            "isolate_seed_neighbor_count": len(isolate_neighbors),
            "isolate_seed_neighbors": ";".join(isolate_neighbors),
            "nonisolate_seed_neighbor_count": len(nonisolate_neighbors),
            "nonisolate_seed_neighbors": ";".join(nonisolate_neighbors),
            "nonisolate_module_count": len(nonisolate_modules),
            "nonisolate_modules": ";".join(nonisolate_modules),
            "baseline_seed_component_count": len(baseline_components),
            "baseline_seed_components": ";".join(
                str(x) for x in baseline_components
            ),
            "component_pair_count": len(pairs),
            "_pairs": pairs,
        }

    # --------------------------------------------------------
    # Pair redundancy
    # --------------------------------------------------------

    pair_rows = []

    for pair, genes in sorted(pair_to_connectors.items()):
        pair_rows.append({
            "component_A": pair[0],
            "component_B": pair[1],
            "connector_support_count": len(genes),
            "supporting_connectors": ";".join(sorted(genes)),
            "uniquely_supported": len(genes) == 1,
        })

    pair_df = pd.DataFrame(pair_rows)

    if not pair_df.empty:
        pair_df.sort_values(
            [
                "connector_support_count",
                "component_A",
                "component_B",
            ],
            ascending=[
                False,
                True,
                True,
            ],
            inplace=True,
        )

    pair_df.to_csv(
        output_dir / "component_pair_support.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Marginal removal + unique pair support
    # --------------------------------------------------------

    audit_rows = []

    for gene in sorted(accepted_genes):
        info = connector_info[gene]

        H = G.copy()

        if gene in H:
            H.remove_node(gene)

        removed_stats = seed_component_stats(
            H,
            seed_genes,
        )

        marginal_component_increase = (
            removed_stats["component_count"]
            - full_stats["component_count"]
        )

        largest_component_drop = (
            full_stats["largest_seed_component"]
            - removed_stats["largest_seed_component"]
        )

        unique_pairs = [
            pair
            for pair in info["_pairs"]
            if len(
                pair_to_connectors[pair]
            ) == 1
        ]

        shared_pairs = [
            pair
            for pair in info["_pairs"]
            if len(
                pair_to_connectors[pair]
            ) > 1
        ]

        connector_class = classify_connector(
            info["robust_causal_core"],
            info["seed_neighbor_count"],
            info["isolate_seed_neighbor_count"],
            info["nonisolate_seed_neighbor_count"],
            info["nonisolate_module_count"],
            info["baseline_seed_component_count"],
        )

        # A deliberately descriptive promotion flag:
        # this is NOT a biological rank.
        structurally_distinct = (
            connector_class
            in {
                "cross_module_bridge",
                "mixed_module_isolate_bridge",
            }
            or len(unique_pairs) > 0
            or marginal_component_increase > 0
            or info["robust_causal_core"]
        )

        row = {
            key: value
            for key, value in info.items()
            if key != "_pairs"
        }

        row.update({
            "connector_class": connector_class,
            "unique_component_pair_count": len(unique_pairs),
            "shared_component_pair_count": len(shared_pairs),
            "unique_component_pairs": ";".join(
                f"{a}-{b}"
                for a, b in sorted(unique_pairs)
            ),
            "marginal_seed_component_increase_on_removal":
                marginal_component_increase,
            "largest_seed_component_drop_on_removal":
                largest_component_drop,
            "structurally_distinct":
                structurally_distinct,
        })

        # Preserve useful Phase D provenance.
        if gene in accepted_lookup.index:
            source = accepted_lookup.loc[gene]

            for col in [
                "acceptance_reason",
                "query_support_count",
                "supporting_queries",
                "supporting_query_types",
                "max_string_score",
                "max_non_text_evidence",
                "text_mining_dominant_edge_fraction",
                "settings_union_nonanchor_bottleneck",
                "max_union_nonanchor_pairs_lost",
            ]:
                if col in source.index:
                    row[col] = source[col]

        audit_rows.append(row)

    audit = pd.DataFrame(audit_rows)

    audit.sort_values(
        [
            "structurally_distinct",
            "nonisolate_module_count",
            "unique_component_pair_count",
            "marginal_seed_component_increase_on_removal",
            "baseline_seed_component_count",
            "seed_neighbor_count",
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

    audit.to_csv(
        output_dir / "connector_audit.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Class summary
    # --------------------------------------------------------

    class_summary = (
        audit.groupby(
            "connector_class",
            dropna=False,
        )
        .agg(
            connector_count=("Gene", "count"),
            structurally_distinct_count=(
                "structurally_distinct",
                "sum",
            ),
            mean_seed_neighbors=(
                "seed_neighbor_count",
                "mean",
            ),
            mean_components_touched=(
                "baseline_seed_component_count",
                "mean",
            ),
            total_unique_component_pairs=(
                "unique_component_pair_count",
                "sum",
            ),
        )
        .reset_index()
    )

    class_summary.to_csv(
        output_dir / "connector_class_summary.csv",
        index=False,
    )

    # --------------------------------------------------------
    # Figures
    # --------------------------------------------------------

    if not audit.empty:
        plot_df = audit.copy()

        plot_df.sort_values(
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

        plot_df = plot_df.head(25).iloc[::-1]

        plt.figure(
            figsize=(10, max(6, 0.34 * len(plot_df) + 2))
        )

        y = range(len(plot_df))

        plt.barh(
            y,
            plot_df["baseline_seed_component_count"],
        )

        plt.yticks(
            y,
            plot_df["Gene"],
        )

        plt.xlabel(
            "Original seed components touched"
        )

        plt.title(
            "Accepted connector reach after final fixed-network query"
        )

        plt.tight_layout()

        plt.savefig(
            figures_dir / "connector_component_reach.png",
            dpi=300,
            bbox_inches="tight",
        )

        plt.close()

        marginal_df = audit.sort_values(
            [
                "marginal_seed_component_increase_on_removal",
                "largest_seed_component_drop_on_removal",
            ],
            ascending=[
                False,
                False,
            ],
        ).head(25).iloc[::-1]

        plt.figure(
            figsize=(10, max(6, 0.34 * len(marginal_df) + 2))
        )

        y = range(len(marginal_df))

        plt.barh(
            y,
            marginal_df[
                "marginal_seed_component_increase_on_removal"
            ],
        )

        plt.yticks(
            y,
            marginal_df["Gene"],
        )

        plt.xlabel(
            "Increase in seed components after connector removal"
        )

        plt.title(
            "Marginal structural contribution of accepted connectors"
        )

        plt.tight_layout()

        plt.savefig(
            figures_dir / "connector_marginal_compression.png",
            dpi=300,
            bbox_inches="tight",
        )

        plt.close()

    diagnostics = {
        "full_expanded_seed_component_count":
            full_stats["component_count"],
        "full_expanded_largest_seed_component":
            full_stats["largest_seed_component"],
        "accepted_connector_count":
            len(accepted_genes),
        "connector_classes":
            (
                class_summary
                .set_index("connector_class")[
                    "connector_count"
                ]
                .to_dict()
                if not class_summary.empty
                else {}
            ),
        "component_pairs_supported":
            len(pair_df),
        "component_pairs_with_single_connector":
            (
                int(
                    pair_df["uniquely_supported"].sum()
                )
                if not pair_df.empty
                else 0
            ),
    }

    write_json(
        output_dir / "connector_audit_diagnostics.json",
        diagnostics,
    )

    # --------------------------------------------------------
    # Console
    # --------------------------------------------------------

    print()
    print("R-Noel Phase D connector audit")
    print("==============================")

    print(
        f"Accepted connectors: {len(accepted_genes)}"
    )

    print(
        "Expanded seed components:",
        full_stats["component_count"],
    )

    print(
        "Largest expanded seed component:",
        full_stats["largest_seed_component"],
    )

    print()
    print("Connector classes")
    print("-----------------")

    if class_summary.empty:
        print("None.")
    else:
        print(
            class_summary.to_string(
                index=False
            )
        )

    print()
    print("Connector audit")
    print("---------------")

    display_cols = [
        "Gene",
        "connector_class",
        "seed_neighbor_count",
        "isolate_seed_neighbor_count",
        "nonisolate_seed_neighbor_count",
        "nonisolate_module_count",
        "baseline_seed_component_count",
        "unique_component_pair_count",
        "marginal_seed_component_increase_on_removal",
        "largest_seed_component_drop_on_removal",
        "structurally_distinct",
    ]

    print(
        audit[
            display_cols
        ]
        .head(30)
        .to_string(
            index=False
        )
    )

    print()
    print(
        "Audit complete."
    )

    print(
        f"Processed: {output_dir}"
    )

    print(
        f"Figures: {figures_dir}"
    )

    print()
    print(
        "Interpret 'isolate_rescue' separately from a true "
        "cross-module bridge: multiple alternative proteins can "
        "connect the same pair of previously isolated seed genes."
    )


if __name__ == "__main__":
    main()
