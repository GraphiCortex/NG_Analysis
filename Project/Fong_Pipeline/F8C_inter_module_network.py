
#!/usr/bin/env python3
"""
F8C — Functional Intermodule Connectivity

Reconstruct the original F3 STRING networks, project gene-gene
associations onto the overlapping F8B biological modules,
and identify CORE-INTERFACE structural bridges.

All STRING edges are UNDIRECTED functional associations.

Requirements: Python standard library only.
Previous pipeline results are never modified.
"""

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


# ============================================================
# 1. PATHS
# ============================================================

PROJECT = Path(__file__).resolve().parents[1]

F3 = (
    PROJECT / "Data/Processed/FongDownstream/F3_STRING"
)

F8 = (
    PROJECT
    / "Data/Processed/FongDownstream/F8_FunctionalAtlas"
)

F8B = F8 / "F8B_BiologicalModules"
OUT = F8 / "F8C_IntermoduleConnectivity"

# Display placeholder, NOT a twentieth biological mechanism.
UNRESOLVED = "UNRESOLVED"


# ============================================================
# 2. UTILITIES
# ============================================================

def read_csv(path):
    with path.open(
        newline="", encoding="utf-8-sig"
    ) as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows, fields):
    with path.open(
        "w", newline="", encoding="utf-8-sig"
    ) as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=fields,
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def boolean(value):
    return str(value).strip().lower() in {
        "true", "1", "yes"
    }


def gene_pair(row):
    a = row["Gene_A"].strip().upper()
    b = row["Gene_B"].strip().upper()

    if a == b:
        raise ValueError(f"Self-edge: {a}")

    return tuple(sorted((a, b)))


def collect_edges(rows, allowed_genes):
    """
    Normalize undirected STRING associations.

    Enforce unique gene pairs and required confidence >= 700.
    """
    unique = {}

    for row in rows:
        key = gene_pair(row)

        if not set(key).issubset(allowed_genes):
            raise ValueError(
                f"Unmapped STRING endpoints: {key}"
            )

        if key in unique:
            raise ValueError(
                f"Duplicate undirected edge: {key}"
            )

        score = float(row["score_1000"])

        if not (700 <= score <= 1000):
            raise ValueError(
                f"Unexpected STRING score: {key} = {score}"
            )

        unique[key] = row

    return unique


def component_sizes(node_set, edge_keys):
    """
    Reconstruct connected components without NetworkX.
    Include isolated nodes.
    """
    neighbors = {
        gene: set() for gene in node_set
    }

    for a, b in edge_keys:
        neighbors[a].add(b)
        neighbors[b].add(a)

    seen = set()
    sizes = []

    for gene in sorted(node_set):

        if gene in seen:
            continue

        pending = [gene]
        seen.add(gene)
        size = 0

        while pending:
            current = pending.pop()
            size += 1

            for neighbor in neighbors[current] - seen:
                seen.add(neighbor)
                pending.append(neighbor)

        sizes.append(size)

    return sizes, neighbors


# ============================================================
# 3. LOAD AND VALIDATE
# ============================================================

def main():

    atlas_rows = read_csv(
        F8B / "F8B_gene_role_atlas.csv"
    )

    summary_rows = read_csv(
        F8B / "F8B_module_summary.csv"
    )

    core_nodes = read_csv(
        F3 / "CORE/nodes.csv"
    )

    extended_nodes = read_csv(
        F3 / "EXTENDED/nodes.csv"
    )

    core_edge_rows = read_csv(
        F3 / "CORE/edges.csv"
    )

    extended_edge_rows = read_csv(
        F3 / "EXTENDED/edges.csv"
    )

    metadata = json.loads(
        (F3 / "F3_metadata.json").read_text(
            encoding="utf-8"
        )
    )

    atlas = {
        r["Gene"].strip().upper(): r
        for r in atlas_rows
    }

    core_node_map = {
        r["Gene"].strip().upper(): r
        for r in core_nodes
    }

    extended_node_map = {
        r["Gene"].strip().upper(): r
        for r in extended_nodes
    }

    if (
        len(atlas) != 321
        or len(extended_node_map) != 321
        or len(core_node_map) != 194
    ):
        raise ValueError(
            "Incorrect node counts or duplicate gene IDs"
        )

    if set(atlas) != set(extended_node_map):
        raise ValueError(
            "F8B and F3 gene universes differ"
        )

    tier = {
        gene: row["Tier"]
        for gene, row in atlas.items()
    }

    assert Counter(tier.values()) == {
        "CORE": 194,
        "INTERFACE": 127
    }

    assert set(core_node_map) == {
        g for g in atlas if tier[g] == "CORE"
    }

    assert sum(
        boolean(r["Shared37"]) for r in atlas_rows
    ) == 37

    catalog = {
        r["Module_id"]: r
        for r in summary_rows
    }

    assert len(catalog) == 19

    catalog[UNRESOLVED] = {
        "Module_id": UNRESOLVED,
        "Module_label": "Role unassigned from F2/F1 evidence",
        "Family": "UNASSIGNED",
    }

    modules = {}

    for gene, row in atlas.items():

        assigned = tuple(sorted(
            x.strip()
            for x in row["Biological_modules"].split(";")
            if x.strip()
        ))

        if (
            len(assigned) != len(set(assigned))
            or any(x not in catalog for x in assigned)
        ):
            raise ValueError(
                f"Invalid module membership: {gene}"
            )

        modules[gene] = assigned or (UNRESOLVED,)

    assert sum(
        modules[g] == (UNRESOLVED,)
        for g in modules
    ) == 29


    # ========================================================
    # 4. RECONSTRUCT THE ORIGINAL STRING NETWORKS
    # ========================================================

    original_core_isolates = {
        gene
        for gene, node in core_node_map.items()
        if boolean(node["Isolate"])
    }

    assert len(original_core_isolates) == 70

    core_edges = collect_edges(
        core_edge_rows,
        set(core_node_map)
    )

    extended_edges = collect_edges(
        extended_edge_rows,
        set(atlas)
    )

    assert len(core_edges) == 359
    assert len(extended_edges) == 1139

    induced_core = {
        key
        for key in extended_edges
        if tier[key[0]] == tier[key[1]] == "CORE"
    }

    if induced_core != set(core_edges):
        raise ValueError(
            "CORE edges are not preserved in EXTENDED"
        )

    for key in core_edges:

        old_score = float(
            core_edges[key]["score_1000"]
        )

        new_score = float(
            extended_edges[key]["score_1000"]
        )

        if abs(old_score - new_score) > 1e-8:
            raise ValueError(
                f"Changed STRING score: {key}"
            )

    core_sizes, _ = component_sizes(
        set(core_node_map),
        core_edges
    )

    ext_sizes, neighbors = component_sizes(
        set(atlas),
        extended_edges
    )

    assert max(core_sizes) == 117
    assert max(ext_sizes) == 219

    assert len(core_sizes) == 74
    assert len(ext_sizes) == 98

    for gene, node in extended_node_map.items():

        if int(node["Degree"]) != len(neighbors[gene]):
            raise ValueError(
                f"F3 degree mismatch: {gene}"
            )

    print("STRING network reconstruction passed.")


    # ========================================================
    # 5. RECONSTRUCT CORE-INTERFACE CONNECTIONS
    # ========================================================

    edge_partition = Counter()

    interface_core_neighbors = defaultdict(set)

    interface_original_isolate_neighbors = defaultdict(set)

    interface_core_communities = defaultdict(set)

    interface_core_biological_roles = defaultdict(set)

    gene_edges = []

    projection_rows = []

    projected_counts = Counter()
    projected_weight = defaultdict(float)

    projected_single_label_counts = Counter()
    projected_gene_edges = defaultdict(set)


    # ========================================================
    # 6. PROJECT GENE EDGES ONTO BIOLOGICAL MODULES
    # ========================================================

    for (a, b), source_row in sorted(
        extended_edges.items()
    ):

        pair_class = "-".join(sorted((
            tier[a], tier[b]
        )))

        edge_partition[pair_class] += 1

        if pair_class == "CORE-INTERFACE":

            core = a if tier[a] == "CORE" else b
            inter = b if core == a else a

            interface_core_neighbors[inter].add(core)

            if core in original_core_isolates:

                interface_original_isolate_neighbors[
                    inter
                ].add(core)

            else:

                interface_core_communities[inter].add(
                    atlas[core]["Original_CORE_community"]
                )

            interface_core_biological_roles[inter].update(
                modules[core]
            )

        raw = {

            "Gene_A": a,
            "Gene_B": b,

            "STRING_ID_A":
                extended_node_map[a]["STRING_ID"],

            "STRING_ID_B":
                extended_node_map[b]["STRING_ID"],

            "Edge_class": pair_class,

            "STRING_score_1000":
                int(float(source_row["score_1000"])),

            "STRING_combined_score":
                source_row["combined_score"],

            "F8B_modules_A": "; ".join(modules[a]),
            "F8B_modules_B": "; ".join(modules[b]),

            "Shared37_A":
                boolean(atlas[a]["Shared37"]),

            "Shared37_B":
                boolean(atlas[b]["Shared37"]),

            "STRING_evidence_type":
                "UNDIRECTED_FUNCTIONAL_ASSOCIATION",
        }

        gene_edges.append(raw)


        # ----------------------------------------------------
        # Fractional projection
        # ----------------------------------------------------
        #
        # A gene may belong to multiple biological modules.
        #
        # Each STRING edge contributes TOTAL mass 1.
        # Its mass is distributed equally over endpoint-role
        # combinations.
        #
        # This avoids treating one edge as multiple independent
        # STRING observations.

        allocation = defaultdict(float)

        amount = 1.0 / (
            len(modules[a]) * len(modules[b])
        )

        for left in modules[a]:
            for right in modules[b]:

                key = tuple(sorted((left, right)))

                allocation[key] += amount

        if abs(sum(allocation.values()) - 1.0) > 1e-10:
            raise ValueError(
                f"Edge mass not conserved: {a}-{b}"
            )

        for (m1, m2), weight in sorted(
            allocation.items()
        ):

            key = (m1, m2)

            projected_counts[key] += 1
            projected_weight[key] += weight

            projected_gene_edges[key].add((a, b))

            if (
                len(modules[a]) == 1
                and len(modules[b]) == 1
            ):
                projected_single_label_counts[key] += 1

            projection_rows.append({

                "Gene_A": a,
                "Gene_B": b,

                "Module_A": m1,
                "Module_B": m2,

                "Module_pair_type":
                    "WITHIN" if m1 == m2 else "BETWEEN",

                "Fractional_edge_weight":
                    f"{weight:.12f}",

                "STRING_score_1000":
                    raw["STRING_score_1000"],

                "Edge_class": pair_class,

                "Interpretation":
                    "Undirected STRING edge projected onto F8B roles",
            })


    # ========================================================
    # 7. VALIDATE EDGE PARTITION AND CONSERVATION
    # ========================================================

    expected = {
        "CORE-CORE": 359,
        "CORE-INTERFACE": 405,
        "INTERFACE-INTERFACE": 375,
    }

    if dict(edge_partition) != expected:
        raise ValueError(
            f"Unexpected edge partition: {edge_partition}"
        )

    total_projected_weight = sum(
        projected_weight.values()
    )

    assert abs(
        total_projected_weight - 1139
    ) < 1e-7


    # ========================================================
    # 8. INTERFACE BRIDGE AUDIT
    # ========================================================

    bridges = []
    reconciliations = []

    for gene in sorted(
        g for g in atlas
        if tier[g] == "INTERFACE"
    ):

        observed = len(
            interface_core_neighbors[gene]
        )

        isolate_links = len(
            interface_original_isolate_neighbors[gene]
        )

        nonisolate_count = observed - isolate_links

        f4_neighbors = int(float(
            atlas[gene]["CORE_neighbors"] or 0
        ))

        f4_modules = int(float(
            atlas[gene]["CORE_modules_touched"] or 0
        ))

        # F4 counts exactly reproduce F3 counts after removing
        # CORE neighbors that were isolates in the CORE network.

        if nonisolate_count != f4_neighbors:
            raise ValueError(
                f"F3/F4 CORE neighbor discrepancy: {gene}"
            )

        if (
            len(interface_core_communities[gene])
            != f4_modules
        ):
            raise ValueError(
                f"F3/F4 CORE community discrepancy: {gene}"
            )

        delta = observed - f4_neighbors

        if delta:

            reconciliations.append({

                "Gene": gene,
                "F3": observed,
                "F4": f4_neighbors,
                "Difference": delta,

                "Linked_original_CORE_isolates":
                    "; ".join(sorted(
                        interface_original_isolate_neighbors[
                            gene
                        ]
                    )),
            })

        bridges.append({

            "Gene": gene,

            "F3_direct_CORE_neighbors": observed,

            "F3_CORE_neighbor_genes":
                "; ".join(sorted(
                    interface_core_neighbors[gene]
                )),

            "F3_original_CORE_isolate_neighbors":
                isolate_links,

            "F3_original_CORE_isolate_genes":
                "; ".join(sorted(
                    interface_original_isolate_neighbors[
                        gene
                    ]
                )),

            "F3_nonisolate_CORE_neighbors":
                nonisolate_count,

            "F3_distinct_CORE_community_labels":
                len(interface_core_communities[gene]),

            "F3_CORE_community_labels":
                "; ".join(sorted(
                    interface_core_communities[gene]
                )),

            "F3_distinct_CORE_biological_modules":
                len(interface_core_biological_roles[gene]),

            "F3_CORE_biological_modules":
                "; ".join(sorted(
                    interface_core_biological_roles[gene]
                )),

            "F3_total_degree":
                len(neighbors[gene]),

            "F4_reported_CORE_neighbors":
                f4_neighbors,

            "F4_reported_CORE_modules_touched":
                f4_modules,

            "F3_minus_F4_CORE_neighbors":
                delta,

            "F8B_flagged_F4_cross_module_bridge":
                boolean(
                    atlas[gene]["Cross_module_bridge"]
                ),

            "F8B_roles":
                "; ".join(modules[gene]),

            "Direction":
                atlas[gene]["Direction"],

            "log2FC":
                atlas[gene]["log2FC"],
        })

    bridges.sort(
        key=lambda r: (
            -r["F3_direct_CORE_neighbors"],
            -r["F3_distinct_CORE_community_labels"],
            r["Gene"]
        )
    )


    # ========================================================
    # 9. BUILD VISUALIZATION-READY GENE NODES
    # ========================================================

    network_nodes = []

    for gene in sorted(atlas):

        row = atlas[gene]
        node = extended_node_map[gene]

        network_nodes.append({

            "id": gene,
            "label": gene,

            "STRING_ID": node["STRING_ID"],

            "tier": tier[gene],

            "modules": list(modules[gene]),

            "direction": row["Direction"],

            "log2FC": float(row["log2FC"]),
            "padj": float(row["padj"]),

            "degree": len(neighbors[gene]),

            "shared37":
                boolean(row["Shared37"]),

            "promoter_evidence":
                boolean(row["Promoter_evidence"]),

            "f4_cross_module_bridge":
                boolean(row["Cross_module_bridge"]),

            "f3_core_neighbor_count": (
                len(interface_core_neighbors[gene])
                if tier[gene] == "INTERFACE"
                else 0
            ),

            "f3_original_core_isolate_neighbors": (
                len(
                    interface_original_isolate_neighbors[gene]
                )
                if tier[gene] == "INTERFACE"
                else 0
            ),

            "f4_core_neighbor_count":
                int(float(row["CORE_neighbors"] or 0)),

            "original_extended_community":
                row["Original_EXT_community"],
        })


    # ========================================================
    # 10. AGGREGATE MODULE-LEVEL CONNECTIONS
    # ========================================================

    module_edges = []

    for (m1, m2) in sorted(projected_counts):

        module_edges.append({

            "Module_A": m1,
            "Module_B": m2,

            "Module_label_A":
                catalog[m1]["Module_label"],

            "Module_label_B":
                catalog[m2]["Module_label"],

            "Module_pair_type":
                "WITHIN" if m1 == m2 else "BETWEEN",

            "Unique_STRING_gene_edges":
                len(projected_gene_edges[(m1, m2)]),

            "Fractional_STRING_edge_weight":
                round(projected_weight[(m1, m2)], 9),

            "Unambiguous_one_role_per_endpoint_edges":
                projected_single_label_counts[(m1, m2)],
        })

    module_edges.sort(
        key=lambda r: (
            -r["Fractional_STRING_edge_weight"],
            r["Module_A"],
            r["Module_B"]
        )
    )


    # ========================================================
    # 11. SUMMARY METADATA
    # ========================================================

    report = {

        "phase":
            "F8C_inter_module_connectivity",

        "STRING_version":
            metadata["string_version"],

        "confidence_threshold_score_1000":
            metadata["required_score"],

        "core_genes": 194,
        "extended_genes": 321,

        "core_edges": len(core_edges),
        "extended_edges": len(extended_edges),

        "edge_partition": dict(edge_partition),

        "CORE_network_edges_exactly_retained_in_EXTENDED":
            True,

        "biological_modules": 19,

        "unresolved_pseudo_module":
            UNRESOLVED,

        "unresolved_genes": sum(
            modules[g] == (UNRESOLVED,)
            for g in modules
        ),

        "F3_interface_genes_with_CORE_neighbor":
            sum(
                bool(interface_core_neighbors[g])
                for g in tier
                if tier[g] == "INTERFACE"
            ),

        "F4_interface_genes_with_CORE_neighbor":
            sum(
                int(float(
                    atlas[g]["CORE_neighbors"] or 0
                )) > 0
                for g in tier
                if tier[g] == "INTERFACE"
            ),

        "Original_CORE_isolates":
            len(original_core_isolates),

        "F3_to_original_CORE_isolate_edges":
            sum(
                len(v)
                for v in
                interface_original_isolate_neighbors.values()
            ),

        "F3_F4_isolate_reconciliations":
            reconciliations,

        "F8B_shared37": sum(
            boolean(r["Shared37"])
            for r in atlas_rows
        ),

        "projected_module_pairs":
            len(module_edges),

        "projected_fractional_edge_weight_sum":
            round(total_projected_weight, 7),

        "warnings": [

            "STRING associations are undirected.",

            "STRING connectivity is not causal evidence.",

            "F8B biological roles overlap.",

            "Fractional allocation prevents edge-count inflation.",

            "UNRESOLVED is not a biological mechanism.",

            "Former CORE isolates explain F3/F4 count differences.",

            "External directed-prior X nodes are not added in F8C.",
        ],
    }


    # ========================================================
    # 12. EXPORT RESULTS
    # ========================================================

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    write_csv(
        OUT / "F8C_STRING_gene_edges.csv",
        gene_edges,
        list(gene_edges[0])
    )

    write_csv(
        OUT / "F8C_module_edge_contributions.csv",
        projection_rows,
        list(projection_rows[0])
    )

    write_csv(
        OUT / "F8C_module_connectivity.csv",
        module_edges,
        list(module_edges[0])
    )

    write_csv(
        OUT / "F8C_interface_bridge_audit.csv",
        bridges,
        list(bridges[0])
    )

    write_csv(
        OUT / "F8C_F3_F4_isolate_reconciliation.csv",
        reconciliations,
        [
            "Gene",
            "F3",
            "F4",
            "Difference",
            "Linked_original_CORE_isolates"
        ]
    )

    (
        OUT / "F8C_summary.json"
    ).write_text(
        json.dumps(report, indent=2),
        encoding="utf-8"
    )


    # ========================================================
    # 13. INTERACTIVE-NETWORK DATA EXPORT
    # ========================================================

    graph = {

        "metadata": {

            "source":
                "F3 STRING v12.5; F8B functional annotations",

            "edge_interpretation":
                "undirected functional association",
        },

        "nodes": network_nodes,

        "edges": [

            {
                "id":
                    f"{r['Gene_A']}--{r['Gene_B']}",

                "source": r["Gene_A"],
                "target": r["Gene_B"],

                "score_1000":
                    r["STRING_score_1000"],

                "tier_edge_class":
                    r["Edge_class"],
            }

            for r in gene_edges
        ],

        "modules": [

            {
                "id": r["Module_id"],
                "label": r["Module_label"],
                "family": r["Family"],

                "pseudo":
                    r["Module_id"] == UNRESOLVED,
            }

            for r in catalog.values()
        ],

        "module_edges": module_edges,
    }

    (
        OUT / "F8C_network.json"
    ).write_text(
        json.dumps(graph, indent=2),
        encoding="utf-8"
    )


    # ========================================================
    # 14. FINAL TERMINAL REPORT
    # ========================================================

    print()
    print("F8C completed successfully.")

    print(
        "CORE/EXTENDED nodes:",
        len(core_node_map),
        len(extended_node_map)
    )

    print(
        "CORE/EXTENDED edges:",
        len(core_edges),
        len(extended_edges)
    )

    print(
        "Edge partition:",
        dict(edge_partition)
    )

    print(
        "F3/F4 INTERFACE genes contacting CORE:",
        report["F3_interface_genes_with_CORE_neighbor"],
        report["F4_interface_genes_with_CORE_neighbor"]
    )

    print(
        "CORE-isolate reconciliations:",
        len(reconciliations)
    )

    for row in reconciliations:

        print(
            f"  {row['Gene']}: "
            f"F3={row['F3']} "
            f"F4={row['F4']} "
            f"former isolate="
            f"{row['Linked_original_CORE_isolates']}"
        )

    print(
        "Projected association weight sum:",
        report["projected_fractional_edge_weight_sum"]
    )

    print(
        "Module pairs (including within-module):",
        len(module_edges)
    )

    print("Outputs:", OUT)


if __name__ == "__main__":
    main()
