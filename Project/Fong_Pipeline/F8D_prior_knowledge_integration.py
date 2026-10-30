
#!/usr/bin/env python3
"""
F8D: Current F5 directed prior-knowledge integration.

Integrates:
  - F8C STRING network
  - Corrected F5 directed edges
  - Exact reachability under E4_X0, E4_X1, E4_X2, E5_X2
  - Representative routes
  - Exact-required intermediates
  - Historical connector comparison from F6

Python standard library only.
Does not modify earlier pipeline outputs.
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


# ============================================================
# 1. PATHS AND DEFAULTS
# ============================================================

PROJECT = Path(__file__).resolve().parents[1]

BASE = PROJECT / "Data/Processed/FongDownstream"

F5 = BASE / "F5_Mechanistic/Corrected"

F8 = BASE / "F8_FunctionalAtlas"

OUT = F8 / "F8D_PriorKnowledge"

# Maximum 4 edges, maximum 1 non-F2 intermediate.
DEFAULT_SETTING = "E4_X1"


# ============================================================
# 2. UTILITIES
# ============================================================

def load_csv(path: Path) -> list[dict]:

    with path.open(
        newline="",
        encoding="utf-8-sig"
    ) as f:
        return list(csv.DictReader(f))


def write_csv(
    path: Path,
    rows: list[dict],
    columns: list[str]
) -> None:

    with path.open(
        "w",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=columns,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def dump(path: Path, data: object) -> None:

    # Compact the large viewer JSON for easier distribution.
    compact = path.name == "F8D_viewer_data.json"

    path.write_text(
        json.dumps(
            data,
            ensure_ascii=False,
            indent=None if compact else 2,
            separators=(",", ":") if compact else None
        ),
        encoding="utf-8"
    )


def yes(value: object) -> bool:

    return str(value).strip().lower() in {
        "true", "yes", "1"
    }


def semi(value: str) -> list[str]:

    return [
        x.strip()
        for x in (value or "").split(";")
        if x.strip()
    ]


def num_or_none(value: str):

    return None if not value else float(value)


# ============================================================
# 3. MAIN
# ============================================================

def main() -> None:

    # --------------------------------------------------------
    # 3A. Read previously generated evidence
    # --------------------------------------------------------

    f8c_path = (
        F8
        / "F8C_IntermoduleConnectivity/F8C_network.json"
    )

    f8c = json.loads(
        f8c_path.read_text(encoding="utf-8")
    )

    metadata = json.loads(
        (F5 / "F5_corrected_metadata.json").read_text(
            encoding="utf-8"
        )
    )

    union = load_csv(
        F5 / "representative_union_edges.csv"
    )

    paths = load_csv(
        F5 / "representative_paths.csv"
    )

    pairs = load_csv(
        F5 / "pair_setting_summary.csv"
    )

    required = load_csv(
        F5 / "exact_required_nodes.csv"
    )

    historical = load_csv(
        F8 / "prior_node_inventory.csv"
    )

    assert (
        metadata["phase"]
        == "F5_corrected_omnipath_downstream"
    )

    settings = metadata["settings"]

    assert set(settings) == {
        "E4_X0",
        "E4_X1",
        "E4_X2",
        "E5_X2",
    }

    assert DEFAULT_SETTING in settings

    assert (
        settings[DEFAULT_SETTING]["max_edges"] == 4
    )

    assert (
        settings[DEFAULT_SETTING]["max_external"] == 1
    )


    # ========================================================
    # 4. DOCUMENTED IDENTIFIER NORMALIZATION
    # ========================================================

    # Do not invent aliases.
    # Normalize only mappings supplied in corrected F5 metadata.

    aliases = {}

    for canonical, synonyms in metadata[
        "anchor_aliases"
    ].items():

        for synonym in synonyms:

            if (
                synonym in aliases
                and aliases[synonym] != canonical
            ):
                raise ValueError(
                    f"Conflicting alias: {synonym}"
                )

            aliases[synonym] = canonical


    def normalize(name: str) -> str:

        name = name.strip()

        return aliases.get(name, name)


    # ========================================================
    # 5. VERIFY EXPERIMENTAL NETWORK
    # ========================================================

    observed = {
        n["id"]: n
        for n in f8c["nodes"]
    }

    assert len(observed) == 321

    assert len(f8c["edges"]) == 1139

    assert sum(
        bool(n["shared37"])
        for n in observed.values()
    ) == 37

    assert len({
        edge["id"]
        for edge in f8c["edges"]
    }) == 1139

    assert len(historical) == 39

    historical_by_gene = {
        r["Gene"]: r
        for r in historical
    }

    assert len(historical_by_gene) == 39


    # ========================================================
    # 6. EXACT SOURCE-TARGET REACHABILITY
    # ========================================================

    # Reachability comes from the exact bounded-state search.
    # It is NOT inferred from the representative path sample.

    pair_by_key = {}

    pair_rows = []

    for r in pairs:

        setting = r["Setting"]

        key = (
            setting,
            r["Source"],
            r["Target"]
        )

        if (
            setting not in settings
            or key in pair_by_key
        ):
            raise ValueError(
                f"Invalid/duplicate setting-pair: {key}"
            )

        if (
            int(r["Max_edges"])
            != settings[setting]["max_edges"]
            or int(r["Max_external"])
            != settings[setting]["max_external"]
        ):
            raise ValueError(
                f"Setting mismatch: {key}"
            )

        record = {

            "setting": setting,

            "source": normalize(
                r["Resolved_source"] or r["Source"]
            ),

            "source_anchor": r["Source"],

            "source_family": r["Source_family"],

            "target": normalize(
                r["Resolved_target"] or r["Target"]
            ),

            "target_anchor": r["Target"],

            "target_family": r["Target_family"],

            "testable": yes(r["Testable"]),

            "reachable_exact":
                yes(r["Reachable_exact"]),

            "min_external_exact":
                num_or_none(
                    r["Min_external_exact"]
                ),

            "shortest_edges_any_exact":
                num_or_none(
                    r["Shortest_edges_any_exact"]
                ),

            "shortest_edges_at_min_external":
                num_or_none(
                    r["Shortest_edges_at_min_external"]
                ),

            "representative_path_count":
                int(r["Representative_paths"]),

            "representative_state_cap_hit":
                yes(r["Representative_state_cap_hit"]),
        }

        if (
            record["reachable_exact"]
            and not record["testable"]
        ):
            raise ValueError(
                f"Reachable but not testable: {key}"
            )

        pair_by_key[key] = record

        pair_rows.append(record)


    # ========================================================
    # 7. DIRECTED PRIOR-KNOWLEDGE EDGE UNION
    # ========================================================

    # One edge can recur under multiple settings.
    # We retain one edge identity plus per-setting provenance.

    directed_by_key = {}

    setting_edges = set()

    for r in union:

        setting = r["Setting"]

        source = normalize(r["Source"])
        target = normalize(r["Target"])

        edge_key = (source, target)

        setting_key = (
            setting,
            source,
            target
        )

        if (
            setting not in settings
            or setting_key in setting_edges
        ):
            raise ValueError(
                f"Duplicate/invalid directed edge: "
                f"{setting_key}"
            )

        setting_edges.add(setting_key)

        if (
            yes(r["Source_F2_seed"])
            != (source in observed)
        ):
            raise ValueError(
                f"Source F2 flag mismatch: {setting_key}"
            )

        if (
            yes(r["Target_F2_seed"])
            != (target in observed)
        ):
            raise ValueError(
                f"Target F2 flag mismatch: {setting_key}"
            )

        if source == target:
            raise ValueError(
                f"Self-directed edge: {setting_key}"
            )

        effect = int(r["Effect"])

        label = r["Effect_label"]

        if (
            effect not in {-1, 0, 1}
            or (
                (effect == 0)
                != (label in {"unsigned", "ambiguous"})
            )
        ):
            raise ValueError(
                f"Invalid effect/label: {setting_key}"
            )

        if (
            (effect == -1 and label != "inhibition")
            or (
                effect == 1
                and label != "stimulation"
            )
        ):
            raise ValueError(
                f"Sign-label mismatch: {setting_key}"
            )

        if edge_key not in directed_by_key:

            directed_by_key[edge_key] = {

                "id": f"DIR::{source}::{target}",

                "source": source,
                "target": target,

                "source_original": r["Source"],
                "target_original": r["Target"],

                "effect_numeric": effect,
                "effect_label": label,

                "source_F2": source in observed,
                "target_F2": target in observed,

                "layers": set(),

                "settings": {},

                "kind":
                    "DIRECTED_PRIOR_KNOWLEDGE",
            }

        d = directed_by_key[edge_key]

        if (
            d["effect_numeric"],
            d["effect_label"]
        ) != (
            effect,
            label
        ):
            raise ValueError(
                f"Conflicting signs across settings: "
                f"{edge_key}"
            )

        d["layers"].update(
            semi(r["Layers"])
        )

        d["settings"][setting] = {

            "layers": semi(r["Layers"]),

            "supporting_records":
                int(r["Supporting_records"]),

            "max_curation_effort":
                float(r["Max_curation_effort"]),

            "representative_path_usage":
                int(r["Representative_path_usage"]),

            "representative_pair_usage":
                int(r["Representative_pair_usage"]),
        }

    assert len(directed_by_key) == 393

    assert len(union) == 1438

    external_from_edges = {

        x

        for (a, b) in directed_by_key

        for x in (a, b)

        if x not in observed
    }


    # ========================================================
    # 8. REPRESENTATIVE PATH RECORDS
    # ========================================================

    # Retain all capped representative paths.
    #
    # These are not an exhaustive set of all possible routes.
    # Their recurrence is not independent mechanistic evidence.

    route_rows = []

    route_keys = set()

    external_from_paths = set()

    for r in paths:

        setting = r["Setting"]

        key = (
            setting,
            r["Source"],
            r["Target"]
        )

        if (
            key not in pair_by_key
            or not pair_by_key[key]["reachable_exact"]
        ):
            raise ValueError(
                f"Route absent from reachable pairs: {key}"
            )

        route_id = (
            setting,
            r["Source"],
            r["Target"],
            int(r["Path_index"])
        )

        if route_id in route_keys:
            raise ValueError(
                f"Duplicate route: {route_id}"
            )

        route_keys.add(route_id)

        nodes_original = [
            x.strip()
            for x in r["Path"].split("->")
        ]

        nodes = [
            normalize(x)
            for x in nodes_original
        ]

        if (
            nodes[0] != pair_by_key[key]["source"]
            or nodes[-1] != pair_by_key[key]["target"]
        ):
            raise ValueError(
                f"Route endpoint mismatch: {route_id}"
            )

        if (
            len(nodes) - 1
            != int(r["Edges"])
        ):
            raise ValueError(
                f"Path-length mismatch: {route_id}"
            )

        intermediates = nodes[1:-1]

        external = [
            x for x in intermediates
            if x not in observed
        ]

        internal = [
            x for x in intermediates
            if x in observed
        ]

        reported_external = [
            normalize(x)
            for x in semi(
                r["External_intermediates"]
            )
        ]

        reported_internal = [
            normalize(x)
            for x in semi(
                r["F2_internal_intermediates"]
            )
        ]

        if (
            external != reported_external
            or internal != reported_internal
        ):
            raise ValueError(
                f"Intermediate classification mismatch: "
                f"{route_id}"
            )

        if (
            len(external)
            != int(r["External_intermediate_count"])
            or len(internal)
            != int(r["F2_internal_intermediate_count"])
        ):
            raise ValueError(
                f"Intermediate count mismatch: {route_id}"
            )

        if (
            len(nodes) - 1
            > settings[setting]["max_edges"]
            or len(external)
            > settings[setting]["max_external"]
        ):
            raise ValueError(
                f"Route violates setting: {route_id}"
            )

        if (
            int(r["Path_index"])
            > metadata["representative_path_limit"]
        ):
            raise ValueError(
                f"Representative path exceeds cap: {route_id}"
            )

        route_edge_ids = []

        for a, b in zip(
            nodes[:-1],
            nodes[1:]
        ):

            if (
                setting, a, b
            ) not in setting_edges:

                raise ValueError(
                    f"Unreported directed route edge: "
                    f"{setting}, {a}, {b}"
                )

            route_edge_ids.append(
                directed_by_key[(a, b)]["id"]
            )

        external_from_paths.update(external)

        route_rows.append({

            "id": (
                f"PATH::{setting}::"
                f"{r['Source']}::"
                f"{r['Target']}::"
                f"{r['Path_index']}"
            ),

            "setting": setting,

            "source": nodes[0],
            "target": nodes[-1],

            "source_family":
                r["Source_family"],

            "target_family":
                r["Target_family"],

            "path_index":
                int(r["Path_index"]),

            "nodes": nodes,

            "nodes_original": nodes_original,

            "edge_ids": route_edge_ids,

            "edges": len(nodes) - 1,

            "external_intermediates":
                external,

            "internal_F2_intermediates":
                internal,

            "effect_reported":
                r["Path_effect"],

            "effect_numeric_reported":
                num_or_none(
                    r["Path_effect_numeric"]
                ),

            "evidence_layers_reported":
                semi(r["Path_layers"]),

            "representative_only": True,
        })

    assert len(route_rows) == 8039

    assert (
        external_from_edges
        == external_from_paths
    )

    assert len(external_from_paths) == 22


    # ========================================================
    # 9. CHECK REPRESENTATIVE COUNTS AGAINST EXACT PAIRS
    # ========================================================

    count_by_pair = Counter(

        (
            r["setting"],
            r["source"],
            r["target"]
        )

        for r in route_rows
    )

    for (
        setting,
        source_original,
        target_original
    ), p in pair_by_key.items():

        count = count_by_pair[(
            setting,
            p["source"],
            p["target"]
        )]

        if (
            count
            != p["representative_path_count"]
        ):
            raise ValueError(
                "Representative path total mismatch: "
                f"{setting}, "
                f"{source_original}, "
                f"{target_original}"
            )

    expected_reachable = {

        "E4_X0": 98,
        "E4_X1": 108,
        "E4_X2": 108,
        "E5_X2": 108,
    }

    for setting, count in expected_reachable.items():

        assert sum(
            r["reachable_exact"]
            for r in pair_rows
            if r["setting"] == setting
        ) == count

        assert sum(
            r["testable"]
            for r in pair_rows
            if r["setting"] == setting
        ) == 117


    # ========================================================
    # 10. EXACT-REQUIRED INTERMEDIATES
    # ========================================================

    # Keep exact-required tests separate from ordinary
    # representative-path occurrence counts.

    required_rows = []

    req_by_pair = defaultdict(list)

    for r in required:

        key = (
            r["Setting"],
            r["Source"],
            r["Target"]
        )

        if (
            key not in pair_by_key
            or not pair_by_key[key]["reachable_exact"]
        ):
            raise ValueError(
                f"Unrecognized required-node test: {key}"
            )

        entry = {

            "setting": r["Setting"],

            "source": normalize(r["Source"]),
            "target": normalize(r["Target"]),

            "node": normalize(r["Node"]),

            "node_class": r["Node_class"],

            "exact_required":
                yes(r["Exact_required"]),

            "representative_occurrences":
                int(r["Representative_occurrences"]),

            "representative_path_count":
                int(r["Representative_path_count"]),

            "definition":
                metadata["required_node_definition"],
        }

        required_rows.append(entry)

        req_by_pair[key].append(entry)

    assert len(required_rows) == 12

    for route in route_rows:

        route["exact_required_intermediates"] = [

            x["node"]

            for x in req_by_pair[(
                route["setting"],
                route["source"],
                route["target"]
            )]

            if (
                x["exact_required"]
                and x["node"] in route["nodes"][1:-1]
            )
        ]


    # ========================================================
    # 11. HISTORICAL-CONNECTOR COMPARISON
    # ========================================================

    # This comparison does NOT reconstruct historical edges.
    #
    # "Not recovered" means absent from current selected
    # representative intermediates, not biologically disproven.

    history_rows = []

    for gene, r in sorted(
        historical_by_gene.items()
    ):

        old = yes(
            r["Old_external_connector"]
        )

        current = (
            gene in external_from_paths
        )

        if (
            current
            != yes(r["New_F5_representative_external"])
        ):
            raise ValueError(
                f"F6/F5 node mismatch: {gene}"
            )

        if old and current:
            status = "BOTH"

        elif current:
            status = "CURRENT_ONLY"

        elif old:
            status = "HISTORICAL_ONLY"

        else:
            status = "NEITHER"

        history_rows.append({

            "gene": gene,

            "comparison_status": status,

            "historical_connector": old,

            "current_representative_external":
                current,

            "strict_Fong_DEG":
                yes(
                    r["Observed_in_current_strict_Fong"]
                ),

            "F2_member":
                yes(r["Observed_in_F2"]),

            "old_master_member":
                yes(r["Observed_in_old_master"]),

            "Julian_promoter_evidence":
                yes(r["Julian_sig_promoter"]),

            "historical_edges_available": False,

            "notes": (
                "Historical connector links "
                "are not reconstructed here."
            ),
        })

    status_counts = Counter(
        r["comparison_status"]
        for r in history_rows
    )

    assert status_counts["BOTH"] == 2

    assert status_counts["CURRENT_ONLY"] == 20

    assert status_counts["HISTORICAL_ONLY"] == 17


    # ========================================================
    # 12. BUILD THE COMBINED NODE INVENTORY
    # ========================================================

    # An external node is outside F2, not necessarily
    # absent from all Fong experimental expression data.

    usages = Counter(

        node

        for r in route_rows

        if r["setting"] == DEFAULT_SETTING

        for node in r["external_intermediates"]
    )

    gene_nodes = []

    for n in f8c["nodes"]:

        enriched = dict(n)

        enriched.update({

            "origin": "F2_EXPERIMENTAL",

            "has_Fong_expression": True,
        })

        gene_nodes.append(enriched)


    # Add the 22 external intermediates from CURRENT F5.
    # Do not add historical-only nodes to the main graph.

    for name in sorted(external_from_paths):

        r = historical_by_gene[name]

        if yes(r["Observed_in_F2"]):
            raise ValueError(
                "F5 external intermediate unexpectedly "
                f"belongs to F2: {name}"
            )

        gene_nodes.append({

            "id": name,

            "label": name,

            "tier": "F5_EXTERNAL",

            "origin": "F5_EXTERNAL_INTERMEDIATE",

            "modules": [],

            "direction": None,
            "log2FC": None,
            "padj": None,

            "shared37": False,

            "promoter_evidence":
                yes(r["Julian_sig_promoter"]),

            "has_Fong_expression":
                yes(
                    r["Observed_in_current_strict_Fong"]
                ),

            "strict_Fong_DEG":
                yes(
                    r["Observed_in_current_strict_Fong"]
                ),

            "old_connector":
                yes(r["Old_external_connector"]),

            "default_setting_representative_occurrences":
                usages[name],
        })

    assert len(gene_nodes) == 343


    # ========================================================
    # 13. PREPARE DIRECTED EDGE DATA
    # ========================================================

    directed_edges = []

    edge_audit = []

    for key, d in sorted(
        directed_by_key.items()
    ):

        row = dict(d)

        row["layers"] = sorted(
            row["layers"]
        )

        row["available_settings"] = [

            s for s in settings

            if s in row["settings"]
        ]

        directed_edges.append(row)

        edge_audit.append({

            "source": row["source"],
            "target": row["target"],

            "original_source":
                row["source_original"],

            "original_target":
                row["target_original"],

            "effect_label":
                row["effect_label"],

            "effect_numeric":
                row["effect_numeric"],

            "layers":
                "; ".join(row["layers"]),

            "settings":
                "; ".join(
                    row["available_settings"]
                ),

            "source_F2":
                row["source_F2"],

            "target_F2":
                row["target_F2"],

            "edge_class":
                "PRIOR_DIRECTED_NOT_EXPERIMENTAL_CAUSATION",
        })


    # ========================================================
    # 14. VIEWER-READY JSON
    # ========================================================

    # This JSON will be consumed by our React/Cytoscape viewer.
    #
    # STRING and directed prior knowledge remain separate layers.
    #
    # Historical-only connectors remain comparison metadata,
    # without fabricated gene-gene edges.

    bundle = {

        "metadata": {

            "phase": "F8D_PriorKnowledge",

            "default_setting": DEFAULT_SETTING,

            "settings": settings,

            "species":
                metadata["organism"],

            "node_count_F2": 321,

            "current_F5_external_node_count": 22,

            "source_network": (
                "F8C undirected STRING "
                "functional associations"
            ),

            "directed_source": (
                "F5 corrected OmniPath-based "
                "prior knowledge"
            ),

            "reachability_definition": (
                "Exact within specified bounded setting"
            ),

            "representative_warning":
                metadata["important_method_note"],

            "exact_required_definition":
                metadata["required_node_definition"],

            "historical_policy": (
                "Historical-only nodes inventoried, "
                "not inserted into the main graph "
                "without old edge provenance."
            ),

            "alias_normalization": aliases,
        },

        "nodes": gene_nodes,

        "string_edges": f8c["edges"],

        "directed_edges": directed_edges,

        "modules": f8c["modules"],

        "module_edges": f8c["module_edges"],

        "exact_reachability": pair_rows,

        "representative_routes": route_rows,

        "exact_required_nodes": required_rows,

        "historical_comparison": history_rows,
    }


    # ========================================================
    # 15. EXPORT
    # ========================================================

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    dump(
        OUT / "F8D_viewer_data.json",
        bundle
    )

    write_csv(
        OUT / "F8D_directed_edges.csv",
        edge_audit,
        list(edge_audit[0])
    )

    write_csv(
        OUT / "F8D_external_node_inventory.csv",
        history_rows,
        list(history_rows[0])
    )

    write_csv(
        OUT / "F8D_exact_reachability.csv",
        pair_rows,
        list(pair_rows[0])
    )

    write_csv(
        OUT / "F8D_exact_required_nodes.csv",
        required_rows,
        list(required_rows[0])
    )


    # ========================================================
    # 16. VALIDATION REPORT
    # ========================================================

    report = {

        "F2_nodes": 321,

        "current_external_nodes":
            len(external_from_paths),

        "strict_Fong_DEGs_among_external": sorted(

            gene

            for gene in external_from_paths

            if yes(
                historical_by_gene[gene][
                    "Observed_in_current_strict_Fong"
                ]
            )
        ),

        "viewer_nodes_including_current_external":
            len(gene_nodes),

        "original_STRING_edges":
            len(f8c["edges"]),

        "unique_prior_directed_edges":
            len(directed_edges),

        "setting_specific_prior_edge_rows":
            len(union),

        "representative_routes_all_settings":
            len(route_rows),

        "historical_comparison_status":
            dict(status_counts),

        "reachable_pairs_by_setting":
            expected_reachable,

        "setting_directed_edges":
            dict(
                Counter(
                    r["Setting"]
                    for r in union
                )
            ),

        "setting_representative_routes":
            dict(
                Counter(
                    r["Setting"]
                    for r in paths
                )
            ),

        "true_exact_required_nodes": [

            x for x in required_rows

            if x["exact_required"]
        ],

        "notes": [

            (
                "Prior interactions are directed, but do "
                "not prove signaling activation in Fong cells."
            ),

            (
                "Path effects are prior-derived and "
                "not measured causal effects."
            ),

            "No historical-only edges are fabricated.",

            (
                "Path recurrence is affected by the "
                "representative-path cap."
            ),

            (
                "Experimental expression is not inferred "
                "from prior-network membership."
            ),
        ],
    }

    dump(
        OUT / "F8D_validation_report.json",
        report
    )


    # ========================================================
    # 17. TERMINAL SUMMARY
    # ========================================================

    print(
        "F8D PASS — viewer-ready prior-knowledge layer"
    )

    print(
        "F2 genes:",
        report["F2_nodes"]
    )

    print(
        "Current F5 external intermediates:",
        report["current_external_nodes"]
    )

    print(
        "External intermediates also strict Fong DEGs:",
        ", ".join(
            report["strict_Fong_DEGs_among_external"]
        ) or "none"
    )

    print(
        "Total nodes in viewer dataset:",
        report[
            "viewer_nodes_including_current_external"
        ]
    )

    print(
        "STRING edges / unique directed edges:",
        report["original_STRING_edges"],
        report["unique_prior_directed_edges"]
    )

    print(
        "Representative paths:",
        report["representative_routes_all_settings"]
    )

    print(
        "Historical comparison:",
        report["historical_comparison_status"]
    )

    print(
        "Exact reachable pairs:",
        report["reachable_pairs_by_setting"]
    )

    print(
        "Exact-required nodes:",
        [

            (
                r["setting"],
                r["source"],
                r["target"],
                r["node"]
            )

            for r in report[
                "true_exact_required_nodes"
            ]
        ]
    )

    print("Outputs:", OUT)


if __name__ == "__main__":
    main()
