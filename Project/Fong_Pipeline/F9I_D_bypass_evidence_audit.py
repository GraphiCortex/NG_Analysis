
#!/usr/bin/env python3
r"""F9I-D — CHEK1-bypass evidence audit of Fong-derived candidate routes.

Run from Project/:
    py .\Fong_Pipeline\F9I_D_bypass_evidence_audit.py --self-test
    py .\Fong_Pipeline\F9I_D_bypass_evidence_audit.py

Requires adjacent (unchanged) scripts:
    F9F_fongwide_directed_connectivity.py
    F9I_mechanistic_route_audit.py
and original cached F5 inputs, corrected F5 metadata, F9F gene table,
F9I-B representative routes, and F9I-C deletion settings.

Scientific intent:
    * Independently reconstruct SHA-verified original mouse prior network.
    * Reproduce F9I-C exact deletion statuses (9 selected pairs, 4 settings).
    * Find a small set of simple, constrained bypass routes after blocking
      CHEK1 for RAD51/TIMELESS and E2F1 for TOPBP1 -> PYCARD.
    * Prefer routes without accession-like identifiers, but explicitly audit
      the impact of excluding them rather than silently treating them as invalid.
    * Collect database source tokens and citations for every displayed edge.
    * Produce manual *study-context* review template: no assertion that a
      cited paper used neurons, mouse, or the Fong system until reviewed.

We do NOT infer causal activity, tissue relevance, or statistical significance.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "Data/Processed/FongDownstream"
F5 = BASE / "F5_Mechanistic/Corrected"
F9F = BASE / "F9F_FongWideDirected"
F9IB = BASE / "F9I_MechanisticRoutes"
F9IC = BASE / "F9I_C_DependencySensitivity"
OUT = BASE / "F9I_D_EvidenceQuality"
SETTING = ("E4_X1", 4, 1)
SETTINGS = (("E4_X0", 4, 0), ("E4_X1", 4, 1),
            ("E4_X2", 4, 2), ("E5_X2", 5, 2))
DEATH = ("FAS", "TNFRSF10B", "PYCARD", "MLKL")
PAIRS = tuple((s, t, "CHEK1") for s in ("RAD51", "TIMELESS") for t in DEATH) + (
    ("TOPBP1", "PYCARD", "E2F1"),
)
ACCESSION_LIKE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
PMID_LIKE = re.compile(r"(?<!\d)[1-9]\d{6,8}(?!\d)")


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, records, fields):
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(records)


def truth(x):
    return str(x).strip().lower() in {"true", "1", "yes"}


def exists_within(states, max_edges, max_external):
    return any(e <= max_edges and x <= max_external for e, x in states)


def best_state(states, max_edges, max_external):
    ok = [(e, x) for e, x in states if e <= max_edges and x <= max_external]
    if not ok:
        return None
    min_x = min(x for _, x in ok)
    return min_x, min(e for e, x in ok if x == min_x)


def remove_nodes(outgoing, excluded):
    """Remove each prohibited node and all its incident edges."""
    return {
        u: tuple(v for v in vs if v not in excluded)
        for u, vs in outgoing.items() if u not in excluded
    }


def reverse(outgoing):
    incoming = defaultdict(list)
    for u, vs in outgoing.items():
        for v in vs:
            incoming[v].append(u)
    return incoming


def labeled_sign(mask):
    return {0: "unsigned", 1: "stimulation", 2: "inhibition", 3: "ambiguous"}[mask]


def select_references(text):
    """Potential PubMed numeric tokens, NOT independently verified PMIDs."""
    return sorted(set(PMID_LIKE.findall(text or "")))


def self_test():
    from F9I_mechanistic_route_audit import (
        exact_states, rev_distances, choose_paths, path_effect
    )
    graph = {
        "A": ("B", "X"),
        "B": ("C",),
        "X": ("C",),
        "C": ("D",),
        "D": (),
    }
    seeds = {"A", "B", "C", "D"}
    records = {(u, v): {"mask": 1, "effort": 1.0} for u, vs in graph.items() for v in vs}
    baseline = exact_states(graph, "A", "D", seeds, rev_distances(reverse(graph), "D"), 4, 1)
    after = remove_nodes(graph, {"B"})
    after_states = exact_states(after, "A", "D", seeds, rev_distances(reverse(after), "D"), 4, 1)
    assert exists_within(baseline, 4, 0)
    assert not exists_within(after_states, 4, 0)
    assert exists_within(after_states, 4, 1)
    min_x, shortest_e = best_state(after_states, 4, 1)
    assert (min_x, shortest_e) == (1, 3)
    chosen, capped, _ = choose_paths(after, records, "A", "D", seeds,
                                    rev_distances(reverse(after), "D"), 4, min_x, shortest_e)
    assert chosen == [("A", "X", "C", "D")] and not capped
    assert path_effect(chosen[0], records) == "net_positive"
    assert select_references("A:12345678; B:12345678; X:32145678") == ["12345678", "32145678"]
    assert not select_references("no PMID")
    print("F9I-D SELF-TEST PASSED: blocked-route discovery, constraints and provenance parsing.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--max-paths", type=int, default=8)
    args = parser.parse_args()
    if args.self_test:
        return self_test()
    if not 1 <= args.max_paths <= 30:
        raise ValueError("--max-paths must be 1–30")

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from F9F_fongwide_directed_connectivity import load_graph
    from F9I_mechanistic_route_audit import (
        exact_states, rev_distances, choose_paths, get_provenance, path_effect
    )

    meta = json.loads((F5 / "F5_corrected_metadata.json").read_text(encoding="utf-8"))
    f9f = read_csv(F9F / "F9F_gene_connectivity.csv")
    assert len(f9f) == 321
    mapping = {r["Gene"]: r["Graph_identifier"] for r in f9f if truth(r["Mapped_in_prior"])}
    assert len(mapping) == 257
    seeds = {r["Gene"] for r in f9f}
    for gene, graphid in mapping.items():
        if graphid != gene:
            seeds.add(graphid)
    assert len(seeds) == meta["graph_seed_gene_count_after_alias_resolution"]

    print("F9I-D: reconstructing SHA-verified F5 network...")
    records, nodes, layers = load_graph(meta, set(mapping.values()))
    outgoing = defaultdict(list)
    for u, v in records:
        outgoing[u].append(v)
    outgoing = {u: tuple(sorted(vs)) for u, vs in outgoing.items()}
    incoming = reverse(outgoing)
    accession_nodes = {g for g in nodes if g not in seeds and ACCESSION_LIKE.fullmatch(g)}

    crows = read_csv(F9IC / "F9I_C_dependency_settings.csv")
    cdict = {
        (r["Source"], r["Target"], r["Blocked_intermediate"], r["Setting"]): r
        for r in crows
    }
    assert len(cdict) == 360
    previous_paths = {
        (r["Source"], r["Target"], r["Route"])
        for r in read_csv(F9IB / "F9I_representative_routes.csv")
    }

    pair_rows, route_rows, route_edges = [], [], set()
    checked = 0
    route_caps = 0

    for source, target, blocker in PAIRS:
        src, dst = mapping[source], mapping[target]
        d0 = rev_distances(incoming, dst, 5)
        original_states = exact_states(outgoing, src, dst, seeds, d0, 5, 2)
        after = remove_nodes(outgoing, {blocker})
        dr = rev_distances(reverse(after), dst, 5)
        after_states = exact_states(after, src, dst, seeds, dr, 5, 2)
        after_clean = remove_nodes(outgoing, {blocker} | accession_nodes)
        dc = rev_distances(reverse(after_clean), dst, 5)
        clean_states = exact_states(after_clean, src, dst, seeds, dc, 5, 2)

        for label, max_e, max_x in SETTINGS:
            baseline = exists_within(original_states, max_e, max_x)
            bypass = exists_within(after_states, max_e, max_x)
            clean_bypass = exists_within(clean_states, max_e, max_x)
            old = cdict[source, target, blocker, label]
            assert baseline == truth(old["Baseline_reachable"]), (source, target, label, "baseline")
            assert bypass == truth(old["After_deletion_reachable"]), (source, target, label, "blocked")
            checked += 1
            pair_rows.append({
                "Source": source, "Target": target, "Removed_node": blocker,
                "Setting": label, "Baseline_reachable": baseline,
                "Bypass_after_deletion": bypass,
                "Bypass_without_accession_like": clean_bypass,
                "Bypass_min_external": (best_state(after_states, max_e, max_x) or ("", ""))[0],
                "Bypass_shortest_at_min_external": (best_state(after_states, max_e, max_x) or ("", ""))[1],
                "Interpretation": "BYPASS_EXISTS" if bypass else "NO_BYPASS_UNDER_CONSTRAINT"
            })

        # The earliest relaxed setting E4_X1 is the specific hypothesis.
        # Prefer accession-free routes while retaining an exact audit of
        # whether pruning changed reachability.
        state = best_state(clean_states, 4, 1)
        graph = after_clean
        dist = dc
        pool = "accession_filtered"
        if state is None:
            state = best_state(after_states, 4, 1)
            graph = after
            dist = dr
            pool = "full_prior_with_accession_like"
        if state is None:
            print("No E4_X1 bypass:", source, target, blocker)
            continue

        min_external, shortest_edges = state
        reps, capped, expansions = choose_paths(
            graph, records, src, dst, seeds, dist,
            max_edges=4, min_x=min_external,
            shortest_at_min_x=shortest_edges,
            max_paths=args.max_paths, state_cap=200_000
        )
        if capped:
            route_caps += 1
        if not reps:
            raise RuntimeError(f"Route extraction failed despite exact reachability: {source}->{target}")

        for ix, route in enumerate(reps, 1):
            canonical = [source if g == src else target if g == dst else g for g in route]
            intermediates = route[1:-1]
            external = [g for g in intermediates if g not in seeds]
            bad = [g for g in intermediates if g in accession_nodes]
            for a, b in zip(route, route[1:]):
                route_edges.add((a, b))
            effects = [records[a, b]["mask"] for a, b in zip(route, route[1:])]
            route_rows.append({
                "Route_ID": f"F9ID_{source}_{target}_{ix:02d}",
                "Source": source, "Target": target,
                "Blocked_intermediate": blocker,
                "Setting": "E4_X1",
                "Route": " -> ".join(canonical),
                "Edges": len(route) - 1,
                "Outside_F2_count": len(external),
                "Outside_F2_symbols": ";".join(external),
                "Accession_like_symbols": ";".join(bad),
                "Unknown_or_ambiguous_edge_count": sum(m in (0, 3) for m in effects),
                "Net_database_sign": path_effect(route, records),
                "Route_pool": pool,
                "Already_in_F9I_B_display": (source, target, " -> ".join(route)) in previous_paths,
                "Provenance_status": "AWAITING_PER_EDGE_REVIEW",
                "Biological_cell_type_validated": "NOT_ASSESSED",
                "Route_search_state_cap_hit": capped,
                "Display_only": True,
            })
        print(f"  {source} -> {target}, remove {blocker}: {len(reps)} example bypasses "
              f"(minX={min_external}; {pool})")

    # Get EXACT source references from original cached rows for the
    # displayed edges, rather than inferring papers from general biology.
    print(f"Recovering original source/reference records for {len(route_edges)} unique edges...")
    provenance = get_provenance(route_edges)
    edge_rows, literature_rows, mentioned = [], [], set()
    for a, b in sorted(route_edges):
        info = records[a, b]
        prov = provenance[a, b]
        resources = sorted(prov["sources"])
        refs = sorted(prov["references"])
        token_text = ";".join(refs)
        pmids = select_references(token_text)
        edge_rows.append({
            "Source": a, "Target": b,
            "Source_in_F2": a in seeds,
            "Target_in_F2": b in seeds,
            "Prior_effect": labeled_sign(info["mask"]),
            "Layers": ";".join(sorted(info["layers"])),
            "Supporting_records": info["count"],
            "Max_curation_effort_NOT_confidence_score": info["effort"],
            "Source_databases": ";".join(resources),
            "Reference_tokens": token_text,
            "Candidate_PubMed_IDs_unverified": ";".join(pmids),
            "Has_reference_tokens": bool(refs),
            "Neuron_context_verified": "NOT_ASSESSED",
            "Fong_cell_context_verified": "NOT_ASSESSED",
            "Protein_or_transcript_activity_verified": "NOT_ASSESSED",
            "Manual_review_required": True,
        })
        if not refs:
            refs = [""]
        for ref in refs:
            k = (a, b, ref)
            if k in mentioned:
                continue
            mentioned.add(k)
            literature_rows.append({
                "Source": a, "Target": b,
                "Reference_token": ref,
                "Candidate_PMID_unverified": ";".join(select_references(ref)),
                "Source_databases": ";".join(resources),
                "Publication_verified": "",
                "Organism_in_original_study": "",
                "Cell_type_or_tissue": "",
                "Experimental_system": "",
                "Perturbation_or_stimulus": "",
                "Measured_outcome": "",
                "Supports_direction": "",
                "Supports_effect_sign": "",
                "Evidence_directness": "",
                "Relevant_to_neural_cells": "",
                "Relevant_to_Rb_E2F_loss": "",
                "Supports_death_endpoint_mechanism": "",
                "Review_notes": "",
            })

    # Add per-route reference coverage after provenance has been collected.
    ec = {(r["Source"], r["Target"]): r for r in edge_rows}
    for route in route_rows:
        nodes_in_route = route["Route"].split(" -> ")
        attached = [ec[(a, b)] for a, b in zip(nodes_in_route, nodes_in_route[1:])]
        route["Edges_with_reference_tokens"] = sum(r["Has_reference_tokens"] for r in attached)
        route["All_edges_have_reference_tokens"] = all(r["Has_reference_tokens"] for r in attached)
        route["Database_sign_interpretable"] = route["Unknown_or_ambiguous_edge_count"] == 0
        route["References_DO_NOT_validate_cell_type"] = True

    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / "F9ID_bypass_pair_sensitivity.csv", pair_rows, list(pair_rows[0]))
    write_csv(OUT / "F9ID_bypass_representative_routes.csv", route_rows, list(route_rows[0]))
    write_csv(OUT / "F9ID_route_edge_provenance.csv", edge_rows, list(edge_rows[0]))
    write_csv(OUT / "F9ID_manual_reference_review_template.csv", literature_rows, list(literature_rows[0]))
    report = {
        "analysis": "F9I-D CHEK1/E2F1 bypass evidence triage",
        "full_F5_nodes": len(nodes), "full_F5_directed_edges": len(records),
        "selected_source_endpoint_blocker_triples": len(PAIRS),
        "F9I_C_pair_setting_checks_reproduced": checked,
        "selected_setting": "E4_X1",
        "bypass_triples_E4_X1": sum(r["Setting"] == "E4_X1" and r["Bypass_after_deletion"] for r in pair_rows),
        "bypass_triples_E4_X1_without_accession_like": sum(r["Setting"] == "E4_X1" and r["Bypass_without_accession_like"] for r in pair_rows),
        "representative_routes": len(route_rows),
        "unique_route_edges": len(edge_rows),
        "distinct_reference_tokens_to_review": len(literature_rows),
        "state_cap_hits": route_caps,
        "prior_layers": layers,
        "limitations": [
            "All route results conditional on original SHA-verified OmniPath F5 mouse graph.",
            "Routes are non-exhaustive examples; E/X reachability and deletions are exact.",
            "Literature review fields intentionally blank: resource citations are not evidence of neuronal or Fong-context relevance.",
            "OmniPath signed edges and net signs are database annotations, not experimental response predictions.",
            "Degree/path significance, phosphoproteomics, protein activity, and direct causal attribution were not tested.",
        ]
    }
    (OUT / "F9ID_validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("F9I-D PASS — exact F9I-C statuses reproduced; provenance exported")
    print("  Checked:", checked, "pair-setting deletion comparisons")
    print("  E4_X1 bypasses:", report["bypass_triples_E4_X1"],
          "accession-free:", report["bypass_triples_E4_X1_without_accession_like"])
    print("  Illustrative bypass routes:", len(route_rows),
          "edges:", len(edge_rows), "literature rows:", len(literature_rows))
    print("  Output:", OUT)


if __name__ == "__main__":
    main()
