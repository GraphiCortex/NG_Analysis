#!/usr/bin/env python3
r"""F9I-C: exact node-deletion robustness across the original F5 E/X settings.

Run from Project/:
    py .\Fong_Pipeline\F9I_C_dependency_sensitivity.py --self-test
    py .\Fong_Pipeline\F9I_C_dependency_sensitivity.py

Place beside F9F_fongwide_directed_connectivity.py and
F9I_mechanistic_route_audit.py. Reads original SHA-verified F5 prior data;
compares baseline reachability against F9I-B; does not edit prior outputs.

In these tests, every source and target is a mapped F2 gene. Therefore an
outside-F2 node can only occur as an intermediate, so multi-target bounded
reachability is equivalent to the original F5 pairwise E/X definition.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict, Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "Data/Processed/FongDownstream"
F5 = BASE / "F5_Mechanistic/Corrected"
F9F = BASE / "F9F_FongWideDirected"
F9IB = BASE / "F9I_MechanisticRoutes"
OUT = BASE / "F9I_C_DependencySensitivity"

SETTINGS = (("E4_X0", 4, 0), ("E4_X1", 4, 1),
            ("E4_X2", 4, 2), ("E5_X2", 5, 2))
# These are diagnostic perturbations of the prior graph, NOT interventions
# performed in cells. Their selection is exploratory, not a p-value test.
BLOCKERS = {
    "TOPBP1": ("E2F1", "CHEK1"),
    "RAD51": ("CHEK1", "BRCA2"),
    "TIMELESS": ("CHEK1", "ATM"),
    "CLSPN": ("CHEK1", "CHEK2"),
    "RPA2": ("CHEK1", "BLM"),
}
# Older 04b robustness script used this *heuristic*; it may misclassify
# unusual legitimate symbols, so every matched identifier is exported.
ACCESSION_LIKE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|"
    r"[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)


def rows(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def save_csv(path, data, names):
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=names, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(data)


def truth(v):
    return str(v).strip().lower() in ("1", "yes", "true")


def bounded_states(adj, seed_graph, source, targets, excluded=frozenset()):
    """Exact E/X-bounded reachability for F2 endpoints simultaneously.

    States track (vertex, length, external-intermediate-count). A walk
    witnessing bounded reachability can be cycle-erased, yielding a simple
    path meeting the same or stronger bounds. Results are thresholded (<=).
    """
    target_set = set(targets)
    if source in excluded:
        return {target: set() for target in targets}
    frontier = [set([source]), set(), set()]
    found = {target: set() for target in targets}
    for length in range(1, 6):
        new = [set(), set(), set()]
        for x, layer in enumerate(frontier):
            for u in layer:
                for v in adj.get(u, ()):
                    if v in excluded:
                        continue
                    # All queried targets are F2, so none is outside-F2.
                    nx = x + (v not in seed_graph)
                    if nx <= 2:
                        new[nx].add(v)
        for x, layer in enumerate(new):
            for target in target_set & layer:
                if target != source:
                    found[target].add((length, x))
        frontier = new
    return found


def reachable(states, setting):
    _, max_e, max_x = setting
    return any(e <= max_e and x <= max_x for e, x in states)


def self_test():
    # A->B->C and A->X->C, X is outside F2. B is required for E4_X0
    # but not E4_X1. The third route checks an additional source.
    adj = {"A": ["B", "X"], "B": ["C"], "X": ["C"], "C": []}
    seeds = {"A", "B", "C"}
    before = bounded_states(adj, seeds, "A", ["B", "C"])
    after = bounded_states(adj, seeds, "A", ["B", "C"], {"B"})
    assert reachable(before["C"], SETTINGS[0])
    assert not reachable(after["C"], SETTINGS[0])
    assert reachable(after["C"], SETTINGS[1])
    assert reachable(before["B"], SETTINGS[0])
    assert not reachable(after["B"], SETTINGS[3])
    # No outside-F2 intermediates: excluding X does not alter E4_X0.
    no_accession = bounded_states(adj, seeds, "A", ["C"], {"X"})
    assert reachable(no_accession["C"], SETTINGS[0])
    print("F9I-C SELF-TEST PASSED: bounded multi-target deletion and E/X sensitivity.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        return self_test()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from F9F_fongwide_directed_connectivity import load_graph

    meta = json.loads((F5 / "F5_corrected_metadata.json").read_text(encoding="utf-8"))
    gene_rows = rows(F9F / "F9F_gene_connectivity.csv")
    assert len(gene_rows) == 321 and len({r["Gene"] for r in gene_rows}) == 321
    mapped = {r["Gene"]: r["Graph_identifier"]
              for r in gene_rows if truth(r["Mapped_in_prior"])}
    assert len(mapped) == 257

    graph_seed = {r["Gene"] for r in gene_rows}
    for original, synonyms in meta["anchor_aliases"].items():
        if original in mapped and mapped[original] in synonyms:
            graph_seed.add(mapped[original])
    assert len(graph_seed) == meta["graph_seed_gene_count_after_alias_resolution"]

    print("Loading original F5 network (SHA256 and dimensions validated)...")
    edges, all_nodes, _ = load_graph(meta, set(mapped.values()))
    adj = defaultdict(list)
    for a, b in edges:
        adj[a].append(b)

    families = meta["target_families"]
    endpoints = [(t, family) for family, targets in families.items() for t in targets]
    mapped_endpoints = [t for t, _ in endpoints if t in mapped]
    assert len(endpoints) == 13 and len(mapped_endpoints) == 9
    assert all(s in mapped for s in BLOCKERS)

    accession_like = sorted(n for n in all_nodes if n not in graph_seed and ACCESSION_LIKE.fullmatch(n))
    accession_set = set(accession_like)
    baseline = {}
    no_accession = {}
    deleted = {}
    for source in BLOCKERS:
        src = mapped[source]
        target_ids = [mapped[t] for t in mapped_endpoints]
        baseline[source] = bounded_states(adj, graph_seed, src, target_ids)
        no_accession[source] = bounded_states(adj, graph_seed, src, target_ids, accession_set)
        for blocker in BLOCKERS[source]:
            if blocker not in all_nodes:
                raise ValueError(f"Blocked gene missing from graph: {blocker}")
            deleted[source, blocker] = bounded_states(
                adj, graph_seed, src, target_ids, {blocker}
            )

    # Reproduce EVERY F9I-B source-endpoint-setting baseline flag,
    # including the four F2 endpoints absent from the original graph.
    originals = rows(F9IB / "F9I_endpoint_exact_reachability.csv")
    assert len(originals) == len(BLOCKERS) * len(endpoints)
    checked = 0
    for record in originals:
        source, target = record["Source"], record["Target"]
        assert source in BLOCKERS
        for setting in SETTINGS:
            actual = (reachable(baseline[source][mapped[target]], setting)
                      if target in mapped else False)
            assert actual == truth(record["Reachable_" + setting[0]]), (
                "F9I-B reproduction failed", source, target, setting[0]
            )
            checked += 1

    result = []
    for source, blockers in BLOCKERS.items():
        for target, family in endpoints:
            target_id = mapped.get(target)
            for setting in SETTINGS:
                if target_id is None:
                    continue
                original = reachable(baseline[source][target_id], setting)
                after_accession = reachable(no_accession[source][target_id], setting)
                for blocker in blockers:
                    after = reachable(deleted[source, blocker][target_id], setting)
                    result.append({
                        "Source": source,
                        "Target": target,
                        "Target_family": family,
                        "Setting": setting[0],
                        "Blocked_intermediate": blocker,
                        "Baseline_reachable": original,
                        "After_deletion_reachable": after,
                        "Exact_required_at_setting": original and not after,
                        "After_accession_filter_reachable": after_accession,
                        "Accession_filter_changes_reachability": original and not after_accession,
                        "Interpretation": ("REQUIRED" if original and not after
                            else "ALTERNATIVE_ROUTE_REMAINS" if original
                            else "BASELINE_UNREACHABLE")
                    })

    # Verify the relevant original F9I-B deletion tests, including its
    # exact positives and the negatives for our selected blocker genes.
    originals_deletions = rows(F9IB / "F9I_node_deletion_tests.csv")
    relevant = [r for r in originals_deletions if
                r["Candidate_intermediate"] in BLOCKERS.get(r["Source"], ())]
    lookup = {(r["Source"], r["Target"], r["Setting"], r["Blocked_intermediate"]): r
              for r in result}
    for old in relevant:
        key = (old["Source"], old["Target"], old["Setting"], old["Candidate_intermediate"])
        assert key in lookup and truth(old["Exact_required"]) == lookup[key]["Exact_required_at_setting"], key

    OUT.mkdir(parents=True, exist_ok=True)
    save_csv(OUT / "F9I_C_dependency_settings.csv", result, list(result[0]))
    save_csv(OUT / "F9I_C_accession_like_identifiers.csv",
             [{"Identifier": g, "F2_seed": False} for g in accession_like],
             ["Identifier", "F2_seed"])

    unique_reqs = sorted({(r["Source"], r["Blocked_intermediate"]) for r in result
                          if r["Exact_required_at_setting"]})
    changes = sorted({(r["Source"], r["Target"], r["Setting"]) for r in result
                      if r["Accession_filter_changes_reachability"]})
    report = {
        "analysis": "F9I-C exact deletion sensitivity (no statistical inference)",
        "prior_nodes": len(all_nodes), "prior_edges": len(edges),
        "F9I_B_baseline_setting_checks_passed": checked,
        "F9I_B_intermediate_deletion_checks_passed": len(relevant),
        "sources": list(BLOCKERS), "blockers": BLOCKERS,
        "mapped_F2_endpoints": mapped_endpoints,
        "setting_rows": len(result),
        "exact_required_source_blocker_pairs": [list(x) for x in unique_reqs],
        "accession_like_identifier_count": len(accession_like),
        "accession_filter_reachability_changes": [list(x) for x in changes],
        "notes": [
            "All deletion conclusions are specific to directed E/X-bounded paths in the original mouse OmniPath graph.",
            "All tested endpoints are F2 genes, making the multi-target exact E/X state search equivalent to pairwise reachability.",
            "Accession-like identifiers are a regex-based heuristic, not validated UniProt matches; inspect the CSV.",
            "No interventions were performed in neuronal cells; these are in-silico node removals.",
            "Tests are conditional sensitivity diagnostics, not p-values or experimental causal identification."
        ]
    }
    (OUT / "F9I_C_validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("F9I-C PASS: original graph and F9I-B checks reproduced")
    print("  Baseline comparisons:", checked, "old deletion comparisons:", len(relevant))
    print("  Required source/blocker combinations:", unique_reqs)
    print("  Accession-like symbols:", len(accession_like),
          "reachability changes:", len(changes))
    print("  Output:", OUT)


if __name__ == "__main__":
    main()
