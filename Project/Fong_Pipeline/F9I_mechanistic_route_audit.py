#!/usr/bin/env python3
"""F9I-B: exact, provenance-preserving directed route audit for five
experimentally annotated INTERFACE candidates.

Save beside F9F_fongwide_directed_connectivity.py in Project/Fong_Pipeline.
Run from Project/: py .\\Fong_Pipeline\\F9I_mechanistic_route_audit.py
Optional:         py .\\Fong_Pipeline\\F9I_mechanistic_route_audit.py --self-test

Uses the ORIGINAL corrected F5 cached four mouse OmniPath tables via the
F9F loader (including SHA256 validation). Does not download anything and
never edits F5/F9F. Tests exact bounded reachability; capped routes are
EXAMPLES only and do NOT determine whether a pair is reachable.
"""
from __future__ import annotations

import argparse
import csv
import random
from functools import lru_cache
import json
import math
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "Data/Processed/FongDownstream"
F5 = BASE / "F5_Mechanistic/Corrected"
F9F = BASE / "F9F_FongWideDirected"
F9I = BASE / "F9I_CandidateEvidence"
OUT = BASE / "F9I_MechanisticRoutes"
RAW = ROOT / "Data/Raw/Network/FongDownstream/OmniPath"
CANDIDATES = ("TOPBP1", "CLSPN", "RAD51", "TIMELESS", "RPA2")
SETTINGS = (("E4_X0", 4, 0), ("E4_X1", 4, 1),
            ("E4_X2", 4, 2), ("E5_X2", 5, 2))
TIERS = {"E4_X0": "A", "E4_X1": "B", "E4_X2": "C", "E5_X2": "D"}
RAW_FILES = {
    "SIGNALING": "omnipath_mouse_signaling.tsv",
    "COLLECTRI": "omnipath_mouse_collectri.tsv",
    "DOROTHEA_AB": "omnipath_mouse_dorothea_AB.tsv",
    "TF_TARGET": "omnipath_mouse_tf_target.tsv",
}


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        out = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        out.writeheader()
        out.writerows(rows)


def yes(value):
    return str(value).strip().upper() in {"TRUE", "YES", "1"}


def normalized(v):
    return str(v or "").strip().upper()


def rev_distances(incoming, target, max_depth=5):
    dist = {target: 0}
    queue = deque([target])
    while queue:
        v = queue.popleft()
        if dist[v] >= max_depth:
            continue
        for u in incoming.get(v, ()):
            if u not in dist:
                dist[u] = dist[v] + 1
                queue.append(u)
    return dist


def exact_states(outgoing, source, target, seed_graph, reverse_dist,
                 max_edges=5, max_external=2, blocked_node=None):
    """Same bounded (node, length, outside-count) state formulation as F5.
    A reachable bounded walk implies a simple path satisfying no worse bounds.
    """
    if (source == target or source not in outgoing or
            source not in reverse_dist or blocked_node in {source, target}):
        return set()
    queue = deque([(source, 0, 0)])
    seen = {(source, 0, 0)}
    found = set()
    while queue:
        node, e, x = queue.popleft()
        if e >= max_edges:
            continue
        for v in outgoing.get(node, ()):
            if v == blocked_node:
                continue
            ne = e + 1
            if v != target and reverse_dist.get(v, 999) > max_edges - ne:
                continue
            nx = x + int(v != target and v not in seed_graph)
            if nx > max_external:
                continue
            if v == target:
                found.add((ne, nx))
                continue
            state = (v, ne, nx)
            if state not in seen:
                seen.add(state)
                queue.append(state)
    return found


def sign_label(info):
    mask = info["mask"]
    return {0: "unsigned", 1: "stimulation", 2: "inhibition",
            3: "ambiguous"}[mask]


def path_effect(path, records):
    effects = [records[a, b]["mask"] for a, b in zip(path, path[1:])]
    if any(m in (0, 3) for m in effects):
        return "undetermined"
    return "net_positive" if sum(m == 2 for m in effects) % 2 == 0 else "net_negative"


def choose_paths(outgoing, records, source, target, seed_graph,
                 dist, max_edges, min_x, shortest_at_min_x,
                 max_paths=5, state_cap=200_000):
    """Find a few *simple* routes at exact minimum external support.

    Reverse bounded feasibility prunes the search to paths capable of
    satisfying the exact (edge count, outside-F2 count) constraint. Routes
    are displayed shortest first; signed/curated edge preference within
    each depth is only a display heuristic, not a confidence score.
    """
    max_e = min(max_edges, shortest_at_min_x + 1)
    chosen = []
    chosen_set = set()
    per_first_hop = Counter()
    expansions = 0
    capped = False

    def neighbor_order(u):
        return sorted(outgoing.get(u, ()), key=lambda v: (
            int(records[u, v]["mask"] in (0, 3)),
            -float(records[u, v]["effort"]),
            v,
        ))

    ordered = {u: neighbor_order(u) for u in outgoing}

    for length in range(shortest_at_min_x, max_e + 1):
        @lru_cache(None)
        def can_finish(u, rem_e, rem_x):
            if rem_e < 0 or rem_x < 0:
                return False
            if u == target:
                return rem_e == 0 and rem_x == 0
            if rem_e == 0 or dist.get(u, 999) > rem_e:
                return False
            for v in ordered.get(u, ()):
                cost = int(v != target and v not in seed_graph)
                if cost <= rem_x and can_finish(v, rem_e - 1, rem_x - cost):
                    return True
            return False

        if not can_finish(source, length, min_x):
            continue

        def dfs(path, rem_e, rem_x):
            nonlocal expansions, capped
            if capped or len(chosen) >= max_paths:
                return
            expansions += 1
            if expansions > state_cap:
                capped = True
                return
            u = path[-1]
            if u == target:
                if rem_e == 0 and rem_x == 0 and path not in chosen_set:
                    first_hop = path[1]
                    if per_first_hop[first_hop] < 2:
                        chosen.append(path)
                        chosen_set.add(path)
                        per_first_hop[first_hop] += 1
                return
            if rem_e == 0:
                return
            for v in ordered.get(u, ()):
                if v in path or dist.get(v, 999) > rem_e - 1:
                    continue
                cost = int(v != target and v not in seed_graph)
                if cost > rem_x:
                    continue
                if can_finish(v, rem_e - 1, rem_x - cost):
                    dfs(path + (v,), rem_e - 1, rem_x - cost)
                    if capped or len(chosen) >= max_paths:
                        return

        dfs((source,), length, min_x)
        if capped or len(chosen) >= max_paths:
            break

    return chosen, capped, expansions


def get_provenance(edge_keys):
    """Second pass retains citations only for edges on chosen routes."""
    provenance = {key: {"sources": set(), "references": set()} for key in edge_keys}
    for layer, filename in RAW_FILES.items():
        with (RAW / filename).open(newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fields = {name.lower(): name for name in (reader.fieldnames or [])}
            def find(*opts):
                return next((fields[o] for o in opts if o in fields), None)
            sc = find("source_genesymbol", "source_genesymbols", "source_genesymb")
            tc = find("target_genesymbol", "target_genesymbols", "target_genesymb")
            dc = find("is_directed", "consensus_direction")
            if not sc or not tc:
                raise ValueError(f"Unexpected gene column schema: {filename}")
            for row in reader:
                if dc and not yes(row.get(dc)):
                    continue
                key = (normalized(row.get(sc)), normalized(row.get(tc)))
                if key not in provenance:
                    continue
                for k in ("sources", "references"):
                    provenance[key][k].update(
                        term.strip()
                        for term in (row.get(k) or "").split(";")
                        if term.strip()
                    )
    return provenance


def self_test():
    # Distinct checks: no-intermediate and one-outside-intermediate paths;
    # evidence that removing a required intermediate removes all routes.
    adj = {"A": ("B", "X"), "B": ("C",), "X": ("C",), "C": ("D",),
           "D": (), "Y": ()}
    inc = defaultdict(list)
    for a, dests in adj.items():
        for b in dests:
            inc[b].append(a)
    seeds = {"A", "B", "C", "D"}
    d = rev_distances(inc, "D")
    states = exact_states(adj, "A", "D", seeds, d, 4, 1)
    assert (3, 0) in states and (3, 1) in states
    assert exact_states(adj, "A", "D", seeds, d, 4, 0) == {(3, 0)}
    assert exact_states(adj, "A", "D", seeds, d, 4, 0, blocked_node="B") == set()
    records = {(a, b): {"mask": 1, "effort": 2.0} for a, ds in adj.items() for b in ds}
    paths, hit, _ = choose_paths(adj, records, "A", "D", seeds, d, 4, 0, 3)
    assert ("A", "B", "C", "D") in paths and not hit
    assert path_effect(paths[0], records) == "net_positive"
    # Independent exhaustive simple-path checks on small random graphs.
    for seed in range(12):
        rng = random.Random(seed)
        names = tuple("ABCDEFX")
        sources = set(names[:4])
        random_out = {u: tuple(v for v in names if v != u and
                       rng.random() < 0.26) for u in names}
        reverse = defaultdict(list)
        for u, vs in random_out.items():
            for v in vs:
                reverse[v].append(u)
        for a in names[:3]:
            for b in names[3:6]:
                distance = rev_distances(reverse, b)
                exact = exact_states(random_out, a, b, sources, distance, 5, 2)
                brute = set()
                stack = [(a, (a,), 0)]
                while stack:
                    u, path, x = stack.pop()
                    if u == b:
                        brute.add((len(path) - 1, x))
                        continue
                    if len(path) > 5:
                        continue
                    for v in random_out[u]:
                        if v not in path:
                            nx = x + int(v != b and v not in sources)
                            if nx <= 2:
                                stack.append((v, path + (v,), nx))
                for _, max_e, max_x in SETTINGS:
                    assert (any(e <= max_e and x <= max_x for e, x in exact)
                            == any(e <= max_e and x <= max_x for e, x in brute)), (seed, a, b)
    print("F9I-B SELF-TEST PASSED: random exhaustive reachability, deletion, routes.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--max-routes", type=int, default=5)
    parser.add_argument("--state-cap", type=int, default=200_000)
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not (1 <= args.max_routes <= 25):
        parser.error("--max-routes must be between 1 and 25")

    # The F9F loader is deliberately reused to enforce the SAME original
    # SHA256, directed-row, graph-node and graph-edge checks as F9F.
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from F9F_fongwide_directed_connectivity import load_graph

    meta = json.loads((F5 / "F5_corrected_metadata.json").read_text(encoding="utf-8"))
    genes = {row["Gene"]: row for row in read_csv(F9F / "F9F_gene_connectivity.csv")}
    assert len(genes) == 321
    mapped = {name: row["Graph_identifier"] for name, row in genes.items()
              if yes(row["Mapped_in_prior"])}
    assert len(mapped) == 257
    seeds = set(genes)
    # Identical to F5's recorded H2AFX alias-to-experimental-seed treatment.
    for canonical, aliases in meta["anchor_aliases"].items():
        if canonical in mapped and mapped[canonical] in aliases:
            seeds.add(mapped[canonical])
    assert len(seeds) == meta["graph_seed_gene_count_after_alias_resolution"]

    print("Reconstructing original F5 prior graph with hash validation...")
    records, all_nodes, layer_counts = load_graph(meta, set(mapped.values()))
    outgoing = defaultdict(list)
    incoming = defaultdict(list)
    for a, b in records:
        outgoing[a].append(b)
        incoming[b].append(a)
    for u in outgoing:
        outgoing[u].sort()
    print(f"Validated {len(all_nodes):,} nodes; {len(records):,} edges.")

    target_families = meta["target_families"]
    targets = [(name, family) for family, names in target_families.items() for name in names]
    target_dist = {t: rev_distances(incoming, mapped[t])
                   for t, _ in targets if t in mapped}
    f9f_pairs = {(r["Source"], r["Target"]): r
                 for r in read_csv(F9F / "F9F_pair_reachability.csv")
                 if r["Source"] in CANDIDATES and r["Target"] in {t for t, _ in targets}}
    assert len(f9f_pairs) == len(CANDIDATES) * len(targets)
    originals = {(r["Setting"], r["Source"], r["Target"]): r
                 for r in read_csv(F5 / "pair_setting_summary.csv")
                 if r["Source"] in CANDIDATES}

    outcome = []
    route_rows = []
    edge_rows = []
    deletion_rows = []
    selected_edges = set()
    validated_f9f = 0
    validated_f5 = 0
    caps_hit = 0
    route_count = 0
    for source in CANDIDATES:
        for target, family in targets:
            old = f9f_pairs[source, target]
            source_id = mapped.get(source)
            target_id = mapped.get(target)
            if source_id and target_id:
                states = exact_states(outgoing, source_id, target_id, seeds,
                                      target_dist[target], 5, 2)
            else:
                states = set()
            flags = {setting: any(e <= max_e and x <= max_x for e, x in states)
                     for setting, max_e, max_x in SETTINGS}
            for setting, _, _ in SETTINGS:
                expected = yes(old["Reachable_" + setting])
                if flags[setting] != expected:
                    raise ValueError(f"F9F exact reachability mismatch: {source}->{target} {setting}: {flags[setting]} != {expected}")
                validated_f9f += 1
                prev = originals.get((setting, source, target))
                if prev is not None:
                    if flags[setting] != yes(prev["Reachable_exact"]):
                        raise ValueError(f"Original F5 discrepancy: {setting}: {source}->{target}")
                    validated_f5 += 1
            tier_setting = next((name for name, _, _ in SETTINGS if flags[name]), None)
            base = dict(Source=source, Target=target, Target_family=family,
                        Source_mapped=bool(source_id), Target_mapped=bool(target_id),
                        Support_tier=TIERS[tier_setting] if tier_setting else
                        "UNREACHABLE" if source_id and target_id else "UNTESTABLE",
                        Strictest_setting=tier_setting or "",
                        **{"Reachable_" + setting: flags[setting] for setting, _, _ in SETTINGS})
            selected = []
            if tier_setting:
                setting, ecap, xcap = next(item for item in SETTINGS if item[0] == tier_setting)
                min_x = min(x for e, x in states if e <= ecap and x <= xcap)
                shortest_at_min_x = min(e for e, x in states if e <= ecap and x == min_x)
                selected, capped, explored = choose_paths(
                    outgoing, records, source_id, target_id, seeds, target_dist[target],
                    ecap, min_x, shortest_at_min_x,
                    max_paths=args.max_routes, state_cap=args.state_cap,
                )
                caps_hit += int(capped)
                if not selected:
                    raise ValueError(f"No representative found despite exact reachability: {source}->{target} {setting}")
                base.update(Min_external_exact=min_x,
                            Shortest_edges_at_min_external_exact=shortest_at_min_x,
                            Representative_routes=len(selected), Search_state_cap_hit=capped)
                # F5's exact required-node logic: a truly required internal
                # node must be present in every retained representative route.
                # Use that only to nominate candidates; deletion is tested
                # against the COMPLETE graph, not the representative list.
                shared_intermediates = set.intersection(
                    *(set(path[1:-1]) for path in selected)
                )
                for node in sorted(shared_intermediates):
                    after_deletion = exact_states(
                        outgoing, source_id, target_id, seeds, target_dist[target],
                        ecap, xcap, blocked_node=node
                    )
                    deletion_rows.append(dict(
                        Source=source, Target=target, Target_family=family,
                        Setting=setting, Candidate_intermediate=node,
                        Intermediate_class="F2" if node in seeds else "outside_F2",
                        Exact_required=(not bool(after_deletion)),
                        Test_definition="Deleting this node removes all routes satisfying E/X"
                    ))

                for i, path in enumerate(selected, 1):
                    route_count += 1
                    path_id = f"F9IB_{source}_{target}_{i:02d}"
                    outside = [node for node in path[1:-1] if node not in seeds]
                    unknown_edges = sum(records[a, b]["mask"] in (0, 3)
                                        for a, b in zip(path, path[1:]))
                    used_layers = sorted({layer for a, b in zip(path, path[1:])
                                          for layer in records[a, b]["layers"]})
                    route_rows.append(dict(Path_ID=path_id, Source=source,
                        Target=target, Target_family=family, Support_tier=TIERS[setting],
                        Strictest_setting=setting, Route=" -> ".join(path),
                        Edges=len(path)-1, External_intermediates=";".join(outside),
                        External_count=len(outside), Unknown_or_ambiguous_edges=unknown_edges,
                        Net_prior_sign=path_effect(path, records), Layers=";".join(used_layers),
                        Display_only_not_exhaustive=True))
                    for step, (a, b) in enumerate(zip(path, path[1:]), 1):
                        info = records[a, b]
                        selected_edges.add((a, b))
                        edge_rows.append(dict(Path_ID=path_id, Step=step, Source=a, Target=b,
                            Effect_label=sign_label(info), Layers=";".join(sorted(info["layers"])),
                            Supporting_records=info["count"], Max_curation_effort=info["effort"]))
            else:
                base.update(Min_external_exact="", Shortest_edges_at_min_external_exact="",
                            Representative_routes=0, Search_state_cap_hit=False)
            outcome.append(base)
    # Regression control from corrected F5: TOPBP1 -> PYCARD under E4_X0
    # requires E2F1. Recheck its original required-node evidence on the full
    # graph; it should not depend on which example routes were retained.
    original_required = [r for r in read_csv(F5 / "exact_required_nodes.csv")
                         if yes(r.get("Exact_required")) and
                         r["Source"] in CANDIDATES]
    validated_required = 0
    for r in original_required:
        source, target, setting = r["Source"], r["Target"], r["Setting"]
        if target not in mapped:
            continue
        _, ecap, xcap = next(x for x in SETTINGS if x[0] == setting)
        after = exact_states(outgoing, mapped[source], mapped[target], seeds,
                             target_dist[target], ecap, xcap,
                             blocked_node=r["Node"])
        if after:
            raise ValueError(f"F5 exact-required control failed: {setting} {source}->{target} without {r['Node']}")
        validated_required += 1

    # Add original evidence sources and publication references for EACH route edge.
    print(f"Reconstructing provenance for {len(selected_edges)} distinct route edges...")
    provenance = get_provenance(selected_edges)
    for edge in edge_rows:
        info = provenance[edge["Source"], edge["Target"]]
        edge["Sources"] = ";".join(sorted(info["sources"]))
        edge["References"] = ";".join(sorted(info["references"]))
    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / "F9I_endpoint_exact_reachability.csv", outcome,
              ["Source", "Target", "Target_family", "Source_mapped", "Target_mapped",
               "Support_tier", "Strictest_setting"] +
              ["Reachable_" + s for s, _, _ in SETTINGS] +
              ["Min_external_exact", "Shortest_edges_at_min_external_exact",
               "Representative_routes", "Search_state_cap_hit"])
    write_csv(OUT / "F9I_representative_routes.csv", route_rows,
              ["Path_ID", "Source", "Target", "Target_family", "Support_tier", "Strictest_setting",
               "Route", "Edges", "External_intermediates", "External_count",
               "Unknown_or_ambiguous_edges", "Net_prior_sign", "Layers",
               "Display_only_not_exhaustive"])
    write_csv(OUT / "F9I_node_deletion_tests.csv", deletion_rows,
              ["Source", "Target", "Target_family", "Setting", "Candidate_intermediate",
               "Intermediate_class", "Exact_required", "Test_definition"])
    write_csv(OUT / "F9I_route_edge_provenance.csv", edge_rows,
              ["Path_ID", "Step", "Source", "Target", "Effect_label", "Layers",
               "Supporting_records", "Max_curation_effort", "Sources", "References"])
    summary = dict(analysis="F9I-B directed candidate to F5 endpoint audit",
                   original_prior_nodes=len(all_nodes), original_prior_edges=len(records),
                   original_raw_files_sha256_validated=True, candidates=list(CANDIDATES),
                   endpoints={family: members for family, members in target_families.items()},
                   source_target_pairs=len(outcome),
                   testable_pairs=sum(bool(r["Source_mapped"] and r["Target_mapped"]) for r in outcome),
                   tier_counts=dict(Counter(r["Support_tier"] for r in outcome)),
                   cross_checks_F9F=validated_f9f, cross_checks_original_F5=validated_f5,
                   representative_routes=len(route_rows), route_unique_edges=len(selected_edges),
                   proposed_nodes_exactly_tested=len(deletion_rows),
                   exact_required_intermediates=sum(yes(r["Exact_required"]) for r in deletion_rows),
                   original_F5_exact_required_controls_passed=validated_required,
                   path_state_cap_hit_pairs=caps_hit, raw_layer_directed_counts=layer_counts,
                   limitations=["Exact route existence is limited to bounded F5 E/X constraints and original mouse OmniPath graph.",
                     "Absence of a graph node is coverage failure, not biological absence.",
                     "Representative routes are capped illustrations, selected using display heuristics, not confidence-ranked mechanisms.",
                     "Edge sign is prior-database annotation; net sign is undetermined for unsigned or ambiguous edges.",
                     "Neither directed prior paths nor Fong expression establish causal activity in neuronal tissue.",
                     "Source gene prioritization was post-hoc and requires direct experimental validation."])
    (OUT / "F9I_route_validation.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("F9I-B PASS: F9F and original F5 reachability reproduced")
    print("Exact pair checks:", validated_f9f, "F9F and", validated_f5, "original F5")
    print("Support tiers:", summary["tier_counts"])
    print("Illustrative routes:", len(route_rows), "/", len(outcome), "candidate-endpoint pairs")
    print("Exact required-node tests:", len(deletion_rows),
          "found:", summary["exact_required_intermediates"],
          "original F5 controls:", validated_required)
    print("Capped route searches:", caps_hit)
    print("Outputs:", OUT)


if __name__ == "__main__":
    main()
