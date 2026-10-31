
#!/usr/bin/env python3
"""
F9F — Fong-wide directed functional connectivity.

Exact bounded reachability between all ordered pairs of
321 selected Fong F2 genes, using the ORIGINAL corrected
F5 mouse OmniPath network.

Settings:
    E4_X0, E4_X1, E4_X2, E5_X2

Algorithm:
    Multi-source bitset dynamic programming.

Validation:
    - Original F5 raw-cache SHA256 hashes
    - Original F5 graph dimensions and edge counts
    - Original F5 seed/alias conventions
    - Original F5 reachability results for eligible pairs

Standard library only. No downloads.

Run from Project:
    py .\\Fong_Pipeline\\F9F_fongwide_directed_connectivity.py

Optional:
    py .\\Fong_Pipeline\\F9F_fongwide_directed_connectivity.py --self-test
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random

from collections import Counter, defaultdict
from pathlib import Path


# ============================================================
# PATHS
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

BASE = ROOT / "Data/Processed/FongDownstream"

RAW = (
    ROOT
    / "Data/Raw/Network/FongDownstream/OmniPath"
)

F5 = BASE / "F5_Mechanistic/Corrected"

F8B = (
    BASE
    / "F8_FunctionalAtlas"
    / "F8B_BiologicalModules"
    / "F8B_gene_role_atlas.csv"
)

OUT = BASE / "F9F_FongWideDirected"


INPUTS = {
    "SIGNALING":
        "omnipath_mouse_signaling.tsv",
    "COLLECTRI":
        "omnipath_mouse_collectri.tsv",
    "DOROTHEA_AB":
        "omnipath_mouse_dorothea_AB.tsv",
    "TF_TARGET":
        "omnipath_mouse_tf_target.tsv",
}

HASH_KEYS = {
    "SIGNALING": "signaling_sha256",
    "COLLECTRI": "collectri_sha256",
    "DOROTHEA_AB": "dorothea_sha256",
    "TF_TARGET": "tf_target_sha256",
}

SETTINGS = (
    "E4_X0",
    "E4_X1",
    "E4_X2",
    "E5_X2",
)

NO_PATH = 255


# ============================================================
# UTILITIES
# ============================================================

def rows(path, sep=","):

    with path.open(
        "r",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        yield from csv.DictReader(
            f,
            delimiter=sep
        )


def write_csv(path, data, fields):

    with path.open(
        "w",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fields,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(data)


def write_json(path, obj):

    path.write_text(
        json.dumps(
            obj,
            indent=2,
            ensure_ascii=False
        ),
        encoding="utf-8"
    )


def true(v):

    return str(v).strip().upper() in {
        "TRUE", "1", "YES", "Y"
    }


def number(v):

    try:
        x = float(v)

        return (
            x if math.isfinite(x)
            else 0.0
        )

    except (ValueError, TypeError):
        return 0.0


def symbol(v):

    s = (v or "").strip().upper()

    return (
        ""
        if s in {
            "", "NA", "N/A", "NONE", "NULL"
        }
        else s
    )


def sha256(path):

    h = hashlib.sha256()

    with path.open("rb") as f:

        for block in iter(
            lambda: f.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def effect(row):
    """
    Match the original F5 regulatory-sign rules.

    Consensus flags take precedence when present.
    """

    a = true(
        row.get("consensus_stimulation")
    )

    b = true(
        row.get("consensus_inhibition")
    )

    if a or b:
        return (
            1 if a and not b
            else -1 if b and not a
            else 0
        )

    a = true(
        row.get("is_stimulation")
    )

    b = true(
        row.get("is_inhibition")
    )

    return (
        1 if a and not b
        else -1 if b and not a
        else 0
    )


# ============================================================
# RECONSTRUCT THE ORIGINAL F5 PRIOR NETWORK
# ============================================================

def load_graph(metadata, keep_provenance_for):

    records = {}
    directed_rows = Counter()

    original_meta = metadata["raw_metadata"]

    for layer, filename in INPUTS.items():

        path = RAW / filename

        if not path.exists():
            raise FileNotFoundError(
                f"Missing original F5 cache: {path}"
            )

        observed_sha = sha256(path)

        expected_sha = original_meta[
            HASH_KEYS[layer]
        ]

        if observed_sha != expected_sha:

            raise ValueError(
                f"F5 cache SHA mismatch: {filename}\n"
                f"Expected {expected_sha}\n"
                f"Observed {observed_sha}"
            )

        it = iter(
            rows(path, sep="\t")
        )

        try:
            first = next(it)

        except StopIteration:
            raise ValueError(
                f"Empty source: {path}"
            )

        lower = {
            k.lower(): k
            for k in first
        }

        def col(names):

            return next(
                (
                    lower[x]
                    for x in names
                    if x in lower
                ),
                None
            )

        sc = col((
            "source_genesymbol",
            "source_genesymbols",
            "source_genesymb"
        ))

        tc = col((
            "target_genesymbol",
            "target_genesymbols",
            "target_genesymb"
        ))

        dc = col((
            "is_directed",
            "consensus_direction"
        ))

        if not sc or not tc:

            raise ValueError(
                f"Gene-symbol columns missing "
                f"in {filename}: {list(first)}"
            )

        import itertools

        for row in itertools.chain(
            (first,),
            it
        ):

            if (
                dc is not None
                and not true(row.get(dc))
            ):
                continue

            a = symbol(row.get(sc))
            b = symbol(row.get(tc))

            if (
                not a
                or not b
                or a == b
            ):
                continue

            directed_rows[layer] += 1

            key = (a, b)

            rec = records.get(key)

            if rec is None:

                rec = {
                    "mask": 0,
                    "layers": set(),
                    "count": 0,
                    "effort": 0.0,
                    "sources": set(),
                    "references": set()
                }

                records[key] = rec

            e = effect(row)

            if e == 1:
                rec["mask"] |= 1

            elif e == -1:
                rec["mask"] |= 2

            rec["layers"].add(layer)
            rec["count"] += 1

            rec["effort"] = max(
                rec["effort"],
                number(
                    row.get("curation_effort")
                )
            )

            # Preserve detailed provenance where
            # both endpoints may map to F2.
            # Avoid storing full citations for
            # all 175k edges unnecessarily.

            if (
                a in keep_provenance_for
                and b in keep_provenance_for
            ):

                for field in (
                    "sources",
                    "references"
                ):

                    value = row.get(field) or ""

                    for item in value.split(";"):

                        if item.strip():
                            rec[field].add(
                                item.strip()
                            )

        print(
            f"  {layer}: "
            f"{directed_rows[layer]:,} "
            f"directed records"
        )

    nodes = set()

    for a, b in records:
        nodes.add(a)
        nodes.add(b)

    expected = metadata["network"]

    if (
        len(nodes) != expected["nodes"]
        or len(records) != expected["edges"]
    ):

        raise ValueError(
            "Reconstructed F5 network differs "
            "from metadata: "
            f"{len(nodes)} nodes / "
            f"{len(records)} edges; expected "
            f"{expected['nodes']} / "
            f"{expected['edges']}."
        )

    by_layer_meta = {
        "SIGNALING":
            "directed_signaling_rows",
        "COLLECTRI":
            "directed_collectri_rows",
        "DOROTHEA_AB":
            "directed_dorothea_AB_rows",
        "TF_TARGET":
            "directed_tf_target_rows"
    }

    for layer, key in by_layer_meta.items():

        if directed_rows[layer] != expected[key]:

            raise ValueError(
                f"F5 directed record count mismatch "
                f"for {layer}: "
                f"{directed_rows[layer]} "
                f"versus {expected[key]}"
            )

    return (
        records,
        nodes,
        dict(directed_rows)
    )


# ============================================================
# EXACT MULTI-SOURCE BOUNDED REACHABILITY
# ============================================================

def compute_reachability(
    edges,
    node_ids,
    selected_sources,
    seed_graph
):
    """
    Exact bounded reachability under the four E/X settings.

    Each state tracks:
        - graph node
        - edges used
        - non-F2 intermediate count

    A Python integer bitmask stores every F2 source
    capable of reaching that state.

    The algorithm propagates all sources simultaneously.

    Returns:
        minlen[x][source_index * N + target_index]

    for exact external counts x = 0, 1, 2.

    We use these values only to answer bounded
    <=E and <=X reachability questions.

    Walk cycles do not create false positive bounded
    reachability: removing cycles cannot increase E or X.
    """

    gnodes = {
        g: i
        for i, g in enumerate(
            sorted(node_ids)
        )
    }

    graph_edges = [
        (
            gnodes[a],
            gnodes[b],
            int(b not in seed_graph)
        )
        for a, b in edges
    ]

    sources = [
        (
            gnodes[g]
            if g is not None
            else None
        )
        for g in selected_sources
    ]

    n = len(sources)

    minlen = [
        bytearray([NO_PATH]) * (n * n)
        for _ in range(3)
    ]

    frontier = [
        {},
        {},
        {}
    ]

    for si, g in enumerate(sources):

        if g is not None:

            frontier[0][g] = (
                frontier[0].get(g, 0)
                | (1 << si)
            )

    for step in range(1, 6):

        following = [
            {},
            {},
            {}
        ]

        for u, v, external in graph_edges:

            if external:

                for x in (0, 1):

                    bits = frontier[x].get(
                        u, 0
                    )

                    if bits:

                        following[x + 1][v] = (
                            following[x + 1].get(
                                v, 0
                            )
                            | bits
                        )

            else:

                for x in (0, 1, 2):

                    bits = frontier[x].get(
                        u, 0
                    )

                    if bits:

                        following[x][v] = (
                            following[x].get(
                                v, 0
                            )
                            | bits
                        )

        # Extract first reachable length for
        # every F2 source -> F2 target.

        for target_idx, target_graph_id in enumerate(
            sources
        ):

            if target_graph_id is None:
                continue

            for x in (0, 1, 2):

                bits = following[x].get(
                    target_graph_id,
                    0
                )

                arr = minlen[x]

                while bits:

                    bit = bits & -bits

                    si = bit.bit_length() - 1

                    pos = (
                        si * n
                        + target_idx
                    )

                    if arr[pos] == NO_PATH:
                        arr[pos] = step

                    bits ^= bit

        frontier = following

        print(
            "  Exact propagation: <=",
            step,
            "edge(s); reachable states:",
            sum(
                len(x)
                for x in following
            )
        )

    return minlen


# ============================================================
# CONSTRAINED REACHABILITY METRICS
# ============================================================

def shortest(
    values,
    edge_cap,
    external_cap
):

    valid = [
        values[x]
        for x in range(external_cap + 1)
        if values[x] <= edge_cap
    ]

    if not valid:
        return None, None, None

    shortest_any = min(valid)

    min_external = next(
        x
        for x in range(external_cap + 1)
        if values[x] <= edge_cap
    )

    return (
        min_external,
        shortest_any,
        values[min_external]
    )


def tier_of(v):

    if v[0] <= 4:
        return "A"

    if min(v[0], v[1]) <= 4:
        return "B"

    if min(v) <= 4:
        return "C"

    if min(v) <= 5:
        return "D"

    return "UNREACHABLE"


# ============================================================
# INDEPENDENT ALGORITHM SELF-TEST
# ============================================================

def self_test():
    """
    Compare bounded DP reachability with exhaustive
    simple-path enumeration on small random graphs.
    """

    for seed in range(8):

        rng = random.Random(seed)

        names = [
            str(i)
            for i in range(8)
        ]

        chosen = names[:4]

        edges = {
            (a, b): {}
            for a in names
            for b in names
            if (
                a != b
                and rng.random() < 0.30
            )
        }

        obs = compute_reachability(
            edges,
            names,
            chosen,
            set(chosen)
        )

        n = len(chosen)

        adj = defaultdict(list)

        for a, b in edges:
            adj[a].append(b)

        for i, a in enumerate(chosen):

            for j, b in enumerate(chosen):

                if i == j:
                    continue

                found = [
                    255,
                    255,
                    255
                ]

                stack = [
                    (a, [a], 0)
                ]

                while stack:

                    u, path, x = stack.pop()

                    if u == b:

                        found[x] = min(
                            found[x],
                            len(path) - 1
                        )

                        continue

                    if len(path) == 6:
                        continue

                    for v in adj[u]:

                        if v in path:
                            continue

                        nx = (
                            x
                            + int(
                                v not in chosen
                                and v != b
                            )
                        )

                        if nx <= 2:

                            stack.append(
                                (
                                    v,
                                    path + [v],
                                    nx
                                )
                            )

                actual = [
                    obs[x][i * n + j]
                    for x in range(3)
                ]

                # Compare bounded reachability,
                # not arbitrary exact-X walk states.

                for emax, xmax in (
                    (4, 0),
                    (4, 1),
                    (4, 2),
                    (5, 2)
                ):

                    got = any(
                        actual[k] <= emax
                        for k in range(xmax + 1)
                    )

                    exp = any(
                        found[k] <= emax
                        for k in range(xmax + 1)
                    )

                    if got != exp:

                        raise AssertionError(
                            (
                                seed,
                                a,
                                b,
                                emax,
                                xmax,
                                actual,
                                found
                            )
                        )

    print(
        "SELF TEST PASSED: bitset DP matches "
        "exhaustive simple-path reachability."
    )


# ============================================================
# MAIN ANALYSIS
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--self-test",
        action="store_true"
    )

    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    meta = json.loads(
        (
            F5 / "F5_corrected_metadata.json"
        ).read_text(
            encoding="utf-8"
        )
    )

    assert (
        meta["phase"]
        == "F5_corrected_omnipath_downstream"
    )

    assert tuple(meta["settings"]) == SETTINGS

    f2 = list(
        rows(F8B)
    )

    assert len(f2) == 321

    lookup = {
        symbol(r["Gene"]): r
        for r in f2
    }

    assert len(lookup) == 321

    assert Counter(
        r["Tier"]
        for r in f2
    ) == {
        "CORE": 194,
        "INTERFACE": 127
    }

    print(
        "Reconstructing and validating "
        "original F5 prior network..."
    )

    original_set = set(lookup)

    potential_aliases = {
        symbol(v)
        for values in meta[
            "anchor_aliases"
        ].values()
        for v in values
    }

    edges, nodes, layer_counts = load_graph(
        meta,
        original_set | potential_aliases
    )

    # Match F5 anchor alias resolution.
    # No generic fuzzy gene-symbol substitution.

    additional_aliases = {}

    for anchor in [
        g
        for group in (
            "source_families",
            "target_families"
        )
        for items in meta[group].values()
        for g in items
    ]:

        g = symbol(anchor)

        if g in nodes:
            continue

        matches = [
            symbol(v)
            for v in meta[
                "anchor_aliases"
            ].get(g, [])
            if symbol(v) in nodes
        ]

        if matches and g in original_set:

            additional_aliases[g] = matches[0]

    graph_seeds = (
        original_set
        | set(
            additional_aliases.values()
        )
    )

    if len(graph_seeds) != meta[
        "graph_seed_gene_count_after_alias_resolution"
    ]:

        raise ValueError(
            "Seed/alias count mismatch: "
            "confirm original F5 alias audit."
        )

    canonical = list(lookup)

    resolved = [
        (
            g
            if g in nodes
            else additional_aliases.get(g)
        )
        for g in canonical
    ]

    mapped = {
        a: b
        for a, b in zip(
            canonical,
            resolved
        )
        if b is not None
    }

    if len(set(mapped.values())) != len(mapped):

        raise ValueError(
            "Two distinct F2 genes map to "
            "the same graph identifier."
        )

    if not mapped:

        raise ValueError(
            "No F2 genes found in original "
            "prior graph."
        )

    print(
        f"Full prior: {len(nodes):,} nodes; "
        f"{len(edges):,} directed edges"
    )

    print(
        f"F2 gene coverage: {len(mapped)}/321; "
        f"aliases: {additional_aliases}"
    )

    # --------------------------------------------------------
    # ALL-F2 MULTI-SOURCE COMPUTATION
    # --------------------------------------------------------

    print(
        "Computing exact bounded reachability "
        "for ALL sources simultaneously..."
    )

    mins = compute_reachability(
        edges,
        nodes,
        resolved,
        graph_seeds
    )

    n = len(canonical)

    pair_rows = []
    pair_flags = {}

    summary = Counter()
    summary_directions = Counter()

    summaries = {
        g: Counter()
        for g in canonical
    }

    present = {
        g: bool(mapped.get(g))
        for g in canonical
    }

    # Selected F8B biological modules may overlap.

    memberships = {
        g: [
            m.strip()
            for m in lookup[g][
                "Biological_modules"
            ].split(";")
            if (
                m.strip()
                and m.strip() != "UNRESOLVED"
            )
        ]
        for g in canonical
    }

    module_summary = defaultdict(Counter)

    for i, source in enumerate(canonical):

        for j, target in enumerate(canonical):

            if i == j:
                continue

            testable = (
                present[source]
                and present[target]
            )

            v = tuple(
                mins[x][i * n + j]
                for x in range(3)
            )

            flags = [
                any(
                    v[k] <= e
                    for k in range(x + 1)
                )
                for e, x in (
                    (4, 0),
                    (4, 1),
                    (4, 2),
                    (5, 2)
                )
            ]

            if not testable:
                assert not any(flags)

            tier = (
                tier_of(v)
                if testable
                else "UNTESTABLE"
            )

            (
                min_x,
                shortest_any,
                shortest_at_min_x
            ) = shortest(
                v,
                5,
                2
            )

            pair_rows.append({

                "Source": source,
                "Target": target,

                "Source_tier":
                    lookup[source]["Tier"],

                "Target_tier":
                    lookup[target]["Tier"],

                "Source_direction":
                    lookup[source]["Direction"],

                "Target_direction":
                    lookup[target]["Direction"],

                "Source_log2FC":
                    lookup[source]["log2FC"],

                "Target_log2FC":
                    lookup[target]["log2FC"],

                "Testable": testable,

                "Reachable_E4_X0": flags[0],
                "Reachable_E4_X1": flags[1],
                "Reachable_E4_X2": flags[2],
                "Reachable_E5_X2": flags[3],

                "Support_tier": tier,

                "Min_external_E5_X2":
                    min_x,

                "Shortest_edges_E5_X2":
                    shortest_any,

                "Shortest_edges_at_min_external_E5_X2":
                    shortest_at_min_x,
            })

            pair_flags[
                (source, target)
            ] = (
                testable,
                flags
            )

            summary[tier] += 1

            if testable:

                summary_directions[
                    (
                        lookup[source]["Direction"],
                        lookup[target]["Direction"],
                        tier
                    )
                ] += 1

                summaries[source][
                    "out_testable"
                ] += 1

                summaries[target][
                    "in_testable"
                ] += 1

                for k, s in enumerate(SETTINGS):

                    if flags[k]:

                        summaries[source][
                            "out_" + s
                        ] += 1

                        summaries[target][
                            "in_" + s
                        ] += 1

            # Descriptive program connectivity.
            # Overlapping modules are not
            # independent observations.

            for sm in memberships[source]:

                for tm in memberships[target]:

                    c = module_summary[
                        (sm, tm)
                    ]

                    c["total_pairs"] += 1

                    if testable:

                        c["testable_pairs"] += 1

                        for k, s in enumerate(
                            SETTINGS
                        ):

                            if flags[k]:

                                c[
                                    "reachable_" + s
                                ] += 1

    assert len(pair_rows) == 321 * 320

    # --------------------------------------------------------
    # VALIDATE AGAINST ORIGINAL CORRECTED F5 RESULTS
    # --------------------------------------------------------

    print(
        "Historical corrected F5 cross-check "
        "(pairs with BOTH endpoints in F2)..."
    )

    mismatches = []
    compared = 0

    for row in rows(
        F5 / "pair_setting_summary.csv"
    ):

        a = symbol(row["Source"])
        b = symbol(row["Target"])

        if (
            a not in lookup
            or b not in lookup
            or a == b
        ):
            continue

        got_testable, flags = pair_flags[
            (a, b)
        ]

        s = row["Setting"]

        j = SETTINGS.index(s)

        expected_testable = true(
            row["Testable"]
        )

        expected_reachable = true(
            row["Reachable_exact"]
        )

        compared += 1

        if (
            expected_testable != got_testable
            or expected_reachable != flags[j]
        ):

            mismatches.append(
                (
                    s, a, b,
                    expected_testable,
                    got_testable,
                    expected_reachable,
                    flags[j]
                )
            )

    if mismatches:

        raise ValueError(
            "F5 reproduction failures: "
            f"{mismatches[:12]}"
        )

    if not compared:

        raise ValueError(
            "No eligible original F5 pairs "
            "for cross-validation."
        )

    print(
        f"  Original pair-setting checks "
        f"passed: {compared}"
    )

    # --------------------------------------------------------
    # DIRECT F2-F2 PRIOR INTERACTIONS
    # --------------------------------------------------------

    direct_f2 = []

    inverse = {
        v: k
        for k, v in mapped.items()
    }

    for (a, b), info in edges.items():

        if (
            a not in inverse
            or b not in inverse
            or a == b
        ):
            continue

        effect_mask = info["mask"]

        label = (
            "ambiguous"
            if effect_mask == 3
            else "stimulation"
            if effect_mask == 1
            else "inhibition"
            if effect_mask == 2
            else "unsigned"
        )

        direct_f2.append({

            "Source": inverse[a],
            "Target": inverse[b],

            "Effect_label": label,

            "Layers": "; ".join(
                sorted(
                    info["layers"]
                )
            ),

            "Supporting_records":
                info["count"],

            "Max_curation_effort":
                info["effort"],

            "Sources": "; ".join(
                sorted(
                    info["sources"]
                )
            ),

            "References": "; ".join(
                sorted(
                    info["references"]
                )
            ),
        })

    direct_f2.sort(
        key=lambda r: (
            r["Source"],
            r["Target"]
        )
    )

    # --------------------------------------------------------
    # GENE-LEVEL SUMMARY
    # --------------------------------------------------------

    per_gene = []

    for g in canonical:

        rec = {

            "Gene": g,

            "Tier":
                lookup[g]["Tier"],

            "Direction":
                lookup[g]["Direction"],

            "log2FC":
                lookup[g]["log2FC"],

            "padj":
                lookup[g]["padj"],

            "Graph_identifier":
                mapped.get(g) or "",

            "Mapped_in_prior":
                present[g],

            "Biological_modules":
                lookup[g]["Biological_modules"],

            "Original_CORE_community":
                lookup[g]["Original_CORE_community"],

            "F4_structural_bridge":
                lookup[g]["Cross_module_bridge"],

            "F4_CORE_neighbors":
                lookup[g]["CORE_neighbors"],

            "In_testable":
                summaries[g]["in_testable"],

            "Out_testable":
                summaries[g]["out_testable"],
        }

        for s in SETTINGS:

            rec[
                "Incoming_" + s
            ] = summaries[g][
                "in_" + s
            ]

            rec[
                "Outgoing_" + s
            ] = summaries[g][
                "out_" + s
            ]

        per_gene.append(rec)

    # --------------------------------------------------------
    # BIOLOGICAL MODULE-LEVEL SUMMARY
    # --------------------------------------------------------

    modrows = []

    for (
        sm, tm
    ), counts in sorted(
        module_summary.items()
    ):

        rec = {

            "Source_module": sm,
            "Target_module": tm,

            "All_ordered_gene_pairs":
                counts["total_pairs"],

            "Testable_pairs":
                counts["testable_pairs"],
        }

        for s in SETTINGS:

            count = counts[
                "reachable_" + s
            ]

            rec[
                "Reachable_" + s
            ] = count

            rec[
                "Fraction_reachable_" + s
            ] = (
                round(
                    count
                    / counts["testable_pairs"],
                    6
                )
                if counts["testable_pairs"]
                else None
            )

        modrows.append(rec)

    # --------------------------------------------------------
    # EXPORT
    # --------------------------------------------------------

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    write_csv(
        OUT / "F9F_pair_reachability.csv",
        pair_rows,
        list(pair_rows[0])
    )

    write_csv(
        OUT / "F9F_gene_connectivity.csv",
        per_gene,
        list(per_gene[0])
    )

    write_csv(
        OUT / "F9F_module_connectivity.csv",
        modrows,
        list(modrows[0])
    )

    write_csv(
        OUT / "F9F_direct_F2_interactions.csv",
        direct_f2,
        [
            "Source",
            "Target",
            "Effect_label",
            "Layers",
            "Supporting_records",
            "Max_curation_effort",
            "Sources",
            "References"
        ]
    )

    viewer = {

        "description":
            "F9F bounded directed prior-knowledge "
            "reachability; not experimental "
            "pathway activity",

        "settings":
            meta["settings"],

        "gene_connectivity":
            per_gene,

        "module_connectivity":
            modrows,

        "direct_F2_interactions":
            direct_f2,
    }

    write_json(
        OUT / "F9F_viewer_overlay.json",
        viewer
    )

    validation = {

        "phase":
            "F9F_FongWideDirected",

        "source":
            "Original F5 corrected mouse "
            "OmniPath caches "
            "(SHA256 validated)",

        "nodes": len(nodes),

        "directed_edges":
            len(edges),

        "input_directed_rows":
            layer_counts,

        "F2_genes": n,

        "F2_mapped_prior_graph":
            len(mapped),

        "F2_unmapped":
            sorted(
                set(canonical) - set(mapped)
            ),

        "F2_aliases":
            additional_aliases,

        "pairs_total":
            len(pair_rows),

        "pairs_testable":
            sum(
                1
                for r in pair_rows
                if r["Testable"]
            ),

        "support_tiers":
            dict(summary),

        "reachable_pairs": {

            s: sum(
                1
                for r in pair_rows
                if r[
                    "Reachable_" + s
                ]
            )
            for s in SETTINGS
        },

        "direct_F2_edges":
            len(direct_f2),

        "original_F5_pair_setting_checks_passed":
            compared,

        "Fong_direction_pair_tier_counts": [

            {
                "source_direction": a,
                "target_direction": b,
                "tier": t,
                "count": cnt
            }

            for (
                a, b, t
            ), cnt in sorted(
                summary_directions.items()
            )
        ],

        "interpretation_limitations": [

            "Bounded reachability is exact for "
            "F2-to-F2 pairs in this prior "
            "network only.",

            "Directed paths are not evidence "
            "that the relationships were "
            "active in Fong cells.",

            "Biological modules overlap; "
            "the reported fractions are "
            "descriptive, not p-values.",

            "No statistically significant "
            "directed drivers or causal "
            "effects are inferred.",

            "Path counts are not computed. "
            "Representative paths require "
            "a separate evidence-ranked audit.",

            "No new databases or downloads "
            "were used.",
        ],
    }

    write_json(
        OUT / "F9F_validation.json",
        validation
    )

    print()
    print(
        "F9F PASS — exact all-F2-pair "
        "directed prior reachability"
    )

    print(
        "Pair count:",
        len(pair_rows),
        "testable:",
        validation["pairs_testable"]
    )

    print(
        "Support tiers:",
        dict(summary)
    )

    print(
        "Reachability:",
        validation["reachable_pairs"]
    )

    print(
        "Original F5 pair-setting comparisons:",
        compared
    )

    print(
        "Direct F2→F2 interactions:",
        len(direct_f2)
    )

    print(
        "Outputs:",
        OUT
    )


if __name__ == "__main__":
    main()
