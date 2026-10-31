
#!/usr/bin/env python3
"""
F9H — Overlap-aware sensitivity analysis.

Tests directed reachability between genes exclusively assigned
to opposite modules, excluding their shared members.

Uses the full F9F pair matrix and permutes complete annotation
profiles within the same connectivity strata used in F9G.

Standard library only. No new database queries.
"""

import argparse
import csv
import json
import random
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "Data/Processed/FongDownstream"

F9F = BASE / "F9F_FongWideDirected"
F9G = BASE / "F9G_FunctionalSpecificity"
OUT = BASE / "F9H_OverlapSensitivity"

SETTINGS = ("E4_X0", "E4_X1")


def read(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def export(path, data):
    if not data:
        return

    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def truth(x):
    return str(x).strip().lower() in ("true", "1", "yes")


def indices(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def masks_from(profiles, donors):
    result = defaultdict(int)

    for i, donor in enumerate(donors):
        for module in profiles[donor]:
            result[module] |= 1 << i

    return result


def exclusive(masks, a, b):
    if a == b:
        return 0, 0, 0

    common = masks[a] & masks[b]

    return (
        masks[a] & ~common,
        masks[b] & ~common,
        common.bit_count(),
    )


def count_paths(matrix, source_mask, target_mask):
    return sum(
        (matrix[i] & target_mask).bit_count()
        for i in indices(source_mask)
    )


def strata(row):
    incoming = int(row["Incoming_E4_X0"])
    outgoing = int(row["Outgoing_E4_X0"])

    ob = (
        "O0" if outgoing == 0
        else "O1_200" if outgoing <= 200
        else "O201_PLUS"
    )

    ib = (
        "I0" if incoming == 0
        else "I1_85" if incoming <= 85
        else "I86_PLUS"
    )

    return row["Tier"], ob, ib


def bh(pvalues):
    q = [None] * len(pvalues)

    order = sorted(
        (i for i, p in enumerate(pvalues) if p is not None),
        key=lambda i: pvalues[i],
    )

    best = 1.0

    for rank in range(len(order), 0, -1):
        i = order[rank - 1]
        best = min(best, pvalues[i] * len(order) / rank)
        q[i] = best

    return q


def self_test():
    profiles = [
        ("A",),
        ("A", "B"),
        ("B",),
        ("A",),
        ("B",),
    ]

    masks = masks_from(profiles, range(5))
    a, b, shared = exclusive(masks, "A", "B")

    assert (
        a.bit_count(),
        b.bit_count(),
        shared
    ) == (2, 2, 1)

    matrix = [
        1 << 2,
        0,
        0,
        (1 << 2) | (1 << 4),
        0,
    ]

    assert count_paths(matrix, a, b) == 3

    print("F9H SELF TEST PASSED.")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--self-test",
        action="store_true",
    )

    parser.add_argument(
        "--permutations",
        type=int,
        default=2000,
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20261107,
    )

    args = parser.parse_args()

    if args.self_test:
        return self_test()

    if args.permutations < 99:
        raise ValueError("At least 99 permutations required.")

    # ========================================================
    # 1. F9F GENE UNIVERSE
    # ========================================================

    all_genes = read(
        F9F / "F9F_gene_connectivity.csv"
    )

    mapped = [
        r for r in all_genes
        if truth(r["Mapped_in_prior"])
    ]

    if len(all_genes) != 321 or len(mapped) != 257:
        raise ValueError("F9F mapping changed.")

    ids = {
        r["Gene"]: i
        for i, r in enumerate(mapped)
    }

    profiles = [
        tuple(
            x.strip()
            for x in r["Biological_modules"].split(";")
            if x.strip() and x.strip() != "UNRESOLVED"
        )
        for r in mapped
    ]

    names = sorted({
        m
        for profile in profiles
        for m in profile
    })

    if len(names) != 19:
        raise ValueError("Expected 19 biological modules.")

    # ========================================================
    # 2. LOAD EXACT F9F REACHABILITY
    # ========================================================

    matrix = {
        s: [0] * 257
        for s in SETTINGS
    }

    nrows = 0
    ntest = 0

    with (
        F9F / "F9F_pair_reachability.csv"
    ).open(newline="", encoding="utf-8-sig") as f:

        for pair in csv.DictReader(f):
            nrows += 1

            i = ids.get(pair["Source"])
            j = ids.get(pair["Target"])

            testable = i is not None and j is not None

            if testable != truth(pair["Testable"]):
                raise ValueError("F9F pair coverage mismatch.")

            if not testable:
                continue

            ntest += 1

            for s in SETTINGS:
                if truth(pair["Reachable_" + s]):
                    matrix[s][i] |= 1 << j

    if (nrows, ntest) != (102720, 65792):
        raise ValueError("Incorrect F9F pair universe.")

    totals = [
        sum(mask.bit_count() for mask in matrix[s])
        for s in SETTINGS
    ]

    if totals != [19190, 35669]:
        raise ValueError("F9F reachability totals changed.")

    # ========================================================
    # 3. LOAD ORIGINAL F9G RESULTS
    # ========================================================

    info = json.loads(
        (F9G / "F9G_validation.json").read_text(
            encoding="utf-8"
        )
    )

    screen = read(
        F9G / "F9G_module_specificity.csv"
    )

    scores = {
        (
            r["Setting"],
            r["Source_module"],
            r["Target_module"]
        ): r
        for r in screen
    }

    if len(scores) != 722:
        raise ValueError("F9G module results incomplete.")

    original_hits = [
        r for r in screen
        if (
            r["Setting"] == "E4_X0"
            and r["Status"] == "TESTED"
            and float(r["BH_q_two_sided"]) < 0.05
        )
    ]

    print(
        f"Original F9G: {info['permutations']} permutations, "
        f"{len(original_hits)} significant E4_X0 pairs."
    )

    if info["permutations"] < 20000:
        print("WARNING: F9G source is an exploratory run.")

    # ========================================================
    # 4. VALIDATE MODULE COUNTS
    # ========================================================

    original = masks_from(profiles, range(257))

    module_summary = {
        (
            r["Source_module"],
            r["Target_module"]
        ): r
        for r in read(
            F9F / "F9F_module_connectivity.csv"
        )
    }

    for a in names:
        for b in names:
            source_mask = original[a]
            target_mask = original[b]

            denominator = (
                source_mask.bit_count()
                * target_mask.bit_count()
                - (source_mask & target_mask).bit_count()
            )

            if denominator != int(
                module_summary[a, b]["Testable_pairs"]
            ):
                raise ValueError("Module denominator mismatch.")

            for s in SETTINGS:
                actual = count_paths(
                    matrix[s],
                    source_mask,
                    target_mask,
                )

                expected = int(
                    module_summary[a, b]["Reachable_" + s]
                )

                if actual != expected:
                    raise ValueError(
                        f"Module reachability mismatch: {s}, {a}, {b}"
                    )

    print("F9F graph/module validation passed.")

    # ========================================================
    # 5. EXCLUSIVE MODULE PAIRS
    # ========================================================

    candidates = {}
    eligible = []

    for a in names:
        for b in names:
            sm, tm, overlap = exclusive(original, a, b)

            ns = sm.bit_count()
            nt = tm.bit_count()

            valid = (
                a != b
                and min(ns, nt) >= 5
                and ns * nt >= 100
            )

            candidates[a, b] = {
                "overlap": overlap,
                "source_only": ns,
                "target_only": nt,
                "possible": ns * nt,
                "valid": valid,
            }

            if valid:
                eligible.append((a, b))

    print(
        "Eligible exclusive comparisons per setting:",
        len(eligible)
    )

    # ========================================================
    # 6. F9G-COMPATIBLE STRATIFICATION
    # ========================================================

    groups = defaultdict(list)

    for i, row in enumerate(mapped):
        groups[strata(row)].append(i)

    movable = sum(
        len(v) for v in groups.values()
        if len(v) > 1
    )

    if movable < 220:
        raise ValueError("Too many singleton strata.")

    # ========================================================
    # 7. OBSERVED EXCLUSIVE REACHABILITY
    # ========================================================

    observed = {
        (s, a, b): count_paths(
            matrix[s],
            *exclusive(original, a, b)[:2]
        )
        for s in SETTINGS
        for a, b in eligible
    }

    # Sum, upper-tail hits, lower-tail hits.
    stats = {
        (s, a, b): [0.0, 0, 0]
        for s in SETTINGS
        for a, b in eligible
    }

    rng = random.Random(args.seed)

    # ========================================================
    # 8. COMPLETE-PROFILE PERMUTATIONS
    # ========================================================

    for rep in range(args.permutations):
        donors = list(range(257))

        for group in groups.values():
            shuffled = group[:]
            rng.shuffle(shuffled)

            for i, j in zip(group, shuffled):
                donors[i] = j

        shuffled_masks = masks_from(profiles, donors)

        for a, b in eligible:
            sm, tm, _ = exclusive(shuffled_masks, a, b)

            if (
                sm.bit_count() * tm.bit_count()
                != candidates[a, b]["possible"]
            ):
                raise ValueError(
                    "Permutation changed exclusive pair denominator."
                )

            for s in SETTINGS:
                key = (s, a, b)

                value = count_paths(
                    matrix[s], sm, tm
                )

                st = stats[key]

                st[0] += value
                st[1] += int(value >= observed[key])
                st[2] += int(value <= observed[key])

        if (rep + 1) % 500 == 0:
            print(
                f"Permutations: {rep + 1}/{args.permutations}"
            )

    # ========================================================
    # 9. EMPIRICAL P-VALUES AND BH CORRECTION
    # ========================================================

    output = []

    for s in SETTINGS:
        for a in names:
            for b in names:
                item = candidates[a, b]
                valid = item["valid"]
                count = observed.get((s, a, b))

                row = {
                    "Setting": s,
                    "Source_module": a,
                    "Target_module": b,

                    "Shared_members": item["overlap"],
                    "Source_only": item["source_only"],
                    "Target_only": item["target_only"],
                    "Exclusive_pairs": item["possible"],

                    "Exclusive_reachable":
                        count if valid else "",

                    "Exclusive_fraction":
                        round(count / item["possible"], 7)
                        if valid else "",

                    "Status":
                        "TESTED" if valid
                        else "WITHIN_MODULE" if a == b
                        else "LOW_EXCLUSIVE_COVERAGE",

                    "Original_F9G_BH_q":
                        scores[s, a, b]["BH_q_two_sided"],

                    "Null_mean": "",
                    "Excess_fraction": "",
                    "Empirical_two_sided_p": "",
                    "BH_q_two_sided": "",
                }

                if valid:
                    total, ge, le = stats[s, a, b]

                    average = total / args.permutations

                    p_high = (ge + 1) / (
                        args.permutations + 1
                    )

                    p_low = (le + 1) / (
                        args.permutations + 1
                    )

                    p_two = min(
                        1.0,
                        2 * min(p_high, p_low)
                    )

                    row.update({
                        "Null_mean": round(average, 4),

                        "Excess_fraction":
                            round(
                                (count - average)
                                / item["possible"],
                                7
                            ),

                        "Empirical_two_sided_p":
                            round(p_two, 7),
                    })

                output.append(row)

    for s in SETTINGS:
        subset = [
            r for r in output
            if r["Setting"] == s
        ]

        pvalues = [
            float(r["Empirical_two_sided_p"])
            if r["Status"] == "TESTED"
            else None
            for r in subset
        ]

        for row, q in zip(subset, bh(pvalues)):
            if q is not None:
                row["BH_q_two_sided"] = round(q, 7)

    # ========================================================
    # 10. AUDIT ORIGINAL SIGNIFICANT PAIRS
    # ========================================================

    prior_edges = read(
        F9F / "F9F_direct_F2_interactions.csv"
    )

    if len(prior_edges) != 963:
        raise ValueError("Direct interaction count changed.")

    membership = {
        r["Gene"]: set(profiles[i])
        for i, r in enumerate(mapped)
    }

    audit = []
    exclusive_edges = []

    for hit in original_hits:
        a = hit["Source_module"]
        b = hit["Target_module"]

        src = {
            g for g, p in membership.items()
            if a in p
        }

        dst = {
            g for g, p in membership.items()
            if b in p
        }

        shared = src & dst

        matched = [
            e for e in prior_edges
            if (
                e["Source"] in src
                and e["Target"] in dst
            )
        ]

        distinct = [
            e for e in matched
            if (
                e["Source"] not in shared
                and e["Target"] not in shared
            )
        ]

        stat = next(
            r for r in output
            if (
                r["Setting"] == "E4_X0"
                and r["Source_module"] == a
                and r["Target_module"] == b
            )
        )

        audit.append({
            "Source_module": a,
            "Target_module": b,
            "Original_BH_q": hit["BH_q_two_sided"],

            "Shared_genes": len(shared),
            "Direct_interactions": len(matched),
            "Exclusive_direct_interactions": len(distinct),

            "Source_only": stat["Source_only"],
            "Target_only": stat["Target_only"],

            "Exclusive_fraction":
                stat["Exclusive_fraction"],

            "Exclusive_excess_fraction":
                stat["Excess_fraction"],

            "Exclusive_BH_q":
                stat["BH_q_two_sided"],

            "Status": stat["Status"],
            "Shared_gene_names":
                "; ".join(sorted(shared)),
        })

        for edge in distinct:
            exclusive_edges.append({
                "Source_module": a,
                "Target_module": b,
                "Source_gene": edge["Source"],
                "Target_gene": edge["Target"],
                "Effect_label": edge["Effect_label"],
                "Layers": edge["Layers"],
            })

    # ========================================================
    # 11. EXPORT
    # ========================================================

    OUT.mkdir(parents=True, exist_ok=True)

    export(
        OUT / "F9H_exclusive_module_results.csv",
        output,
    )

    export(
        OUT / "F9H_original_hits_audit.csv",
        audit,
    )

    if exclusive_edges:
        export(
            OUT / "F9H_exclusive_direct_edges.csv",
            exclusive_edges,
        )

    summary = {
        "phase": "F9H_OverlapSensitivity",
        "permutations": args.permutations,
        "seed": args.seed,

        "F9G_input_permutations":
            info["permutations"],

        "eligible_exclusive_pairs_per_setting":
            len(eligible),

        "original_E4_X0_hits":
            len(original_hits),

        "movable_profiles": movable,

        "exclusive_E4_X0_BH_under_005": sum(
            r["Status"] == "TESTED"
            and r["BH_q_two_sided"] != ""
            and float(r["BH_q_two_sided"]) < 0.05
            for r in output
            if r["Setting"] == "E4_X0"
        ),

        "limitations": [
            "This is a sensitivity analysis, not independent replication.",
            "The original F9G network and degree strata are retained.",
            "Shared module members are excluded from opposite-module comparisons.",
            "Within-module connectivity is evaluated separately.",
            "Database reachability does not establish biological activation."
        ],
    }

    (
        OUT / "F9H_validation.json"
    ).write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8"
    )

    print()
    print("F9H PASS — overlap-aware analysis")

    print(
        "E4_X0 exclusive comparisons with BH q < 0.05:",
        summary["exclusive_E4_X0_BH_under_005"]
    )

    for r in audit:
        print(
            f"{r['Source_module']} -> {r['Target_module']}: "
            f"shared={r['Shared_genes']}, "
            f"direct={r['Direct_interactions']}, "
            f"exclusive_direct={r['Exclusive_direct_interactions']}, "
            f"exclusive_q={r['Exclusive_BH_q'] or 'not tested'}"
        )

    print("Outputs:", OUT)


if __name__ == "__main__":
    main()
