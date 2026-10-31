
#!/usr/bin/env python3
"""
F9G — Directed functional-module specificity.

Tests whether selected biological annotations concentrate
around F9F directed reachability beyond coarse connectivity
and CORE/INTERFACE strata.

This is NOT a causal or pathway-activation test.

Inputs: local F9F outputs. No downloads.
Python standard library only.
"""

import argparse
import csv
import json
import math
import random

from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

F9F = (
    ROOT
    / "Data/Processed/FongDownstream/F9F_FongWideDirected"
)

OUT = (
    ROOT
    / "Data/Processed/FongDownstream/F9G_FunctionalSpecificity"
)

SEED = 20261101
N_PERM = 2000

SETTINGS = ("E4_X0", "E4_X1")

ALL_SETTINGS = (
    "E4_X0",
    "E4_X1",
    "E4_X2",
    "E5_X2",
)

MIN_MAPPED = 5
MIN_PAIRS = 100


# ============================================================
# UTILITIES
# ============================================================

def read_csv(path):

    with path.open(
        newline="",
        encoding="utf-8-sig"
    ) as f:

        return list(csv.DictReader(f))


def write_csv(path, records, headers):

    with path.open(
        "w",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=headers,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(records)


def yes(value):

    return str(value).strip().lower() in (
        "true", "yes", "1"
    )


def bits_to_indices(bits):

    while bits:

        bit = bits & -bits
        yield bit.bit_length() - 1

        bits ^= bit


def membership_masks(profiles, donor):
    """
    Shuffle complete joint annotation profiles,
    preserving GO/Reactome annotation overlaps.
    """

    masks = defaultdict(int)

    for i, j in enumerate(donor):

        for module in profiles[j]:

            masks[module] |= 1 << i

    return masks


def counts_for_modules(rowbits, masks, mods):
    """
    Count reachable directed gene pairs for each
    ordered biological-module pair.

    Exclude self-pairs.
    """

    indices = {

        m: list(
            bits_to_indices(
                masks.get(m, 0)
            )
        )

        for m in mods
    }

    results = {}

    for source_module in mods:

        sources = indices[source_module]

        for target_module in mods:

            target_mask = masks.get(
                target_module, 0
            )

            count = sum(

                (
                    rowbits[i]
                    & target_mask
                ).bit_count()

                for i in sources
            )

            possible = (

                len(sources)
                * len(indices[target_module])

                - (
                    masks.get(source_module, 0)
                    & target_mask
                ).bit_count()
            )

            results[
                (source_module, target_module)
            ] = (
                count,
                possible
            )

    return results


def bh(pvals):
    """
    Benjamini-Hochberg FDR adjustment.
    None means a comparison was not tested.
    """

    result = [None] * len(pvals)

    order = sorted(

        (
            i
            for i, p in enumerate(pvals)
            if p is not None
        ),

        key=lambda i: pvals[i]
    )

    m = len(order)
    best = 1.0

    for rank in range(m, 0, -1):

        i = order[rank - 1]

        best = min(
            best,
            pvals[i] * m / rank
        )

        result[i] = best

    return result


def strata_key(row):
    """
    Predeclared broad connectivity classes.

    Preserve CORE/INTERFACE and approximate
    source/target reachability propensity.
    """

    outgoing = int(
        row["Outgoing_E4_X0"]
    )

    incoming = int(
        row["Incoming_E4_X0"]
    )

    ob = (
        "O0"
        if outgoing == 0

        else "O1_200"
        if outgoing <= 200

        else "O201_PLUS"
    )

    ib = (
        "I0"
        if incoming == 0

        else "I1_85"
        if incoming <= 85

        else "I86_PLUS"
    )

    return (
        row["Tier"],
        ob,
        ib
    )


# ============================================================
# INDEPENDENT SELF-TEST
# ============================================================

def self_test():

    profiles = [
        ("A",),
        ("A", "B"),
        ("B",),
        ("C",),
        ("A", "C"),
    ]

    edges = [
        (0, 1),
        (0, 2),
        (1, 2),
        (1, 3),
        (3, 0),
        (4, 1),
    ]

    rows = [0] * len(profiles)

    for a, b in edges:
        rows[a] |= 1 << b

    mods = ["A", "B", "C"]

    donor = list(range(5))

    got = counts_for_modules(
        rows,
        membership_masks(profiles, donor),
        mods
    )

    for sa in mods:
        for tb in mods:

            expected = sum(

                sa in profiles[i]
                and tb in profiles[j]

                for i, j in edges
            )

            denominator = sum(

                sa in profiles[i]
                and tb in profiles[j]
                and i != j

                for i in range(5)
                for j in range(5)
            )

            assert got[(sa, tb)] == (
                expected,
                denominator
            )

    assert bh(
        [0.01, 0.04, None, 0.03]
    ) == [
        0.03, 0.04, None, 0.04
    ]

    rng = random.Random(SEED)

    for _ in range(100):

        shuffled = donor[:]
        rng.shuffle(shuffled)

        masks = membership_masks(
            profiles,
            shuffled
        )

        for mod in mods:

            assert masks[mod].bit_count() == sum(
                mod in profile
                for profile in profiles
            )

    print(
        "SELF TEST PASSED: module counts, "
        "self-pair exclusion, and profile shuffling."
    )


# ============================================================
# MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--self-test",
        action="store_true"
    )

    parser.add_argument(
        "--permutations",
        type=int,
        default=N_PERM
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=SEED
    )

    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    if args.permutations < 99:

        raise ValueError(
            "Use at least 99 permutations; "
            "prefer 2000 or more."
        )

    genes_all_rows = read_csv(
        F9F / "F9F_gene_connectivity.csv"
    )

    old_module_rows = read_csv(
        F9F / "F9F_module_connectivity.csv"
    )

    pairs_path = (
        F9F / "F9F_pair_reachability.csv"
    )

    genes_all = {
        r["Gene"]: r
        for r in genes_all_rows
    }

    if (
        len(genes_all) != 321
        or len(genes_all_rows) != 321
    ):

        raise ValueError(
            "Expected 321 unique F2 genes."
        )

    mapped = [
        r for r in genes_all_rows
        if yes(r["Mapped_in_prior"])
    ]

    if len(mapped) != 257:

        raise ValueError(
            f"Expected 257 mapped genes; "
            f"found {len(mapped)}."
        )

    names = [
        r["Gene"]
        for r in mapped
    ]

    index = {
        g: i
        for i, g in enumerate(names)
    }

    n = len(names)

    profiles = [

        tuple(

            s.strip()

            for s in r[
                "Biological_modules"
            ].split(";")

            if (
                s.strip()
                and s.strip() != "UNRESOLVED"
            )
        )

        for r in mapped
    ]

    mods = sorted({

        module

        for profile in profiles
        for module in profile
    })

    assert len(mods) == 19

    old = {

        (
            r["Source_module"],
            r["Target_module"]
        ): r

        for r in old_module_rows
    }

    assert len(old) == 361

    # ========================================================
    # LOAD FULL F9F PAIRWISE REACHABILITY
    # ========================================================

    matrix = {

        setting: [0] * n

        for setting in ALL_SETTINGS
    }

    read_count = 0
    valid_count = 0

    with pairs_path.open(
        newline="",
        encoding="utf-8-sig"
    ) as stream:

        for pair in csv.DictReader(stream):

            read_count += 1

            a = pair["Source"]
            b = pair["Target"]

            if (
                a == b
                or a not in genes_all
                or b not in genes_all
            ):

                raise ValueError(
                    f"Unexpected F9F pair: {a}->{b}"
                )

            pair_mapped = (
                a in index
                and b in index
            )

            if (
                yes(pair["Testable"])
                != pair_mapped
            ):

                raise ValueError(
                    f"Testable flag mismatch: {a}->{b}"
                )

            if pair_mapped:

                valid_count += 1

                i = index[a]
                j = index[b]

                for setting in ALL_SETTINGS:

                    if yes(
                        pair["Reachable_" + setting]
                    ):

                        matrix[setting][i] |= 1 << j

            elif any(

                yes(pair["Reachable_" + s])

                for s in ALL_SETTINGS
            ):

                raise ValueError(
                    f"Unmapped pair marked reachable: "
                    f"{a}->{b}"
                )

    if (
        read_count != 321 * 320
        or valid_count != 257 * 256
    ):

        raise ValueError(
            "F9F pair-file dimensions differ "
            f"from expected values: "
            f"{read_count}, {valid_count}"
        )

    # ========================================================
    # REPRODUCE ORIGINAL F9F COUNTS
    # ========================================================

    for i, row in enumerate(mapped):

        for setting in ALL_SETTINGS:

            outgoing = (
                matrix[setting][i].bit_count()
            )

            incoming = sum(

                (bits >> i) & 1

                for bits in matrix[setting]
            )

            if (
                outgoing != int(
                    row["Outgoing_" + setting]
                )

                or incoming != int(
                    row["Incoming_" + setting]
                )
            ):

                raise ValueError(
                    f"F9F gene summary mismatch: "
                    f"{row['Gene']} {setting}"
                )

    original_masks = membership_masks(
        profiles,
        list(range(n))
    )

    baseline = {

        setting: counts_for_modules(
            matrix[setting],
            original_masks,
            mods
        )

        for setting in ALL_SETTINGS
    }

    for setting in ALL_SETTINGS:

        for key, (
            reachable,
            possible
        ) in baseline[setting].items():

            previous = old.get(key)

            if (

                previous is None

                or reachable != int(
                    previous[
                        "Reachable_" + setting
                    ]
                )

                or possible != int(
                    previous["Testable_pairs"]
                )
            ):

                raise ValueError(
                    "F9F module mismatch: "
                    f"{setting} {key}"
                )

    print(
        "F9F cross-validation passed: "
        f"{valid_count:,} testable pairs, "
        f"{len(mods)} modules."
    )

    # ========================================================
    # ELIGIBLE MODULE COMPARISONS
    # ========================================================

    sizes = {

        mod: original_masks[mod].bit_count()

        for mod in mods
    }

    eligible = {

        key

        for key, (_, denominator)
        in baseline["E4_X0"].items()

        if (
            sizes[key[0]] >= MIN_MAPPED
            and sizes[key[1]] >= MIN_MAPPED
            and denominator >= MIN_PAIRS
        )
    }

    print(
        f"Comparisons eligible: {len(eligible)}/361"
    )

    # ========================================================
    # DEGREE/TIER STRATA
    # ========================================================

    groups = defaultdict(list)

    for i, row in enumerate(mapped):

        groups[strata_key(row)].append(i)

    movable = sum(

        len(group)

        for group in groups.values()

        if len(group) >= 2
    )

    print(
        f"Permutation strata: {len(groups)}; "
        f"movable assignments: {movable}/{n}"
    )

    if movable < 220:

        raise ValueError(
            "Too many singleton strata."
        )

    # ========================================================
    # PERMUTATION TEST
    # ========================================================

    keys = sorted(eligible)

    observed = {

        (setting, key):
            baseline[setting][key][0]

        for setting in SETTINGS
        for key in keys
    }

    # sum, sum of squares, upper-tail hits,
    # lower-tail hits, number of draws

    stats = {

        (setting, key): [0.0, 0.0, 0, 0, 0]

        for setting in SETTINGS
        for key in keys
    }

    rng = random.Random(args.seed)

    members = list(range(n))

    for rep in range(args.permutations):

        shuffled = members[:]

        for group in groups.values():

            donors = group[:]
            rng.shuffle(donors)

            for i, j in zip(group, donors):

                shuffled[i] = j

        shuffled_masks = membership_masks(
            profiles,
            shuffled
        )

        for setting in SETTINGS:

            counts = counts_for_modules(
                matrix[setting],
                shuffled_masks,
                mods
            )

            for key in keys:

                value, possible = counts[key]

                if (
                    possible
                    != baseline[setting][key][1]
                ):

                    raise ValueError(
                        "Permutation altered the "
                        "module-pair denominator."
                    )

                st = stats[(setting, key)]

                st[0] += value
                st[1] += value * value

                st[2] += int(
                    value >= observed[(setting, key)]
                )

                st[3] += int(
                    value <= observed[(setting, key)]
                )

                st[4] += 1

        if (rep + 1) % 500 == 0:

            print(
                f"  Permutations: "
                f"{rep + 1}/{args.permutations}"
            )

    # ========================================================
    # SUMMARIZE INFERENCE
    # ========================================================

    output = []

    for source_mod in mods:

        for target_mod in mods:

            key = (
                source_mod,
                target_mod
            )

            denominator = (
                baseline["E4_X0"][key][1]
            )

            keep = key in eligible

            for setting in SETTINGS:

                actual = (
                    baseline[setting][key][0]
                )

                row = {

                    "Source_module": source_mod,
                    "Target_module": target_mod,
                    "Setting": setting,

                    "Source_mapped_genes":
                        sizes[source_mod],

                    "Target_mapped_genes":
                        sizes[target_mod],

                    "Testable_pairs": denominator,

                    "Observed_reachable": actual,

                    "Observed_fraction":
                        round(actual / denominator, 7)
                        if denominator else "",

                    "Status":
                        "TESTED"
                        if keep
                        else "LOW_COVERAGE_NOT_TESTED",

                    "Null_mean_count": "",
                    "Null_sd_count": "",
                    "Null_mean_fraction": "",

                    "Observed_minus_null_fraction": "",

                    "Empirical_p_enriched": "",
                    "Empirical_p_depleted": "",
                    "Empirical_p_two_sided": "",
                    "BH_q_two_sided": "",

                    "Direction_of_deviation": "",
                }

                if keep:

                    (
                        total,
                        squares,
                        ge,
                        le,
                        reps
                    ) = stats[(setting, key)]

                    average = total / reps

                    sd = math.sqrt(
                        max(
                            0,
                            squares / reps
                            - average * average
                        )
                    )

                    p_high = (
                        ge + 1
                    ) / (reps + 1)

                    p_low = (
                        le + 1
                    ) / (reps + 1)

                    p_two = min(
                        1.0,
                        2 * min(p_high, p_low)
                    )

                    row.update({

                        "Null_mean_count":
                            round(average, 4),

                        "Null_sd_count":
                            round(sd, 4),

                        "Null_mean_fraction":
                            round(average / denominator, 7),

                        "Observed_minus_null_fraction":
                            round(
                                (actual - average)
                                / denominator,
                                7
                            ),

                        "Empirical_p_enriched":
                            round(p_high, 7),

                        "Empirical_p_depleted":
                            round(p_low, 7),

                        "Empirical_p_two_sided":
                            round(p_two, 7),

                        "Direction_of_deviation":
                            "ENRICHED"
                            if actual > average
                            else "DEPLETED"
                            if actual < average
                            else "EQUAL",
                    })

                output.append(row)

    # ========================================================
    # MULTIPLE TESTING
    # ========================================================

    # Each setting constitutes its own testing family.
    # Apply BH to two-sided empirical p-values.

    for setting in SETTINGS:

        positions = [

            i for i, row in enumerate(output)

            if row["Setting"] == setting
        ]

        pvalues = [

            float(
                output[i]["Empirical_p_two_sided"]
            )

            if output[i]["Status"] == "TESTED"
            else None

            for i in positions
        ]

        qvalues = bh(pvalues)

        for i, q in zip(positions, qvalues):

            if q is not None:

                output[i][
                    "BH_q_two_sided"
                ] = round(q, 7)

    # ========================================================
    # EXPORT
    # ========================================================

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    write_csv(
        OUT / "F9G_module_specificity.csv",
        output,
        list(output[0])
    )

    priorities = sorted(

        [
            r for r in output

            if (
                r["Setting"] == "E4_X0"
                and r["Status"] == "TESTED"
            )
        ],

        key=lambda r: (

            float(r["BH_q_two_sided"]),

            -abs(float(
                r["Observed_minus_null_fraction"]
            ))
        )
    )

    write_csv(
        OUT / "F9G_strict_setting_review.csv",
        priorities,
        list(output[0])
    )

    coverage = []

    for mod in mods:

        related = [

            r for r in mapped

            if mod in r[
                "Biological_modules"
            ].split("; ")
        ]

        coverage.append({

            "Module": mod,

            "Mapped_F2_genes": sizes[mod],

            "Original_F2_genes": sum(

                mod in r[
                    "Biological_modules"
                ].split("; ")

                for r in genes_all_rows
            ),

            "Outgoing_E4_X0_per_gene_mean":
                round(
                    sum(
                        int(r["Outgoing_E4_X0"])
                        for r in related
                    ) / sizes[mod],
                    3
                ),

            "Incoming_E4_X0_per_gene_mean":
                round(
                    sum(
                        int(r["Incoming_E4_X0"])
                        for r in related
                    ) / sizes[mod],
                    3
                ),
        })

    write_csv(
        OUT / "F9G_module_coverage.csv",
        coverage,
        list(coverage[0])
    )

    report = {

        "analysis":
            "F9G_conditioned_module_annotation_permutation",

        "input":
            "F9F exact F2-pair reachability",

        "settings_tested":
            list(SETTINGS),

        "permutations":
            args.permutations,

        "seed":
            args.seed,

        "mapped_F2_genes":
            n,

        "module_count":
            len(mods),

        "comparisons_tested_per_setting":
            len(eligible),

        "stratum_definition":
            "Tier x outgoing_E4_X0 "
            "(0, 1-200, 201+) x incoming_E4_X0 "
            "(0, 1-85, 86+)",

        "stratum_sizes": {

            "/".join(k): len(v)

            for k, v in sorted(
                groups.items()
            )
        },

        "movable_annotation_profiles":
            movable,

        "preserved": (
            "Full directed reachability network, "
            "joint annotation profiles, "
            "CORE/INTERFACE tier, and broad "
            "incoming/outgoing connectivity classes."
        ),

        "not_preserved": (
            "Exact degree within bins, expression "
            "direction, and pathway-edge signs."
        ),

        "warnings": [

            "This tests annotation association, "
            "not biological causation.",

            "Module annotations overlap.",

            "Results depend on the chosen "
            "connectivity strata.",

            "Permutations do not preserve exact "
            "degree within each stratum.",

            "Missing OmniPath genes are not "
            "biologically disconnected.",

            "No pathway activation is inferred.",
        ],
    }

    (
        OUT / "F9G_validation.json"
    ).write_text(
        json.dumps(
            report,
            indent=2
        ),
        encoding="utf-8"
    )

    significant = [

        r for r in priorities

        if float(
            r["BH_q_two_sided"]
        ) < 0.05
    ]

    print()
    print(
        "F9G PASS — conditional annotation "
        "permutation analysis"
    )

    print(
        "E4_X0 comparisons with BH q < 0.05:",
        len(significant)
    )

    for row in significant[:12]:

        print(
            row["Source_module"],
            "->",
            row["Target_module"],
            row["Direction_of_deviation"],
            "delta=",
            row["Observed_minus_null_fraction"],
            "q=",
            row["BH_q_two_sided"]
        )

    print(
        "Outputs:",
        OUT
    )


if __name__ == "__main__":
    main()
