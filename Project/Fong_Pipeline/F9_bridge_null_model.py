
#!/usr/bin/env python3
"""
F9E: Degree-preserving permutation analysis of
CORE-INTERFACE community connectivity.

Preserves:
  - 194 CORE and 127 INTERFACE gene identities
  - the 405 CORE-INTERFACE STRING edges in number
  - each gene's CORE-INTERFACE degree
  - original CORE community assignments

Randomizes ONLY the pairing of CORE-INTERFACE edges.

Standard library only.
"""

import csv
import json
import random

from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, pstdev


ROOT = Path(__file__).resolve().parents[1]

BASE = (
    ROOT
    / "Data/Processed/FongDownstream/F8_FunctionalAtlas"
)

F8D = (
    BASE
    / "F8D_PriorKnowledge/F8D_viewer_data.json"
)

F8B = (
    BASE
    / "F8B_BiologicalModules/F8B_gene_role_atlas.csv"
)

OUT = BASE / "F9E_BridgeNullModel"

SEED = 20261030
PERMUTATIONS = 2000
SWAPS_PER_DRAW = 4000


def csv_write(path, rows, columns):

    with path.open(
        "w",
        newline="",
        encoding="utf-8-sig"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=columns
        )

        writer.writeheader()
        writer.writerows(rows)


def benjamini_hochberg(values):

    q = [1.0] * len(values)
    m = len(values)

    ordered = sorted(
        range(m),
        key=lambda i: values[i]
    )

    previous = 1.0

    for rank in range(m, 0, -1):

        index = ordered[rank - 1]

        previous = min(
            previous,
            values[index] * m / rank
        )

        q[index] = previous

    return q


def main():

    data = json.loads(
        F8D.read_text(encoding="utf-8")
    )

    with F8B.open(
        newline="",
        encoding="utf-8-sig"
    ) as f:

        rowmap = {
            r["Gene"]: r
            for r in csv.DictReader(f)
        }

    genes = {
        n["id"]: n
        for n in data["nodes"]
    }

    core = {
        g for g in rowmap
        if rowmap[g]["Tier"] == "CORE"
    }

    interface = sorted(
        g for g in rowmap
        if rowmap[g]["Tier"] == "INTERFACE"
    )

    assert len(core) == 194
    assert len(interface) == 127

    assert sum(
        bool(n.get("f4_cross_module_bridge"))
        for n in genes.values()
    ) == 25

    labels = {
        g: rowmap[g]["Original_CORE_community"]
        for g in core
    }

    assert sum(
        value == "STRING isolate"
        for value in labels.values()
    ) == 70

    # --------------------------------------------------------
    # EXTRACT OBSERVED CORE-INTERFACE EDGES
    # --------------------------------------------------------

    bipartite = []

    for edge in data["string_edges"]:

        a = edge["source"]
        b = edge["target"]

        if a in interface and b in core:
            bipartite.append((a, b))

        elif b in interface and a in core:
            bipartite.append((b, a))

    assert len(bipartite) == 405
    assert len(set(bipartite)) == 405

    # --------------------------------------------------------
    # COUNT DISTINCT ORIGINAL CORE COMMUNITIES
    # --------------------------------------------------------

    def group_counts(edge_list):

        groups = defaultdict(set)

        for inter, core_gene in edge_list:

            community = labels[core_gene]

            # Preserve the F4 convention:
            # original CORE isolates do not contribute
            # an original non-isolate community label.
            if community != "STRING isolate":
                groups[inter].add(community)

        return {
            gene: len(groups[gene])
            for gene in interface
        }

    observed = group_counts(bipartite)

    for gene in interface:

        expected = int(
            rowmap[gene]["CORE_modules_touched"] or 0
        )

        if observed[gene] != expected:
            raise ValueError(
                f"F4 community mismatch: {gene}"
            )

    observed_total = sum(
        count >= 2
        for count in observed.values()
    )

    assert observed_total == 25

    degree = Counter(
        gene for gene, _ in bipartite
    )

    # --------------------------------------------------------
    # DEGREE-PRESERVING BIPARTITE EDGE SWAPS
    # --------------------------------------------------------

    rng = random.Random(SEED)

    edges = bipartite.copy()
    current = set(edges)

    def shuffle_edges(accepted_target):

        accepted = 0
        attempts = 0

        while (
            accepted < accepted_target
            and attempts < accepted_target * 30
        ):

            attempts += 1

            a = rng.randrange(len(edges))
            b = rng.randrange(len(edges))

            if a == b:
                continue

            i, c = edges[a]
            j, d = edges[b]

            if i == j or c == d:
                continue

            # Reject duplicate edges.
            if (
                (i, d) in current
                or (j, c) in current
            ):
                continue

            current.remove((i, c))
            current.remove((j, d))

            current.add((i, d))
            current.add((j, c))

            edges[a] = (i, d)
            edges[b] = (j, c)

            accepted += 1

        if accepted < accepted_target:
            raise RuntimeError(
                "Insufficient accepted edge swaps"
            )

    # Initial mixing.
    shuffle_edges(SWAPS_PER_DRAW * 2)

    null_totals = []

    null_per_gene = {
        gene: []
        for gene in interface
    }

    for _ in range(PERMUTATIONS):

        shuffle_edges(SWAPS_PER_DRAW)

        counts = group_counts(edges)

        null_totals.append(
            sum(
                count >= 2
                for count in counts.values()
            )
        )

        for gene, count in counts.items():
            null_per_gene[gene].append(count)

    # --------------------------------------------------------
    # GLOBAL EMPIRICAL P-VALUES
    # --------------------------------------------------------

    global_low_p = (
        1 + sum(
            value <= observed_total
            for value in null_totals
        )
    ) / (PERMUTATIONS + 1)

    global_high_p = (
        1 + sum(
            value >= observed_total
            for value in null_totals
        )
    ) / (PERMUTATIONS + 1)

    # --------------------------------------------------------
    # GENE-LEVEL NULL RESULTS
    # --------------------------------------------------------

    results = []

    for gene in interface:

        samples = null_per_gene[gene]
        observed_count = observed[gene]

        p_fewer = (
            1 + sum(
                value <= observed_count
                for value in samples
            )
        ) / (PERMUTATIONS + 1)

        p_more = (
            1 + sum(
                value >= observed_count
                for value in samples
            )
        ) / (PERMUTATIONS + 1)

        results.append({

            "Gene": gene,

            "CORE_links": degree[gene],

            "Observed_distinct_CORE_communities":
                observed_count,

            "F4_structural_bridge":
                observed_count >= 2,

            "Null_mean_distinct_communities":
                round(mean(samples), 4),

            "Null_sd_distinct_communities":
                round(pstdev(samples), 4),

            "P_fewer_or_equal_communities":
                p_fewer,

            "P_more_or_equal_communities":
                p_more,
        })

    # Correct separate lower/upper-tail gene-level
    # tests for multiple comparisons.

    for side in ("fewer", "more"):

        pfield = (
            "P_fewer_or_equal_communities"
            if side == "fewer"
            else "P_more_or_equal_communities"
        )

        qvals = benjamini_hochberg([
            row[pfield]
            for row in results
        ])

        for row, q in zip(results, qvals):
            row["BH_q_" + side] = round(q, 7)

    results.sort(
        key=lambda r: (
            -r["Observed_distinct_CORE_communities"],
            -r["CORE_links"],
            r["Gene"]
        )
    )

    # --------------------------------------------------------
    # EXPORT
    # --------------------------------------------------------

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    csv_write(
        OUT / "F9E_gene_bridge_null.csv",
        results,
        list(results[0])
    )

    csv_write(
        OUT / "F9E_null_distribution.csv",
        [
            {
                "Permutation": j + 1,
                "Bridge_count": count,
            }
            for j, count in enumerate(null_totals)
        ],
        ["Permutation", "Bridge_count"]
    )

    summary = {

        "analysis":
            "Degree-preserving bipartite CORE-INTERFACE null",

        "seed": SEED,

        "permutations": PERMUTATIONS,

        "accepted_swaps_between_draws":
            SWAPS_PER_DRAW,

        "CORE_genes": 194,
        "INTERFACE_genes": 127,

        "CORE_INTERFACE_edges": 405,

        "original_CORE_isolates_not_counted_as_communities":
            70,

        "observed_cross_community_bridges":
            observed_total,

        "null_bridge_count_mean":
            mean(null_totals),

        "null_bridge_count_sd":
            pstdev(null_totals),

        "null_range": [
            min(null_totals),
            max(null_totals)
        ],

        "one_sided_empirical_p_observed_fewer":
            global_low_p,

        "one_sided_empirical_p_observed_more":
            global_high_p,

        "limitations": [

            "Only CORE-INTERFACE edges are rewired.",

            "Original CORE communities and all bipartite degrees are fixed.",

            "The null does not preserve biological edge propensities.",

            "Results concern network structure, not causal effects.",

            "MCMC draws may be correlated; test robustness across seeds and mixing lengths.",
        ],
    }

    (
        OUT / "F9E_summary.json"
    ).write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8"
    )

    print(
        "F9E PASS — CORE/INTERFACE degree-preserving null"
    )

    print(
        "Observed bridges:",
        observed_total
    )

    print(
        "Null mean:",
        round(summary["null_bridge_count_mean"], 3)
    )

    print(
        "Null range:",
        summary["null_range"]
    )

    print(
        "Empirical one-sided p (fewer):",
        round(global_low_p, 6)
    )

    print("Outputs:", OUT)


if __name__ == "__main__":
    main()
