
#!/usr/bin/env python3
"""
F8B — Biological role hierarchy from audited F8A evidence.

No new pathway claims. Preserves source membership and overlaps.
Standard library only. Reads F8A; never modifies earlier results.
"""

import csv
from collections import Counter
from itertools import combinations
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]

SRC = (
    PROJECT
    / "Data/Processed/FongDownstream/F8_FunctionalAtlas"
    / "F8A_FunctionalRoles/F8A_gene_role_atlas.csv"
)

OUT = (
    PROJECT
    / "Data/Processed/FongDownstream/F8_FunctionalAtlas"
    / "F8B_BiologicalModules"
)


# ============================================================
# 1. BIOLOGICAL MODULE DEFINITIONS
# ============================================================

# Each entry:
# (ID, label, family, evidence semantics, GO IDs,
#  named F8A Reactome groups)
#
# Gene Ontology annotations indicate membership or involvement.
# They do not demonstrate pathway activation in Rb-TKO cells.
#
# Reactome annotations come from the F1-significant pathways,
# not from a complete Reactome annotation universe.

MODULES = [

    (
        "INTRINSIC",
        "Intrinsic apoptotic signaling",
        "DEATH",
        "pathway_annotation",
        ["GO:0097193"],
        []
    ),

    (
        "EXTRINSIC",
        "Extrinsic / death-receptor apoptotic signaling",
        "DEATH",
        "pathway_annotation",
        ["GO:0097191"],
        []
    ),

    (
        "APOP_EXEC",
        "Apoptotic execution-phase annotation",
        "DEATH",
        "phase_annotation_review_required",
        ["GO:0097194"],
        []
    ),

    (
        "NECROPTOSIS",
        "Necroptosis / necroptotic signaling",
        "DEATH",
        "pathway_annotation",
        ["GO:0070266", "GO:0097527"],
        []
    ),

    (
        "PROG_NECROSIS",
        "Programmed necrotic death (broader category)",
        "DEATH",
        "broader_death_annotation",
        ["GO:0097300"],
        []
    ),

    (
        "AUTOPHAGY",
        "Autophagic process",
        "AUTOPHAGY",
        "process_annotation",
        ["GO:0061919"],
        []
    ),

    (
        "AUTOPHAGY_REG",
        "Regulation of autophagy",
        "AUTOPHAGY",
        "regulator_annotation",
        ["GO:0010506"],
        []
    ),

    (
        "MITOPHAGY",
        "Mitophagy",
        "AUTOPHAGY",
        "process_annotation",
        ["GO:0000422"],
        []
    ),

    (
        "PRO_DEATH_REG",
        "Positive regulation of apoptosis / cell death",
        "REGULATION",
        "regulator_annotation",
        ["GO:0043065", "GO:0043068"],
        []
    ),

    (
        "ANTI_DEATH_REG",
        "Negative regulation of apoptosis / cell death",
        "REGULATION",
        "regulator_annotation",
        ["GO:0043066", "GO:0043069", "GO:0060548"],
        []
    ),

    (
        "ATR_STRESS",
        "ATR-associated replication stress",
        "GENOME",
        "Reactome_pathway_membership",
        [],
        ["ATR / replication-stress pathway"]
    ),

    (
        "REPLICATION",
        "DNA replication and replication initiation",
        "GENOME",
        "Reactome_pathway_membership",
        [],
        ["Replication initiation / DNA synthesis"]
    ),

    (
        "CHECKPOINT",
        "DNA-damage / cell-cycle checkpoint",
        "GENOME",
        "pathway_annotation",
        ["GO:0000077"],
        ["DNA-damage / cell-cycle checkpoint pathway"]
    ),

    (
        "DNA_REPAIR",
        "DNA repair and lesion bypass",
        "GENOME",
        "pathway_annotation",
        ["GO:0006281"],
        ["DNA repair / lesion bypass pathway"]
    ),

    (
        "HR_DSB",
        "Homologous recombination / DSB repair",
        "GENOME",
        "Reactome_pathway_membership",
        [],
        ["Homologous recombination / DSB repair"]
    ),

    (
        "G1S",
        "G1/S and cell-cycle regulation",
        "CELL_CYCLE",
        "Reactome_pathway_membership",
        [],
        ["G1/S / cell-cycle regulation pathway"]
    ),

    (
        "MITOSIS",
        "Mitosis and chromosome segregation",
        "CELL_CYCLE",
        "Reactome_pathway_membership",
        [],
        ["Mitosis / chromosome segregation pathway"]
    ),

    (
        "ER_STRESS",
        "ER-stress response",
        "STRESS",
        "response_annotation",
        ["GO:0034976"],
        []
    ),

    (
        "OXIDATIVE",
        "Oxidative-stress response",
        "STRESS",
        "response_annotation",
        ["GO:0006979"],
        []
    ),
]


# ============================================================
# 2. REVIEW ANCHORS FROM THE OCTOBER REPORT
# ============================================================

# These are biologically motivated review suggestions.
# They are NOT automatically treated as verified GO memberships.

REVIEW_ANCHORS = {

    "PMAIP1":
        "NOXA / BH3-only mitochondrial apoptosis candidate",

    "BBC3":
        "PUMA / BH3-only mitochondrial apoptosis candidate",

    "FAS":
        "Death receptor; apoptosis versus alternative signaling needs context",

    "TNFRSF10B":
        "TRAIL-related death receptor candidate",

    "ZBP1":
        "Innate sensing and regulated-death candidate",

    "PYCARD":
        "ASC inflammasome adaptor; not a necroptosis executor",

    "MLKL":
        "Necroptosis effector; requires activation-level confirmation",

    "BDNF":
        "Neurotrophic signaling candidate; distinct from generic anti-death GO",

    "KDR":
        "VEGF receptor / trophic signaling candidate",

    "CLSPN":
        "Replication-checkpoint mediator",

    "TOPBP1":
        "DNA replication and ATR checkpoint scaffold",

    "CHEK1":
        "Checkpoint kinase; function may be protective or damaging",

    "H2AFX":
        "DNA damage-associated histone; gamma-H2AX is the phosphorylated marker",

    "E2F1":
        "E2F-dependent cycling / potential fate branch point",

    "MYC":
        "Pleiotropic transcription factor, fate-interface hypothesis",

    "DDIT3":
        "CHOP stress-response transcription factor",
}


# ============================================================
# 3. UTILITIES
# ============================================================

def split_semicolon(value):
    return {
        s.strip()
        for s in (value or "").split(";")
        if s.strip()
    }


def flag(value):
    return str(value).strip().lower() in {
        "true", "1", "yes"
    }


def write_csv(filename, rows, fields):

    with (
        OUT / filename
    ).open("w", newline="", encoding="utf-8-sig") as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


# ============================================================
# 4. MAIN ANALYSIS
# ============================================================

def main():

    if not SRC.exists():
        raise FileNotFoundError(
            f"Missing F8A input: {SRC}"
        )

    with SRC.open(
        newline="",
        encoding="utf-8-sig"
    ) as handle:

        genes = list(csv.DictReader(handle))


    # --------------------------------------------------------
    # 4A. Preserve the existing experimental universe
    # --------------------------------------------------------

    if (
        len(genes) != 321
        or len({r["Gene"] for r in genes}) != 321
    ):
        raise ValueError(
            "F8A must contain 321 unique genes"
        )

    if Counter(
        r["Tier"] for r in genes
    ) != {
        "CORE": 194,
        "INTERFACE": 127
    }:
        raise ValueError(
            "F8A CORE/INTERFACE membership changed"
        )

    if sum(
        flag(r["Shared37"])
        for r in genes
    ) != 37:
        raise ValueError(
            "Expected exactly 37 shared genes"
        )

    if sum(
        flag(r["Cross_module_bridge"])
        for r in genes
    ) != 25:
        raise ValueError(
            "Expected 25 cross-module INTERFACE bridges"
        )


    # --------------------------------------------------------
    # 4B. Construct biological module catalogue
    # --------------------------------------------------------

    catalogue = []

    for mid, label, family, kind, go, reactome in MODULES:

        catalogue.append({

            "Module_id": mid,
            "Module_label": label,
            "Family": family,
            "Evidence_semantics": kind,
            "GO_IDs": "; ".join(go),
            "F8A_Reactome_groups": "; ".join(reactome),

        })


    # --------------------------------------------------------
    # 4C. Assign evidence-supported module memberships
    # --------------------------------------------------------

    membership = []
    gene_atlas = []

    module_sets = {
        entry[0]: set()
        for entry in MODULES
    }

    priority_map = {
        "HIGH": 0,
        "MEDIUM": 1,
        "ROUTINE": 2
    }

    for row in genes:

        gene = row["Gene"]

        go_set = split_semicolon(
            row["GO_ids"]
        )

        reactome_set = split_semicolon(
            row["Reactome_pathway_modules"]
        )

        active_ids = []
        families = set()

        for (
            mid, label, family, kind, go, reactome
        ) in MODULES:

            matched_go = sorted(
                go_set.intersection(go)
            )

            matched_reactome = sorted(
                reactome_set.intersection(reactome)
            )

            if not (
                matched_go or matched_reactome
            ):
                continue

            active_ids.append(mid)
            families.add(family)
            module_sets[mid].add(gene)

            membership.append({

                "Gene": gene,
                "Tier": row["Tier"],

                "Module_id": mid,
                "Module_label": label,
                "Family": family,

                "Evidence_semantics": kind,

                "GO_evidence": "; ".join(matched_go),

                "F1_Reactome_evidence":
                    "; ".join(matched_reactome),

                "Shared37": flag(row["Shared37"]),

                "Cross_module_bridge":
                    flag(row["Cross_module_bridge"]),

                "Direction": row["Direction"],
                "log2FC": row["log2FC"],

            })


        # ----------------------------------------------------
        # 4D. Flag ambiguous biological interpretations
        # ----------------------------------------------------

        death_mechanisms = {
            "INTRINSIC",
            "EXTRINSIC",
            "NECROPTOSIS"
        }

        conflicting_death = (
            len(
                death_mechanisms.intersection(active_ids)
            ) >= 2
        )

        no_refined = not active_ids

        shared = flag(
            row["Shared37"]
        )

        cross_bridge = flag(
            row["Cross_module_bridge"]
        )

        review_reasons = []

        if no_refined:
            review_reasons.append(
                "No specific mechanism from selected F2/F1 evidence"
            )

        if conflicting_death:
            review_reasons.append(
                "Multiple death-modality annotations; inspect gene-specific role"
            )

        if shared:
            review_reasons.append(
                "37-gene old/new/promoter convergence"
            )

        if cross_bridge:
            review_reasons.append(
                "Cross-module INTERFACE bridge"
            )

        if gene in REVIEW_ANCHORS:
            review_reasons.append(
                "Hypothesis anchor in October report"
            )

        if "APOP_EXEC" in active_ids:
            review_reasons.append(
                "Execution-phase GO does not imply direct executor protein"
            )

        if (
            "ANTI_DEATH_REG" in active_ids
            and "PRO_DEATH_REG" in active_ids
        ):
            review_reasons.append(
                "Both positive and negative cell-death regulation annotations"
            )


        # Review priority is a workflow decision.
        # It is NOT a numerical measure of biological importance.

        if (
            shared
            or cross_bridge
            or conflicting_death
            or no_refined
        ):
            priority = "HIGH"

        elif (
            gene in REVIEW_ANCHORS
            or len(families) >= 3
        ):
            priority = "MEDIUM"

        else:
            priority = "ROUTINE"


        # ----------------------------------------------------
        # 4E. Build the consolidated gene atlas
        # ----------------------------------------------------

        gene_atlas.append({

            "Gene": gene,
            "Ensembl_ID": row["Ensembl_ID"],

            "Tier": row["Tier"],
            "Direction": row["Direction"],

            "log2FC": row["log2FC"],
            "padj": row["padj"],

            "Biological_modules":
                "; ".join(active_ids),

            "Biological_families":
                "; ".join(sorted(families)),

            "Number_of_modules": len(active_ids),

            # We deliberately do not force a single biological
            # mechanism before gene-level curation.
            "Unresolved_primary_mechanism": "YES",

            "Report_anchor_hypothesis_UNVALIDATED":
                REVIEW_ANCHORS.get(gene, ""),

            "Shared37": shared,

            "Promoter_evidence":
                flag(row["Promoter_evidence"]),

            "CORE_neighbors":
                row["CORE_neighbors"],

            "CORE_modules_touched":
                row["CORE_modules_touched"],

            "Cross_module_bridge": cross_bridge,

            "Original_CORE_community":
                row["Original_CORE_community"],

            "Original_EXT_community":
                row["Original_EXT_community"],

            "Curation_priority": priority,

            "Review_reasons":
                "; ".join(review_reasons),

            "Context_broad_GO_ids":
                "; ".join(sorted(
                    go_set.intersection({

                        "GO:0008219",
                        "GO:0012501",
                        "GO:0006915",
                        "GO:0097190",
                        "GO:0006974"

                    })
                )),

        })


    # ========================================================
    # 5. EXPORT
    # ========================================================

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    write_csv(
        "F8B_module_catalog.csv",
        catalogue,
        list(catalogue[0])
    )

    write_csv(
        "F8B_gene_module_membership.csv",
        membership,
        list(membership[0])
    )

    write_csv(
        "F8B_gene_role_atlas.csv",
        gene_atlas,
        list(gene_atlas[0])
    )


    # ========================================================
    # 6. MODULE-LEVEL SUMMARY
    # ========================================================

    index = {
        row["Gene"]: row
        for row in gene_atlas
    }

    module_summary = []

    for (
        mid, label, family, kind, go, reactome
    ) in MODULES:

        subset = [
            index[g]
            for g in sorted(module_sets[mid])
        ]

        module_summary.append({

            "Module_id": mid,
            "Module_label": label,
            "Family": family,

            "Evidence_semantics": kind,

            "Genes": len(subset),

            "CORE": sum(
                r["Tier"] == "CORE"
                for r in subset
            ),

            "INTERFACE": sum(
                r["Tier"] == "INTERFACE"
                for r in subset
            ),

            "Up": sum(
                r["Direction"] == "Up"
                for r in subset
            ),

            "Down": sum(
                r["Direction"] == "Down"
                for r in subset
            ),

            "Shared37": sum(
                r["Shared37"]
                for r in subset
            ),

            "INTERFACE_cross_module_bridges": sum(
                r["Cross_module_bridge"]
                for r in subset
            ),

        })

    write_csv(
        "F8B_module_summary.csv",
        module_summary,
        list(module_summary[0])
    )


    # ========================================================
    # 7. BIOLOGICAL ROLE OVERLAP
    # ========================================================

    # These overlaps are shared-gene annotations.
    #
    # They are NOT STRING gene-gene edges, directed signaling
    # interactions, or evidence of mechanistic transmission.

    overlaps = []

    module_pairs = [
        (m[0], m[1])
        for m in MODULES
    ]

    for (a, x), (b, y) in combinations(
        module_pairs, 2
    ):

        shared_genes = (
            module_sets[a] & module_sets[b]
        )

        if not shared_genes:
            continue

        union = (
            module_sets[a] | module_sets[b]
        )

        overlaps.append({

            "Module_A": a,
            "Module_B": b,

            "Shared_gene_count":
                len(shared_genes),

            "Jaccard_overlap":
                round(
                    len(shared_genes) / len(union),
                    5
                ),

            "Shared_genes":
                "; ".join(sorted(shared_genes)),

            "Warning":
                "Co-annotation only; NOT STRING edges or causal links",

        })

    overlaps.sort(
        key=lambda r: (
            -r["Shared_gene_count"],
            r["Module_A"],
            r["Module_B"]
        )
    )

    overlap_columns = [
        "Module_A",
        "Module_B",
        "Shared_gene_count",
        "Jaccard_overlap",
        "Shared_genes",
        "Warning"
    ]

    write_csv(
        "F8B_role_coannotation_overlaps.csv",
        overlaps,
        overlap_columns
    )


    # ========================================================
    # 8. GENE-LEVEL CURATION QUEUE
    # ========================================================

    queue = sorted(

        gene_atlas,

        key=lambda r: (

            priority_map[r["Curation_priority"]],

            not r["Shared37"],

            not r["Cross_module_bridge"],

            r["Gene"]

        )
    )

    write_csv(
        "F8B_curation_queue.csv",
        queue,
        list(queue[0])
    )


    # ========================================================
    # 9. FINAL QC
    # ========================================================

    no_role = [
        row["Gene"]
        for row in gene_atlas
        if row["Number_of_modules"] == 0
    ]

    print("\nF8B biological role hierarchy completed")

    print(
        f"Genes: {len(gene_atlas)} "
        "(194 CORE, 127 INTERFACE)"
    )

    print(
        "Shared37:",
        sum(r["Shared37"] for r in gene_atlas)
    )

    print(
        "Selected modules:",
        len(MODULES)
    )

    print(
        "Gene-module memberships:",
        len(membership),
        "(overlapping)"
    )

    print(
        "Genes without refined F2/F1 membership:",
        len(no_role)
    )

    print(
        "High-priority review genes:",
        sum(
            r["Curation_priority"] == "HIGH"
            for r in gene_atlas
        )
    )

    print(
        "No-role genes:",
        ", ".join(sorted(no_role))
    )

    print("Outputs:", OUT)

    print(
        "REMINDER: no module-module STRING "
        "connectivity has been inferred in this step."
    )


if __name__ == "__main__":
    main()
