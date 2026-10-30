
#!/usr/bin/env python3
"""
F8A — Evidence-based functional decomposition.

Inputs:
    F2_CellFate/cell_fate_master.csv
    F2_CellFate/reactome_annotations.csv
    F2_CellFate/GO_term_audit.csv
    F2_CellFate/F2_metadata.json
    F8_FunctionalAtlas/gene_role_draft.csv

Outputs:
    F8_FunctionalAtlas/F8A_FunctionalRoles/

Uses the Python standard library only.
Does not modify previous analysis outputs.
"""

import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


# ============================================================
# 1. PATH CONFIGURATION
# ============================================================

PROJECT = Path(__file__).resolve().parents[1]

F2 = (
    PROJECT
    / "Data/Processed/FongDownstream/F2_CellFate"
)

F8 = (
    PROJECT
    / "Data/Processed/FongDownstream/F8_FunctionalAtlas"
)

OUT = F8 / "F8A_FunctionalRoles"


# ============================================================
# 2. SPECIFIC GO-BASED BIOLOGICAL ROLES
# ============================================================

# Each GO ID corresponds to a term actually present in the
# original F2 analysis.
#
# The Boolean indicates whether the annotation is sufficiently
# specific to contribute to our preliminary role decomposition.
#
# Specific GO membership is NOT proof that the associated
# biological process is functionally active in the Fong model.

GO_RULES = {
    "GO:0097193": (
        "Intrinsic apoptosis signaling", True
    ),
    "GO:0097191": (
        "Extrinsic apoptosis signaling", True
    ),
    "GO:0097194": (
        "Apoptotic execution phase", True
    ),
    "GO:0097190": (
        "Apoptotic signaling, unspecified", False
    ),
    "GO:0006915": (
        "Apoptotic process, broad", False
    ),
    "GO:0008219": (
        "Cell death, broad", False
    ),
    "GO:0012501": (
        "Programmed cell death, broad", False
    ),
    "GO:0043065": (
        "Positive regulation of apoptosis", True
    ),
    "GO:0043068": (
        "Positive regulation of programmed death", True
    ),
    "GO:0043066": (
        "Negative regulation of apoptosis", True
    ),
    "GO:0043069": (
        "Negative regulation of programmed death", True
    ),
    "GO:0060548": (
        "Negative regulation of cell death", True
    ),
    "GO:0070266": (
        "Necroptotic process", True
    ),
    "GO:0097527": (
        "Necroptotic signaling", True
    ),
    "GO:0097300": (
        "Programmed necrotic cell death", True
    ),
    "GO:0000422": (
        "Mitophagy", True
    ),
    "GO:0010506": (
        "Regulation of autophagy", True
    ),
    "GO:0061919": (
        "Autophagic process", True
    ),
    "GO:0000077": (
        "DNA-damage checkpoint", True
    ),
    "GO:0006281": (
        "DNA repair", True
    ),
    "GO:0006974": (
        "DNA-damage response, broad", False
    ),
    "GO:0034976": (
        "ER-stress response", True
    ),
    "GO:0006979": (
        "Oxidative-stress response", True
    ),
}


# ============================================================
# 3. REACTOME PATHWAY CLASSIFICATION
# ============================================================

# These rules classify membership in F1-significant Reactome
# pathways by their reported pathway names.
#
# This is a supplementary annotation layer, not a new
# enrichment analysis.
#
# Pathway FDR is a pathway-level statistic. It must not be
# interpreted as a gene-specific significance value.

REACTOME_RULES = {

    "ATR / replication-stress pathway": (
        r"ATR in response to replication stress"
    ),

    "Replication initiation / DNA synthesis": (
        r"DNA Replication|"
        r"DNA replication|"
        r"S Phase|"
        r"Synthesis of DNA|"
        r"DNA strand elongation|"
        r"pre-replicative|"
        r"Leading Strand|"
        r"Lagging Strand|"
        r"replication initiation|"
        r"processive synthesis"
    ),

    "DNA-damage / cell-cycle checkpoint pathway": (
        r"checkpoint|"
        r"Chk1/Chk2|"
        r"G1 Cell Cycle Arrest"
    ),

    "Homologous recombination / DSB repair": (
        r"Homologous Recombination|"
        r"Homology Directed Repair|"
        r"Homologous DNA|"
        r"HDR through|"
        r"D-Loop|"
        r"Holliday Junction|"
        r"Single Strand Annealing|"
        r"double-strand break ends"
    ),

    "DNA repair / lesion bypass pathway": (
        r"DNA Repair|"
        r"DNA Double-Strand Break Repair|"
        r"Translesion|"
        r"DNA Damage Bypass|"
        r"base excision|"
        r"AP sites|"
        r"Abasic Sites|"
        r"GG-NER|"
        r"Flap Intermediate|"
        r"Polymerase switching"
    ),

    "Mitosis / chromosome segregation pathway": (
        r"^Cell Cycle, Mitotic$|"
        r"Mitotic (?:Prometaphase|Metaphase|Anaphase|Spindle)|"
        r"^M Phase$|"
        r"Metaphase|"
        r"Anaphase|"
        r"chromatid|"
        r"Prometaphase Chromosomes|"
        r"spindle|"
        r"kinetochore|"
        r"centrosome|"
        r"PLK1|"
        r"Polo-like kinase|"
        r"AURKA|"
        r"APC/C|"
        r"G2/M Transition"
    ),

    "G1/S / cell-cycle regulation pathway": (
        r"^Cell Cycle$|"
        r"^G1/S Transition$|"
        r"^G0 and Early G1$|"
        r"^G2 Phase$|"
        r"Cell Cycle Genes|"
        r"Mitotic G1 phase"
    ),
}


# ============================================================
# 4. FILE UTILITIES
# ============================================================

def load_csv(path):
    with path.open(
        newline="",
        encoding="utf-8-sig"
    ) as handle:
        return list(csv.DictReader(handle))


def save_csv(path, rows, columns):
    with path.open(
        "w",
        newline="",
        encoding="utf-8-sig"
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def tokens(value):
    return {
        part.strip()
        for part in (value or "").split(";")
        if part.strip()
    }


def yes(value):
    return str(value).strip().lower() in {
        "true", "yes", "1"
    }


# ============================================================
# 5. LOAD DATA AND VERIFY PROVENANCE
# ============================================================

def main():

    f2 = load_csv(
        F2 / "cell_fate_master.csv"
    )

    reactome = load_csv(
        F2 / "reactome_annotations.csv"
    )

    audit = load_csv(
        F2 / "GO_term_audit.csv"
    )

    metadata = json.loads(
        (F2 / "F2_metadata.json").read_text(
            encoding="utf-8"
        )
    )

    previous = load_csv(
        F8 / "gene_role_draft.csv"
    )

    assert len(f2) == metadata["selected_gene_count"] == 321

    assert Counter(
        row["F2_tier"] for row in f2
    ) == {
        "CORE": 194,
        "INTERFACE": 127
    }

    assert len({
        row["Gene"] for row in f2
    }) == 321

    assert {
        row["go_id"]
        for row in metadata["go_panel"]
    } == set(GO_RULES)

    previous_by_gene = {
        row["Gene"].upper(): row
        for row in previous
    }

    assert len(previous_by_gene) == 321

    go_label_by_id = {
        row["go_id"]: row["label"]
        for row in metadata["go_panel"]
    }

    # Independently verify each GO term's reported gene count.
    for item in audit:

        go_id = item["GO_ID"]

        actual = sum(
            go_id in tokens(row["F2_GO_terms"])
            for row in f2
        )

        assert actual == int(item["Strict_Fong_hits"]), (
            go_id, actual
        )

    print("F2 source and GO-count validation passed.")


    # ========================================================
    # 6. INDEX REACTOME PATHWAYS BY ENSEMBL ID
    # ========================================================

    reac_by_ensembl = defaultdict(list)

    for pathway in reactome:

        for ens in tokens(pathway["F2_Ensembl_genes"]):

            reac_by_ensembl[ens].append(pathway)


    # ========================================================
    # 7. BUILD THE FUNCTIONAL ATLAS
    # ========================================================

    atlas = []
    evidence = []
    matrix = []

    go_programs = {
        name
        for name, specific in GO_RULES.values()
        if specific
    }

    reactome_programs = set(REACTOME_RULES)

    program_columns = sorted(
        go_programs | reactome_programs
    )

    for gene_row in f2:

        gene = gene_row["Gene"].upper()
        old = previous_by_gene[gene]

        go_ids = tokens(
            gene_row["F2_GO_terms"]
        )

        unexpected = go_ids - set(GO_RULES)

        assert not unexpected, (
            gene, unexpected
        )

        direct = set()
        broad = set()
        reac = set()
        matched_pathways = set()

        matching_reactome = reac_by_ensembl.get(
            gene_row["Ensembl_ID"], []
        )


        # ----------------------------------------------------
        # 7A. GO-derived evidence
        # ----------------------------------------------------

        for go_id in sorted(go_ids):

            program, specific = GO_RULES[go_id]

            if specific:
                direct.add(program)
            else:
                broad.add(program)

            evidence.append({
                "Gene": gene,
                "Tier": gene_row["F2_tier"],
                "Program": program,
                "Source": "GO_selected_F2",
                "Specific_evidence": specific,
                "Evidence_ID": go_id,
                "Evidence_label": go_label_by_id[go_id],
                "Reactome_analysis": "",
                "Reactome_pathway_FDR": "",
            })


        # ----------------------------------------------------
        # 7B. Reactome-derived evidence
        # ----------------------------------------------------

        for pathway in matching_reactome:

            name = pathway["Reactome_name"]

            matched_pathways.add((
                pathway["Analysis"],
                pathway["Reactome_ID"],
                name
            ))

            for program, pattern in REACTOME_RULES.items():

                if not re.search(
                    pattern,
                    name,
                    flags=re.IGNORECASE
                ):
                    continue

                reac.add(program)

                evidence.append({
                    "Gene": gene,
                    "Tier": gene_row["F2_tier"],
                    "Program": program,
                    "Source": "F1_enriched_Reactome_intersection",
                    "Specific_evidence": False,
                    "Evidence_ID": pathway["Reactome_ID"],
                    "Evidence_label": name,
                    "Reactome_analysis": pathway["Analysis"],
                    "Reactome_pathway_FDR": pathway["FDR"],
                })


        # ----------------------------------------------------
        # 7C. Identify ambiguous and high-priority annotations
        # ----------------------------------------------------

        precise = direct | reac

        review_reasons = []

        if not precise:
            review_reasons.append(
                "No mechanism-specific GO or selected Reactome support"
            )

        if len(precise) > 1:
            review_reasons.append(
                "Multiple potentially overlapping functional programs"
            )

        if (
            "Intrinsic apoptosis signaling" in direct
            and "Extrinsic apoptosis signaling" in direct
        ):
            review_reasons.append(
                "Both intrinsic and extrinsic GO annotations"
            )

        if (
            any("Necropt" in x for x in direct)
            and any("apoptos" in x.lower() for x in direct)
        ):
            review_reasons.append(
                "Death modalities overlap; contextual validation required"
            )

        if yes(old["Shared37"]):
            review_reasons.append(
                "Shared 37: prioritize biological review"
            )

        if yes(old["Cross_module_bridge"]):
            review_reasons.append(
                "INTERFACE connects multiple CORE modules"
            )

        if old["Proposed_primary_role"] != "UNREVIEWED":
            review_reasons.append(
                "Earlier seed annotation requires literature review"
            )


        # ----------------------------------------------------
        # 7D. Master annotation record
        # ----------------------------------------------------

        atlas.append({

            "Gene": gene,
            "Ensembl_ID": gene_row["Ensembl_ID"],
            "Tier": gene_row["F2_tier"],

            "Direction": gene_row["Direction"],
            "log2FC": gene_row["log2FoldChange"],
            "padj": gene_row["padj"],

            "GO_specific_roles": "; ".join(sorted(direct)),

            "GO_broad_annotations": "; ".join(sorted(broad)),

            "Reactome_pathway_modules": "; ".join(sorted(reac)),

            "GO_ids": "; ".join(sorted(go_ids)),

            "Reactome_unique_pathways": len({
                p[1] for p in matched_pathways
            }),

            "Reactome_pathway_details": "; ".join(
                f"{p[0]}:{p[1]} | {p[2]}"
                for p in sorted(matched_pathways)
            ),

            "Shared37": yes(old["Shared37"]),

            "Promoter_evidence": yes(
                old["Promoter_evidence"]
            ),

            "CORE_neighbors": old["CORE_neighbors"],

            "CORE_modules_touched": old["CORE_modules_touched"],

            "Cross_module_bridge": yes(
                old["Cross_module_bridge"]
            ),

            "Original_CORE_community": old["F4_CORE_module"],

            "Original_EXT_community": old["F4_EXT_module"],

            "Prior_seed_role_unverified": (
                ""
                if old["Proposed_primary_role"] == "UNREVIEWED"
                else old["Proposed_primary_role"]
            ),

            # We do not impose a single biological role yet.
            "Curated_primary_role": "UNASSIGNED",

            "Review_reasons": "; ".join(review_reasons),

            "Review_needed": bool(review_reasons),
        })


        # ----------------------------------------------------
        # 7E. Gene-by-program incidence matrix
        # ----------------------------------------------------

        matrix.append({
            "Gene": gene,
            "Tier": gene_row["F2_tier"],

            **{
                program: int(program in precise)
                for program in program_columns
            }
        })


    # ========================================================
    # 8. FINAL VALIDATION
    # ========================================================

    assert len(atlas) == 321

    assert sum(
        row["Shared37"] for row in atlas
    ) == 37

    assert sum(
        row["Cross_module_bridge"] for row in atlas
    ) == 25

    f2_by_gene = {
        row["Gene"].upper(): row
        for row in f2
    }

    assert all(
        row["Shared37"]
        == yes(
            f2_by_gene[row["Gene"]]["Old_126_gene_master"]
        )
        for row in atlas
    )

    print("F8/F2 integration validation passed.")


    # ========================================================
    # 9. EXPORT GENE-LEVEL TABLES
    # ========================================================

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )

    save_csv(
        OUT / "F8A_gene_role_atlas.csv",
        atlas,
        list(atlas[0])
    )

    save_csv(
        OUT / "F8A_gene_program_matrix.csv",
        matrix,
        list(matrix[0])
    )

    save_csv(
        OUT / "F8A_role_evidence_long.csv",
        evidence,
        list(evidence[0])
    )


    # ========================================================
    # 10. PROGRAM SUMMARY
    # ========================================================

    summaries = []

    for program in program_columns:

        members = [
            row for row in matrix
            if row[program] == 1
        ]

        summaries.append({

            "Program": program,

            "Evidence_class": (
                "specific_GO"
                if program in go_programs
                else "Reactome_F1_intersection"
            ),

            "Genes": len(members),

            "CORE": sum(
                row["Tier"] == "CORE"
                for row in members
            ),

            "INTERFACE": sum(
                row["Tier"] == "INTERFACE"
                for row in members
            ),

            "Note": (
                "Overlapping membership; "
                "not proof of pathway activation"
            )
        })

    save_csv(
        OUT / "F8A_program_summary.csv",
        summaries,
        list(summaries[0])
    )


    # ========================================================
    # 11. BIOLOGICAL REVIEW QUEUE
    # ========================================================

    review = [
        row for row in atlas
        if row["Review_needed"]
    ]

    # Shared genes and cross-module bridges get priority.
    # This is an ordering, not a biological importance score.

    review.sort(
        key=lambda row: (
            not row["Shared37"],
            not row["Cross_module_bridge"],
            -int(row["CORE_neighbors"] or 0),
            row["Gene"]
        )
    )

    save_csv(
        OUT / "F8A_review_queue.csv",
        review,
        list(atlas[0])
    )


    # ========================================================
    # 12. AUDIT REPORT
    # ========================================================

    report = [

        "F8A evidence-based role audit",

        f"Genes: {len(atlas)} (194 CORE + 127 INTERFACE)",

        f"Old/new/promoter shared: "
        f"{sum(row['Shared37'] for row in atlas)}",

        f"Genes with specific selected GO roles: "
        f"{sum(bool(row['GO_specific_roles']) for row in atlas)}",

        f"Genes with mapped F1 Reactome roles: "
        f"{sum(bool(row['Reactome_pathway_modules']) for row in atlas)}",

        f"Genes without a refined role: "
        f"{sum(not row['GO_specific_roles'] and not row['Reactome_pathway_modules'] for row in atlas)}",

        f"Gene-program evidence rows: {len(evidence)}",

        "WARNING: F2 covers only 23 selected GO terms.",

        "WARNING: Reactome results cover F1-significant pathways only.",

        "WARNING: pathway FDR is not gene-specific.",

        "WARNING: membership does not prove pathway activation.",

        "WARNING: STRING communities and biological functions differ.",

        "PRIMARY ROLES REMAIN UNASSIGNED pending curation.",

        "",
        "Program membership counts:"
    ]

    for row in summaries:
        report.append(
            f"  {row['Program']}: {row['Genes']}"
        )

    (OUT / "F8A_quality_report.txt").write_text(
        "\n".join(report) + "\n",
        encoding="utf-8"
    )

    print()
    print("\n".join(report))
    print()
    print("Outputs:", OUT)


if __name__ == "__main__":
    main()
