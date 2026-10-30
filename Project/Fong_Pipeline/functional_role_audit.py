
#!/usr/bin/env python3
"""F8: Non-destructive first pass for the Fong functional atlas.

Requires Python standard library only. Reads F6 as .xlsx without altering it.
"""
from __future__ import annotations

import csv
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "Data/Processed/FongDownstream/F6_Integrated/F6_integrated_evidence.xlsx"
OUTPUT = ROOT / "Data/Processed/FongDownstream/F8_FunctionalAtlas"

S = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
P = "{http://schemas.openxmlformats.org/package/2006/relationships}"

# A small, source-informed starting set, NOT an exhaustive biological annotation.
# Each entry is a hypothesis for manual curation, not an experimental pathway claim.
SEEDS = {
    "PMAIP1": ("Intrinsic apoptosis", "BH3-only NOXA; mitochondrial apoptotic regulation"),
    "BBC3": ("Intrinsic apoptosis", "BH3-only PUMA; mitochondrial apoptotic regulation"),
    "FAS": ("Extrinsic apoptosis", "Death receptor; context-dependent outcome"),
    "TNFRSF10B": ("Extrinsic apoptosis", "TRAIL death-receptor signaling"),
    "MLKL": ("Necroptosis", "Execution machinery; activation requires phosphorylation"),
    "ZBP1": ("Regulated inflammatory death", "Innate sensor associated with regulated death pathways"),
    "PYCARD": ("Inflammasome signaling", "ASC adaptor; do not equate with MLKL-mediated necroptosis"),
    "BDNF": ("Trophic / survival signaling", "Neurotrophic ligand; context-dependent transcriptional direction"),
    "KDR": ("Trophic / survival signaling", "Receptor tyrosine kinase; VEGF-related signaling"),
    "IRS2": ("Trophic / survival signaling", "Insulin/IGF signal transduction adaptor"),
    "PCNA": ("DNA replication / checkpoint", "Replication clamp and DNA-repair support"),
    "CLSPN": ("DNA replication / checkpoint", "Replication-stress checkpoint mediator"),
    "TOPBP1": ("DNA replication / checkpoint", "ATR checkpoint and replication initiation scaffold"),
    "CHEK1": ("DNA replication / checkpoint", "DNA-damage and replication-stress checkpoint kinase"),
    "CHEK2": ("DNA damage response", "DNA-damage checkpoint kinase"),
    "H2AFX": ("DNA damage response", "Histone H2AX; DNA-damage signaling marker when phosphorylated"),
    "EXO1": ("DNA repair", "DNA-end resection and mismatch-repair-associated nuclease"),
    "RAD51": ("DNA repair", "Homologous-recombination recombinase"),
    "BRCA1": ("DNA repair", "DNA-damage response and homologous recombination"),
    "CDK1": ("Cell cycle / mitosis", "Mitotic cyclin-dependent kinase"),
    "CDK2": ("Cell cycle / replication", "G1/S and S-phase cyclin-dependent kinase"),
    "AURKB": ("Cell cycle / mitosis", "Chromosome passenger complex kinase"),
    "BIRC5": ("Cell cycle / mitosis", "Survivin; mitotic chromosome passenger and antiapoptotic association"),
    "E2F1": ("Cell-cycle transcriptional control", "E2F activator with context-dependent fate roles"),
    "MYC": ("Cell-cycle transcriptional control", "Pleiotropic regulator with fate-linked target evidence"),
    "DDIT3": ("ER stress / UPR", "CHOP stress-induced transcription factor"),
    "TRIB3": ("ER stress / UPR", "Stress-response pseudokinase"),
}


def col_number(a1: str) -> int:
    letters = re.match(r"[A-Z]+", a1.upper()).group()
    n = 0
    for char in letters:
        n = n * 26 + ord(char) - ord("A") + 1
    return n - 1


def read_sheet(path: Path, title: str) -> list[dict[str, str]]:
    """Read simple worksheet values (including shared strings and booleans)."""
    with zipfile.ZipFile(path) as z:
        book = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        targets = {
            r.attrib["Id"]: r.attrib["Target"]
            for r in rels.findall(f"{P}Relationship")
        }

        for sheet in book.findall(f"{S}sheets/{S}sheet"):
            if sheet.attrib["name"] == title:
                target = targets[sheet.attrib[f"{R}id"]]
                member = (
                    target.lstrip("/")
                    if target.startswith("/")
                    else "xl/" + target
                )
                break
        else:
            raise KeyError(f"Sheet not found: {title}")

        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            shared = ET.fromstring(z.read("xl/sharedStrings.xml"))
            strings = [
                "".join(t.text or "" for t in si.iter(f"{S}t"))
                for si in shared.findall(f"{S}si")
            ]

        page = ET.fromstring(z.read(member))
        matrix = []

        for row in page.findall(f"{S}sheetData/{S}row"):
            values: dict[int, str] = {}

            for cell in row.findall(f"{S}c"):
                v = cell.find(f"{S}v")
                typ = cell.attrib.get("t")

                if typ == "inlineStr":
                    ins = cell.find(f"{S}is")
                    val = (
                        "".join(t.text or "" for t in ins.iter(f"{S}t"))
                        if ins is not None
                        else ""
                    )
                elif v is None:
                    val = ""
                elif typ == "s":
                    val = strings[int(v.text)]
                elif typ == "b":
                    val = "TRUE" if v.text == "1" else "FALSE"
                else:
                    val = v.text or ""

                values[col_number(cell.attrib["r"])] = val

            matrix.append(values)

    if not matrix:
        return []

    columns = [
        matrix[0].get(i, "")
        for i in range(max(matrix[0]) + 1)
    ]

    return [
        {
            name: row.get(i, "")
            for i, name in enumerate(columns)
            if name
        }
        for row in matrix[1:]
    ]


def true(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def number(value: str) -> float:
    try:
        return float(value)
    except (ValueError, TypeError):
        return 0.0


def write_csv(filename: str, rows: list[dict], fields: list[str]) -> None:
    with (
        OUTPUT / filename
    ).open("w", newline="", encoding="utf-8-sig") as out:
        writer = csv.DictWriter(
            out,
            fieldnames=fields,
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    if not INPUT.is_file():
        sys.exit(f"Missing F6 workbook: {INPUT}")

    evidence = read_sheet(INPUT, "Experimental evidence")
    prior = read_sheet(INPUT, "External prior nodes")

    assert len({r["Gene"] for r in evidence}) == len(evidence), (
        "Duplicate evidence gene IDs; audit required"
    )

    genes = [
        r for r in evidence
        if r["F2_tier"] in {"CORE", "INTERFACE"}
    ]

    assert len(genes) == 321, (
        f"Expected 321 F2 genes, got {len(genes)}"
    )

    assert Counter(r["F2_tier"] for r in genes) == {
        "CORE": 194,
        "INTERFACE": 127
    }

    rows = []

    for r in genes:
        gene = r["Gene"].upper().strip()
        tier = r["F2_tier"]

        n_core = (
            int(number(r["F4_interface_CORE_neighbor_count"]))
            if tier == "INTERFACE"
            else 0
        )

        n_modules = (
            int(number(r["F4_interface_CORE_modules_touched_count"]))
            if tier == "INTERFACE"
            else 0
        )

        shared37 = (
            r["Integration_class"]
            == "Old + new + promoter convergence"
        )

        seed = SEEDS.get(gene)

        rows.append({
            "Gene": gene,
            "Tier": tier,
            "Fong_direction": r["Fong_direction"],
            "Fong_log2FC": r["Fong_log2FC"],
            "Fong_padj": r["Fong_padj"],
            "Original_F2_categories": r["F2_categories"],
            "F4_CORE_module": r["F4_CORE_module_label"],
            "F4_EXT_module": r["F4_EXT_module_label"],
            "Proposed_primary_role": (
                seed[0] if seed else "UNREVIEWED"
            ),
            "Proposed_role_note": seed[1] if seed else "",
            "Biological_review_status": (
                "SEED_REQUIRES_REVIEW"
                if seed
                else "UNASSIGNED"
            ),
            "Secondary_roles": "",
            "Shared37": shared37,
            "Promoter_evidence": true(
                r["Julian_E2F4_vs_E2F3_sig_promoter"]
            ),
            "Old_master": true(r["Old_master"]),
            "CORE_neighbors": n_core,
            "CORE_modules_touched": n_modules,
            "Structural_bridge_to_CORE": (
                tier == "INTERFACE" and n_core >= 1
            ),
            "Cross_module_bridge": (
                tier == "INTERFACE" and n_modules >= 2
            ),
            "F4_interface_betweenness": (
                r["F4_interface_betweenness"]
                if tier == "INTERFACE"
                else ""
            ),
            "Prior_knowledge_directed_support": (
                int(number(r["F5_TierA_outgoing_count"])) > 0
                or int(number(r["F5_TierA_incoming_count"])) > 0
            ),
            "Integration_class": r["Integration_class"],
            "Candidate_panels": r["Candidate_panels"],
        })

    shared = [r for r in rows if r["Shared37"]]

    contact = [
        r for r in rows
        if r["Structural_bridge_to_CORE"]
    ]

    cross = [
        r for r in rows
        if r["Cross_module_bridge"]
    ]

    assert len(shared) == 37, (
        f"Expected shared 37, got {len(shared)}"
    )

    assert all(r["Promoter_evidence"] for r in shared), (
        "Shared37/promoter mismatch"
    )

    assert len(contact) == 86, (
        f"Expected 86 CORE-contact bridges, got {len(contact)}"
    )

    assert len(cross) == 25, (
        f"Expected 25 multi-module bridges, got {len(cross)}"
    )

    OUTPUT.mkdir(parents=True, exist_ok=True)

    columns = list(rows[0])

    write_csv(
        "gene_role_draft.csv",
        sorted(rows, key=lambda x: (x["Tier"], x["Gene"])),
        columns
    )

    write_csv(
        "shared_37_genes.csv",
        sorted(shared, key=lambda x: x["Gene"]),
        columns
    )

    write_csv(
        "interface_bridge_candidates.csv",
        sorted(
            contact,
            key=lambda x: (
                -x["CORE_modules_touched"],
                -x["CORE_neighbors"],
                -number(x["F4_interface_betweenness"])
            )
        ),
        columns
    )

    prior_cols = [
        "Gene",
        "Old_external_connector",
        "New_F5_representative_external",
        "New_F5_path_occurrences",
        "New_F5_distinct_pairs",
        "New_F5_settings",
        "Observed_in_current_strict_Fong",
        "Observed_in_F2",
        "Observed_in_old_master",
        "Julian_sig_promoter",
    ]

    write_csv(
        "prior_node_inventory.csv",
        prior,
        prior_cols
    )

    with (
        OUTPUT / "audit_summary.txt"
    ).open("w", encoding="utf-8") as out:
        out.write(
            "F8 Functional atlas: first-pass audit "
            "(NOT curated GO/pathway roles)\n"
        )

        out.write(f"Total F2 genes: {len(rows)}\n")

        out.write(
            f"CORE: {sum(r['Tier'] == 'CORE' for r in rows)}\n"
        )

        out.write(
            f"INTERFACE: {sum(r['Tier'] == 'INTERFACE' for r in rows)}\n"
        )

        out.write(f"Shared 37: {len(shared)}\n")

        out.write(
            f"INTERFACE with CORE contact: {len(contact)}\n"
        )

        out.write(
            "INTERFACE touching 2+ CORE modules: "
            f"{len(cross)}\n"
        )

        out.write(
            "Source-informed seed roles (unreviewed): "
            f"{sum(r['Biological_review_status'] == 'SEED_REQUIRES_REVIEW' for r in rows)}\n"
        )

        out.write(
            f"External-prior sheet rows: {len(prior)} "
            "(NOT all necessarily external-only)\n"
        )

        out.write(
            "IMPORTANT: STRING links are functional associations, "
            "not demonstrated causal regulation.\n"
        )

        out.write(
            "IMPORTANT: directed prior nodes can overlap observed genes; "
            "do not simply label every row external.\n"
        )

    print(
        (OUTPUT / "audit_summary.txt").read_text(
            encoding="utf-8"
        )
    )

    print(f"Outputs: {OUTPUT}")


if __name__ == "__main__":
    main()
