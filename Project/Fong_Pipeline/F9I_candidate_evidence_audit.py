#!/usr/bin/env python3
"""F9I-A: audit nine predeclared F4 INTERFACE connectors with F6/F8B/F9F evidence.

Place in Project/Fong_Pipeline/F9I_candidate_evidence_audit.py.
Run from Project/: py .\\Fong_Pipeline\\F9I_candidate_evidence_audit.py

No network requests. Standard library only. This is an evidence audit,
not an independent significance test or causal-driver ranking.
"""
from __future__ import annotations
import csv
import json
import re
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "Data" / "Processed" / "FongDownstream"
F6 = BASE / "F6_Integrated" / "F6_integrated_evidence.xlsx"
F8B = BASE / "F8_FunctionalAtlas" / "F8B_BiologicalModules" / "F8B_gene_role_atlas.csv"
F9F = BASE / "F9F_FongWideDirected"
OUT = BASE / "F9I_CandidateEvidenceAudit"

CANDIDATES = ("EXO1", "H2AFX", "TOPBP1", "CLSPN", "RAD51", "RAD54L", "TIMELESS", "RPA2", "LIG1")
S = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
R = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
PKG = "{http://schemas.openxmlformats.org/package/2006/relationships}"


def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def to_bool(v):
    return str(v).strip().lower() in ("1", "true", "yes")


def col_index(address):
    letters = re.match(r"[A-Z]+", address).group(0)
    n = 0
    for c in letters:
        n = n * 26 + ord(c) - 64
    return n - 1


def read_xlsx_sheet(path, sheet_name):
    """Read cell values from xlsx without changing the source workbook."""
    with zipfile.ZipFile(path) as z:
        workbook = ET.fromstring(z.read("xl/workbook.xml"))
        relmap = {r.attrib["Id"]: r.attrib["Target"] for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
        sheet_ref = next((s for s in workbook.iter(S + "sheet") if s.attrib.get("name") == sheet_name), None)
        if sheet_ref is None:
            raise ValueError(f"Missing F6 worksheet: {sheet_name}")
        target = relmap[sheet_ref.attrib[R + "id"]].lstrip("/")
        if not target.startswith("xl/"):
            target = "xl/" + target
        target = re.sub(r"/\./", "/", target)
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            ss = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in ss.iter(S + "si"):
                shared.append("".join(t.text or "" for t in si.iter(S + "t")))
        xml = ET.fromstring(z.read(target))
        table = []
        for row in xml.iter(S + "row"):
            cells = {}
            for c in row.findall(S + "c"):
                typ = c.attrib.get("t")
                v = c.find(S + "v")
                if typ == "inlineStr":
                    value = "".join(t.text or "" for t in c.iter(S + "t"))
                elif v is None or v.text is None:
                    continue
                elif typ == "s":
                    value = shared[int(v.text)]
                elif typ == "b":
                    value = v.text == "1"
                elif typ in ("str", "e"):
                    value = v.text
                else:
                    try:
                        value = float(v.text)
                    except ValueError:
                        value = v.text
                cells[col_index(c.attrib["r"])] = value
            if cells:
                table.append([cells.get(i) for i in range(max(cells) + 1)])
    if not table:
        raise ValueError("Empty F6 worksheet")
    headers = [str(x) if x is not None else "" for x in table[0]]
    return [{h: (row[i] if i < len(row) else None) for i, h in enumerate(headers)} for row in table[1:]]


def write_csv(path, records):
    if not records:
        raise ValueError(f"Refusing to write an empty table: {path}")
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    paths = [F6, F8B, F9F / "F9F_gene_connectivity.csv", F9F / "F9F_direct_F2_interactions.csv"]
    for p in paths:
        if not p.is_file():
            raise FileNotFoundError(f"Required input missing: {p}")

    f6rows = read_xlsx_sheet(F6, "Experimental evidence")
    f6 = {str(x["Gene"]).upper(): x for x in f6rows if x.get("Gene")}
    atlas = {x["Gene"]: x for x in read_csv(F8B)}
    f9f = {x["Gene"]: x for x in read_csv(F9F / "F9F_gene_connectivity.csv")}
    prior = read_csv(F9F / "F9F_direct_F2_interactions.csv")

    assert len(atlas) == 321 and len(f9f) == 321 and len(prior) == 963
    audit = []
    outgoing_rows = []
    for g in CANDIDATES:
        r = atlas[g]
        h = f6[g]
        p = f9f[g]
        assert r["Tier"] == "INTERFACE" and to_bool(r["Cross_module_bridge"])
        assert r["Direction"] == "Up" and to_bool(r["Promoter_evidence"])
        assert to_bool(h["Julian_E2F4_vs_E2F3_sig_promoter"])
        assert to_bool(h["Fong_F2"]) and to_bool(h["Fong_strict_DE"])
        out = [x for x in prior if x["Source"] == g]
        inc = [x for x in prior if x["Target"] == g]
        mapped = to_bool(p["Mapped_in_prior"])
        if not mapped:
            assert not out and not inc
        audit.append({
            "Gene": g,
            "Fong_direction": r["Direction"],
            "Fong_log2FC": r["log2FC"],
            "Fong_adjusted_p": r["padj"],
            "F4_CORE_neighbors": r["CORE_neighbors"],
            "F4_CORE_communities": r["CORE_modules_touched"],
            "Julian_E2F4_vs_E2F3_sig_promoter": bool(h["Julian_E2F4_vs_E2F3_sig_promoter"]),
            "F6_Julian_binding_expression_overlap": bool(h["Julian_binding_expression_overlap"]),
            "F6_Julian_significant_binding_and_expression": bool(h["Julian_significant_binding_and_expression"]),
            "F6_old_master": bool(h["Old_master"]),
            "F9F_prior_graph_mapped": mapped,
            "F9F_E4_X0_outgoing": p["Outgoing_E4_X0"],
            "F9F_E4_X1_outgoing": p["Outgoing_E4_X1"],
            "F9F_E5_X2_outgoing": p["Outgoing_E5_X2"],
            "F9F_direct_outgoing_F2": len(out),
            "F9F_direct_incoming_F2": len(inc),
            "Direct_F2_targets": "; ".join(f'{x["Target"]} ({x["Effect_label"]})' for x in out),
            "Biological_modules": r["Biological_modules"],
            "Interpretation_note": (
                "No node in selected prior graph" if not mapped else
                "No outgoing path in tested F9F settings" if int(p["Outgoing_E5_X2"]) == 0 else
                "Outside-F2 intermediate needed to reach F2 genes under E4_X1" if int(p["Outgoing_E4_X0"]) == 0 else
                "Outgoing directed relationships available; evaluate pathway provenance individually"
            ),
        })
        for x in out:
            outgoing_rows.append({"Candidate": g, **x})

    OUT.mkdir(parents=True, exist_ok=True)
    write_csv(OUT / "F9I_nine_candidate_evidence_audit.csv", audit)
    if outgoing_rows:
        write_csv(OUT / "F9I_direct_candidate_interactions.csv", outgoing_rows)
    report = {
        "analysis": "F9I-A descriptive nine-gene provenance audit",
        "candidates": list(CANDIDATES),
        "all_Fong_up": sum(x["Fong_direction"] == "Up" for x in audit),
        "F6_Julian_E2F4_vs_E2F3_promoter_flag": sum(x["Julian_E2F4_vs_E2F3_sig_promoter"] for x in audit),
        "F6_Julian_binding_expression_overlap": sum(x["F6_Julian_binding_expression_overlap"] for x in audit),
        "F6_Julian_significant_binding_and_expression": sum(x["F6_Julian_significant_binding_and_expression"] for x in audit),
        "mapped_in_F9F": sum(x["F9F_prior_graph_mapped"] for x in audit),
        "no_outgoing_E5_X2": [x["Gene"] for x in audit if x["F9F_prior_graph_mapped"] and int(x["F9F_E5_X2_outgoing"]) == 0],
        "direct_outgoing_edge_count": len(outgoing_rows),
        "limitations": [
            "Selection is post-hoc: nine preselected F4 bridges with a Julien promoter flag.",
            "F6 promoter differential flag, binding-expression overlap, and significant binding-expression are three distinct fields.",
            "Fong gene expression, Julien promoter analyses, and network priors are from different measurements/contexts.",
            "No direct target validation, causal effect, or gene-level enrichment is established.",
            "Prior reachability is limited by graph coverage; unmapped does not mean biologically disconnected."
        ],
    }
    (OUT / "F9I_validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("F9I-A PASS: audited", len(audit), "preselected genes")
    print("Julien promoter flags / binding-expression overlap / jointly significant:",
          report["F6_Julian_E2F4_vs_E2F3_promoter_flag"],
          report["F6_Julian_binding_expression_overlap"],
          report["F6_Julian_significant_binding_and_expression"])
    print("Mapped:", report["mapped_in_F9F"], "direct outgoing F2 edges:", len(outgoing_rows))
    print("No outgoing F9F routes:", report["no_outgoing_E5_X2"])
    print("Output:", OUT)


if __name__ == "__main__":
    main()
