
#!/usr/bin/env python3
r"""Export validated F9 results for the React 3D website.

Run from Project:
    py .\Fong_Pipeline\F9_website_export.py

Then copy the generated F9_overlay.json into:
    website/public/data/F9_overlay.json
"""

import csv
import json
from pathlib import Path

PROJECT = Path(__file__).resolve().parents[1]
BASE = PROJECT / "Data/Processed/FongDownstream"
OUT = BASE / "F9_Website/F9_overlay.json"


def load(folder, filename):
    path = BASE / folder / filename

    with path.open(
        encoding="utf-8-sig",
        newline=""
    ) as handle:
        return list(csv.DictReader(handle))


def yes(value):
    return str(value).strip().lower() in (
        "true", "yes", "1"
    )


def main():
    genes = load(
        "F9F_FongWideDirected",
        "F9F_gene_connectivity.csv"
    )

    interactions = load(
        "F9F_FongWideDirected",
        "F9F_direct_F2_interactions.csv"
    )

    routes = load(
        "F9I_D_EvidenceQuality",
        "F9ID_bypass_representative_routes.csv"
    )

    provenance = load(
        "F9I_D_EvidenceQuality",
        "F9ID_route_edge_provenance.csv"
    )

    candidate_folder = (
        BASE / "F9I_CandidateEvidenceAudit"
    )

    candidate_file = (
        candidate_folder
        / "F9I_nine_candidate_evidence_audit.csv"
    )

    candidates = (
        load(
            "F9I_CandidateEvidenceAudit",
            "F9I_nine_candidate_evidence_audit.csv"
        )
        if candidate_file.exists()
        else []
    )

    expected = (321, 963, 25, 42)
    actual = (
        len(genes),
        len(interactions),
        len(routes),
        len(provenance)
    )

    if actual != expected:
        raise ValueError(
            f"Unexpected source dimensions: "
            f"{actual}; expected {expected}"
        )

    gene_lookup = {}

    for row in genes:
        gene_lookup[row["Gene"]] = {
            "mapped": yes(row["Mapped_in_prior"]),
            "inE4X0": int(row["Incoming_E4_X0"]),
            "outE4X0": int(row["Outgoing_E4_X0"]),
            "inE4X1": int(row["Incoming_E4_X1"]),
            "outE4X1": int(row["Outgoing_E4_X1"]),
            "inE5X2": int(row["Incoming_E5_X2"]),
            "outE5X2": int(row["Outgoing_E5_X2"])
        }

    mapped = sum(
        item["mapped"]
        for item in gene_lookup.values()
    )

    if mapped != 257:
        raise ValueError(
            f"Expected 257 mapped genes; got {mapped}"
        )

    direct_edges = []

    for row in interactions:
        direct_edges.append({
            "source": row["Source"],
            "target": row["Target"],
            "effect": row["Effect_label"],
            "layers": row["Layers"].split("; "),
            "records": int(row["Supporting_records"]),
            "sources": row["Sources"],
            "references": row["References"]
        })

    bypass_routes = []

    for row in routes:
        bypass_routes.append({
            "id": row["Route_ID"],
            "source": row["Source"],
            "target": row["Target"],
            "blocker": row["Blocked_intermediate"],
            "setting": row["Setting"],
            "nodes": [
                value.strip()
                for value in row["Route"].split(" -> ")
            ],
            "outside": int(row["Outside_F2_count"]),
            "netSign": row["Net_database_sign"],
            "unknownEdges": int(
                row["Unknown_or_ambiguous_edge_count"]
            ),
            "hasReferences": yes(
                row["All_edges_have_reference_tokens"]
            ),
            "context": "NOT_ASSESSED"
        })

    edge_evidence = {}

    for row in provenance:
        edge_key = (
            f"{row['Source']}|{row['Target']}"
        )

        edge_evidence[edge_key] = {
            "effect": row["Prior_effect"],
            "layers": row["Layers"].split(";"),
            "sources": row["Source_databases"],
            "references": row["Reference_tokens"],
            "records": int(row["Supporting_records"]),
            "neuralValidated": False
        }

    candidate_ids = [
        "EXO1",
        "H2AFX",
        "TOPBP1",
        "CLSPN",
        "RAD51",
        "RAD54L",
        "TIMELESS",
        "RPA2",
        "LIG1"
    ]

    candidate_evidence = {}

    for row in candidates:
        candidate_evidence[row["Gene"]] = {
            "promoterFlag": yes(
                row["Julian_E2F4_vs_E2F3_sig_promoter"]
            ),
            "f6Overlap": yes(
                row["F6_Julian_binding_expression_overlap"]
            ),
            "significantBindingExpression": yes(
                row[
                    "F6_Julian_significant_binding_and_expression"
                ]
            )
        }

    output = {
        "metadata": {
            "phase": "F9_WEB_OVERLAY",
            "directInteractions": 963,
            "mappedF2Genes": 257,
            "candidateCount": 9,
            "exampleBypassRoutes": 25,
            "edgeProvenanceRecords": 42,
            "note": (
                "Database associations and illustrative paths, "
                "not demonstrated neuronal mechanisms."
            )
        },

        "genes": gene_lookup,
        "directEdges": direct_edges,
        "candidateIds": candidate_ids,
        "candidateEvidence": candidate_evidence,
        "bypassRoutes": bypass_routes,
        "edgeEvidence": edge_evidence,

        "findings": {
            "F9E": {
                "observedBridges": 25,
                "nullMean": 53.40,
                "empiricalPApprox": 0.00050,
                "scope": (
                    "Bipartite degree-preserving STRING null; "
                    "not individual-gene significance"
                )
            },
            "F9G_F9H": {
                "initialSignificantPairs": 6,
                "exclusiveSignificantPairs": 0,
                "scope": (
                    "Module annotations overlap; "
                    "do not claim distinct pathway communication"
                )
            },
            "F9I": {
                "chek1RequiredE4X0": 16,
                "otherRequiredE4X0": 1,
                "requiredRelaxed": 0,
                "scope": (
                    "Constrained graph dependency, "
                    "not a cellular perturbation result"
                )
            }
        }
    }

    OUT.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    OUT.write_text(
        json.dumps(
            output,
            ensure_ascii=False,
            separators=(",", ":")
        ),
        encoding="utf-8"
    )

    print("F9 website overlay PASS:", OUT)
    print(
        f"  F9F: {len(direct_edges)} direct interactions; "
        f"{len(gene_lookup)} F2 genes"
    )
    print(
        f"  F9I-D: {len(bypass_routes)} example routes; "
        f"{len(edge_evidence)} distinct route edges"
    )
    print(
        "  Copy F9_overlay.json to "
        "website/public/data/F9_overlay.json"
    )


if __name__ == "__main__":
    main()
