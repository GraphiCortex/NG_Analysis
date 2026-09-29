#!/usr/bin/env python3
"""
07_final_synthesis.py

Final synthesis/reporting stage for the R-Noel / BIOL_295 network project.

This script DOES NOT perform new discovery analysis.
It summarizes Phases 01-06 into:
    1. an experimental-candidate table,
    2. a separate external-connector table,
    3. a compact mechanistic hypothesis map,
    4. a Markdown synthesis suitable for discussion with Dr. Noel,
    5. publication-style summary figures.

Run from Project/
-----------------
    py .\Scripts\07_final_synthesis.py

INPUTS
------
Data/Processed/Network/Integrated/
    experimental_seed_evidence_matrix.csv
    external_connector_evidence_matrix.csv

Data/Processed/Network/Modules/
    module_summary.csv

OPTIONAL
--------
Data/Processed/Network/Mechanistic/Robustness/
    robustness_pair_metrics.csv

OUTPUTS
-------
Data/Processed/Network/FinalSynthesis/
    final_experimental_candidates.csv
    final_external_connectors.csv
    experimental_hypotheses.csv
    mechanistic_model_edges.csv
    final_synthesis.md
    final_synthesis_diagnostics.json

Figures/Network/FinalSynthesis/
    experimental_candidate_evidence.png
    external_connector_evidence.png
    mechanistic_model.png

IMPORTANT
---------
This stage distinguishes:
    - experimentally observed candidate genes
    - network-inferred external connectors

It does not merge them into one ranked list.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


# ============================================================
# Helpers
# ============================================================

def clean_gene(value):
    if pd.isna(value):
        return None
    value = str(value).strip().upper()
    return value or None


def write_json(path, obj):
    path.write_text(
        json.dumps(
            obj,
            indent=2,
            ensure_ascii=False,
            default=str,
        ),
        encoding="utf-8",
    )


def bool_value(value):
    if pd.isna(value):
        return False
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return value != 0
    return str(value).strip().upper() in {
        "TRUE", "YES", "Y", "1"
    }


def safe_num(value):
    try:
        if pd.isna(value):
            return 0.0
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def semicolon(values):
    return "; ".join(
        str(v)
        for v in values
        if v is not None
        and str(v).strip()
    )


# ============================================================
# Biological interpretation maps
# ============================================================

SEED_THEMES = {
    # Upstream / cell cycle
    "RB1": "Rb/E2F-cell-cycle control",
    "CCND1": "Rb/E2F-cell-cycle control",
    "CDKN2A": "Rb/E2F-cell-cycle control",
    "CDKN2D": "Rb/E2F-cell-cycle control",

    # Mitotic/checkpoint program
    "MYBL2": "aberrant mitotic/checkpoint program",
    "MELK": "aberrant mitotic/checkpoint program",
    "KIF14": "aberrant mitotic/checkpoint program",
    "BIRC5": "aberrant mitotic/checkpoint program",
    "AURKB": "aberrant mitotic/checkpoint program",
    "BUB1": "aberrant mitotic/checkpoint program",

    # Extrinsic death
    "FAS": "FAS/FADD extrinsic-death signaling",
    "FADD": "FAS/FADD extrinsic-death signaling",
    "MLKL": "regulated necrosis / necroptosis",

    # Intrinsic / mitochondrial death
    "BCL2L11": "mitochondrial apoptosis",
    "BOK": "mitochondrial apoptosis",

    # Mitophagy / autophagy
    "BNIP3L": "mitophagy / mitochondrial quality control",
    "ATG7": "autophagy / mitophagy",
    "GABARAPL1": "autophagy / mitophagy",

    # Trophic / signaling context
    "CXCL12": "trophic / survival signaling context",
    "BDNF": "trophic / survival signaling context",
    "KDR": "trophic / survival signaling context",

    # Other important seed connectors
    "HMGB1": "stress / extracellular signaling context",
    "HMGB2": "chromatin / stress context",
    "CAV1": "membrane-signaling bridge",
    "YAP1": "growth / stress signaling",
}

EXTERNAL_THEMES = {
    "ESR1": "causal + structural signaling connector",
    "ABL1": "causal + structural stress/signaling connector",
    "TP53": "causal transcriptional death-response connector",
    "MAPK9": "nonredundant structural JNK-family bridge",
    "PRKCA": "broad PKC-family structural bridge",
    "PRKCB": "redundant PKC-family structural bridge",
    "PRKCG": "mixed module/isolate PKC-family bridge",
    "FADS2": "nonredundant isolate-rescue lipid-metabolism connector",
    "PANK2": "nonredundant isolate-rescue mitochondrial/lipid connector",
}


def seed_hypothesis(gene, row):
    theme = SEED_THEMES.get(
        gene,
        "module-associated candidate",
    )

    if gene == "FAS":
        return (
            "Rb-family loss may engage FAS/FADD-dependent death "
            "in a developmental-stage-dependent manner."
        )

    if gene == "BNIP3L":
        return (
            "Mitophagy / mitochondrial quality-control signaling may "
            "connect Rb/E2F dysregulation to neuronal death."
        )

    if gene in {
        "MYBL2", "MELK", "KIF14",
        "BIRC5", "AURKB", "BUB1",
    }:
        return (
            "A specific aberrant mitotic/checkpoint program may lie "
            "upstream of neuronal death rather than generic cell-cycle "
            "re-entry alone."
        )

    if gene in {"CXCL12", "BDNF", "KDR"}:
        return (
            "Altered trophic/survival signaling may modify how neural "
            "cells respond to Rb-family loss."
        )

    if gene in {"BCL2L11", "BOK"}:
        return (
            "Mitochondrial apoptosis may contribute to the death program "
            "downstream of cell-cycle dysregulation."
        )

    if gene in {"RB1", "CCND1", "CDKN2D"}:
        return (
            "Cell-cycle control may provide the upstream state from which "
            "downstream stress and death pathways emerge."
        )

    return (
        f"{theme} is supported by one or more experimental/network "
        "evidence layers and merits context-specific validation."
    )


def seed_experiment(gene):
    if gene in {"FAS", "FADD"}:
        return (
            "Perturb FAS/FADD after pocket-protein loss; quantify neuronal "
            "survival/death and compare embryonic versus adult neural contexts."
        )

    if gene == "BNIP3L":
        return (
            "Perturb BNIP3L after pocket-protein loss; measure survival, "
            "mitochondrial membrane potential, and mitophagy/autophagic flux."
        )

    if gene in {
        "MYBL2", "MELK", "KIF14",
        "BIRC5", "AURKB", "BUB1",
    }:
        return (
            "Perturb the candidate after pocket-protein loss; measure "
            "cell-cycle/mitotic markers, DNA-damage/checkpoint activation, "
            "and neuronal survival."
        )

    if gene in {"CXCL12", "BDNF", "KDR"}:
        return (
            "Modulate the signaling axis after pocket-protein loss; measure "
            "survival and downstream MAPK/PI3K-AKT signaling."
        )

    if gene in {"BCL2L11", "BOK"}:
        return (
            "Perturb the mitochondrial-death effector and assess rescue of "
            "neuronal death together with mitochondrial apoptotic readouts."
        )

    if gene in {"RB1", "CCND1", "CDKN2D"}:
        return (
            "Measure temporal activation of the cell-cycle state relative to "
            "damage and death markers after pocket-protein loss."
        )

    return (
        "Validate expression/protein change in the relevant neural model, "
        "then perturb the gene and assay survival plus pathway-specific readouts."
    )


def external_hypothesis(gene, row):
    if gene == "ESR1":
        return (
            "ESR1 may connect otherwise distinct signaling modules while also "
            "participating in robust directed causal paths."
        )

    if gene == "ABL1":
        return (
            "ABL1 may provide a stress/signaling route linking Rb-associated "
            "states to downstream death programs."
        )

    if gene == "TP53":
        return (
            "TP53 may act as a causal transcriptional connector even though it "
            "is not a high-confidence direct STRING bridge to the seed set."
        )

    if gene == "MAPK9":
        return (
            "MAPK9 may be a nonredundant structural bridge linking several "
            "otherwise separated experimental programs."
        )

    if gene == "PRKCA":
        return (
            "PRKCA may broadly connect multiple seed modules, although its "
            "structural contribution is redundant with alternative routes."
        )

    if gene in {"FADS2", "PANK2"}:
        return (
            "The connector may uniquely rescue a previously isolated "
            "metabolic/lipid-associated seed neighborhood."
        )

    return (
        "The protein is a network-inferred connector and should be validated "
        "before being treated as part of the Rb-loss mechanism."
    )


def external_experiment(gene):
    if gene in {"ESR1", "ABL1", "TP53", "MAPK9", "PRKCA"}:
        return (
            "First confirm pathway/protein activation after pocket-protein loss "
            "(including phospho-state where relevant), then perturb the node and "
            "test whether death and module-specific readouts change."
        )

    if gene in {"FADS2", "PANK2"}:
        return (
            "Confirm expression/activity in the neural model, perturb the node, "
            "and test whether the isolated lipid/mitochondrial seed program and "
            "cell survival are altered."
        )

    return (
        "Validate presence/activity in the model before functional perturbation; "
        "treat as hypothesis-generating rather than experimentally established."
    )


# ============================================================
# Final candidate tables
# ============================================================

def prepare_seed_table(seed_df):
    df = seed_df.copy()

    df["Gene"] = df["Gene"].map(clean_gene)

    if "independent_evidence_flag_count" in df.columns:
        df.rename(
            columns={
                "independent_evidence_flag_count":
                    "evidence_feature_count"
            },
            inplace=True,
        )

    df["biological_theme"] = (
        df["Gene"]
        .map(SEED_THEMES)
        .fillna(
            "module-associated candidate"
        )
    )

    df["working_hypothesis"] = df.apply(
        lambda row:
            seed_hypothesis(
                row["Gene"],
                row,
            ),
        axis=1,
    )

    df["proposed_experiment"] = (
        df["Gene"]
        .map(seed_experiment)
    )

    # Descriptive follow-up tier, not a numerical rank.
    def followup_class(row):
        evidence_class = str(
            row.get(
                "integrated_evidence_class",
                "",
            )
        )

        gene = row["Gene"]

        if evidence_class == "multi_layer_strong":
            return "cross-layer experimental follow-up"

        if gene == "BNIP3L":
            return "mechanistic follow-up"

        if gene in {
            "MYBL2", "MELK", "KIF14",
            "BIRC5", "AURKB", "BUB1",
        }:
            return "mitotic-program follow-up"

        if bool_value(
            row.get(
                "adult_embryo_shared",
                False,
            )
        ):
            return "cross-dataset follow-up"

        return "context-dependent follow-up"

    df["followup_category"] = (
        df.apply(
            followup_class,
            axis=1,
        )
    )

    sort_cols = [
        col
        for col in [
            "evidence_feature_count",
            "phaseC_union_nonanchor_bottleneck_settings",
            "phaseC_required_internal_settings",
            "adult_embryo_shared",
            "andrusiak_overlap",
        ]
        if col in df.columns
    ]

    if sort_cols:
        df.sort_values(
            sort_cols,
            ascending=[False] * len(sort_cols),
            inplace=True,
        )

    return df


def prepare_external_table(external_df):
    df = external_df.copy()

    df["Gene"] = df["Gene"].map(clean_gene)

    if "independent_evidence_flag_count" in df.columns:
        df.rename(
            columns={
                "independent_evidence_flag_count":
                    "evidence_feature_count"
            },
            inplace=True,
        )

    df["biological_theme"] = (
        df["Gene"]
        .map(EXTERNAL_THEMES)
        .fillna(
            "network-inferred external connector"
        )
    )

    df["working_hypothesis"] = df.apply(
        lambda row:
            external_hypothesis(
                row["Gene"],
                row,
            ),
        axis=1,
    )

    df["proposed_experiment"] = (
        df["Gene"]
        .map(external_experiment)
    )

    sort_cols = [
        col
        for col in [
            "evidence_feature_count",
            "causal_union_nonanchor_bottleneck_settings",
            "causal_required_internal_settings",
            "unique_component_pair_count",
            "marginal_seed_component_increase_on_removal",
        ]
        if col in df.columns
    ]

    if sort_cols:
        df.sort_values(
            sort_cols,
            ascending=[False] * len(sort_cols),
            inplace=True,
        )

    return df


# ============================================================
# Experimental hypotheses table
# ============================================================

def make_hypothesis_table(seed_df, external_df):
    rows = [
        {
            "hypothesis_id": "H1",
            "hypothesis": (
                "FAS/FADD-dependent extrinsic-death signaling contributes "
                "to neuronal death after Rb-family loss and may differ "
                "between developmental stages."
            ),
            "primary_genes": "FAS; FADD",
            "evidence_origin": "experimental seeds + Phase B + Phase C",
            "critical_test": (
                "Perturb FAS/FADD after pocket-protein loss and compare "
                "survival/death phenotypes across embryonic and adult "
                "neural contexts."
            ),
        },
        {
            "hypothesis_id": "H2",
            "hypothesis": (
                "A specific aberrant mitotic/checkpoint program is causally "
                "upstream of neuronal death after pocket-protein loss."
            ),
            "primary_genes": (
                "MYBL2; MELK; KIF14; BIRC5; AURKB; BUB1"
            ),
            "evidence_origin": (
                "cross-dataset recurrence + Andrusiak overlap + "
                "Phase B module structure"
            ),
            "critical_test": (
                "Perturb selected mitotic/checkpoint genes after Rb-family "
                "loss and assay checkpoint/DNA-damage markers together with "
                "neuronal survival."
            ),
        },
        {
            "hypothesis_id": "H3",
            "hypothesis": (
                "Mitophagy / mitochondrial quality control contributes to "
                "the transition from Rb/E2F dysregulation to neuronal death."
            ),
            "primary_genes": "BNIP3L; ATG7; GABARAPL1",
            "evidence_origin": (
                "Phase B autophagy module + robust Phase C path analysis"
            ),
            "critical_test": (
                "Perturb BNIP3L and related autophagy machinery; measure "
                "mitochondrial integrity, autophagic flux, and survival."
            ),
        },
        {
            "hypothesis_id": "H4",
            "hypothesis": (
                "Altered trophic/survival signaling modulates the neuronal "
                "response to Rb-family loss."
            ),
            "primary_genes": "CXCL12; BDNF; KDR",
            "evidence_origin": (
                "cross-dataset recurrence/reversal + Phase B signaling module"
            ),
            "critical_test": (
                "Modulate the trophic signaling axis and assay survival plus "
                "MAPK/PI3K-AKT activity after pocket-protein loss."
            ),
        },
        {
            "hypothesis_id": "H5",
            "hypothesis": (
                "ESR1, ABL1, TP53, and MAPK9 are candidate missing connectors "
                "between experimentally observed modules."
            ),
            "primary_genes": "ESR1; ABL1; TP53; MAPK9",
            "evidence_origin": (
                "Phase C causal robustness + Phase D controlled expansion"
            ),
            "critical_test": (
                "Confirm activation/state in the neural model, then perturb "
                "each connector and test whether multiple downstream modules "
                "and cell-death phenotypes change."
            ),
        },
    ]

    return pd.DataFrame(rows)


# ============================================================
# Mechanistic model
# ============================================================

def make_model_edges():
    rows = [
        {
            "source": "Rb-family loss",
            "target": "Rb/E2F-cell-cycle dysregulation",
            "edge_type": "working model",
            "evidence_layer": "project premise + integrated analysis",
        },
        {
            "source": "Rb/E2F-cell-cycle dysregulation",
            "target": "aberrant mitotic/checkpoint program",
            "edge_type": "working model",
            "evidence_layer": "cross-dataset + module evidence",
        },
        {
            "source": "aberrant mitotic/checkpoint program",
            "target": "stress / damage-response state",
            "edge_type": "working model",
            "evidence_layer": "network/pathway synthesis",
        },
        {
            "source": "stress / damage-response state",
            "target": "FAS/FADD extrinsic death",
            "edge_type": "hypothesized branch",
            "evidence_layer": "experimental + structural + mechanistic",
        },
        {
            "source": "stress / damage-response state",
            "target": "mitochondrial apoptosis",
            "edge_type": "hypothesized branch",
            "evidence_layer": "experimental + pathway evidence",
        },
        {
            "source": "stress / damage-response state",
            "target": "mitophagy / mitochondrial quality control",
            "edge_type": "hypothesized branch",
            "evidence_layer": "module + mechanistic evidence",
        },
        {
            "source": "trophic / survival signaling",
            "target": "cell-death susceptibility",
            "edge_type": "context-modifying branch",
            "evidence_layer": "cross-dataset + structural evidence",
        },
        {
            "source": "external connectors",
            "target": "stress / signaling transitions",
            "edge_type": "network-inferred",
            "evidence_layer": "OmniPath + STRING expansion",
        },
    ]

    return pd.DataFrame(rows)


# ============================================================
# Figures
# ============================================================

def evidence_heatmap_seed(df, output_path):
    cols = [
        col
        for col in [
            "adult_embryo_shared",
            "reversed",
            "andrusiak_overlap",
            "multi_mechanism",
            "phaseB_cross_module_degree",
            "phaseC_required_internal_settings",
            "phaseC_union_nonanchor_bottleneck_settings",
        ]
        if col in df.columns
    ]

    if not cols:
        return

    top = df.head(20).copy()

    matrix = []

    for _, row in top.iterrows():
        vals = []

        for col in cols:
            value = row[col]

            if col in {
                "phaseB_cross_module_degree",
                "phaseC_required_internal_settings",
                "phaseC_union_nonanchor_bottleneck_settings",
            }:
                vals.append(
                    1 if safe_num(value) > 0 else 0
                )
            else:
                vals.append(
                    1 if bool_value(value) else 0
                )

        matrix.append(vals)

    plt.figure(
        figsize=(max(9, 1.2 * len(cols)), max(7, 0.35 * len(top) + 2))
    )

    plt.imshow(
        matrix,
        aspect="auto",
        vmin=0,
        vmax=1,
    )

    plt.xticks(
        range(len(cols)),
        cols,
        rotation=45,
        ha="right",
    )

    plt.yticks(
        range(len(top)),
        top["Gene"],
    )

    plt.title(
        "Experimental candidates: evidence-layer presence"
    )

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def evidence_heatmap_external(df, output_path):
    cols = [
        col
        for col in [
            "causal_required_internal_settings",
            "causal_union_nonanchor_bottleneck_settings",
            "nonisolate_module_count",
            "unique_component_pair_count",
            "marginal_seed_component_increase_on_removal",
        ]
        if col in df.columns
    ]

    if not cols:
        return

    top = df.head(19).copy()

    matrix = []

    for _, row in top.iterrows():
        vals = [
            1 if safe_num(row[col]) > 0 else 0
            for col in cols
        ]
        matrix.append(vals)

    plt.figure(
        figsize=(max(8, 1.4 * len(cols)), max(6, 0.35 * len(top) + 2))
    )

    plt.imshow(
        matrix,
        aspect="auto",
        vmin=0,
        vmax=1,
    )

    plt.xticks(
        range(len(cols)),
        cols,
        rotation=45,
        ha="right",
    )

    plt.yticks(
        range(len(top)),
        top["Gene"],
    )

    plt.title(
        "External connectors: evidence-layer presence"
    )

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def mechanistic_model_figure(edges, output_path):
    G = nx.DiGraph()

    for _, row in edges.iterrows():
        G.add_edge(
            row["source"],
            row["target"],
        )

    # Fixed positions for readability.
    pos = {
        "Rb-family loss": (0, 4),
        "Rb/E2F-cell-cycle dysregulation": (0, 3),
        "aberrant mitotic/checkpoint program": (0, 2),
        "stress / damage-response state": (0, 1),
        "FAS/FADD extrinsic death": (-2, 0),
        "mitochondrial apoptosis": (0, 0),
        "mitophagy / mitochondrial quality control": (2, 0),
        "trophic / survival signaling": (-3, 2),
        "cell-death susceptibility": (-3, 1),
        "external connectors": (3, 2),
        "stress / signaling transitions": (3, 1),
    }

    plt.figure(
        figsize=(12, 8)
    )

    nx.draw_networkx_nodes(
        G,
        pos,
        node_size=2600,
    )

    nx.draw_networkx_edges(
        G,
        pos,
        arrows=True,
        arrowsize=20,
        width=1.5,
    )

    nx.draw_networkx_labels(
        G,
        pos,
        font_size=8,
    )

    plt.title(
        "Working mechanistic synthesis of the R-Noel analysis"
    )

    plt.axis("off")

    plt.tight_layout()

    plt.savefig(
        output_path,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


# ============================================================
# Markdown report
# ============================================================

def write_markdown(
    seed_df,
    external_df,
    hypotheses,
    module_summary,
    output_path,
):
    lines = []

    lines.append(
        "# R-Noel Final Computational Synthesis"
    )
    lines.append("")
    lines.append(
        "This document summarizes the results of Phases 01-06. "
        "It does not introduce a new discovery analysis."
    )
    lines.append("")

    lines.append(
        "## Central biological interpretation"
    )
    lines.append("")
    lines.append(
        "The integrated analysis does not support a single generic "
        "cell-death program. Instead, it supports a working model in "
        "which Rb/E2F-associated cell-cycle dysregulation converges with "
        "an aberrant mitotic/checkpoint program and then branches into "
        "distinct death/stress pathways, including FAS/FADD-associated "
        "extrinsic death, mitochondrial apoptosis, and mitophagy/"
        "mitochondrial quality-control programs."
    )
    lines.append("")
    lines.append(
        "FAS is the clearest cross-layer experimental candidate in the "
        "current analysis. BNIP3L is a strong mechanistically inferred "
        "seed connector. MYBL2/MELK/KIF14/BIRC5/AURKB/BUB1 define a "
        "recurrent mitotic/checkpoint axis. CXCL12/BDNF/KDR define a "
        "parallel trophic/signaling context."
    )
    lines.append("")

    lines.append(
        "## Experimental candidates"
    )
    lines.append("")

    top_seed = seed_df.head(15)

    for _, row in top_seed.iterrows():
        lines.append(
            f"### {row['Gene']}"
        )
        lines.append("")
        lines.append(
            f"- Evidence class: `{row.get('integrated_evidence_class', '')}`"
        )
        lines.append(
            f"- Biological theme: {row.get('biological_theme', '')}"
        )
        lines.append(
            f"- Working hypothesis: {row.get('working_hypothesis', '')}"
        )
        lines.append(
            f"- Proposed experiment: {row.get('proposed_experiment', '')}"
        )
        lines.append("")

    lines.append(
        "## External network-inferred connectors"
    )
    lines.append("")
    lines.append(
        "These genes were not treated as equivalent to experimental seed "
        "genes. They entered through prior-knowledge networks and therefore "
        "remain hypothesis-generating until validated in the neural model."
    )
    lines.append("")

    for _, row in external_df.head(12).iterrows():
        lines.append(
            f"### {row['Gene']}"
        )
        lines.append("")
        lines.append(
            f"- Evidence class: `{row.get('integrated_evidence_class', '')}`"
        )
        lines.append(
            f"- Connector class: `{row.get('connector_class', '')}`"
        )
        lines.append(
            f"- Biological theme: {row.get('biological_theme', '')}"
        )
        lines.append(
            f"- Working hypothesis: {row.get('working_hypothesis', '')}"
        )
        lines.append(
            f"- Proposed experiment: {row.get('proposed_experiment', '')}"
        )
        lines.append("")

    lines.append(
        "## Experimental hypotheses"
    )
    lines.append("")

    for _, row in hypotheses.iterrows():
        lines.append(
            f"### {row['hypothesis_id']}"
        )
        lines.append("")
        lines.append(
            row["hypothesis"]
        )
        lines.append("")
        lines.append(
            f"Primary genes: {row['primary_genes']}"
        )
        lines.append("")
        lines.append(
            f"Evidence origin: {row['evidence_origin']}"
        )
        lines.append("")
        lines.append(
            f"Critical test: {row['critical_test']}"
        )
        lines.append("")

    lines.append(
        "## Interpretation boundaries"
    )
    lines.append("")
    lines.append(
        "- STRING associations are functional associations, not causal edges."
    )
    lines.append(
        "- OmniPath-directed paths are prior-knowledge-supported hypotheses "
        "and include mouse interactions that can derive from homology-translated "
        "knowledge."
    )
    lines.append(
        "- Small graph modules should be interpreted as network neighborhoods, "
        "not automatically as complete biological pathways."
    )
    lines.append(
        "- External connectors must be validated experimentally before being "
        "incorporated into the Rb-loss mechanism."
    )
    lines.append(
        "- Evidence-feature counts summarize detected features; they are not "
        "statistically independent evidence scores."
    )
    lines.append(
        "- The final mechanistic diagram is a working synthesis to be tested, "
        "not a demonstrated causal chain."
    )
    lines.append("")

    output_path.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


# ============================================================
# Main
# ============================================================

def main():
    project_root = (
        Path(__file__)
        .resolve()
        .parent
        .parent
    )

    processed_network = (
        project_root
        / "Data"
        / "Processed"
        / "Network"
    )

    integrated_dir = (
        processed_network
        / "Integrated"
    )

    modules_dir = (
        processed_network
        / "Modules"
    )

    output_dir = (
        processed_network
        / "FinalSynthesis"
    )

    figures_dir = (
        project_root
        / "Figures"
        / "Network"
        / "FinalSynthesis"
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    figures_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    seed_path = (
        integrated_dir
        / "experimental_seed_evidence_matrix.csv"
    )

    external_path = (
        integrated_dir
        / "external_connector_evidence_matrix.csv"
    )

    module_summary_path = (
        modules_dir
        / "module_summary.csv"
    )

    for path in [
        seed_path,
        external_path,
        module_summary_path,
    ]:
        if not path.exists():
            raise FileNotFoundError(path)

    seed_df = pd.read_csv(
        seed_path
    )

    external_df = pd.read_csv(
        external_path
    )

    module_summary = pd.read_csv(
        module_summary_path
    )

    seed_final = prepare_seed_table(
        seed_df
    )

    external_final = prepare_external_table(
        external_df
    )

    hypotheses = make_hypothesis_table(
        seed_final,
        external_final,
    )

    model_edges = make_model_edges()

    seed_final.to_csv(
        output_dir
        / "final_experimental_candidates.csv",
        index=False,
    )

    external_final.to_csv(
        output_dir
        / "final_external_connectors.csv",
        index=False,
    )

    hypotheses.to_csv(
        output_dir
        / "experimental_hypotheses.csv",
        index=False,
    )

    model_edges.to_csv(
        output_dir
        / "mechanistic_model_edges.csv",
        index=False,
    )

    write_markdown(
        seed_final,
        external_final,
        hypotheses,
        module_summary,
        output_dir
        / "final_synthesis.md",
    )

    evidence_heatmap_seed(
        seed_final,
        figures_dir
        / "experimental_candidate_evidence.png",
    )

    evidence_heatmap_external(
        external_final,
        figures_dir
        / "external_connector_evidence.png",
    )

    mechanistic_model_figure(
        model_edges,
        figures_dir
        / "mechanistic_model.png",
    )

    diagnostics = {
        "experimental_seed_count":
            len(seed_final),
        "external_connector_count":
            len(external_final),
        "hypothesis_count":
            len(hypotheses),
        "top_seed_by_feature_count":
            (
                seed_final[
                    [
                        "Gene",
                        "evidence_feature_count",
                    ]
                ]
                .head(10)
                .to_dict(
                    orient="records"
                )
                if "evidence_feature_count"
                in seed_final.columns
                else []
            ),
        "top_external_by_feature_count":
            (
                external_final[
                    [
                        "Gene",
                        "evidence_feature_count",
                    ]
                ]
                .head(10)
                .to_dict(
                    orient="records"
                )
                if "evidence_feature_count"
                in external_final.columns
                else []
            ),
        "important_note":
            (
                "Evidence-feature counts are descriptive and are not "
                "statistically independent scores."
            ),
    }

    write_json(
        output_dir
        / "final_synthesis_diagnostics.json",
        diagnostics,
    )

    print()
    print(
        "R-Noel final synthesis - Phase 07"
    )
    print(
        "================================"
    )

    print()
    print(
        "Experimental follow-up candidates"
    )
    print(
        "---------------------------------"
    )

    display_cols = [
        col
        for col in [
            "Gene",
            "integrated_evidence_class",
            "evidence_feature_count",
            "biological_theme",
            "followup_category",
        ]
        if col in seed_final.columns
    ]

    print(
        seed_final[
            display_cols
        ]
        .head(20)
        .to_string(
            index=False
        )
    )

    print()
    print(
        "External connector hypotheses"
    )
    print(
        "-----------------------------"
    )

    external_cols = [
        col
        for col in [
            "Gene",
            "integrated_evidence_class",
            "evidence_feature_count",
            "biological_theme",
        ]
        if col in external_final.columns
    ]

    print(
        external_final[
            external_cols
        ]
        .head(15)
        .to_string(
            index=False
        )
    )

    print()
    print(
        "Testable hypotheses"
    )
    print(
        "-------------------"
    )

    print(
        hypotheses[
            [
                "hypothesis_id",
                "primary_genes",
                "hypothesis",
            ]
        ]
        .to_string(
            index=False
        )
    )

    print()
    print(
        "Phase 07 complete."
    )

    print(
        f"Processed: {output_dir}"
    )

    print(
        f"Figures: {figures_dir}"
    )

    print()
    print(
        "The project now ends in testable biological hypotheses, "
        "not another network score."
    )


if __name__ == "__main__":
    main()
