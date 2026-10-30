#!/usr/bin/env python3

"""
Phase F1 — Global functional enrichment of strict Fong DEGs
============================================================

Purpose
-------
Characterize the global biological programs represented among the
strict Fong 2022 differential-expression genes BEFORE selecting
cell-death/survival genes.

Input
-----
Data/Processed/FongDownstream/F0_Audit/
    Fong_downstream_phase_F0_audit.xlsx

Sheet:
    Fong strict DEG master

Queries
-------
1. ALL  : all strict Fong DEGs
2. UP   : log2FoldChange > 0
3. DOWN : log2FoldChange < 0

Functional sources
------------------
- GO Biological Process (GO:BP)
- Reactome (REAC)

Statistics
----------
Over-representation analysis through g:Profiler g:GOSt.

Multiple-testing correction:
    Benjamini-Hochberg FDR

Threshold:
    adjusted p < 0.05

Background
----------
Annotated Mus musculus gene universe.

IMPORTANT:
This is not the ideal experiment-specific background because the currently
available Fong workbook contains only the already-filtered significant genes.
If the complete unfiltered DESeq result table becomes available, the analysis
should be rerun with all tested genes as the custom background.

Outputs
-------
Data/Processed/FongDownstream/F1_Enrichment/

    F1_global_enrichment.xlsx
    all_enrichment.csv
    up_enrichment.csv
    down_enrichment.csv
    enrichment_metadata.json

    figures/
        ALL_GO_BP_top20.png
        ALL_REAC_top20.png
        UP_GO_BP_top20.png
        UP_REAC_top20.png
        DOWN_GO_BP_top20.png
        DOWN_REAC_top20.png

Run
---
From Project/:

    py .\\Scripts\\F1_global_enrichment.py
"""

from __future__ import annotations

from pathlib import Path
import json
import math
import sys
from typing import Any

import pandas as pd
import requests
import matplotlib.pyplot as plt
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F0_Audit"
    / "Fong_downstream_phase_F0_audit.xlsx"
)

INPUT_SHEET = "Fong strict DEG master"

OUTPUT_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F1_Enrichment"
)

FIGURE_DIR = OUTPUT_DIR / "figures"

OUTPUT_XLSX = OUTPUT_DIR / "F1_global_enrichment.xlsx"
OUTPUT_METADATA = OUTPUT_DIR / "enrichment_metadata.json"

GPROFILER_URL = (
    "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
)

GPROFILER_VERSION_URL = (
    "https://biit.cs.ut.ee/gprofiler/api/util/data_versions/"
)

ORGANISM = "mmusculus"

SOURCES = [
    "GO:BP",
    "REAC",
]

FDR_THRESHOLD = 0.05

REQUEST_TIMEOUT = 120

TOP_N_PLOT = 20


# ============================================================
# Helpers
# ============================================================

def normalize_gene(value: Any) -> str | None:
    if pd.isna(value):
        return None

    gene = str(value).strip()

    if not gene:
        return None

    return gene.upper()


def safe_neg_log10(value: float) -> float:
    """
    Convert p-value to -log10(p), avoiding log10(0).
    """
    if value <= 0:
        value = 1e-300

    return -math.log10(value)


def query_gprofiler(
    genes: list[str],
    query_name: str,
) -> tuple[pd.DataFrame, dict]:
    """
    Run one g:Profiler enrichment query.

    Returns
    -------
    results_df
    full_json_response
    """

    payload = {
        "organism": ORGANISM,
        "query": genes,
        "sources": SOURCES,
        "user_threshold": FDR_THRESHOLD,
        "significance_threshold_method": "fdr",
        "domain_scope": "annotated",
        "ordered": False,
        "all_results": False,
        "no_evidences": False,
    }

    print(
        f"[INFO] Querying g:Profiler for {query_name}: "
        f"{len(genes):,} genes..."
    )

    response = requests.post(
        GPROFILER_URL,
        json=payload,
        timeout=REQUEST_TIMEOUT,
        headers={
            "User-Agent": "BIOL_295_FongDownstream/1.0"
        },
    )

    response.raise_for_status()

    data = response.json()

    if "result" not in data:
        raise RuntimeError(
            f"Unexpected g:Profiler response for {query_name}."
        )

    results = pd.DataFrame(data["result"])

    if results.empty:
        print(
            f"[WARNING] No significant enrichment terms for "
            f"{query_name}."
        )

        return results, data

    results.insert(
        0,
        "analysis_group",
        query_name,
    )

    return results, data


def extract_intersection_genes(
    row: pd.Series,
    full_response: dict,
) -> str:
    """
    Recover contributing genes for an enrichment term.

    g:Profiler returns an 'intersections' structure whose entries correspond
    to the Ensembl IDs stored in genes_metadata.

    This function converts that structure into a semicolon-separated gene-ID
    string. In most cases these will be Ensembl IDs.

    Original gene symbols remain available separately in the Fong master.
    """

    intersections = row.get("intersections")

    if not isinstance(intersections, list):
        return ""

    meta = full_response.get(
        "meta",
        {},
    )

    genes_metadata = meta.get(
        "genes_metadata",
        {},
    )

    query_meta = genes_metadata.get(
        "query",
        {},
    )

    if not query_meta:
        return ""

    # Standard single-query API response.
    first_query_key = next(
        iter(query_meta.keys()),
        None,
    )

    if first_query_key is None:
        return ""

    ensgs = query_meta.get(
        first_query_key,
        {},
    ).get(
        "ensgs",
        [],
    )

    if not ensgs:
        return ""

    contributing = []

    for idx, evidence in enumerate(intersections):
        if idx >= len(ensgs):
            break

        if evidence:
            contributing.append(
                str(ensgs[idx])
            )

    return ";".join(contributing)


def add_derived_statistics(
    results: pd.DataFrame,
    query_size_fallback: int,
) -> pd.DataFrame:
    """
    Add interpretable enrichment statistics.
    """

    if results.empty:
        return results

    df = results.copy()

    # g:Profiler returns corrected p-values in p_value.
    df["neg_log10_FDR"] = (
        df["p_value"]
        .astype(float)
        .map(safe_neg_log10)
    )

    if "query_size" not in df.columns:
        df["query_size"] = query_size_fallback

    # Fold enrichment:
    #
    # observed fraction / background fraction
    #
    # (intersection/query_size)
    # -------------------------
    # (term_size/domain_size)
    #
    numerator = (
        df["intersection_size"]
        /
        df["query_size"]
    )

    denominator = (
        df["term_size"]
        /
        df["effective_domain_size"]
    )

    df["fold_enrichment"] = (
        numerator
        /
        denominator
    )

    return df


def build_plot(
    results: pd.DataFrame,
    analysis_group: str,
    source: str,
    output_file: Path,
) -> None:
    """
    Plot the top enriched terms ranked by FDR.
    """

    if results.empty:
        return

    subset = results[
        results["source"] == source
    ].copy()

    if subset.empty:
        print(
            f"[INFO] No {source} terms for {analysis_group}; "
            "skipping plot."
        )
        return

    subset = subset.sort_values(
        by=[
            "p_value",
            "fold_enrichment",
        ],
        ascending=[
            True,
            False,
        ],
    ).head(TOP_N_PLOT)

    subset = subset.sort_values(
        "neg_log10_FDR",
        ascending=True,
    )

    labels = subset["name"].astype(str)

    values = subset["neg_log10_FDR"]

    plt.figure(
        figsize=(10, max(6, len(subset) * 0.38))
    )

    plt.barh(
        labels,
        values,
    )

    plt.xlabel(
        r"$-\log_{10}(\mathrm{FDR})$"
    )

    plt.ylabel("")

    title_source = (
        "GO Biological Process"
        if source == "GO:BP"
        else "Reactome"
    )

    plt.title(
        f"{analysis_group} — {title_source}"
    )

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


def format_excel(path: Path) -> None:
    """
    Apply simple consistent workbook formatting.
    """

    wb = load_workbook(path)

    header_fill = PatternFill(
        fill_type="solid",
        fgColor="D9E1F2",
    )

    green_fill = PatternFill(
        fill_type="solid",
        fgColor="E2F0D9",
    )

    red_fill = PatternFill(
        fill_type="solid",
        fgColor="F4CCCC",
    )

    for ws in wb.worksheets:

        ws.freeze_panes = "A2"

        for cell in ws[1]:
            cell.fill = header_fill
            cell.font = Font(
                bold=True
            )
            cell.alignment = Alignment(
                horizontal="center",
                vertical="center",
            )

        if ws.max_row > 1:
            ws.auto_filter.ref = ws.dimensions

        for column_cells in ws.columns:

            letter = get_column_letter(
                column_cells[0].column
            )

            max_len = 0

            for cell in column_cells[:150]:
                if cell.value is None:
                    continue

                max_len = max(
                    max_len,
                    len(str(cell.value)),
                )

            width = min(
                max(max_len + 2, 12),
                55,
            )

            ws.column_dimensions[
                letter
            ].width = width

        # Highlight UP/DOWN summary sheets slightly.
        if ws.title.startswith("UP"):
            for cell in ws[1]:
                cell.fill = green_fill
                cell.font = Font(bold=True)

        elif ws.title.startswith("DOWN"):
            for cell in ws[1]:
                cell.fill = red_fill
                cell.font = Font(bold=True)

    wb.save(path)


# ============================================================
# Validate input
# ============================================================

print()
print("=" * 72)
print("PHASE F1 — GLOBAL FONG ENRICHMENT")
print("=" * 72)

if not INPUT_FILE.exists():
    raise FileNotFoundError(
        "\nF0 output file not found:\n"
        f"{INPUT_FILE}\n\n"
        "Run F0 first."
    )

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

FIGURE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Load strict Fong DEG master
# ============================================================

print(
    f"\n[INFO] Reading:\n"
    f"       {INPUT_FILE}"
)

fong = pd.read_excel(
    INPUT_FILE,
    sheet_name=INPUT_SHEET,
)

required_columns = {
    "Gene",
    "log2FoldChange",
    "padj",
    "Strict_pass",
}

missing = required_columns.difference(
    fong.columns
)

if missing:
    raise ValueError(
        f"Missing required columns: {sorted(missing)}"
    )

fong["Gene"] = (
    fong["Gene"]
    .map(normalize_gene)
)

fong = fong[
    fong["Gene"].notna()
].copy()

# Excel may import TRUE/FALSE as bool,
# but this also safely handles strings.
strict_mask = (
    fong["Strict_pass"]
    .astype(str)
    .str.upper()
    .isin(
        {
            "TRUE",
            "1",
            "YES",
        }
    )
)

strict = fong[
    strict_mask
].copy()

if strict.empty:
    raise RuntimeError(
        "No strict Fong DEG rows were detected."
    )


# ============================================================
# Construct query groups
# ============================================================

all_genes = sorted(
    strict["Gene"].unique()
)

up_genes = sorted(
    strict.loc[
        strict["log2FoldChange"] > 0,
        "Gene",
    ].unique()
)

down_genes = sorted(
    strict.loc[
        strict["log2FoldChange"] < 0,
        "Gene",
    ].unique()
)

queries = {
    "ALL": all_genes,
    "UP": up_genes,
    "DOWN": down_genes,
}

print()
print(
    f"[INFO] ALL genes  : {len(all_genes):,}"
)
print(
    f"[INFO] UP genes   : {len(up_genes):,}"
)
print(
    f"[INFO] DOWN genes : {len(down_genes):,}"
)


# ============================================================
# Record g:Profiler data versions
# ============================================================

print(
    "\n[INFO] Retrieving g:Profiler data versions..."
)

try:
    version_response = requests.get(
        GPROFILER_VERSION_URL,
        params={
            "organism": ORGANISM,
        },
        timeout=REQUEST_TIMEOUT,
        headers={
            "User-Agent": "BIOL_295_FongDownstream/1.0"
        },
    )

    version_response.raise_for_status()

    data_versions = (
        version_response.json()
    )

except Exception as exc:
    print(
        f"[WARNING] Could not retrieve data versions: {exc}"
    )
    data_versions = {}


# ============================================================
# Run enrichment
# ============================================================

all_results = {}
raw_responses = {}

for query_name, genes in queries.items():

    results, raw = query_gprofiler(
        genes,
        query_name,
    )

    results = add_derived_statistics(
        results,
        query_size_fallback=len(genes),
    )

    if not results.empty:

        results["contributing_ensembl_genes"] = (
            results.apply(
                lambda row: extract_intersection_genes(
                    row,
                    raw,
                ),
                axis=1,
            )
        )

        preferred_columns = [
            "analysis_group",
            "source",
            "native",
            "name",
            "p_value",
            "neg_log10_FDR",
            "fold_enrichment",
            "intersection_size",
            "query_size",
            "term_size",
            "effective_domain_size",
            "precision",
            "recall",
            "contributing_ensembl_genes",
            "description",
            "parents",
        ]

        existing = [
            c
            for c in preferred_columns
            if c in results.columns
        ]

        remaining = [
            c
            for c in results.columns
            if c not in existing
        ]

        results = results[
            existing + remaining
        ]

        results = results.sort_values(
            by=[
                "source",
                "p_value",
                "fold_enrichment",
            ],
            ascending=[
                True,
                True,
                False,
            ],
        )

    all_results[
        query_name
    ] = results

    raw_responses[
        query_name
    ] = raw


# ============================================================
# Export CSV files
# ============================================================

for query_name, results in all_results.items():

    csv_path = (
        OUTPUT_DIR
        / f"{query_name.lower()}_enrichment.csv"
    )

    results.to_csv(
        csv_path,
        index=False,
    )

    print(
        f"[INFO] Wrote {csv_path.name}: "
        f"{len(results):,} terms"
    )


# ============================================================
# Summary
# ============================================================

summary_rows = []

for query_name, genes in queries.items():

    results = all_results[
        query_name
    ]

    go_count = 0
    reac_count = 0

    if not results.empty:

        go_count = int(
            (
                results["source"]
                == "GO:BP"
            ).sum()
        )

        reac_count = int(
            (
                results["source"]
                == "REAC"
            ).sum()
        )

    summary_rows.append(
        {
            "Analysis": query_name,
            "Input genes": len(genes),
            "Significant GO:BP terms": go_count,
            "Significant Reactome terms": reac_count,
            "Background": "annotated mouse genes",
            "Correction": "FDR",
            "Threshold": FDR_THRESHOLD,
        }
    )

summary = pd.DataFrame(
    summary_rows
)


# ============================================================
# Export Excel workbook
# ============================================================

with pd.ExcelWriter(
    OUTPUT_XLSX,
    engine="openpyxl",
) as writer:

    summary.to_excel(
        writer,
        sheet_name="Summary",
        index=False,
    )

    for query_name, results in all_results.items():

        if results.empty:
            pd.DataFrame(
                {
                    "message": [
                        "No significant enrichment terms."
                    ]
                }
            ).to_excel(
                writer,
                sheet_name=query_name,
                index=False,
            )

            continue

        go = results[
            results["source"]
            == "GO:BP"
        ].copy()

        reac = results[
            results["source"]
            == "REAC"
        ].copy()

        results.to_excel(
            writer,
            sheet_name=query_name,
            index=False,
        )

        go.to_excel(
            writer,
            sheet_name=f"{query_name}_GO_BP",
            index=False,
        )

        reac.to_excel(
            writer,
            sheet_name=f"{query_name}_REAC",
            index=False,
        )


format_excel(
    OUTPUT_XLSX
)


# ============================================================
# Figures
# ============================================================

for query_name, results in all_results.items():

    for source in SOURCES:

        source_label = (
            "GO_BP"
            if source == "GO:BP"
            else "REAC"
        )

        output_plot = (
            FIGURE_DIR
            / f"{query_name}_{source_label}_top20.png"
        )

        build_plot(
            results=results,
            analysis_group=query_name,
            source=source,
            output_file=output_plot,
        )


# ============================================================
# Metadata
# ============================================================

metadata = {
    "phase": "F1_global_enrichment",
    "organism": ORGANISM,
    "sources": SOURCES,
    "multiple_testing_correction": "FDR",
    "significance_threshold": FDR_THRESHOLD,
    "domain_scope": "annotated",
    "background_note": (
        "Annotated mouse genes were used as the enrichment background "
        "because the currently available Fong input contains only "
        "already-filtered significant genes. Rerun with a custom "
        "tested-gene background if the complete unfiltered DE table "
        "becomes available."
    ),
    "input_counts": {
        "ALL": len(all_genes),
        "UP": len(up_genes),
        "DOWN": len(down_genes),
    },
    "gprofiler_data_versions": data_versions,
}

with open(
    OUTPUT_METADATA,
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        metadata,
        handle,
        indent=2,
        ensure_ascii=False,
    )


# ============================================================
# Console report
# ============================================================

print()
print("=" * 72)
print("PHASE F1 COMPLETE")
print("=" * 72)

for query_name, genes in queries.items():

    results = all_results[
        query_name
    ]

    n_total = len(results)

    n_go = (
        int(
            (
                results["source"]
                == "GO:BP"
            ).sum()
        )
        if not results.empty
        else 0
    )

    n_reac = (
        int(
            (
                results["source"]
                == "REAC"
            ).sum()
        )
        if not results.empty
        else 0
    )

    print(
        f"{query_name:<5} "
        f"genes={len(genes):>5,} | "
        f"terms={n_total:>5,} | "
        f"GO:BP={n_go:>5,} | "
        f"REAC={n_reac:>5,}"
    )

print()
print(
    "Background:"
)
print(
    "    Annotated Mus musculus gene universe"
)

print()
print(
    "IMPORTANT:"
)
print(
    "    This is an exploratory global enrichment background."
)
print(
    "    If the full Fong tested-gene universe is recovered, "
    "rerun F1 using it as a custom background."
)

print()
print(
    "Outputs:"
)

print(
    f"    {OUTPUT_DIR.relative_to(PROJECT_ROOT)}"
)

print()
print(
    "[NEXT] Inspect the dominant ALL/UP/DOWN programs before "
    "constructing the F2 cell-death/survival master."
)