#!/usr/bin/env python3

"""
Phase F2 — Build the Fong cell-fate / stress master
===================================================

Purpose
-------
Extract a biologically defined subset of the strict Fong DEG universe
that captures:

CORE CELL-FATE BIOLOGY
    - general cell death / programmed cell death
    - apoptosis and apoptotic signaling
    - positive regulation of death
    - negative regulation of death / survival
    - necroptosis / programmed necrotic death
    - autophagy / mitophagy

STRESS-TO-DEATH INTERFACE
    - DNA damage / checkpoint / DNA repair
    - ER stress / oxidative stress

The selection is NOT based on:
    - gene-name keywords
    - manual gene selection
    - whether a term was significant in F1

Instead, predefined GO Biological Process terms are queried through
g:Profiler g:Convert, which returns mouse genes associated with each term.
Those genes are intersected with the 1,891 strict Fong DEGs.

This creates a reproducible ontology-defined candidate universe for
downstream network analysis.

Input
-----
Data/Processed/FongDownstream/F0_Audit/
    Fong_downstream_phase_F0_audit.xlsx

Data/Processed/FongDownstream/F1_Enrichment/
    all_enrichment.csv
    up_enrichment.csv
    down_enrichment.csv

The F1 enrichment files are used ONLY for annotation with significant
Reactome pathways; they do not determine entry into the F2 master.

Outputs
-------
Data/Processed/FongDownstream/F2_CellFate/

    F2_cell_fate_master.xlsx
    cell_fate_master.csv
    core_cell_fate.csv
    stress_interface.csv
    category_summary.csv
    GO_term_audit.csv
    reactome_annotations.csv

    figures/
        F2_category_counts.png
        F2_core_vs_interface.png

Run
---
From Project/:

    py .\\Scripts\\F2_build_cell_fate_master.py
"""

from __future__ import annotations

from pathlib import Path
from collections import defaultdict
from typing import Any
import json
import time

import pandas as pd
import requests
import matplotlib.pyplot as plt

from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


# ============================================================
# Project paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

F0_FILE = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F0_Audit"
    / "Fong_downstream_phase_F0_audit.xlsx"
)

F0_SHEET = "Fong strict DEG master"

F1_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F1_Enrichment"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F2_CellFate"
)

FIGURE_DIR = OUTPUT_DIR / "figures"

OUTPUT_XLSX = (
    OUTPUT_DIR
    / "F2_cell_fate_master.xlsx"
)

OUTPUT_MASTER_CSV = (
    OUTPUT_DIR
    / "cell_fate_master.csv"
)

OUTPUT_CORE_CSV = (
    OUTPUT_DIR
    / "core_cell_fate.csv"
)

OUTPUT_INTERFACE_CSV = (
    OUTPUT_DIR
    / "stress_interface.csv"
)

OUTPUT_CATEGORY_CSV = (
    OUTPUT_DIR
    / "category_summary.csv"
)

OUTPUT_GO_AUDIT_CSV = (
    OUTPUT_DIR
    / "GO_term_audit.csv"
)

OUTPUT_REACTOME_CSV = (
    OUTPUT_DIR
    / "reactome_annotations.csv"
)

OUTPUT_METADATA = (
    OUTPUT_DIR
    / "F2_metadata.json"
)


# ============================================================
# g:Profiler
# ============================================================

GCONVERT_URL = (
    "https://biit.cs.ut.ee/gprofiler/api/convert/convert/"
)

ORGANISM = "mmusculus"

REQUEST_TIMEOUT = 120


# ============================================================
# Predefined ontology panel
# ============================================================

# IMPORTANT:
#
# These terms define the biological scope of F2.
#
# Gene entry into F2 is based on association with one or more of
# these GO terms.
#
# The terms are grouped conceptually so that each gene retains
# its provenance and may belong to several categories.


GO_PANEL = [

    # --------------------------------------------------------
    # CORE: general cell death
    # --------------------------------------------------------

    {
        "category": "Cell_death_general",
        "tier": "CORE",
        "go_id": "GO:0008219",
        "label": "cell death",
    },

    {
        "category": "Cell_death_general",
        "tier": "CORE",
        "go_id": "GO:0012501",
        "label": "programmed cell death",
    },

    # --------------------------------------------------------
    # CORE: apoptosis / apoptotic signaling
    # --------------------------------------------------------

    {
        "category": "Apoptosis_execution",
        "tier": "CORE",
        "go_id": "GO:0006915",
        "label": "apoptotic process",
    },

    {
        "category": "Apoptosis_execution",
        "tier": "CORE",
        "go_id": "GO:0097190",
        "label": "apoptotic signaling pathway",
    },

    {
        "category": "Apoptosis_execution",
        "tier": "CORE",
        "go_id": "GO:0097191",
        "label": "extrinsic apoptotic signaling pathway",
    },

    {
        "category": "Apoptosis_execution",
        "tier": "CORE",
        "go_id": "GO:0097193",
        "label": "intrinsic apoptotic signaling pathway",
    },

    {
        "category": "Apoptosis_execution",
        "tier": "CORE",
        "go_id": "GO:0097194",
        "label": "execution phase of apoptosis",
    },

    # --------------------------------------------------------
    # CORE: positive regulation of death
    # --------------------------------------------------------

    {
        "category": "Pro_death_regulation",
        "tier": "CORE",
        "go_id": "GO:0043065",
        "label": "positive regulation of apoptotic process",
    },

    {
        "category": "Pro_death_regulation",
        "tier": "CORE",
        "go_id": "GO:0043068",
        "label": "positive regulation of programmed cell death",
    },

    # --------------------------------------------------------
    # CORE: survival / anti-death
    # --------------------------------------------------------

    {
        "category": "Survival_anti_death",
        "tier": "CORE",
        "go_id": "GO:0043066",
        "label": "negative regulation of apoptotic process",
    },

    {
        "category": "Survival_anti_death",
        "tier": "CORE",
        "go_id": "GO:0043069",
        "label": "negative regulation of programmed cell death",
    },

    {
        "category": "Survival_anti_death",
        "tier": "CORE",
        "go_id": "GO:0060548",
        "label": "negative regulation of cell death",
    },

    # --------------------------------------------------------
    # CORE: necroptosis / programmed necrotic death
    # --------------------------------------------------------

    {
        "category": "Necroptosis",
        "tier": "CORE",
        "go_id": "GO:0070266",
        "label": "necroptotic process",
    },

    {
        "category": "Necroptosis",
        "tier": "CORE",
        "go_id": "GO:0097527",
        "label": "necroptotic signaling pathway",
    },

    {
        "category": "Necroptosis",
        "tier": "CORE",
        "go_id": "GO:0097300",
        "label": "programmed necrotic cell death",
    },

    # --------------------------------------------------------
    # CORE: autophagy / mitophagy
    # --------------------------------------------------------

    {
        "category": "Autophagy_mitophagy",
        "tier": "CORE",
        "go_id": "GO:0061919",
        "label": "process utilizing autophagic mechanism",
    },

    {
        "category": "Autophagy_mitophagy",
        "tier": "CORE",
        "go_id": "GO:0010506",
        "label": "regulation of autophagy",
    },

    {
        "category": "Autophagy_mitophagy",
        "tier": "CORE",
        "go_id": "GO:0000422",
        "label": "autophagy of mitochondrion",
    },

    # --------------------------------------------------------
    # INTERFACE: DNA damage / checkpoint / repair
    # --------------------------------------------------------

    {
        "category": "DNA_damage_checkpoint_repair",
        "tier": "INTERFACE",
        "go_id": "GO:0006974",
        "label": "DNA damage response",
    },

    {
        "category": "DNA_damage_checkpoint_repair",
        "tier": "INTERFACE",
        "go_id": "GO:0000077",
        "label": "DNA damage checkpoint signaling",
    },

    {
        "category": "DNA_damage_checkpoint_repair",
        "tier": "INTERFACE",
        "go_id": "GO:0006281",
        "label": "DNA repair",
    },

    # --------------------------------------------------------
    # INTERFACE: ER / oxidative stress
    # --------------------------------------------------------

    {
        "category": "ER_oxidative_stress",
        "tier": "INTERFACE",
        "go_id": "GO:0034976",
        "label": "response to endoplasmic reticulum stress",
    },

    {
        "category": "ER_oxidative_stress",
        "tier": "INTERFACE",
        "go_id": "GO:0006979",
        "label": "response to oxidative stress",
    },
]


CORE_CATEGORIES = {
    x["category"]
    for x in GO_PANEL
    if x["tier"] == "CORE"
}

INTERFACE_CATEGORIES = {
    x["category"]
    for x in GO_PANEL
    if x["tier"] == "INTERFACE"
}


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


def normalize_ensembl(value: Any) -> str | None:
    if pd.isna(value):
        return None

    value = str(value).strip()

    if not value:
        return None

    # Remove version number if present:
    #
    # ENSMUSG00000000001.12
    # ->
    # ENSMUSG00000000001

    return value.split(".")[0]


def gconvert_term_to_ensembl(
    go_id: str,
) -> tuple[set[str], list[dict]]:
    """
    Retrieve all mouse genes associated with one GO term.

    g:Profiler g:Convert accepts functional term identifiers and
    returns genes associated with the term.

    Returns
    -------
    ensembl_ids
        Set of normalized Ensembl mouse gene IDs.

    raw_result
        Raw g:Convert result list for audit purposes.
    """

    payload = {
        "organism": ORGANISM,
        "query": [go_id],
        "target": "ENSG",
    }

    response = requests.post(
        GCONVERT_URL,
        json=payload,
        timeout=REQUEST_TIMEOUT,
        headers={
            "User-Agent":
                "BIOL_295_FongDownstream_F2/1.0"
        },
    )

    response.raise_for_status()

    data = response.json()

    result = data.get(
        "result",
        [],
    )

    ensembl_ids = set()

    for item in result:

        converted = item.get(
            "converted"
        )

        if converted is None:
            continue

        converted = normalize_ensembl(
            converted
        )

        if (
            converted
            and converted.upper() != "N/A"
        ):
            ensembl_ids.add(
                converted
            )

    return ensembl_ids, result


def parse_semicolon_genes(
    value: Any,
) -> set[str]:
    """
    Parse semicolon-separated Ensembl IDs from F1 output.
    """

    if pd.isna(value):
        return set()

    value = str(value).strip()

    if not value:
        return set()

    genes = set()

    for item in value.split(";"):

        item = normalize_ensembl(
            item
        )

        if item:
            genes.add(
                item
            )

    return genes


def format_workbook(
    path: Path,
) -> None:
    """
    Apply lightweight workbook formatting.
    """

    wb = load_workbook(
        path
    )

    header_fill = PatternFill(
        fill_type="solid",
        fgColor="D9E1F2",
    )

    core_fill = PatternFill(
        fill_type="solid",
        fgColor="E2F0D9",
    )

    interface_fill = PatternFill(
        fill_type="solid",
        fgColor="FFF2CC",
    )

    up_fill = PatternFill(
        fill_type="solid",
        fgColor="E2F0D9",
    )

    down_fill = PatternFill(
        fill_type="solid",
        fgColor="F4CCCC",
    )

    old_fill = PatternFill(
        fill_type="solid",
        fgColor="D9EAD3",
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
            ws.auto_filter.ref = (
                ws.dimensions
            )

        for column_cells in ws.columns:

            letter = get_column_letter(
                column_cells[0].column
            )

            max_length = 0

            for cell in column_cells[:150]:

                if cell.value is None:
                    continue

                max_length = max(
                    max_length,
                    len(str(cell.value)),
                )

            width = min(
                max(
                    max_length + 2,
                    12,
                ),
                55,
            )

            ws.column_dimensions[
                letter
            ].width = width

    if "Cell fate master" in wb.sheetnames:

        ws = wb[
            "Cell fate master"
        ]

        header_lookup = {
            cell.value: cell.column
            for cell in ws[1]
        }

        tier_col = header_lookup.get(
            "F2_tier"
        )

        dir_col = header_lookup.get(
            "Direction"
        )

        old_col = header_lookup.get(
            "Old_126_gene_master"
        )

        for row in range(
            2,
            ws.max_row + 1,
        ):

            if tier_col:

                value = ws.cell(
                    row=row,
                    column=tier_col,
                ).value

                if value == "CORE":

                    ws.cell(
                        row=row,
                        column=tier_col,
                    ).fill = core_fill

                elif value == "INTERFACE":

                    ws.cell(
                        row=row,
                        column=tier_col,
                    ).fill = interface_fill

            if dir_col:

                value = ws.cell(
                    row=row,
                    column=dir_col,
                ).value

                if value == "Up":

                    ws.cell(
                        row=row,
                        column=dir_col,
                    ).fill = up_fill

                elif value == "Down":

                    ws.cell(
                        row=row,
                        column=dir_col,
                    ).fill = down_fill

            if old_col:

                value = ws.cell(
                    row=row,
                    column=old_col,
                ).value

                if str(value).upper() in {
                    "TRUE",
                    "YES",
                    "1",
                }:

                    ws.cell(
                        row=row,
                        column=old_col,
                    ).fill = old_fill

    wb.save(
        path
    )


# ============================================================
# Validate inputs
# ============================================================

print()
print("=" * 76)
print("PHASE F2 — BUILD FONG CELL-FATE / STRESS MASTER")
print("=" * 76)

if not F0_FILE.exists():

    raise FileNotFoundError(
        "\nF0 workbook not found:\n"
        f"{F0_FILE}"
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
# Load strict Fong master
# ============================================================

print(
    f"\n[INFO] Reading F0 strict master:\n"
    f"       {F0_FILE}"
)

fong = pd.read_excel(
    F0_FILE,
    sheet_name=F0_SHEET,
)

required_columns = {
    "Gene",
    "Row.names",
    "log2FoldChange",
    "padj",
    "Direction",
}

missing_columns = (
    required_columns
    .difference(
        fong.columns
    )
)

if missing_columns:

    raise ValueError(
        "Missing required F0 columns: "
        f"{sorted(missing_columns)}"
    )


fong["Gene"] = (
    fong["Gene"]
    .map(
        normalize_gene
    )
)

fong["Ensembl_ID"] = (
    fong["Row.names"]
    .map(
        normalize_ensembl
    )
)

fong = fong[
    fong["Gene"].notna()
    &
    fong["Ensembl_ID"].notna()
].copy()


# One row per strict DEG was confirmed in F0.

fong_by_ensembl = (
    fong
    .set_index(
        "Ensembl_ID",
        drop=False,
    )
)

fong_ensembl_set = set(
    fong[
        "Ensembl_ID"
    ]
)

print(
    f"[INFO] Strict Fong genes available: "
    f"{len(fong):,}"
)


# ============================================================
# Query GO terms
# ============================================================

gene_categories = defaultdict(
    set
)

gene_go_terms = defaultdict(
    set
)

gene_go_labels = defaultdict(
    set
)

term_audit_rows = []

raw_term_metadata = {}


print()
print(
    "[INFO] Querying predefined GO cell-fate/stress panel..."
)

for i, term in enumerate(
    GO_PANEL,
    start=1,
):

    category = term[
        "category"
    ]

    tier = term[
        "tier"
    ]

    go_id = term[
        "go_id"
    ]

    label = term[
        "label"
    ]

    print(
        f"       [{i:02d}/{len(GO_PANEL):02d}] "
        f"{go_id}  {label}"
    )

    try:

        annotated_genes, raw_result = (
            gconvert_term_to_ensembl(
                go_id
            )
        )

    except Exception as exc:

        print(
            f"       [WARNING] Query failed: {exc}"
        )

        annotated_genes = set()
        raw_result = []

    fong_hits = (
        annotated_genes
        &
        fong_ensembl_set
    )

    for ensembl_id in fong_hits:

        gene_categories[
            ensembl_id
        ].add(
            category
        )

        gene_go_terms[
            ensembl_id
        ].add(
            go_id
        )

        gene_go_labels[
            ensembl_id
        ].add(
            label
        )

    term_audit_rows.append(
        {
            "Category": category,
            "Tier": tier,
            "GO_ID": go_id,
            "GO_label": label,
            "Mouse_genes_in_GO_term":
                len(annotated_genes),
            "Strict_Fong_hits":
                len(fong_hits),
        }
    )

    raw_term_metadata[
        go_id
    ] = {
        "label": label,
        "category": category,
        "tier": tier,
        "mouse_gene_count":
            len(annotated_genes),
        "strict_fong_hits":
            len(fong_hits),
    }

    # Be kind to the public endpoint.
    time.sleep(
        0.15
    )


go_audit = pd.DataFrame(
    term_audit_rows
)


# ============================================================
# Construct F2 master
# ============================================================

selected_ensembl = set(
    gene_categories.keys()
)

print()
print(
    f"[INFO] Fong DEGs matching ≥1 predefined GO term: "
    f"{len(selected_ensembl):,}"
)


master = fong[
    fong["Ensembl_ID"]
    .isin(
        selected_ensembl
    )
].copy()


# ============================================================
# Category flags
# ============================================================

all_categories = sorted(
    CORE_CATEGORIES
    |
    INTERFACE_CATEGORIES
)

for category in all_categories:

    master[
        category
    ] = (
        master[
            "Ensembl_ID"
        ]
        .map(
            lambda x:
                category
                in gene_categories[
                    x
                ]
        )
    )


def category_string(
    ensembl_id: str,
) -> str:

    return "; ".join(
        sorted(
            gene_categories[
                ensembl_id
            ]
        )
    )


def go_id_string(
    ensembl_id: str,
) -> str:

    return "; ".join(
        sorted(
            gene_go_terms[
                ensembl_id
            ]
        )
    )


def go_label_string(
    ensembl_id: str,
) -> str:

    return "; ".join(
        sorted(
            gene_go_labels[
                ensembl_id
            ]
        )
    )


master[
    "F2_categories"
] = (
    master[
        "Ensembl_ID"
    ]
    .map(
        category_string
    )
)

master[
    "F2_GO_terms"
] = (
    master[
        "Ensembl_ID"
    ]
    .map(
        go_id_string
    )
)

master[
    "F2_GO_labels"
] = (
    master[
        "Ensembl_ID"
    ]
    .map(
        go_label_string
    )
)

master[
    "F2_category_count"
] = (
    master[
        "Ensembl_ID"
    ]
    .map(
        lambda x:
            len(
                gene_categories[
                    x
                ]
            )
    )
)


# ============================================================
# CORE vs INTERFACE
# ============================================================

def determine_tier(
    ensembl_id: str,
) -> str:

    categories = gene_categories[
        ensembl_id
    ]

    if categories & CORE_CATEGORIES:
        return "CORE"

    return "INTERFACE"


master[
    "F2_tier"
] = (
    master[
        "Ensembl_ID"
    ]
    .map(
        determine_tier
    )
)


# ============================================================
# Convenience biological flags
# ============================================================

master[
    "Any_death_apoptosis"
] = (
    master[
        [
            "Cell_death_general",
            "Apoptosis_execution",
            "Pro_death_regulation",
            "Necroptosis",
        ]
    ]
    .any(
        axis=1
    )
)

master[
    "Any_survival_anti_death"
] = (
    master[
        "Survival_anti_death"
    ]
)

master[
    "Any_autophagy_mitophagy"
] = (
    master[
        "Autophagy_mitophagy"
    ]
)

master[
    "Any_stress_interface"
] = (
    master[
        [
            "DNA_damage_checkpoint_repair",
            "ER_oxidative_stress",
        ]
    ]
    .any(
        axis=1
    )
)


# ============================================================
# F1 significant Reactome annotations
# ============================================================

reactome_by_gene = {
    "ALL": defaultdict(
        set
    ),
    "UP": defaultdict(
        set
    ),
    "DOWN": defaultdict(
        set
    ),
}

reactome_rows = []


for group in [
    "ALL",
    "UP",
    "DOWN",
]:

    file_path = (
        F1_DIR
        / f"{group.lower()}_enrichment.csv"
    )

    if not file_path.exists():

        print(
            f"[WARNING] Missing F1 file: "
            f"{file_path.name}"
        )

        continue

    enrichment = pd.read_csv(
        file_path
    )

    if (
        "source"
        not in enrichment.columns
        or
        "contributing_ensembl_genes"
        not in enrichment.columns
    ):
        continue

    reactome = enrichment[
        enrichment["source"]
        == "REAC"
    ].copy()

    for _, row in reactome.iterrows():

        ensembl_hits = (
            parse_semicolon_genes(
                row[
                    "contributing_ensembl_genes"
                ]
            )
        )

        ensembl_hits &= (
            selected_ensembl
        )

        if not ensembl_hits:
            continue

        pathway = (
            f"{row.get('native', '')} | "
            f"{row.get('name', '')}"
        )

        for ensembl_id in ensembl_hits:

            reactome_by_gene[
                group
            ][
                ensembl_id
            ].add(
                pathway
            )

        reactome_rows.append(
            {
                "Analysis": group,
                "Reactome_ID":
                    row.get(
                        "native"
                    ),
                "Reactome_name":
                    row.get(
                        "name"
                    ),
                "FDR":
                    row.get(
                        "p_value"
                    ),
                "F2_intersection_size":
                    len(
                        ensembl_hits
                    ),
                "F2_Ensembl_genes":
                    ";".join(
                        sorted(
                            ensembl_hits
                        )
                    ),
            }
        )


for group in [
    "ALL",
    "UP",
    "DOWN",
]:

    column = (
        f"F1_significant_REAC_{group}"
    )

    master[
        column
    ] = (
        master[
            "Ensembl_ID"
        ]
        .map(
            lambda x:
                "; ".join(
                    sorted(
                        reactome_by_gene[
                            group
                        ][x]
                    )
                )
        )
    )


reactome_annotations = pd.DataFrame(
    reactome_rows
)


# ============================================================
# Sort master
# ============================================================

tier_order = {
    "CORE": 0,
    "INTERFACE": 1,
}

master[
    "_tier_order"
] = (
    master[
        "F2_tier"
    ]
    .map(
        tier_order
    )
)

master[
    "_abs_lfc"
] = (
    master[
        "log2FoldChange"
    ]
    .abs()
)

master = master.sort_values(
    by=[
        "_tier_order",
        "F2_category_count",
        "_abs_lfc",
    ],
    ascending=[
        True,
        False,
        False,
    ],
)

master = master.drop(
    columns=[
        "_tier_order",
        "_abs_lfc",
    ]
)


# ============================================================
# Column ordering
# ============================================================

preferred_front = [
    "Gene",
    "Ensembl_ID",
    "baseMean",
    "log2FoldChange",
    "padj",
    "Direction",
    "F2_tier",
    "F2_category_count",
    "F2_categories",
    "F2_GO_terms",
    "F2_GO_labels",
    "Any_death_apoptosis",
    "Any_survival_anti_death",
    "Any_autophagy_mitophagy",
    "Any_stress_interface",
    "Old_126_gene_master",
]

category_columns = (
    all_categories
)

reactome_columns = [
    "F1_significant_REAC_ALL",
    "F1_significant_REAC_UP",
    "F1_significant_REAC_DOWN",
]

front = [
    c
    for c in (
        preferred_front
        +
        category_columns
        +
        reactome_columns
    )
    if c in master.columns
]

remaining = [
    c
    for c in master.columns
    if c not in front
    and c != "Row.names"
]

master = master[
    front
    +
    remaining
]


# ============================================================
# Core/interface subsets
# ============================================================

core = master[
    master[
        "F2_tier"
    ]
    == "CORE"
].copy()

interface = master[
    master[
        "F2_tier"
    ]
    == "INTERFACE"
].copy()


# ============================================================
# Category summary
# ============================================================

category_summary_rows = []

for category in all_categories:

    subset = master[
        master[
            category
        ]
    ]

    tier = (
        "CORE"
        if category in CORE_CATEGORIES
        else "INTERFACE"
    )

    category_summary_rows.append(
        {
            "Category":
                category,
            "Tier":
                tier,
            "Genes":
                len(
                    subset
                ),
            "Up":
                int(
                    (
                        subset[
                            "Direction"
                        ]
                        == "Up"
                    )
                    .sum()
                ),
            "Down":
                int(
                    (
                        subset[
                            "Direction"
                        ]
                        == "Down"
                    )
                    .sum()
                ),
            "Old_master_overlap":
                int(
                    subset[
                        "Old_126_gene_master"
                    ]
                    .astype(str)
                    .str.upper()
                    .isin(
                        {
                            "TRUE",
                            "YES",
                            "1",
                        }
                    )
                    .sum()
                )
                if "Old_126_gene_master"
                in subset.columns
                else 0,
        }
    )


category_summary = pd.DataFrame(
    category_summary_rows
)

category_summary = (
    category_summary
    .sort_values(
        by=[
            "Tier",
            "Genes",
        ],
        ascending=[
            True,
            False,
        ],
    )
)


# ============================================================
# Overall summary
# ============================================================

summary = pd.DataFrame(
    [
        {
            "Metric":
                "Strict Fong DEG universe",
            "Value":
                len(
                    fong
                ),
        },
        {
            "Metric":
                "F2 selected genes",
            "Value":
                len(
                    master
                ),
        },
        {
            "Metric":
                "CORE cell-fate genes",
            "Value":
                len(
                    core
                ),
        },
        {
            "Metric":
                "INTERFACE-only genes",
            "Value":
                len(
                    interface
                ),
        },
        {
            "Metric":
                "F2 upregulated genes",
            "Value":
                int(
                    (
                        master[
                            "Direction"
                        ]
                        == "Up"
                    )
                    .sum()
                ),
        },
        {
            "Metric":
                "F2 downregulated genes",
            "Value":
                int(
                    (
                        master[
                            "Direction"
                        ]
                        == "Down"
                    )
                    .sum()
                ),
        },
        {
            "Metric":
                "F2 genes overlapping old 126-gene master",
            "Value":
                int(
                    master[
                        "Old_126_gene_master"
                    ]
                    .astype(str)
                    .str.upper()
                    .isin(
                        {
                            "TRUE",
                            "YES",
                            "1",
                        }
                    )
                    .sum()
                )
                if "Old_126_gene_master"
                in master.columns
                else 0,
        },
        {
            "Metric":
                "Selection rule",
            "Value":
                (
                    "Strict Fong DEG annotated to >=1 "
                    "predefined GO cell-fate/stress term."
                ),
        },
        {
            "Metric":
                "CORE definition",
            "Value":
                (
                    "Death/apoptosis, pro-death regulation, "
                    "survival/anti-death, necroptosis, "
                    "autophagy/mitophagy."
                ),
        },
        {
            "Metric":
                "INTERFACE definition",
            "Value":
                (
                    "DNA-damage/checkpoint/repair or "
                    "ER/oxidative-stress annotation, "
                    "without CORE annotation."
                ),
        },
    ]
)


# ============================================================
# Save CSVs
# ============================================================

master.to_csv(
    OUTPUT_MASTER_CSV,
    index=False,
)

core.to_csv(
    OUTPUT_CORE_CSV,
    index=False,
)

interface.to_csv(
    OUTPUT_INTERFACE_CSV,
    index=False,
)

category_summary.to_csv(
    OUTPUT_CATEGORY_CSV,
    index=False,
)

go_audit.to_csv(
    OUTPUT_GO_AUDIT_CSV,
    index=False,
)

reactome_annotations.to_csv(
    OUTPUT_REACTOME_CSV,
    index=False,
)


# ============================================================
# Save Excel workbook
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

    master.to_excel(
        writer,
        sheet_name="Cell fate master",
        index=False,
    )

    core.to_excel(
        writer,
        sheet_name="CORE cell fate",
        index=False,
    )

    interface.to_excel(
        writer,
        sheet_name="Stress interface",
        index=False,
    )

    category_summary.to_excel(
        writer,
        sheet_name="Category summary",
        index=False,
    )

    go_audit.to_excel(
        writer,
        sheet_name="GO term audit",
        index=False,
    )

    reactome_annotations.to_excel(
        writer,
        sheet_name="Reactome annotations",
        index=False,
    )


format_workbook(
    OUTPUT_XLSX
)


# ============================================================
# Figures
# ============================================================

plot_df = (
    category_summary
    .sort_values(
        "Genes",
        ascending=True,
    )
)

plt.figure(
    figsize=(
        10,
        max(
            6,
            len(
                plot_df
            )
            * 0.45,
        ),
    )
)

plt.barh(
    plot_df[
        "Category"
    ],
    plot_df[
        "Genes"
    ],
)

plt.xlabel(
    "Number of strict Fong DEGs"
)

plt.ylabel(
    ""
)

plt.title(
    "F2 ontology-defined cell-fate/stress categories"
)

plt.tight_layout()

plt.savefig(
    FIGURE_DIR
    / "F2_category_counts.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close()


tier_counts = (
    master[
        "F2_tier"
    ]
    .value_counts()
    .reindex(
        [
            "CORE",
            "INTERFACE",
        ],
        fill_value=0,
    )
)

plt.figure(
    figsize=(
        6,
        5,
    )
)

plt.bar(
    tier_counts.index,
    tier_counts.values,
)

plt.ylabel(
    "Number of genes"
)

plt.title(
    "F2 core cell-fate vs stress-interface genes"
)

plt.tight_layout()

plt.savefig(
    FIGURE_DIR
    / "F2_core_vs_interface.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close()


# ============================================================
# Metadata
# ============================================================

metadata = {
    "phase":
        "F2_build_cell_fate_master",

    "organism":
        ORGANISM,

    "strict_fong_gene_count":
        len(
            fong
        ),

    "selected_gene_count":
        len(
            master
        ),

    "core_gene_count":
        len(
            core
        ),

    "interface_gene_count":
        len(
            interface
        ),

    "selection_method":
        (
            "Intersection of strict Fong DEGs with "
            "genes associated to predefined GO terms "
            "retrieved using g:Profiler g:Convert."
        ),

    "go_panel":
        GO_PANEL,

    "note":
        (
            "F1 enrichment significance is not an entry "
            "criterion for F2. F1 Reactome results are "
            "used only as annotations."
        ),
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
print("=" * 76)
print("PHASE F2 COMPLETE")
print("=" * 76)

print(
    f"Strict Fong universe        : "
    f"{len(fong):,}"
)

print(
    f"F2 selected genes           : "
    f"{len(master):,}"
)

print(
    f"CORE cell-fate genes        : "
    f"{len(core):,}"
)

print(
    f"INTERFACE-only genes        : "
    f"{len(interface):,}"
)

print()

print(
    f"Upregulated F2 genes        : "
    f"{(master['Direction'] == 'Up').sum():,}"
)

print(
    f"Downregulated F2 genes      : "
    f"{(master['Direction'] == 'Down').sum():,}"
)

if "Old_126_gene_master" in master.columns:

    old_overlap = (
        master[
            "Old_126_gene_master"
        ]
        .astype(str)
        .str.upper()
        .isin(
            {
                "TRUE",
                "YES",
                "1",
            }
        )
        .sum()
    )

    print(
        f"Old-master overlap          : "
        f"{old_overlap:,}"
    )

print()

print(
    "Category counts:"
)

for _, row in (
    category_summary
    .iterrows()
):

    print(
        f"    "
        f"{row['Category']:<32} "
        f"{int(row['Genes']):>4} "
        f"(up={int(row['Up']):>3}, "
        f"down={int(row['Down']):>3})"
    )

print()

print(
    "Output directory:"
)

print(
    f"    "
    f"{OUTPUT_DIR.relative_to(PROJECT_ROOT)}"
)

print()

print(
    "[NEXT] Inspect F2 size/category composition before "
    "building the STRING network. We do not want to send "
    "an unexpectedly huge or biologically diffuse universe "
    "into F3 without checking it first."
)