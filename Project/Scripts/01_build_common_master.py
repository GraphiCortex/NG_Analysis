#!/usr/bin/env python3
"""
01_build_common_master.py

Reproducible preprocessing for the R-Noel pocket-protein / cell-death project.

INPUT
-----
The manually curated Excel workbook containing these five sheets:
    1. fong vs oshikawa-apoptosis
    2. top distinct targets
    3. fong vs oshikawa-necroptosis
    4. fong vs oshikawa-autophagy
    5. andruisak vs fong&oshikawa-apop

WHAT THIS SCRIPT DOES
---------------------
1. Reads the apoptosis, autophagy, and necroptosis candidate sheets.
2. Reconstructs adult-vs-embryo overlap directly from the binary
   Fong / CAG / MAP2 membership columns.
3. Identifies:
       adult shared = Fong == 1 AND (CAG == 1 OR MAP2 == 1)
4. Classifies expression direction across developmental stages as:
       Concordant / REVERSED / MIXED / Not comparable
5. Adds overlap with the Andrusiak dataset.
6. Carries forward:
       - Role
       - Function
       - Pathway
       - top-distinct membership
       - the doctor's yellow manual highlights
7. Creates two new sheets:
       adult-embryo common
       network master
8. Preserves the visual language of the existing workbook rather than
   creating a new Excel-table/dashboard style.
9. Adds empty network-metric columns to "network master" for the next
   analysis script.

The original five sheets are NOT modified.

DEPENDENCY
----------
    pip install openpyxl

USAGE
-----
    python 01_build_common_master.py "cell death project-sorted genes-all datasets- RS - 6-7-26 NG.xlsx"

Optional:
    python 01_build_common_master.py input.xlsx -o output.xlsx
"""

from __future__ import annotations

import argparse
import re
import pandas as pd
from copy import copy
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from openpyxl import load_workbook
from openpyxl.cell.cell import Cell
from openpyxl.styles import PatternFill, Font, Alignment
from openpyxl.utils import get_column_letter


# ---------------------------------------------------------------------
# Workbook configuration
# ---------------------------------------------------------------------

SHEET_APOPTOSIS = "fong vs oshikawa-apoptosis"
SHEET_TOP_DISTINCT = "top distinct targets"
SHEET_NECROPTOSIS = "fong vs oshikawa-necroptosis"
SHEET_AUTOPHAGY = "fong vs oshikawa-autophagy"
SHEET_ANDRUSIAK = "andruisak vs fong&oshikawa-apop"

OUT_COMMON = "adult-embryo common"
OUT_MASTER = "network master"

MECHANISM_PRIORITY = {
    "Apoptosis": 0,
    "Autophagy": 1,
    "Necroptosis": 2,
}

# Annotation columns differ slightly among the source sheets.
ANNOTATION_COLS = {
    SHEET_APOPTOSIS: (9, 10, 11),     # I:J:K
    SHEET_AUTOPHAGY: (8, 9, 10),      # H:I:J
    SHEET_NECROPTOSIS: (8, 9, 10),    # H:I:J
}

# Common yellow shades used for manual highlighting.
YELLOW_RGB = {
    "FFC000",
    "FFFF00",
    "FFF2CC",
    "FFE699",
    "FFD966",
}

# Used only for newly created metadata cells.
MANUAL_YELLOW = PatternFill(fill_type="solid", fgColor="FFC000")


# ---------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------

@dataclass
class GeneRecord:
    gene: str
    mechanisms: set[str] = field(default_factory=set)

    fong: int = 0
    cag: int = 0
    map2: int = 0

    expr_fong: Optional[float] = None
    expr_cag: Optional[float] = None
    expr_map2: Optional[float] = None

    role: Optional[str] = None
    function: Optional[str] = None
    pathway: Optional[str] = None

    # Where the preferred row/style came from.
    source_sheet: Optional[str] = None
    source_row: Optional[int] = None
    source_priority: int = 999

    # Earliest ordering position within a preferred mechanism.
    sort_key: tuple = (999, 999)

    doctor_highlight: bool = False
    yellow_fill: Optional[PatternFill] = None

    andrusiak: bool = False
    top_distinct: bool = False


# ---------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------

def clean_gene(x):
    if pd.isna(x):
        return None

    x = str(x).strip().upper()

    if not x:
        return None

    # Canonicalize labels such as:
    # "VNN1 (VANIN-1)" -> "VNN1"
    x = re.sub(r"\s*\([^)]*\)\s*$", "", x)

    # Reject manually written summary labels such as:
    # "ALL THREE DATASETS"
    # "FONG ET AL. ∩ OSHIKAWA ET AL. CAG"
    if any(ch.isspace() for ch in x):
        return None

    return x

def binary(value) -> int:
    """Return 1 only for an explicit binary 1."""
    try:
        return 1 if float(value) == 1 else 0
    except (TypeError, ValueError):
        return 0


def numeric(value) -> Optional[float]:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def sign(x: Optional[float]) -> Optional[int]:
    if x is None:
        return None
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


def rgb_of_cell(cell: Cell) -> Optional[str]:
    """
    Return an RGB color when the cell has a solid RGB fill.
    Theme/indexed colors are intentionally ignored.
    """
    fill = cell.fill
    if not fill or fill.fill_type != "solid":
        return None

    fg = fill.fgColor
    if fg.type != "rgb" or not fg.rgb:
        return None

    rgb = fg.rgb.upper()
    # openpyxl commonly stores alpha + RGB, e.g. FFFFC000.
    return rgb[-6:]


def row_has_yellow(ws, row: int, start_col: int = 1, end_col: int = 12):
    for col in range(start_col, min(end_col, ws.max_column) + 1):
        cell = ws.cell(row=row, column=col)
        rgb = rgb_of_cell(cell)
        if rgb in YELLOW_RGB:
            return True, copy(cell.fill)
    return False, None


def copy_style(src: Cell, dst: Cell):
    """
    Copy a cell's direct formatting.
    We intentionally do not create Excel tables or dashboard styling.
    """
    if src.has_style:
        dst._style = copy(src._style)
    if src.number_format:
        dst.number_format = src.number_format
    dst.alignment = copy(src.alignment)
    dst.protection = copy(src.protection)


def copy_column_widths(src_ws, dst_ws, first_col: int, last_col: int):
    for col in range(first_col, last_col + 1):
        letter = get_column_letter(col)
        width = src_ws.column_dimensions[letter].width
        if width is not None:
            dst_ws.column_dimensions[letter].width = width


def classify_overlap(rec: GeneRecord) -> str:
    labels = []
    if rec.cag == 1:
        labels.append("CAG")
    if rec.map2 == 1:
        labels.append("MAP2")
    return " + ".join(labels)


def classify_direction(rec: GeneRecord) -> str:
    """
    Compare Fong adult expression direction with only those embryo
    datasets in which the gene is actually present.
    """
    adult = sign(rec.expr_fong)
    if adult in (None, 0):
        return "Not comparable"

    embryo_signs = []

    if rec.cag == 1:
        s = sign(rec.expr_cag)
        if s not in (None, 0):
            embryo_signs.append(s)

    if rec.map2 == 1:
        s = sign(rec.expr_map2)
        if s not in (None, 0):
            embryo_signs.append(s)

    if not embryo_signs:
        return "Not comparable"

    same = [s == adult for s in embryo_signs]

    if all(same):
        return "Concordant"
    if not any(same):
        return "REVERSED"
    return "MIXED"


# ---------------------------------------------------------------------
# Reading the three cell-death sheets
# ---------------------------------------------------------------------

def read_mechanism_sheet(ws, mechanism: str) -> list[dict]:
    """
    Rows are accepted only when:
        - column A contains a gene
        - at least one of fong/cag/map2 is explicitly 0 or 1

    This deliberately ignores summary text that appears lower in some
    sheets.
    """
    role_col, function_col, pathway_col = ANNOTATION_COLS[ws.title]

    rows = []

    for r in range(2, ws.max_row + 1):
        gene = clean_gene(ws.cell(r, 1).value)
        if not gene:
            continue

        raw_memberships = [
            ws.cell(r, 2).value,
            ws.cell(r, 3).value,
            ws.cell(r, 4).value,
        ]

        valid_binary_row = any(
            v in (0, 1, 0.0, 1.0, "0", "1") for v in raw_memberships
        )
        if not valid_binary_row:
            continue

        highlighted, yellow_fill = row_has_yellow(ws, r)

        rows.append({
            "gene": gene,
            "mechanism": mechanism,
            "fong": binary(ws.cell(r, 2).value),
            "cag": binary(ws.cell(r, 3).value),
            "map2": binary(ws.cell(r, 4).value),
            "expr_fong": numeric(ws.cell(r, 5).value),
            "expr_cag": numeric(ws.cell(r, 6).value),
            "expr_map2": numeric(ws.cell(r, 7).value),
            "role": ws.cell(r, role_col).value,
            "function": ws.cell(r, function_col).value,
            "pathway": ws.cell(r, pathway_col).value,
            "source_sheet": ws.title,
            "source_row": r,
            "doctor_highlight": highlighted,
            "yellow_fill": yellow_fill,
        })

    return rows


def merge_mechanism_rows(rows: list[dict]) -> dict[str, GeneRecord]:
    merged: dict[str, GeneRecord] = {}

    for row in rows:
        gene = row["gene"]
        mechanism = row["mechanism"]
        priority = MECHANISM_PRIORITY[mechanism]

        if gene not in merged:
            merged[gene] = GeneRecord(gene=gene)

        rec = merged[gene]
        rec.mechanisms.add(mechanism)

        # Dataset membership is a union across sheets.
        rec.fong = max(rec.fong, row["fong"])
        rec.cag = max(rec.cag, row["cag"])
        rec.map2 = max(rec.map2, row["map2"])

        # Expression values should agree across sheets; keep the first
        # available non-null value.
        if rec.expr_fong is None and row["expr_fong"] is not None:
            rec.expr_fong = row["expr_fong"]
        if rec.expr_cag is None and row["expr_cag"] is not None:
            rec.expr_cag = row["expr_cag"]
        if rec.expr_map2 is None and row["expr_map2"] is not None:
            rec.expr_map2 = row["expr_map2"]

        # Prefer annotations/style from apoptosis, then autophagy,
        # then necroptosis.
        if priority < rec.source_priority:
            rec.source_priority = priority
            rec.source_sheet = row["source_sheet"]
            rec.source_row = row["source_row"]
            rec.role = row["role"]
            rec.function = row["function"]
            rec.pathway = row["pathway"]
            rec.sort_key = (priority, row["source_row"])

        # Fill missing annotations from another mechanism when useful.
        if not rec.role and row["role"]:
            rec.role = row["role"]
        if not rec.function and row["function"]:
            rec.function = row["function"]
        if not rec.pathway and row["pathway"]:
            rec.pathway = row["pathway"]

        if row["doctor_highlight"]:
            rec.doctor_highlight = True
            if rec.yellow_fill is None and row["yellow_fill"] is not None:
                rec.yellow_fill = row["yellow_fill"]

    return merged


# ---------------------------------------------------------------------
# Auxiliary workbook sets
# ---------------------------------------------------------------------

def read_andrusiak_genes(ws) -> set[str]:
    genes = set()

    # A = gene, E = andrusiak binary membership in this workbook.
    for r in range(2, ws.max_row + 1):
        gene = clean_gene(ws.cell(r, 1).value)
        if gene and binary(ws.cell(r, 5).value) == 1:
            genes.add(gene)

    return genes


def read_top_distinct_genes(ws) -> set[str]:
    """
    The sheet contains section titles in column A and gene rows with a
    numeric fold-change in column B. We use the latter criterion.
    """
    genes = set()

    for r in range(1, ws.max_row + 1):
        gene = clean_gene(ws.cell(r, 1).value)
        fc = numeric(ws.cell(r, 2).value)

        if gene and fc is not None:
            genes.add(gene)

    return genes


# ---------------------------------------------------------------------
# Writing sheets in the SAME visual language as the source workbook
# ---------------------------------------------------------------------

BASE_HEADERS = [
    "Gene",
    "fong",
    "cag",
    "map2",
    "expr_fong",
    "expr_cag",
    "expr_map2",
    None,
    "Role (cell death)",
    "Function",
    "Pathway",
]

META_HEADERS_COMMON = [
    "Mechanism(s)",
    "Embryo overlap",
    "Direction across stages",
    "Andrusiak overlap",
    "Top-distinct",
    "Doctor highlight",
]

META_HEADERS_MASTER = [
    "Mechanism(s)",
    "Adult-embryo shared?",
    "Embryo overlap",
    "Direction across stages",
    "Andrusiak overlap",
    "Top-distinct",
    "Doctor highlight",
    "Degree",
    "Betweenness",
    "PageRank",
    "Community",
    "Propagation score",
    "Network notes",
]


def remove_sheet_if_present(wb, name: str):
    if name in wb.sheetnames:
        del wb[name]


def style_new_header_cell(reference_cell: Cell, target_cell: Cell):
    copy_style(reference_cell, target_cell)


def write_base_row(
    wb,
    ws_out,
    output_row: int,
    rec: GeneRecord,
):
    """
    Write A:K.

    A:G and row styling are copied from the preferred source row.
    Role/Function/Pathway are copied from the source's annotation cells
    into the common I:J:K layout.
    """
    src_ws = wb[rec.source_sheet]
    src_row = rec.source_row

    values = [
        rec.gene,
        rec.fong,
        rec.cag,
        rec.map2,
        rec.expr_fong,
        rec.expr_cag,
        rec.expr_map2,
        None,
        rec.role,
        rec.function,
        rec.pathway,
    ]

    for col, value in enumerate(values, start=1):
        ws_out.cell(output_row, col).value = value

    # A:G: same physical layout in all source sheets.
    for col in range(1, 8):
        copy_style(
            src_ws.cell(src_row, col),
            ws_out.cell(output_row, col),
        )

    # H is the visual spacer in the apoptosis sheet.
    copy_style(
        wb[SHEET_APOPTOSIS].cell(2, 8),
        ws_out.cell(output_row, 8),
    )

    # Copy annotation styles from the actual annotation columns.
    role_col, function_col, pathway_col = ANNOTATION_COLS[rec.source_sheet]
    for src_col, dst_col in zip(
        (role_col, function_col, pathway_col),
        (9, 10, 11),
    ):
        copy_style(
            src_ws.cell(src_row, src_col),
            ws_out.cell(output_row, dst_col),
        )

    # If the doctor highlighted this candidate in ANY mechanism sheet,
    # preserve that fact visually on the gene cell without changing the
    # role/function color coding.
    if rec.doctor_highlight:
        ws_out.cell(output_row, 1).fill = (
            copy(rec.yellow_fill)
            if rec.yellow_fill is not None
            else copy(MANUAL_YELLOW)
        )


def style_metadata_area(ws, start_col: int, end_col: int, start_row: int, end_row: int):
    """
    Keep metadata intentionally plain: no blue banded Excel tables.
    """
    for row in ws.iter_rows(
        min_row=start_row,
        max_row=end_row,
        min_col=start_col,
        max_col=end_col,
    ):
        for cell in row:
            cell.font = Font(name="Calibri", size=11, color="000000")
            cell.alignment = Alignment(vertical="center")


def write_common_sheet(wb, records: list[GeneRecord]):
    remove_sheet_if_present(wb, OUT_COMMON)
    ws = wb.create_sheet(OUT_COMMON)

    ref = wb[SHEET_APOPTOSIS]

    headers = BASE_HEADERS + META_HEADERS_COMMON

    for col, header in enumerate(headers, start=1):
        ws.cell(1, col).value = header

        if col <= 11:
            # Copy the exact existing source header style.
            source_col = min(col, 11)
            copy_style(ref.cell(1, source_col), ws.cell(1, col))
        else:
            # Metadata headers use the same restrained header style.
            style_new_header_cell(ref.cell(1, 1), ws.cell(1, col))

    row = 2
    for rec in records:
        write_base_row(wb, ws, row, rec)

        meta = [
            ", ".join(sorted(rec.mechanisms, key=lambda x: MECHANISM_PRIORITY[x])),
            classify_overlap(rec),
            classify_direction(rec),
            "YES" if rec.andrusiak else "",
            "YES" if rec.top_distinct else "",
            "YES" if rec.doctor_highlight else "",
        ]

        for j, value in enumerate(meta, start=12):
            ws.cell(row, j).value = value

        row += 1

    style_metadata_area(ws, 12, len(headers), 2, max(2, row - 1))

    # Copy original A:K widths exactly.
    copy_column_widths(ref, ws, 1, 11)

    # New metadata widths.
    widths = {
        12: 22,  # mechanism
        13: 18,  # embryo overlap
        14: 23,  # direction
        15: 18,  # Andrusiak
        16: 16,  # top distinct
        17: 18,  # doctor highlight
    }
    for col, width in widths.items():
        ws.column_dimensions[get_column_letter(col)].width = width

    ws.freeze_panes = "A2"
    return ws


def write_master_sheet(wb, records: list[GeneRecord]):
    remove_sheet_if_present(wb, OUT_MASTER)
    ws = wb.create_sheet(OUT_MASTER)

    ref = wb[SHEET_APOPTOSIS]

    headers = BASE_HEADERS + META_HEADERS_MASTER

    for col, header in enumerate(headers, start=1):
        ws.cell(1, col).value = header

        if col <= 11:
            copy_style(ref.cell(1, col), ws.cell(1, col))
        else:
            style_new_header_cell(ref.cell(1, 1), ws.cell(1, col))

    row = 2
    for rec in records:
        write_base_row(wb, ws, row, rec)

        adult_embryo_shared = rec.fong == 1 and (rec.cag == 1 or rec.map2 == 1)

        meta = [
            ", ".join(sorted(rec.mechanisms, key=lambda x: MECHANISM_PRIORITY[x])),
            "YES" if adult_embryo_shared else "",
            classify_overlap(rec) if adult_embryo_shared else "",
            classify_direction(rec) if adult_embryo_shared else "",
            "YES" if rec.andrusiak else "",
            "YES" if rec.top_distinct else "",
            "YES" if rec.doctor_highlight else "",
            None,   # degree
            None,   # betweenness
            None,   # PageRank
            None,   # community
            None,   # propagation score
            None,   # network notes
        ]

        for j, value in enumerate(meta, start=12):
            ws.cell(row, j).value = value

        row += 1

    style_metadata_area(ws, 12, len(headers), 2, max(2, row - 1))
    copy_column_widths(ref, ws, 1, 11)

    widths = {
        12: 22,
        13: 22,
        14: 18,
        15: 23,
        16: 18,
        17: 16,
        18: 18,
        19: 12,  # degree
        20: 14,  # betweenness
        21: 12,  # PageRank
        22: 12,  # community
        23: 18,  # propagation
        24: 40,  # notes
    }
    for col, width in widths.items():
        ws.column_dimensions[get_column_letter(col)].width = width

    ws.freeze_panes = "A2"
    return ws


# ---------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build adult-embryo common and network-master sheets."
    )
    parser.add_argument(
        "input",
        type=Path,
        help="Path to the original 5-sheet Excel workbook.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output .xlsx path. Default: <input stem>-common-network-master.xlsx",
    )

    args = parser.parse_args()

    input_path = args.input.expanduser().resolve()

    if not input_path.exists():
        raise FileNotFoundError(input_path)

    if args.output is None:
        output_path = input_path.with_name(
            input_path.stem + "-common-network-master.xlsx"
        )
    else:
        output_path = args.output.expanduser().resolve()

    wb = load_workbook(input_path)

    required = [
        SHEET_APOPTOSIS,
        SHEET_TOP_DISTINCT,
        SHEET_NECROPTOSIS,
        SHEET_AUTOPHAGY,
        SHEET_ANDRUSIAK,
    ]
    missing = [s for s in required if s not in wb.sheetnames]

    if missing:
        raise ValueError(
            "Workbook is missing required sheet(s): "
            + ", ".join(missing)
        )

    all_rows = []

    all_rows.extend(
        read_mechanism_sheet(
            wb[SHEET_APOPTOSIS],
            "Apoptosis",
        )
    )
    all_rows.extend(
        read_mechanism_sheet(
            wb[SHEET_AUTOPHAGY],
            "Autophagy",
        )
    )
    all_rows.extend(
        read_mechanism_sheet(
            wb[SHEET_NECROPTOSIS],
            "Necroptosis",
        )
    )

    records = merge_mechanism_rows(all_rows)

    andrusiak = read_andrusiak_genes(wb[SHEET_ANDRUSIAK])
    top_distinct = read_top_distinct_genes(wb[SHEET_TOP_DISTINCT])

    for rec in records.values():
        rec.andrusiak = rec.gene in andrusiak
        rec.top_distinct = rec.gene in top_distinct

    # Preserve the broad ordering of the existing workbook:
    # apoptosis first, then autophagy-only, then necroptosis-only.
    master_records = sorted(
        records.values(),
        key=lambda r: (r.sort_key, r.gene),
    )

    common_records = [
        rec
        for rec in master_records
        if rec.fong == 1 and (rec.cag == 1 or rec.map2 == 1)
    ]

    write_common_sheet(wb, common_records)
    write_master_sheet(wb, master_records)

    wb.save(output_path)

    # -------------------------------------------------------------
    # Transparent console summary for the Methods / lab notebook.
    # -------------------------------------------------------------

    def n_common(mechanism: str):
        return sum(
            mechanism in rec.mechanisms
            for rec in common_records
        )

    direction_counts = {}
    for rec in common_records:
        d = classify_direction(rec)
        direction_counts[d] = direction_counts.get(d, 0) + 1

    print()
    print("R-Noel common-gene preprocessing complete")
    print("----------------------------------------")
    print(f"Input:  {input_path}")
    print(f"Output: {output_path}")
    print()
    print(f"Unique genes in network master: {len(master_records)}")
    print(f"Adult-embryo common genes:       {len(common_records)}")
    print(f"  Apoptosis common:              {n_common('Apoptosis')}")
    print(f"  Autophagy common:              {n_common('Autophagy')}")
    print(f"  Necroptosis common:            {n_common('Necroptosis')}")
    print()
    print("Direction across stages:")
    for key in ["Concordant", "REVERSED", "MIXED", "Not comparable"]:
        print(f"  {key:15s}: {direction_counts.get(key, 0)}")
    print()
    print(
        "Next step: use 'network master' as the node table for "
        "the separate local/global network-analysis script."
    )


if __name__ == "__main__":
    main()
