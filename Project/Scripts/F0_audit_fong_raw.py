#!/usr/bin/env python3

"""
Phase F0 — Audit the Fong 2022 strict DEG dataset
=================================================

Purpose
-------
Establish a clean, reproducible starting point for the downstream
Fong analysis.

This script:

1. Loads the Fong 2022 "Sig L2FC1" differential-expression sheet.
2. Cleans and standardizes gene symbols.
3. Confirms the filtering represented by the sheet:
       padj < 0.05
       |log2FoldChange| >= 1
4. Classifies genes as up- or downregulated.
5. Compares the strict Fong DEG set against the previous
   126-gene "network master".
6. Exports an audit workbook containing:
       - summary
       - Fong strict DEG master
       - overlap with old network master
       - old-master genes absent from strict Fong
       - rows lacking valid gene symbols

Important
---------
The available Fong workbook is already filtered. Therefore this script
does NOT reconstruct the full unfiltered DESeq result universe.

A future |log2FC| >= 0.5 sensitivity analysis will require the original
unfiltered differential-expression table.

Repository layout expected
--------------------------
Project/
├── Data/
│   ├── Raw/
│   │   └── microarray/
│   │       └── Fong2022/
│   └── Processed/
│       └── FongDownstream/
│           └── F0_Audit/
├── Scripts/
│   └── F0_audit_fong_raw.py
└── ...

Run
---
From Project/:

    py .\\Scripts\\F0_audit_fong_raw.py
"""

from __future__ import annotations

from pathlib import Path
from collections import Counter
import sys

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


# ============================================================
# Project paths
# ============================================================

# Script:
# Project/Scripts/F0_audit_fong_raw.py
#
# Therefore parents[1] = Project/
PROJECT_ROOT = Path(__file__).resolve().parents[1]

FONG_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Raw"
    / "microarray"
    / "Fong2022"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F0_Audit"
)

OUTPUT_FILE = OUTPUT_DIR / "Fong_downstream_phase_F0_audit.xlsx"

FONG_SHEET = "Sig L2FC1"
OLD_MASTER_SHEET = "network master"


# ============================================================
# Utility functions
# ============================================================

def normalize_gene_symbol(value) -> str | None:
    """
    Standardize a gene symbol.

    Returns
    -------
    uppercase gene symbol, or None for missing/invalid entries.
    """
    if pd.isna(value):
        return None

    gene = str(value).strip()

    if not gene:
        return None

    if gene.upper() in {"NA", "N/A", "NONE", "NULL", "NAN"}:
        return None

    return gene.upper()


def find_single_file(
    search_dir: Path,
    patterns: list[str],
    description: str,
) -> Path:
    """
    Search recursively for a file matching one of several glob patterns.

    Raises an informative error if no file is found.

    If several files are found, prints them and exits rather than
    silently choosing one.
    """
    matches: list[Path] = []

    for pattern in patterns:
        matches.extend(search_dir.rglob(pattern))

    # Deduplicate resolved paths.
    unique_matches = sorted(set(p.resolve() for p in matches))

    if not unique_matches:
        raise FileNotFoundError(
            f"\nCould not find {description}.\n"
            f"Searched under:\n  {search_dir}\n"
            f"Patterns:\n  {patterns}\n"
        )

    if len(unique_matches) > 1:
        print(f"\n[ERROR] Multiple candidates found for {description}:")
        for path in unique_matches:
            print(f"    {path}")

        print(
            "\nPlease remove/rename ambiguous copies or specify the "
            "desired file directly in the script."
        )
        sys.exit(1)

    return unique_matches[0]


def format_output_workbook(path: Path) -> None:
    """
    Apply lightweight formatting to the exported audit workbook.
    """
    wb = load_workbook(path)

    header_fill = PatternFill(
        fill_type="solid",
        fgColor="D9E1F2",
    )

    gray_fill = PatternFill(
        fill_type="solid",
        fgColor="E7E6E6",
    )

    green_fill = PatternFill(
        fill_type="solid",
        fgColor="E2F0D9",
    )

    red_fill = PatternFill(
        fill_type="solid",
        fgColor="F4CCCC",
    )

    yellow_fill = PatternFill(
        fill_type="solid",
        fgColor="FFF2CC",
    )

    for ws in wb.worksheets:

        # Header formatting
        for cell in ws[1]:
            cell.fill = header_fill
            cell.font = Font(bold=True)
            cell.alignment = Alignment(
                horizontal="center",
                vertical="center",
            )

        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions

        # Sensible column widths
        for column_cells in ws.columns:
            letter = get_column_letter(column_cells[0].column)

            max_length = 0

            for cell in column_cells[:200]:
                if cell.value is None:
                    continue

                length = len(str(cell.value))
                max_length = max(max_length, length)

            width = min(max(max_length + 2, 12), 45)
            ws.column_dimensions[letter].width = width

    # --------------------------------------------------------
    # Master-sheet highlighting
    # --------------------------------------------------------

    if "Fong strict DEG master" in wb.sheetnames:
        ws = wb["Fong strict DEG master"]

        header_lookup = {
            cell.value: cell.column
            for cell in ws[1]
        }

        lfc_col = header_lookup.get("log2FoldChange")
        direction_col = header_lookup.get("Direction")
        old_master_col = header_lookup.get("Old_126_gene_master")

        for row in range(2, ws.max_row + 1):

            if lfc_col:
                value = ws.cell(row=row, column=lfc_col).value

                if isinstance(value, (int, float)):
                    if value > 0:
                        ws.cell(row=row, column=lfc_col).fill = green_fill
                    elif value < 0:
                        ws.cell(row=row, column=lfc_col).fill = red_fill

            if direction_col:
                direction = ws.cell(
                    row=row,
                    column=direction_col,
                ).value

                if direction == "Up":
                    ws.cell(
                        row=row,
                        column=direction_col,
                    ).fill = green_fill

                elif direction == "Down":
                    ws.cell(
                        row=row,
                        column=direction_col,
                    ).fill = red_fill

            if old_master_col:
                value = ws.cell(
                    row=row,
                    column=old_master_col,
                ).value

                if value == "YES":
                    ws.cell(
                        row=row,
                        column=old_master_col,
                    ).fill = yellow_fill

    # --------------------------------------------------------
    # Summary formatting
    # --------------------------------------------------------

    if "F0 Summary" in wb.sheetnames:
        ws = wb["F0 Summary"]

        ws.column_dimensions["A"].width = 46
        ws.column_dimensions["B"].width = 85

        for row in range(2, ws.max_row + 1):
            ws.cell(row=row, column=1).fill = gray_fill

        for row in ws.iter_rows():
            for cell in row:
                cell.alignment = Alignment(
                    vertical="top",
                    wrap_text=True,
                )

    wb.save(path)


# ============================================================
# Locate input files
# ============================================================

print()
print("=" * 72)
print("PHASE F0 — FONG STRICT DEG AUDIT")
print("=" * 72)

print(f"\n[INFO] Project root:")
print(f"       {PROJECT_ROOT}")

if not FONG_DIR.exists():
    raise FileNotFoundError(
        f"Expected Fong directory does not exist:\n{FONG_DIR}"
    )


# Look for the strict Fong workbook.
FONG_FILE = find_single_file(
    FONG_DIR,
    patterns=[
        "*L2FC1*.xlsx",
        "*l2fc1*.xlsx",
    ],
    description="Fong L2FC1 workbook",
)

print(f"\n[INFO] Fong workbook:")
print(f"       {FONG_FILE}")


# ============================================================
# Previous network-master workbook
# ============================================================

OLD_MASTER_FILE = (
    PROJECT_ROOT
    / "cell death project-analysis-ready.xlsx"
)

if not OLD_MASTER_FILE.exists():
    raise FileNotFoundError(
        "Could not find the previous analysis-ready workbook:\n"
        f"{OLD_MASTER_FILE}"
    )

print(f"\n[INFO] Previous analysis-ready workbook:")
print(f"       {OLD_MASTER_FILE}")

# ============================================================
# Load Fong data
# ============================================================

print(f"\n[INFO] Reading sheet '{FONG_SHEET}'...")

fong = pd.read_excel(
    FONG_FILE,
    sheet_name=FONG_SHEET,
)

print(f"[INFO] Rows loaded: {len(fong):,}")
print(f"[INFO] Columns loaded: {len(fong.columns)}")


# ============================================================
# Validate expected columns
# ============================================================

required_columns = {
    "baseMean",
    "log2FoldChange",
    "padj",
}

missing_required = required_columns.difference(
    fong.columns
)

if missing_required:
    raise ValueError(
        "\nMissing required Fong columns:\n"
        f"{sorted(missing_required)}"
    )


# In the current Fong workbook, the third column contains
# the gene symbol.
SOURCE_GENE_COLUMN = fong.columns[2]

print(
    f"[INFO] Gene-symbol source column detected as: "
    f"{repr(SOURCE_GENE_COLUMN)}"
)


# ============================================================
# Clean gene symbols
# ============================================================

fong["Gene"] = fong[SOURCE_GENE_COLUMN].map(
    normalize_gene_symbol
)

missing_symbol = fong[
    fong["Gene"].isna()
].copy()

valid = fong[
    fong["Gene"].notna()
].copy()


# ============================================================
# Check duplicate symbols
# ============================================================

gene_counts = Counter(valid["Gene"])

duplicate_genes = {
    gene: count
    for gene, count in gene_counts.items()
    if count > 1
}

if duplicate_genes:
    print(
        f"[WARNING] {len(duplicate_genes)} duplicated gene symbols "
        "were detected."
    )

    for gene, count in sorted(duplicate_genes.items()):
        print(f"          {gene}: {count}")

else:
    print("[INFO] No duplicated valid gene symbols detected.")


# ============================================================
# Confirm DEG thresholds
# ============================================================

valid["Pass_FDR_0.05"] = (
    valid["padj"] < 0.05
)

valid["Pass_absL2FC_1"] = (
    valid["log2FoldChange"].abs() >= 1.0
)

valid["Strict_pass"] = (
    valid["Pass_FDR_0.05"]
    &
    valid["Pass_absL2FC_1"]
)

strict_failures = valid[
    ~valid["Strict_pass"]
]

if strict_failures.empty:
    print(
        "[INFO] Every valid-symbol row satisfies "
        "padj < 0.05 and |log2FC| >= 1."
    )

else:
    print(
        f"[WARNING] {len(strict_failures):,} genes fail the "
        "expected strict threshold."
    )


# ============================================================
# Direction
# ============================================================

def classify_direction(value: float) -> str:
    if value > 0:
        return "Up"

    if value < 0:
        return "Down"

    return "Zero"


valid["Direction"] = (
    valid["log2FoldChange"]
    .map(classify_direction)
)


# ============================================================
# Load and compare previous network master
# ============================================================

old = None
old_gene_set: set[str] = set()
old_absent = pd.DataFrame()
overlap = pd.DataFrame()

if OLD_MASTER_FILE is not None:

    print(
        f"\n[INFO] Reading previous sheet "
        f"'{OLD_MASTER_SHEET}'..."
    )

    old = pd.read_excel(
        OLD_MASTER_FILE,
        sheet_name=OLD_MASTER_SHEET,
    )

    if "Gene" not in old.columns:
        raise ValueError(
            "Previous network master does not contain a "
            "'Gene' column."
        )

    old["Gene"] = old["Gene"].map(
        normalize_gene_symbol
    )

    old = old[
        old["Gene"].notna()
    ].copy()

    old_gene_set = set(old["Gene"])

    valid["Old_126_gene_master"] = (
        valid["Gene"].isin(old_gene_set)
    )

    # --------------------------------------------------------
    # Annotation columns retained from old network master
    # --------------------------------------------------------

    desired_old_columns = [
        "Gene",
        "Mechanism(s)",
        "Adult-embryo shared?",
        "Embryo overlap",
        "Direction across stages",
        "Andrusiak overlap",
        "Top-distinct",
        "Doctor highlight",
    ]

    available_old_columns = [
        column
        for column in desired_old_columns
        if column in old.columns
    ]

    old_annotations = old[
        available_old_columns
    ].copy()

    overlap = valid[
        valid["Old_126_gene_master"]
    ].copy()

    overlap = overlap.merge(
        old_annotations,
        on="Gene",
        how="left",
    )

    preferred_overlap_columns = [
        "Gene",
        "baseMean",
        "log2FoldChange",
        "padj",
        "Direction",
        "Mechanism(s)",
        "Adult-embryo shared?",
        "Embryo overlap",
        "Direction across stages",
        "Andrusiak overlap",
        "Top-distinct",
        "Doctor highlight",
    ]

    overlap = overlap[
        [
            column
            for column in preferred_overlap_columns
            if column in overlap.columns
        ]
    ]

    fong_gene_set = set(valid["Gene"])

    old_absent = old[
        ~old["Gene"].isin(fong_gene_set)
    ].copy()

else:
    valid["Old_126_gene_master"] = False


# ============================================================
# Summary statistics
# ============================================================

n_total = len(fong)
n_valid = len(valid)
n_missing = len(missing_symbol)
n_unique = valid["Gene"].nunique()
n_duplicate_symbols = len(duplicate_genes)

n_up = int(
    (valid["Direction"] == "Up").sum()
)

n_down = int(
    (valid["Direction"] == "Down").sum()
)

n_strict = int(
    valid["Strict_pass"].sum()
)

n_old_overlap = int(
    valid["Old_126_gene_master"].sum()
)

n_old_absent = (
    len(old_absent)
    if OLD_MASTER_FILE is not None
    else None
)


summary_rows = [
    [
        "Input Fong workbook",
        str(FONG_FILE.relative_to(PROJECT_ROOT)),
    ],
    [
        "Input sheet",
        FONG_SHEET,
    ],
    [
        "Filtering status",
        (
            "The available sheet is already prefiltered. "
            "It should not be interpreted as the full unfiltered "
            "Fong differential-expression table."
        ),
    ],
    [
        "Expected primary threshold",
        "padj < 0.05 and |log2FoldChange| >= 1",
    ],
    [
        "Total rows",
        n_total,
    ],
    [
        "Rows with valid gene symbols",
        n_valid,
    ],
    [
        "Rows with missing/N/A symbols",
        n_missing,
    ],
    [
        "Unique valid genes",
        n_unique,
    ],
    [
        "Duplicated gene symbols",
        n_duplicate_symbols,
    ],
    [
        "Genes passing strict threshold",
        n_strict,
    ],
    [
        "Upregulated genes",
        n_up,
    ],
    [
        "Downregulated genes",
        n_down,
    ],
    [
        "Overlap with old 126-gene network master",
        n_old_overlap,
    ],
    [
        "Old-master genes absent from strict Fong set",
        n_old_absent if n_old_absent is not None else "N/A",
    ],
    [
        "Important limitation",
        (
            "A |log2FC| >= 0.5 sensitivity analysis cannot be "
            "performed from this sheet because genes with "
            "0.5 <= |log2FC| < 1 have already been removed."
        ),
    ],
]


summary = pd.DataFrame(
    summary_rows,
    columns=[
        "Metric",
        "Value",
    ],
)


# ============================================================
# Build clean strict DEG master
# ============================================================

front_columns = [
    "Gene",
    "baseMean",
    "log2FoldChange",
    "padj",
    "Direction",
    "Pass_FDR_0.05",
    "Pass_absL2FC_1",
    "Strict_pass",
    "Old_126_gene_master",
]

remaining_columns = [
    column
    for column in valid.columns
    if column not in front_columns
    and column != SOURCE_GENE_COLUMN
]

master = valid[
    front_columns
    +
    remaining_columns
].copy()


# ============================================================
# Output directory
# ============================================================

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Export workbook
# ============================================================

print(f"\n[INFO] Writing audit workbook:")
print(f"       {OUTPUT_FILE}")

with pd.ExcelWriter(
    OUTPUT_FILE,
    engine="openpyxl",
) as writer:

    summary.to_excel(
        writer,
        sheet_name="F0 Summary",
        index=False,
    )

    master.to_excel(
        writer,
        sheet_name="Fong strict DEG master",
        index=False,
    )

    if OLD_MASTER_FILE is not None:

        overlap.to_excel(
            writer,
            sheet_name="Old-master overlap",
            index=False,
        )

        old_absent.to_excel(
            writer,
            sheet_name="Old-master absent",
            index=False,
        )

    missing_symbol.to_excel(
        writer,
        sheet_name="Missing symbols",
        index=False,
    )


# Apply useful visual formatting.
format_output_workbook(
    OUTPUT_FILE
)


# ============================================================
# Console report
# ============================================================

print()
print("=" * 72)
print("PHASE F0 COMPLETE")
print("=" * 72)

print(
    f"Total rows                     : "
    f"{n_total:,}"
)

print(
    f"Valid gene symbols             : "
    f"{n_valid:,}"
)

print(
    f"Unique valid genes             : "
    f"{n_unique:,}"
)

print(
    f"Missing/N/A symbols            : "
    f"{n_missing:,}"
)

print(
    f"Duplicated symbols             : "
    f"{n_duplicate_symbols:,}"
)

print()

print(
    f"Strict DEG genes               : "
    f"{n_strict:,}"
)

print(
    f"Upregulated                    : "
    f"{n_up:,}"
)

print(
    f"Downregulated                  : "
    f"{n_down:,}"
)

print()

if OLD_MASTER_FILE is not None:

    print(
        f"Overlap with old master        : "
        f"{n_old_overlap:,}"
    )

    print(
        f"Old-master genes absent        : "
        f"{n_old_absent:,}"
    )

else:

    print(
        "Old-master comparison          : "
        "SKIPPED"
    )

print()

print(
    "Output:"
)

print(
    f"    {OUTPUT_FILE.relative_to(PROJECT_ROOT)}"
)

print()
print(
    "[NEXT] Phase F1 will perform global functional enrichment "
    "on the strict Fong DEG universe before any cell-death/"
    "survival preselection."
)