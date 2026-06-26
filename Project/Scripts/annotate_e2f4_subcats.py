#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

# ─── 1) auto‐locate the project’s Processed folder ──────────────────────────────
script_path = Path(__file__).resolve()
project_dir = script_path.parents[1]            # …/Project
processed   = project_dir / "Data" / "Processed"
# ────────────────────────────────────────────────────────────────────────────────

# ─── 2) define paths ────────────────────────────────────────────────────────────
expr_int    = processed / "expression_integration"
limma_dir   = processed / "limma_results"

# your breakdown workbook (from task 1)
in_xlsx     = expr_int / "E2F4_gene_breakdown.xlsx"

# the new ORA symbols CSVs
cell_ora_fp = limma_dir / "E2F4_Cell_cycle_ORA_symbols.csv"
dna_ora_fp  = limma_dir / "E2F4_DNA_repair_ORA_symbols.csv"
apo_ora_fp  = limma_dir / "E2F4_Cell_death_ORA_symbols.csv"

# output annotated workbook
out_xlsx    = expr_int / "E2F4_gene_breakdown_with_subcats.xlsx"
# ────────────────────────────────────────────────────────────────────────────────

# ─── 3) load ORA tables ─────────────────────────────────────────────────────────
cell_ora = pd.read_csv(cell_ora_fp)
dna_ora  = pd.read_csv(dna_ora_fp)
apo_ora  = pd.read_csv(apo_ora_fp)
# ────────────────────────────────────────────────────────────────────────────────

# ─── 4) build GO‐term mapping: Symbol → [term1, term2, …] ────────────────────────
def build_go_map(ora_df, gene_col="geneID", term_col="Description"):
    mapping = {}
    for _, row in ora_df.iterrows():
        term = row[term_col]
        # geneID is now "Symbol1/Symbol2/…"
        genes = [g.strip() for g in str(row[gene_col]).split("/") if g.strip()]
        for g in genes:
            mapping.setdefault(g, []).append(term)
    return mapping

map_cell = build_go_map(cell_ora)
map_dna  = build_go_map(dna_ora)
map_apo  = build_go_map(apo_ora)

# ─── 5) read your three‐sheet breakdown workbook ────────────────────────────────
xls      = pd.ExcelFile(in_xlsx)
cell_df  = xls.parse("Cell_cycle")
dna_df   = xls.parse("DNA_repair")
apo_df   = xls.parse("Apoptosis")

# helper to detect the gene‐symbol column in each sheet
def find_gene_col(df):
    for c in df.columns:
        if "gene" in c.lower():
            return c
    raise KeyError(f"No gene column in {df.columns.tolist()}")

gcol = find_gene_col(cell_df)

# ─── 6) annotate each sheet with its GO subterms ───────────────────────────────
cell_df["GO_subterms"] = cell_df[gcol].map(lambda g: "; ".join(map_cell.get(g, [])))
dna_df["GO_subterms"]  = dna_df[gcol].map(lambda g: "; ".join(map_dna.get(g, [])))
apo_df["GO_subterms"]  = apo_df[gcol].map(lambda g: "; ".join(map_apo.get(g, [])))

# ─── 7) write out the new annotated workbook ───────────────────────────────────
# remove old file if locked (Windows precaution)
if out_xlsx.exists():
    out_xlsx.unlink()

with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
    cell_df.to_excel(writer, sheet_name="Cell_cycle", index=False)
    dna_df.to_excel(writer, sheet_name="DNA_repair", index=False)
    apo_df.to_excel(writer, sheet_name="Apoptosis", index=False)

print(f"✅ Wrote annotated workbook: {out_xlsx}")
