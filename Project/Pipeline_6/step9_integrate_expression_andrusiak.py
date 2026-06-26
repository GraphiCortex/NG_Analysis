#!/usr/bin/env python3
import os
import pandas as pd
from scipy.stats import pearsonr, spearmanr

bind_path = os.path.join("Project", "New_Data", "E2F4_vs_E2F3_sigProbe2Gene.csv")
expr_csv  = os.path.join("Project", "New_Data_6", "Andrusiak_filtered_full.csv")
go_xlsx   = os.path.join("Project", "New_Data", "GO_results.xlsx")
out_xlsx  = os.path.join("Project", "New_Data_6", "binding_vs_Andrusiak_expression.xlsx")

# --- 1) Binding -> gene level ---
bind_df = pd.read_csv(bind_path)
bind_gene = (
    bind_df.groupby("GENE_SYMBOL", as_index=False)
           .agg(mean_logFC_bind=("logFC", "mean"),
                P_Value_bind=("P.Value", "min"),
                adjP_bind=("adj.P.Val", "min"))
           .rename(columns={"GENE_SYMBOL": "Gene"})
)
bind_gene["Gene"] = bind_gene["Gene"].astype(str).str.upper().str.strip()

# --- 2) Andrusiak expression -> gene level ---
expr_df = pd.read_csv(expr_csv)

needed = {"Gene", "log2FoldChange", "baseMean"}
missing = needed - set(expr_df.columns)
if missing:
    raise ValueError(f"Andrusiak file is missing columns: {sorted(missing)}")

expr_df = expr_df.rename(columns={
    "log2FoldChange": "logFC_expr"
})

expr_df["Gene"] = expr_df["Gene"].astype(str).str.upper().str.strip()

expr_gene = (
    expr_df.groupby("Gene", as_index=False)
           .agg(mean_logFC_expr=("logFC_expr", "mean"),
                baseMean_expr=("baseMean", "mean"))
)

# --- 3) Merge binding <-> expression ---
merged = pd.merge(bind_gene, expr_gene, on="Gene", how="inner")

# --- 3.1) Build GO family sets from GO_results.xlsx (Sh_2) ---
go_full = pd.read_excel(go_xlsx, sheet_name="Sh_2")

def build_set(label):
    df = go_full[go_full["Category"] == label]
    if df.empty:
        return set()
    lists = df["geneID"].dropna().astype(str).str.split("/", expand=False)
    return {g.strip().upper() for sub in lists for g in sub}

labels = ["Cell cycle", "DNA repair", "Apoptosis", "Autophagy", "Necroptosis"]
sets = {lab: build_set(lab) for lab in labels}

flag_names = {
    "Cell cycle": "is_CC",
    "DNA repair": "is_DR",
    "Apoptosis": "is_AP",
    "Autophagy": "is_AU",
    "Necroptosis": "is_NE"
}

for lab in labels:
    merged[flag_names[lab]] = merged["Gene"].isin(sets[lab])

def cats(g):
    return ";".join([lab for lab in labels if g in sets[lab]])

merged["Category"] = merged["Gene"].apply(cats)

flag_cols = [flag_names[lab] for lab in labels]
merged = merged[
    ["Gene", "mean_logFC_bind", "P_Value_bind", "adjP_bind",
     "mean_logFC_expr", "baseMean_expr"] + flag_cols + ["Category"]
]

# --- 4) Write Excel with shading ---
with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as writer:
    merged.to_excel(writer, sheet_name="All_genes", index=False)
    wb, ws = writer.book, writer.sheets["All_genes"]

    fills = {
        "Cell cycle": "#FFF2CC",
        "DNA repair": "#D9EAD3",
        "Apoptosis": "#F4CCCC",
        "Autophagy": "#C6EFCE",
        "Necroptosis": "#E6E6FA"
    }
    fmts = {lab: wb.add_format({"bg_color": hexcol}) for lab, hexcol in fills.items()}

    nrows = merged.shape[0] + 1
    last_col_idx = merged.shape[1] - 1

    def idx_to_col(idx: int) -> str:
        s = ""
        i = idx
        while True:
            i, r = divmod(i, 26)
            s = chr(65 + r) + s
            if i == 0:
                break
            i -= 1
        return s

    rng = f"A2:{idx_to_col(last_col_idx)}{nrows}"

    for lab in labels:
        col_idx = merged.columns.get_loc(flag_names[lab])
        col_letter = idx_to_col(col_idx)
        ws.conditional_format(rng, {
            "type": "formula",
            "criteria": f"=${col_letter}2",
            "format": fmts[lab]
        })

    # "Significant_both" here means:
    # binding significant by Julien criteria
    # Andrusiak side already legacy-filtered, so we just require |expr log2FC| >= 0.25
    sig = merged[
        (merged["P_Value_bind"] < 0.05) &
        (merged["adjP_bind"] < 0.05) &
        (merged["mean_logFC_bind"].abs() >= 0.25) &
        (merged["mean_logFC_expr"].abs() >= 0.25)
    ]
    sig.to_excel(writer, sheet_name="Significant_both", index=False)

print(f"→ Wrote integrated workbook to {out_xlsx}")

# --- 5) Overall correlations ---
if merged.shape[0] >= 2:
    r, p = pearsonr(merged["mean_logFC_bind"], merged["mean_logFC_expr"])
    rs, ps = spearmanr(merged["mean_logFC_bind"], merged["mean_logFC_expr"])
    print(f"Overall: Pearson r={r:.3f} (p={p:.2e}); Spearman ρ={rs:.3f} (p={ps:.2e})")

# --- 6) By-category overlap + Pearson ---
print("\nOverlap/correlation by category:")
for lab in labels:
    sub = merged[merged["Gene"].isin(sets[lab])]
    n = sub.shape[0]
    if n == 0:
        print(f" • {lab}: n=0")
    elif n == 1:
        print(f" • {lab}: n=1 (overlap only; correlation not defined)")
    else:
        r, p = pearsonr(sub["mean_logFC_bind"], sub["mean_logFC_expr"])
        print(f" • {lab}: n={n}, r={r:.3f}, p={p:.2e})")