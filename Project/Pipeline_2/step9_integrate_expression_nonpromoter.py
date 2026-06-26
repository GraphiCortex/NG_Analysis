#!/usr/bin/env python3
import pandas as pd
from scipy.stats import pearsonr, spearmanr

# --- Paths (Pipeline_2 / non-promoter) ---
bind_path = "Project/New_Data_2/E2F4_vs_E2F3_sigProbe2Gene_nonpromoter.csv"
expr_path = "Project/Data/Raw/microarray/Fong2022/Fong 2022 L2FC1 microarray data.csv"
go_xlsx   = "Project/New_Data_2/GO_results_nonpromoter.xlsx"
out_xlsx  = "Project/New_Data_2/binding_vs_TKO_expression_nonpromoter.xlsx"

# --- 1) Binding → gene level ---
bind_df = pd.read_csv(bind_path)
bind_gene = (
    bind_df.groupby("GENE_SYMBOL", as_index=False)
    .agg(mean_logFC_bind=("logFC","mean"),
         P_Value_bind   =("P.Value","min"),
         adjP_bind      =("adj.P.Val","min"))
    .rename(columns={"GENE_SYMBOL":"Gene"})
)
bind_gene["Gene"] = bind_gene["Gene"].str.upper().str.strip()
print(f"✓ Loaded {bind_gene.shape[0]} binding genes (non-promoter)")

# --- 2) Expression → gene level ---
expr_df = (
    pd.read_csv(expr_path)
    .rename(columns={"Unnamed: 2":"Gene"})
    [["Gene","log2FoldChange","pvalue","padj"]]
    .rename(columns={"log2FoldChange":"logFC_expr",
                     "pvalue":"P_Value_expr",
                     "padj":"adjP_expr"})
)
expr_df["Gene"] = expr_df["Gene"].str.upper().str.strip()
expr_gene = (
    expr_df.groupby("Gene", as_index=False)
    .agg(mean_logFC_expr=("logFC_expr","mean"),
         P_Value_expr   =("P_Value_expr","min"),
         adjP_expr      =("adjP_expr","min"))
)
print(f"✓ Loaded {expr_gene.shape[0]} expression genes")

# --- 3) Merge ---
merged = pd.merge(bind_gene, expr_gene, on="Gene", how="inner")
print(f"→ Overlapping binding+expression genes: {merged.shape[0]}")

# --- 3.1) Build gene sets from GO workbook (Sh_2) ---
go_full = pd.read_excel(go_xlsx, sheet_name="Sh_2")

def build_set(label):
    df = go_full[go_full["Category"] == label]
    if df.empty:
        return set()
    lists = df["geneID"].dropna().astype(str).str.split("/", expand=False)
    return {g.strip().upper() for sub in lists for g in sub}

# Four categories (no necroptosis here)
labels = ["Cell cycle","DNA repair","Apoptosis","Autophagy"]
sets = {lab: build_set(lab) for lab in labels}

# flags + compact category string
flag_map = {"Cell cycle":"is_CC", "DNA repair":"is_DR", "Apoptosis":"is_AP", "Autophagy":"is_AU"}
for lab, col in flag_map.items():
    merged[col] = merged["Gene"].isin(sets[lab])

def cats(g):
    return ";".join([lab for lab in labels if g in sets[lab]])
merged["Category"] = merged["Gene"].apply(cats)

# Column order
flag_cols = list(flag_map.values())
merged = merged[["Gene","mean_logFC_bind","P_Value_bind","adjP_bind",
                 "mean_logFC_expr","P_Value_expr","adjP_expr"] + flag_cols + ["Category"]]

# --- 4) Excel with shading ---
with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as writer:
    merged.to_excel(writer, sheet_name="All_genes", index=False)
    wb, ws = writer.book, writer.sheets["All_genes"]

    fills = {
        "Cell cycle":"#FFF2CC",
        "DNA repair":"#D9EAD3",
        "Apoptosis":"#F4CCCC",
        "Autophagy":"#C6EFCE"
    }
    fmts  = {lab: wb.add_format({"bg_color":hex}) for lab,hex in fills.items()}
    nrows = merged.shape[0] + 1

    # compute last column letter
    last_col_letter = ""
    x = len(merged.columns)-1
    while True:
        x, r = divmod(x, 26)
        last_col_letter = chr(65+r) + last_col_letter
        if x == 0: break
        x -= 1
    rng = f"A2:{last_col_letter}{nrows}"

    # add conditional formats based on flag columns
    for lab, col in flag_map.items():
        col_idx = merged.columns.get_loc(col)  # 0-based
        idx = col_idx
        col_letter = ""
        while True:
            idx, r = divmod(idx, 26)
            col_letter = chr(65+r) + col_letter
            if idx == 0: break
            idx -= 1
        ws.conditional_format(rng, {"type":"formula",
                                    "criteria":f"=${col_letter}2",
                                    "format":fmts[lab]})

    # Sheet 2: Significant_both (same thresholds as before)
    sig = merged[
        (merged["P_Value_bind"] < 0.05) &
        (merged["adjP_bind"]   < 0.05) &
        (merged["mean_logFC_bind"].abs() >= 0.5) &
        (merged["adjP_expr"]   < 0.05) &
        (merged["mean_logFC_expr"].abs() >= 1.0)
    ]
    sig.to_excel(writer, sheet_name="Significant_both", index=False)

print(f"→ Wrote integrated workbook to {out_xlsx}")

# --- 5) Overall correlations ---
if merged.shape[0] >= 2:
    r, p  = pearsonr(merged["mean_logFC_bind"], merged["mean_logFC_expr"])
    rs, ps = spearmanr(merged["mean_logFC_bind"], merged["mean_logFC_expr"])
    print(f"Overall: Pearson r={r:.3f} (p={p:.2e}); Spearman ρ={rs:.3f} (p={ps:.2e})")

# --- 6) Overlap/correlation by category ---
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
