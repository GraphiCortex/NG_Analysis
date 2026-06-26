#!/usr/bin/env python3
import os
import pandas as pd
from scipy.stats import pearsonr, spearmanr

bind_path = os.path.join("Project","New_Data","E2F4_vs_E2F3_sigProbe2Gene.csv")
expr_sig  = os.path.join("Project","New_Data_5","Oshikawa_CAG_DE_TREAT025_significant.csv")
go_xlsx   = os.path.join("Project","New_Data","GO_results.xlsx")
out_xlsx  = os.path.join("Project","New_Data_5","binding_vs_Oshikawa_CAG_expression.xlsx")

# 1) binding → gene level
bind_df = pd.read_csv(bind_path)
bind_gene = (
    bind_df.groupby("GENE_SYMBOL", as_index=False)
           .agg(mean_logFC_bind=("logFC","mean"),
                P_Value_bind   =("P.Value","min"),
                adjP_bind      =("adj.P.Val","min"))
           .rename(columns={"GENE_SYMBOL":"Gene"})
)
bind_gene["Gene"] = bind_gene["Gene"].str.upper().str.strip()

# 2) expression (CAG TREAT0.25) → gene level
expr_df = (pd.read_csv(expr_sig)
           .rename(columns={"log2FoldChange":"logFC_expr",
                            "pvalue":"P_Value_expr","padj":"adjP_expr"}))
expr_df["Gene"] = expr_df["Gene"].astype(str).str.upper().str.strip()
expr_gene = (expr_df.groupby("Gene", as_index=False)
             .agg(mean_logFC_expr=("logFC_expr","mean"),
                  P_Value_expr   =("P_Value_expr","min"),
                  adjP_expr      =("adjP_expr","min")))

# 3) merge
merged = pd.merge(bind_gene, expr_gene, on="Gene", how="inner")

# 3.1) GO families from Step7 workbook
go_full = pd.read_excel(go_xlsx, sheet_name="Sh_2")
def build_set(label):
    df = go_full[go_full["Category"] == label]
    if df.empty: return set()
    lists = df["geneID"].dropna().astype(str).str.split("/", expand=False)
    return {g.strip().upper() for sub in lists for g in sub}

labels = ["Cell cycle","DNA repair","Apoptosis","Autophagy","Necroptosis"]
sets = {lab: build_set(lab) for lab in labels}
flags = {"Cell cycle":"is_CC","DNA repair":"is_DR","Apoptosis":"is_AP","Autophagy":"is_AU","Necroptosis":"is_NE"}
for lab in labels:
    merged[flags[lab]] = merged["Gene"].isin(sets[lab])

def cats(g): return ";".join([lab for lab in labels if g in sets[lab]])
merged["Category"] = merged["Gene"].apply(cats)
flag_cols = [flags[lab] for lab in labels]
merged = merged[["Gene","mean_logFC_bind","P_Value_bind","adjP_bind",
                 "mean_logFC_expr","P_Value_expr","adjP_expr"] + flag_cols + ["Category"]]

# 4) write Excel + shading
with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as w:
    merged.to_excel(w, sheet_name="All_genes", index=False)
    wb, ws = w.book, w.sheets["All_genes"]
    fills = {"Cell cycle":"#FFF2CC","DNA repair":"#D9EAD3","Apoptosis":"#F4CCCC","Autophagy":"#C6EFCE","Necroptosis":"#E6E6FA"}
    fmts = {lab: wb.add_format({"bg_color":hex}) for lab,hex in fills.items()}
    nrows = merged.shape[0] + 1
    last_col = merged.shape[1]-1
    def idx_to_col(i):
        s=""; i0=i
        while True:
            i0, r = divmod(i0, 26); s = chr(65+r)+s
            if i0==0: break
            i0 -= 1
        return s
    rng = f"A2:{idx_to_col(last_col)}{nrows}"
    for lab in labels:
        c = idx_to_col(merged.columns.get_loc(flags[lab]))
        ws.conditional_format(rng, {"type":"formula","criteria": f"=${c}2","format": fmts[lab]})
    sig = merged[(merged["adjP_bind"]<0.05) &
                 (merged["mean_logFC_bind"].abs()>=0.25) &
                 (merged["adjP_expr"]<0.05) &
                 (merged["mean_logFC_expr"].abs()>=0.25)]
    sig.to_excel(w, sheet_name="Significant_both", index=False)

print(f"→ Wrote integrated workbook to {out_xlsx}")

# 5) quick correlations
if merged.shape[0] >= 2:
    r,p = pearsonr(merged["mean_logFC_bind"], merged["mean_logFC_expr"])
    rs,ps = spearmanr(merged["mean_logFC_bind"], merged["mean_logFC_expr"])
    print(f"Overall: Pearson r={r:.3f} (p={p:.2e}); Spearman ρ={rs:.3f} (p={ps:.2e})")

print("\nOverlap/correlation by category:")
for lab in labels:
    sub = merged[merged["Gene"].isin(sets[lab])]
    n = len(sub)
    if n < 2:
        print(f" • {lab}: n={n}")
    else:
        r,p = pearsonr(sub["mean_logFC_bind"], sub["mean_logFC_expr"])
        print(f" • {lab}: n={n}, r={r:.3f}, p={p:.2e})")
