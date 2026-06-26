#!/usr/bin/env python3
import os
import textwrap
import pandas as pd
import matplotlib.pyplot as plt

IN_DIR  = os.path.join("Project", "New_Data")
OUT_DIR = os.path.join(IN_DIR, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

# Consistent, compact figure sizing for LaTeX triptychs
FIG_W, FIG_H = 7, 6
DEFAULT_TOP_N = 20  # keep panels visually similar

def _load_category(cat: str) -> pd.DataFrame | None:
    """Prefer *_byID.csv; fall back to keyword csv; return None if missing/empty."""
    byid = os.path.join(IN_DIR, f"E2F4_vs_E2F3_GO_BP_{cat}_byID.csv")
    kw   = os.path.join(IN_DIR, f"E2F4_vs_E2F3_GO_BP_{cat}.csv")
    for fp in (byid, kw):
        if os.path.exists(fp):
            df = pd.read_csv(fp)
            if df is not None and len(df) > 0:
                return df
    return None

def _wrap(series, width=50):
    return ["\n".join(textwrap.wrap(str(s), width=width, break_long_words=False)) for s in series]

def plot_category(cat: str, title: str, out_png: str, top_n: int | None = DEFAULT_TOP_N):
    df = _load_category(cat)
    if df is None or df.empty:
        print(f"Skipping {cat}: no enriched terms.")
        return

    if top_n is not None and len(df) > top_n:
        df = df.nlargest(top_n, "Count")

    # sort ascending so smallest is first (largest ends at bottom unless inverted)
    df = df.sort_values("Count", ascending=True).copy()
    df["DescWrapped"] = _wrap(df["Description"], width=50)

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    bars = ax.barh(df["DescWrapped"], df["Count"])

    max_count = df["Count"].max()
    for b in bars:
        w = b.get_width()
        ax.text(w + max_count*0.01, b.get_y()+b.get_height()/2, f"{int(w)}",
                va="center", fontsize=7)

    # grid and styling
    ax.xaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    ax.set_axisbelow(True)
    ax.set_xlabel("Number of E2F4 target promoters", fontsize=9)
    ax.set_title(title, fontsize=11, pad=10)
    ax.tick_params(axis="y", labelsize=7)
    ax.tick_params(axis="x", labelsize=8)

    # give long labels some breathing room
    plt.subplots_adjust(left=0.35, right=0.98, top=0.88, bottom=0.12)

    out_path = os.path.join(OUT_DIR, out_png)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"→ Wrote {out_path}")

def main():
    plot_category("cell_cycle",   "Top 20 Cell-cycle GO terms",          "GO_cell_cycle.png",   top_n=20)
    plot_category("DNA_repair",   "DNA-repair GO terms",                 "GO_DNA_repair.png")
    plot_category("apoptosis",    "Apoptosis GO terms",                  "GO_apoptosis.png")
    plot_category("autophagy",    "Autophagy GO terms",                  "GO_autophagy.png")
    plot_category("necroptosis",  "Necroptosis GO terms",                "GO_necroptosis.png")  # <-- updated

if __name__ == "__main__":
    main()
