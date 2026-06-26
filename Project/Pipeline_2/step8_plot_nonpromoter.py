#!/usr/bin/env python3
import os
import textwrap
import pandas as pd
import matplotlib.pyplot as plt

IN_DIR  = os.path.join("Project", "New_Data_2")
OUT_DIR = os.path.join(IN_DIR, "figures")
os.makedirs(OUT_DIR, exist_ok=True)

def _load_category(cat: str) -> pd.DataFrame | None:
    """
    Prefer *_byID_nonpromoter.csv; fall back to *_nonpromoter.csv.
    Return None if missing or empty.
    """
    paths = [
        os.path.join(IN_DIR, f"E2F4_vs_E2F3_GO_BP_{cat}_byID_nonpromoter.csv"),
        os.path.join(IN_DIR, f"E2F4_vs_E2F3_GO_BP_{cat}_nonpromoter.csv"),
    ]
    for fp in paths:
        if os.path.exists(fp):
            try:
                df = pd.read_csv(fp)
                if df is not None and len(df) > 0:
                    return df
            except Exception:
                pass
    return None

def _wrap(series, width=48):
    return ["\n".join(textwrap.wrap(str(s), width=width, break_long_words=False)) for s in series]

def plot_category(cat: str, title: str, out_png: str, top_n: int | None = None):
    df = _load_category(cat)
    if df is None or df.empty:
        print(f"Skipping {cat}: no enriched terms.")
        return

    if top_n is not None and len(df) > top_n:
        df = df.nlargest(top_n, "Count")

    df = df.sort_values("Count", ascending=True)
    df["DescWrapped"] = _wrap(df["Description"], width=48)

    # Fixed small size so all figures align in the paper
    fig, ax = plt.subplots(figsize=(12, 12))
    bars = ax.barh(df["DescWrapped"], df["Count"])

    # annotate bars
    max_count = df["Count"].max()
    for b in bars:
        w = b.get_width()
        ax.text(w + max_count * 0.01, b.get_y() + b.get_height() / 2, f"{int(w)}",
                va="center", fontsize=8)

    ax.xaxis.grid(True, linestyle="--", linewidth=0.5, alpha=0.7)
    ax.set_axisbelow(True)
    ax.set_xlabel("Number of E2F4 targets (non-promoter)", fontsize=10)
    ax.set_title(title, fontsize=11, pad=10)
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(axis="x", labelsize=9)

    out_path = os.path.join(OUT_DIR, out_png)
    plt.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    print(f"→ Wrote {out_path}")

def main():
    plot_category("cell_cycle", "Top 20 Cell-cycle GO terms (non-promoter)", "GO2_cell_cycle.png", top_n=20)
    plot_category("DNA_repair", "DNA-repair GO terms (non-promoter)",        "GO2_DNA_repair.png")
    plot_category("apoptosis",  "Apoptosis GO terms (non-promoter)",         "GO2_apoptosis.png")
    plot_category("autophagy",  "Autophagy GO terms (non-promoter)",         "GO2_autophagy.png")

if __name__ == "__main__":
    main()
