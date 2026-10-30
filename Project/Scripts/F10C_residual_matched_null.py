#!/usr/bin/env python
"""
F10C_residual_matched_null.py
=============================

Selection- and expression-aware matched-null stress test for the F10 adult
residual programs.

Primary question
----------------
After conditioning on the features that made a gene likely to have a large
adult Rb-TKO signal, do the predefined biological programs still have
unusually large adult-vs-cortical residuals?

Primary null matching:
    1. same adult Rb-TKO direction (Up/Down), exactly;
    2. similar |adult Rb-TKO log2FC|;
    3. similar adult control expression, using mean(RB_THC_1..3) on log1p scale.

The null universe is the F10A complete-context universe, so genes must already
have adult, CAG, and MAP2 measurements. This is important because it avoids
sampling genes absent from the cortical branch.

IMPORTANT LIMITATION
--------------------
This script matches adult baseline expression but does not match cortical
baseline intensity, because that variable is not present in the F10A output
table. Therefore this analysis addresses selection on adult effect size and
adult expression, but it does NOT eliminate all RNA-seq-vs-microarray
platform/detectability effects. A later extension should add normalized
Oshikawa control intensity when the F8A gene-level expression matrix is
available.

Primary inference:
    two-sided matched-null P values + Benjamini-Hochberg FDR.

Directional P values are also reported but are SECONDARY because direction was
not independently preregistered for every program.

Outputs
-------
F10C_matched_null_primary.csv
F10C_matched_null_sensitivity.csv
F10C_matching_diagnostics.csv
F10C_null_distributions.csv.gz
F10C_summary.txt
figures/F10C_matched_null_summary.png
figures/F10C_FATE_DEATH_COMPOSITE_matched_null.png
figures/F10C_INNATE_STRESS_matched_null.png
"""

from __future__ import annotations

import argparse
import gzip
import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DEFAULT_RESIDUALS = (
    "Data/Processed/Dynamics/F10A_AdultResiduals/"
    "F10A_all_adult_residuals_ranked.csv"
)

DEFAULT_PROGRAMS = (
    "Data/Processed/Dynamics/F10B_ResidualPrograms/"
    "F10B_program_gene_contributions.csv"
)

DEFAULT_FONG = (
    "Data/Raw/microarray/Fong2022/"
    "Fong 2022 L2FC1 microarray data.xlsx"
)

DEFAULT_OUTDIR = "Data/Processed/Dynamics/F10C_MatchedNull"


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment without external dependencies."""
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]

    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0.0, 1.0)

    out = np.empty(n, dtype=float)
    out[order] = q
    return out


def robust_z(x: pd.Series) -> pd.Series:
    """Median/MAD robust standardization."""
    arr = pd.to_numeric(x, errors="coerce").to_numpy(dtype=float)
    med = np.nanmedian(arr)
    mad = 1.4826 * np.nanmedian(np.abs(arr - med))

    if not np.isfinite(mad) or mad <= 1e-12:
        mad = np.nanstd(arr, ddof=1)

    if not np.isfinite(mad) or mad <= 1e-12:
        mad = 1.0

    return pd.Series((arr - med) / mad, index=x.index)


def require_columns(df: pd.DataFrame, cols: list[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"{label} is missing required columns:\n"
            + "\n".join(f"  - {c}" for c in missing)
        )


def load_data(
    residual_path: Path,
    program_path: Path,
    fong_path: Path,
    fong_sheet: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    residuals = pd.read_csv(residual_path)
    programs = pd.read_csv(program_path)

    require_columns(
        residuals,
        [
            "Gene",
            "Ensembl_RbTKO",
            "RbTKO_LFC",
            "RbTKO_Direction",
            "CAG_LFC",
            "MAP2_LFC",
            "Adult_residual_z",
        ],
        "F10A residual table",
    )

    require_columns(
        programs,
        ["Program", "Gene", "Adult_residual_z"],
        "F10B program contribution table",
    )

    fong = pd.read_excel(fong_path, sheet_name=fong_sheet)

    if "Row.names" not in fong.columns:
        raise ValueError(
            "Could not find 'Row.names' in the Fong expression table. "
            "This column is expected to contain Ensembl IDs."
        )

    control_cols = ["RB_THC_1", "RB_THC_2", "RB_THC_3"]
    require_columns(fong, control_cols, "Fong expression table")

    fong_small = fong[
        ["Row.names", "baseMean", *control_cols]
    ].copy()

    fong_small = fong_small.rename(
        columns={"Row.names": "Ensembl_RbTKO"}
    )

    # Ensure one row per Ensembl ID before merge.
    fong_small = (
        fong_small
        .dropna(subset=["Ensembl_RbTKO"])
        .drop_duplicates("Ensembl_RbTKO", keep="first")
    )

    merged = residuals.merge(
        fong_small,
        on="Ensembl_RbTKO",
        how="left",
        validate="many_to_one",
    )

    merged["Adult_control_mean"] = merged[control_cols].mean(axis=1)
    merged["Adult_control_log1p"] = np.log1p(merged["Adult_control_mean"])
    merged["Abs_RbTKO_LFC"] = merged["RbTKO_LFC"].abs()

    # Complete-context requirement.
    needed = [
        "Gene",
        "RbTKO_LFC",
        "RbTKO_Direction",
        "CAG_LFC",
        "MAP2_LFC",
        "Adult_residual_z",
        "Adult_control_mean",
        "Adult_control_log1p",
    ]

    before = len(merged)
    merged = merged.dropna(subset=needed).copy()
    after = len(merged)

    if after == 0:
        raise ValueError("No complete rows remained after merging/filtering.")

    if after < before:
        print(
            f"WARNING: {before-after} residual rows were excluded because "
            "expression/matching variables were missing."
        )

    merged["Match_z_abs_lfc"] = robust_z(merged["Abs_RbTKO_LFC"])
    merged["Match_z_log_expr"] = robust_z(
        merged["Adult_control_log1p"]
    )

    merged = merged.reset_index(drop=True)

    return merged, programs


def program_gene_indices(
    universe: pd.DataFrame,
    programs: pd.DataFrame,
    program: str,
) -> np.ndarray:
    requested = (
        programs.loc[programs["Program"] == program, "Gene"]
        .dropna()
        .astype(str)
        .unique()
    )

    idx_by_gene = {
        str(gene): i
        for i, gene in enumerate(universe["Gene"].astype(str))
    }

    available = [
        idx_by_gene[g]
        for g in requested
        if g in idx_by_gene
    ]

    return np.asarray(available, dtype=int)


def build_neighbor_lists(
    universe: pd.DataFrame,
    program_idx: np.ndarray,
    k: int,
    exclude_all_program_genes: bool,
    all_program_gene_names: set[str],
) -> list[np.ndarray]:
    genes = universe["Gene"].astype(str).to_numpy()
    directions = universe["RbTKO_Direction"].astype(str).to_numpy()
    z_lfc = universe["Match_z_abs_lfc"].to_numpy(dtype=float)
    z_expr = universe["Match_z_log_expr"].to_numpy(dtype=float)

    current_program_genes = set(genes[program_idx].tolist())

    if exclude_all_program_genes:
        exclusion = all_program_gene_names
    else:
        exclusion = current_program_genes

    exclude_mask = np.isin(genes, list(exclusion))

    neighbor_lists: list[np.ndarray] = []

    for idx in program_idx:
        eligible = (
            (directions == directions[idx])
            & (~exclude_mask)
        )

        candidate_idx = np.where(eligible)[0]

        if len(candidate_idx) == 0:
            raise RuntimeError(
                f"No matching candidates for gene {genes[idx]} "
                f"with direction {directions[idx]}."
            )

        dist2 = (
            (z_lfc[candidate_idx] - z_lfc[idx]) ** 2
            +
            (z_expr[candidate_idx] - z_expr[idx]) ** 2
        )

        kk = min(k, len(candidate_idx))

        chosen_local = np.argpartition(
            dist2,
            kk - 1,
        )[:kk]

        chosen_local = chosen_local[
            np.argsort(dist2[chosen_local])
        ]

        neighbor_lists.append(candidate_idx[chosen_local])

    return neighbor_lists


def simulate_null(
    universe: pd.DataFrame,
    program_idx: np.ndarray,
    neighbor_lists: list[np.ndarray],
    iterations: int,
    seed: int,
) -> dict[str, np.ndarray]:
    rng = np.random.default_rng(seed)

    residual = universe["Adult_residual_z"].to_numpy(dtype=float)
    abs_lfc = universe["Abs_RbTKO_LFC"].to_numpy(dtype=float)
    log_expr = universe["Adult_control_log1p"].to_numpy(dtype=float)

    n = len(program_idx)

    null_residual = np.empty(iterations, dtype=float)
    null_abs_lfc = np.empty(iterations, dtype=float)
    null_log_expr = np.empty(iterations, dtype=float)

    for b in range(iterations):
        used: set[int] = set()
        selected: list[int] = []

        # Randomize order so one program gene does not always get first
        # access to overlapping matching candidates.
        for local_i in rng.permutation(n):
            candidates = neighbor_lists[local_i]

            # Small program sizes + large candidate universe mean this
            # loop almost always succeeds immediately.
            available = [
                int(c)
                for c in candidates
                if int(c) not in used
            ]

            if not available:
                raise RuntimeError(
                    "Matching failed because no unique candidate remained. "
                    "Increase --k."
                )

            chosen = int(
                available[
                    rng.integers(0, len(available))
                ]
            )

            used.add(chosen)
            selected.append(chosen)

        selected = np.asarray(selected, dtype=int)

        null_residual[b] = residual[selected].mean()
        null_abs_lfc[b] = abs_lfc[selected].mean()
        null_log_expr[b] = log_expr[selected].mean()

    return {
        "residual": null_residual,
        "abs_lfc": null_abs_lfc,
        "log_expr": null_log_expr,
    }


def summarize_test(
    universe: pd.DataFrame,
    program: str,
    program_idx: np.ndarray,
    null: dict[str, np.ndarray],
) -> dict[str, float | int | str]:
    residual = universe["Adult_residual_z"].to_numpy(dtype=float)
    abs_lfc = universe["Abs_RbTKO_LFC"].to_numpy(dtype=float)
    log_expr = universe["Adult_control_log1p"].to_numpy(dtype=float)

    observed = float(residual[program_idx].mean())
    null_values = null["residual"]

    null_mean = float(null_values.mean())
    null_sd = float(null_values.std(ddof=1))

    if null_sd > 0:
        effect_z = (observed - null_mean) / null_sd
    else:
        effect_z = np.nan

    # Primary: two-sided test around matched-null center.
    p_two = (
        1
        +
        np.count_nonzero(
            np.abs(null_values - null_mean)
            >=
            abs(observed - null_mean)
        )
    ) / (len(null_values) + 1)

    # Secondary: one-sided in the observed direction.
    if observed >= null_mean:
        p_directional = (
            1
            +
            np.count_nonzero(null_values >= observed)
        ) / (len(null_values) + 1)
        direction = "Above_matched_null"
    else:
        p_directional = (
            1
            +
            np.count_nonzero(null_values <= observed)
        ) / (len(null_values) + 1)
        direction = "Below_matched_null"

    return {
        "Program": program,
        "Genes_available": int(len(program_idx)),
        "Observed_mean_residual_z": observed,
        "Matched_null_mean": null_mean,
        "Matched_null_SD": null_sd,
        "Matched_effect_z": effect_z,
        "Matched_null_q025": float(
            np.quantile(null_values, 0.025)
        ),
        "Matched_null_q975": float(
            np.quantile(null_values, 0.975)
        ),
        "Matched_percentile": float(
            100.0 * np.mean(null_values < observed)
        ),
        "Primary_P_two_sided": float(p_two),
        "Secondary_P_directional": float(p_directional),
        "Observed_direction_vs_null": direction,
        "Program_mean_abs_RbTKO_LFC": float(
            abs_lfc[program_idx].mean()
        ),
        "Null_mean_abs_RbTKO_LFC": float(
            null["abs_lfc"].mean()
        ),
        "Delta_mean_abs_RbTKO_LFC": float(
            abs_lfc[program_idx].mean()
            -
            null["abs_lfc"].mean()
        ),
        "Program_mean_log1p_control_expr": float(
            log_expr[program_idx].mean()
        ),
        "Null_mean_log1p_control_expr": float(
            null["log_expr"].mean()
        ),
        "Delta_mean_log1p_control_expr": float(
            log_expr[program_idx].mean()
            -
            null["log_expr"].mean()
        ),
    }


def plot_primary_summary(
    summary: pd.DataFrame,
    outpath: Path,
) -> None:
    s = summary.sort_values(
        "Observed_mean_residual_z"
    ).reset_index(drop=True)

    y = np.arange(len(s))

    fig, ax = plt.subplots(figsize=(9, 6))

    low = (
        s["Matched_null_mean"]
        -
        s["Matched_null_q025"]
    )

    high = (
        s["Matched_null_q975"]
        -
        s["Matched_null_mean"]
    )

    ax.errorbar(
        s["Matched_null_mean"],
        y,
        xerr=np.vstack([low, high]),
        fmt="o",
        capsize=4,
        label="Matched-null mean and 95% interval",
    )

    ax.scatter(
        s["Observed_mean_residual_z"],
        y,
        marker="x",
        s=70,
        label="Observed program mean",
    )

    ax.axvline(0, linewidth=1)

    ax.set_yticks(y)
    ax.set_yticklabels(s["Program"])

    ax.set_xlabel("Mean adult residual z")
    ax.set_title(
        "F10C: observed program residuals vs matched null"
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_null_histogram(
    program: str,
    observed: float,
    null_values: np.ndarray,
    outpath: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.2))

    ax.hist(
        null_values,
        bins=45,
        density=True,
        alpha=0.75,
    )

    ax.axvline(
        observed,
        linewidth=2,
        label=f"Observed = {observed:.3f}",
    )

    ax.axvline(
        float(null_values.mean()),
        linestyle="--",
        linewidth=2,
        label=f"Matched-null mean = {null_values.mean():.3f}",
    )

    ax.set_xlabel("Mean adult residual z")
    ax.set_ylabel("Density")
    ax.set_title(
        f"F10C matched null: {program}"
    )
    ax.legend()

    fig.tight_layout()
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def run_primary(
    universe: pd.DataFrame,
    programs: pd.DataFrame,
    iterations: int,
    k: int,
    seed: int,
    outdir: Path,
) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    all_program_genes = set(
        programs["Gene"].dropna().astype(str)
    )

    rows = []
    null_store: dict[str, np.ndarray] = {}
    diagnostic_rows = []

    program_names = list(
        programs["Program"].dropna().astype(str).unique()
    )

    for i, program in enumerate(program_names):
        idx = program_gene_indices(
            universe,
            programs,
            program,
        )

        if len(idx) == 0:
            continue

        neighbors = build_neighbor_lists(
            universe=universe,
            program_idx=idx,
            k=k,
            exclude_all_program_genes=False,
            all_program_gene_names=all_program_genes,
        )

        null = simulate_null(
            universe=universe,
            program_idx=idx,
            neighbor_lists=neighbors,
            iterations=iterations,
            seed=seed + i * 1009,
        )

        row = summarize_test(
            universe=universe,
            program=program,
            program_idx=idx,
            null=null,
        )

        rows.append(row)
        null_store[program] = null["residual"]

        diagnostic_rows.append(
            {
                "Program": program,
                "Genes_available": len(idx),
                "K_nearest_neighbors": k,
                "Iterations": iterations,
                "Delta_mean_abs_RbTKO_LFC":
                    row["Delta_mean_abs_RbTKO_LFC"],
                "Delta_mean_log1p_control_expr":
                    row["Delta_mean_log1p_control_expr"],
            }
        )

    summary = pd.DataFrame(rows)

    summary["Primary_FDR_two_sided"] = bh_fdr(
        summary["Primary_P_two_sided"].to_numpy()
    )

    summary["Secondary_FDR_directional"] = bh_fdr(
        summary["Secondary_P_directional"].to_numpy()
    )

    diagnostics = pd.DataFrame(diagnostic_rows)

    summary.to_csv(
        outdir / "F10C_matched_null_primary.csv",
        index=False,
    )

    diagnostics.to_csv(
        outdir / "F10C_matching_diagnostics.csv",
        index=False,
    )

    return summary, null_store


def run_sensitivity(
    universe: pd.DataFrame,
    programs: pd.DataFrame,
    k_values: list[int],
    iterations: int,
    seed: int,
) -> pd.DataFrame:
    all_program_genes = set(
        programs["Gene"].dropna().astype(str)
    )

    rows = []
    program_names = list(
        programs["Program"].dropna().astype(str).unique()
    )

    setting_id = 0

    for k in k_values:
        for exclude_all in [False, True]:
            current_rows = []

            for i, program in enumerate(program_names):
                idx = program_gene_indices(
                    universe,
                    programs,
                    program,
                )

                if len(idx) == 0:
                    continue

                neighbors = build_neighbor_lists(
                    universe=universe,
                    program_idx=idx,
                    k=k,
                    exclude_all_program_genes=exclude_all,
                    all_program_gene_names=all_program_genes,
                )

                null = simulate_null(
                    universe=universe,
                    program_idx=idx,
                    neighbor_lists=neighbors,
                    iterations=iterations,
                    seed=(
                        seed
                        + setting_id * 100_003
                        + i * 1009
                    ),
                )

                row = summarize_test(
                    universe=universe,
                    program=program,
                    program_idx=idx,
                    null=null,
                )

                row["K"] = k
                row["Exclude_all_program_genes"] = exclude_all
                row["Sensitivity_iterations"] = iterations

                current_rows.append(row)

            current = pd.DataFrame(current_rows)

            current["Primary_FDR_two_sided"] = bh_fdr(
                current["Primary_P_two_sided"].to_numpy()
            )

            current["Secondary_FDR_directional"] = bh_fdr(
                current["Secondary_P_directional"].to_numpy()
            )

            rows.append(current)
            setting_id += 1

    return pd.concat(rows, ignore_index=True)


def write_null_store(
    null_store: dict[str, np.ndarray],
    path: Path,
) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write("Program,Iteration,Matched_null_mean_residual_z\n")

        for program, values in null_store.items():
            safe_program = program.replace(",", "_")

            for i, value in enumerate(values):
                handle.write(
                    f"{safe_program},{i},{value:.10g}\n"
                )


def write_summary_text(
    universe: pd.DataFrame,
    summary: pd.DataFrame,
    sensitivity: pd.DataFrame,
    args: argparse.Namespace,
    outpath: Path,
) -> None:
    lines = []

    lines.append("F10C — MATCHED-NULL ADULT RESIDUAL STRESS TEST")
    lines.append("=" * 52)
    lines.append("")
    lines.append(f"Complete matched-null universe: {len(universe)} genes")
    lines.append(f"Primary iterations: {args.iterations}")
    lines.append(f"Primary k-nearest candidate pool: {args.k}")
    lines.append(f"Random seed: {args.seed}")
    lines.append("")
    lines.append("Primary matching variables:")
    lines.append("  - exact adult Rb-TKO direction")
    lines.append("  - |adult Rb-TKO log2FC|")
    lines.append("  - log1p mean adult control expression (RB_THC_1..3)")
    lines.append("")
    lines.append(
        "Primary inference uses two-sided matched-null P values and BH FDR."
    )
    lines.append(
        "Directional P values are secondary/exploratory."
    )
    lines.append("")
    lines.append("IMPORTANT LIMITATION:")
    lines.append(
        "  Cortical baseline expression intensity is not present in the "
        "F10A residual table, so matching is not yet performed on "
        "Oshikawa control intensity. The universe is nevertheless "
        "restricted to complete-context genes with CAG and MAP2 values."
    )
    lines.append("")
    lines.append("PRIMARY RESULTS")
    lines.append("-" * 52)

    ordered = summary.sort_values(
        "Primary_P_two_sided"
    )

    for _, row in ordered.iterrows():
        lines.append(
            f"{row['Program']}: "
            f"observed={row['Observed_mean_residual_z']:.3f}, "
            f"matched_null={row['Matched_null_mean']:.3f}, "
            f"P2={row['Primary_P_two_sided']:.5g}, "
            f"FDR2={row['Primary_FDR_two_sided']:.5g}, "
            f"Pdir={row['Secondary_P_directional']:.5g}, "
            f"FDRdir={row['Secondary_FDR_directional']:.5g}"
        )

    lines.append("")
    lines.append("SENSITIVITY")
    lines.append("-" * 52)

    for program in summary["Program"]:
        sub = sensitivity[
            sensitivity["Program"] == program
        ]

        if len(sub):
            lines.append(
                f"{program}: two-sided P range "
                f"[{sub['Primary_P_two_sided'].min():.5g}, "
                f"{sub['Primary_P_two_sided'].max():.5g}], "
                f"matched-null mean range "
                f"[{sub['Matched_null_mean'].min():.3f}, "
                f"{sub['Matched_null_mean'].max():.3f}]"
            )

    lines.append("")
    lines.append("INTERPRETATION GUARDRAILS")
    lines.append("-" * 52)
    lines.append(
        "- A program that was significant in F10B but not here is "
        "selection-sensitive; do not call it independently enriched."
    )
    lines.append(
        "- A matched-null signal is evidence of unusual adult-vs-cortical "
        "divergence conditional on adult effect size/expression; it is "
        "not evidence of causality."
    )
    lines.append(
        "- Program tests overlap in gene membership and are not independent "
        "biological replications."
    )
    lines.append(
        "- Cross-platform effects remain possible because Fong is RNA-seq "
        "and Oshikawa is Agilent microarray."
    )

    outpath.write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--residuals",
        default=DEFAULT_RESIDUALS,
        help="F10A_all_adult_residuals_ranked.csv",
    )

    parser.add_argument(
        "--programs",
        default=DEFAULT_PROGRAMS,
        help="F10B_program_gene_contributions.csv",
    )

    parser.add_argument(
        "--fong",
        default=DEFAULT_FONG,
        help="Fong Rb-TKO RNA-seq differential-expression workbook",
    )

    parser.add_argument(
        "--fong-sheet",
        default="Sig L2FC1",
    )

    parser.add_argument(
        "--outdir",
        default=DEFAULT_OUTDIR,
    )

    parser.add_argument(
        "--iterations",
        type=int,
        default=20_000,
    )

    parser.add_argument(
        "--sensitivity-iterations",
        type=int,
        default=2_000,
    )

    parser.add_argument(
        "--k",
        type=int,
        default=25,
    )

    parser.add_argument(
        "--sensitivity-k",
        default="10,25,50,100",
        help="Comma-separated k values",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=20261022,
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    residual_path = Path(args.residuals)
    program_path = Path(args.programs)
    fong_path = Path(args.fong)
    outdir = Path(args.outdir)
    figdir = outdir / "figures"

    outdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    universe, programs = load_data(
        residual_path=residual_path,
        program_path=program_path,
        fong_path=fong_path,
        fong_sheet=args.fong_sheet,
    )

    # Audit the exact analysis universe.
    universe.to_csv(
        outdir / "F10C_analysis_universe.csv",
        index=False,
    )

    print("")
    print("F10C matched-null analysis")
    print("--------------------------")
    print(f"Universe: {len(universe)} genes")
    print(
        "Matching: exact direction + |RbTKO LFC| "
        "+ adult control expression"
    )
    print("")

    primary, null_store = run_primary(
        universe=universe,
        programs=programs,
        iterations=args.iterations,
        k=args.k,
        seed=args.seed,
        outdir=outdir,
    )

    k_values = [
        int(x.strip())
        for x in args.sensitivity_k.split(",")
        if x.strip()
    ]

    sensitivity = run_sensitivity(
        universe=universe,
        programs=programs,
        k_values=k_values,
        iterations=args.sensitivity_iterations,
        seed=args.seed + 9_000_000,
    )

    sensitivity.to_csv(
        outdir / "F10C_matched_null_sensitivity.csv",
        index=False,
    )

    write_null_store(
        null_store,
        outdir / "F10C_null_distributions.csv.gz",
    )

    plot_primary_summary(
        primary,
        figdir / "F10C_matched_null_summary.png",
    )

    primary_idx = primary.set_index("Program")

    for program in [
        "FATE_DEATH_COMPOSITE",
        "INNATE_STRESS",
    ]:
        if program in null_store and program in primary_idx.index:
            plot_null_histogram(
                program=program,
                observed=float(
                    primary_idx.loc[
                        program,
                        "Observed_mean_residual_z",
                    ]
                ),
                null_values=null_store[program],
                outpath=(
                    figdir
                    /
                    f"F10C_{program}_matched_null.png"
                ),
            )

    write_summary_text(
        universe=universe,
        summary=primary,
        sensitivity=sensitivity,
        args=args,
        outpath=outdir / "F10C_summary.txt",
    )

    display_cols = [
        "Program",
        "Genes_available",
        "Observed_mean_residual_z",
        "Matched_null_mean",
        "Matched_effect_z",
        "Primary_P_two_sided",
        "Primary_FDR_two_sided",
        "Secondary_P_directional",
        "Secondary_FDR_directional",
    ]

    print(
        primary[
            display_cols
        ]
        .sort_values("Primary_P_two_sided")
        .to_string(index=False)
    )

    print("")
    print("Outputs written to:")
    print(outdir)
    print("")
    print(
        "IMPORTANT: treat the two-sided matched-null FDR as the "
        "primary result. The directional test is secondary."
    )


if __name__ == "__main__":
    main()
