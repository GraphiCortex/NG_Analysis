from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(".")

INFILE = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F8C_StatePrograms"
    / "F8C_oshikawa_sample_state_scores.csv"
)

OUTDIR = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F8D_StateBalance"
)

OUTDIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Load module scores
# ============================================================

df = pd.read_csv(INFILE)


# ============================================================
# Derived mechanistic balance variables
#
# Positive value means numerator program exceeds denominator
# program within that sample.
# ============================================================

BALANCES = {

    # Can the forced-cycle state proceed into mitotic execution?
    "MITOTIC_MINUS_REPLICATION":
        (
            "MITOTIC_EXECUTION",
            "REPLICATION_ATR",
        ),

    # Is mitotic execution proportional to overall cycle drive?
    "MITOTIC_MINUS_E2F_CYCLE":
        (
            "MITOTIC_EXECUTION",
            "E2F_CELL_CYCLE",
        ),

    # Is checkpoint/repair capacity keeping pace with
    # replication burden?
    "REPAIR_MINUS_REPLICATION":
        (
            "DNA_REPAIR_CHECKPOINT",
            "REPLICATION_ATR",
        ),

    # Is intrinsic apoptosis increasing relative to
    # replication-state activation?
    "APOPTOSIS_MINUS_REPLICATION":
        (
            "INTRINSIC_APOPTOSIS",
            "REPLICATION_ATR",
        ),
}


for name, (a, b) in BALANCES.items():

    df[name] = (
        df[a]
        -
        df[b]
    )


# ============================================================
# Design coding
# ============================================================

df["TKO"] = (
    df["Group"]
    .str.endswith("TKO")
    .astype(int)
)

df["MAP2"] = (
    df["Group"]
    .str.startswith("MAP2")
    .astype(int)
)

df["Interaction"] = (
    df["TKO"]
    *
    df["MAP2"]
)


X = np.column_stack([
    np.ones(len(df)),
    df["TKO"].values,
    df["MAP2"].values,
    df["Interaction"].values,
])


# ============================================================
# Contrast fitting
# ============================================================

def fit_model(y):

    y = np.asarray(y, dtype=float)

    beta = np.linalg.lstsq(
        X,
        y,
        rcond=None,
    )[0]

    resid = y - X @ beta

    df_resid = (
        len(y)
        -
        X.shape[1]
    )

    sigma2 = (
        np.sum(resid ** 2)
        /
        df_resid
    )

    cov = (
        sigma2
        *
        np.linalg.inv(
            X.T @ X
        )
    )


    def contrast(c):

        c = np.asarray(
            c,
            dtype=float,
        )

        est = float(
            c @ beta
        )

        se = float(
            np.sqrt(
                c @ cov @ c
            )
        )

        t = (
            est / se
            if se > 0
            else np.nan
        )

        p = (
            2
            *
            stats.t.sf(
                abs(t),
                df=df_resid,
            )
            if np.isfinite(t)
            else np.nan
        )

        return est, se, t, p


    return {
        "CAG": contrast(
            [0, 1, 0, 0]
        ),

        "MAP2": contrast(
            [0, 1, 0, 1]
        ),

        "Interaction": contrast(
            [0, 0, 0, 1]
        ),
    }


# ============================================================
# BH correction
# ============================================================

def bh_adjust(values):

    p = np.asarray(
        values,
        dtype=float,
    )

    order = np.argsort(p)

    ranked = p[order]

    n = len(p)

    adj = (
        ranked
        *
        n
        /
        np.arange(
            1,
            n + 1,
        )
    )

    adj = np.minimum.accumulate(
        adj[::-1]
    )[::-1]

    adj = np.minimum(
        adj,
        1,
    )

    result = np.empty(n)

    result[order] = adj

    return result


# ============================================================
# Fit all balance variables
# ============================================================

rows = []


for name in BALANCES:

    fit = fit_model(
        df[name]
    )

    cag = fit["CAG"]
    map2 = fit["MAP2"]
    interaction = fit["Interaction"]


    means = (
        df.groupby(
            "Group"
        )[name]
        .mean()
    )


    rows.append({

        "State_balance":
            name,

        "CAG_Control_mean":
            means.get(
                "CAG_Control",
                np.nan,
            ),

        "CAG_TKO_mean":
            means.get(
                "CAG_TKO",
                np.nan,
            ),

        "MAP2_Control_mean":
            means.get(
                "MAP2_Control",
                np.nan,
            ),

        "MAP2_TKO_mean":
            means.get(
                "MAP2_TKO",
                np.nan,
            ),

        "CAG_TKO_effect":
            cag[0],

        "CAG_P":
            cag[3],

        "MAP2_TKO_effect":
            map2[0],

        "MAP2_P":
            map2[3],

        "MAP2_minus_CAG_interaction":
            interaction[0],

        "Interaction_P":
            interaction[3],
    })


results = pd.DataFrame(
    rows
)


results["Interaction_FDR"] = (
    bh_adjust(
        results[
            "Interaction_P"
        ].values
    )
)


# ============================================================
# Export sample scores + results
# ============================================================

df.to_csv(
    OUTDIR
    / "F8D_sample_balance_scores.csv",
    index=False,
)

results.to_csv(
    OUTDIR
    / "F8D_state_balance_results.csv",
    index=False,
)


# ============================================================
# Console summary
# ============================================================

print()
print("=" * 105)
print("F8D — DERIVED STATE-BALANCE ANALYSIS")
print("=" * 105)

print()

print(
    results[
        [
            "State_balance",
            "CAG_TKO_effect",
            "CAG_P",
            "MAP2_TKO_effect",
            "MAP2_P",
            "MAP2_minus_CAG_interaction",
            "Interaction_P",
            "Interaction_FDR",
        ]
    ]
    .sort_values(
        "Interaction_P"
    )
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.4g}",
    )
)


print()
print("GROUP MEANS")
print("-" * 105)

print(
    results[
        [
            "State_balance",
            "CAG_Control_mean",
            "CAG_TKO_mean",
            "MAP2_Control_mean",
            "MAP2_TKO_mean",
        ]
    ]
    .to_string(
        index=False,
        float_format=lambda x: f"{x:+.3f}",
    )
)


print()
print(
    "Interpretation:"
)

print(
    "  MITOTIC_MINUS_REPLICATION < 0 means replication-state "
    "activation exceeds mitotic-execution activation."
)

print(
    "  A negative MAP2-vs-CAG interaction means this imbalance "
    "is stronger after MAP2-stage Rb-family loss."
)

print()
print(
    "Outputs written to:"
)
print(OUTDIR)