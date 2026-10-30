from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize


ROOT = Path(".")

INPUT = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F8C_StatePrograms"
    / "F8C_oshikawa_sample_state_scores.csv"
)

OUT = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F9B_HeldoutMitosis"
)

OUT.mkdir(parents=True, exist_ok=True)


MODULES = {
    "E": "E2F_CELL_CYCLE",
    "R": "REPLICATION_ATR",
    "Q": "DNA_REPAIR_CHECKPOINT",
    "M": "MITOTIC_EXECUTION",
}


df = pd.read_csv(INPUT)


# ============================================================
# Empirical effect ratios
# ============================================================

def effect(data, module, context):

    ctrl = data.loc[
        data["Group"] == f"{context}_Control",
        module,
    ].astype(float)

    tko = data.loc[
        data["Group"] == f"{context}_TKO",
        module,
    ].astype(float)

    return tko.mean() - ctrl.mean()


def get_ratios(data):

    out = {}

    for short, module in MODULES.items():

        cag = effect(
            data,
            module,
            "CAG",
        )

        map2 = effect(
            data,
            module,
            "MAP2",
        )

        if abs(cag) < 1e-8:
            return None

        out[short] = map2 / cag

    return out


obs = get_ratios(df)


# ============================================================
# Hill model
# ============================================================

def hill(x, K, n):

    x = max(float(x), 1e-12)

    return (
        x ** n
        /
        (K ** n + x ** n)
    )


def normalized_hill(x, K, n):

    return (
        hill(x, K, n)
        /
        hill(1.0, K, n)
    )


def predict(E, K, n):

    R = normalized_hill(
        E,
        K,
        n,
    )

    Q = normalized_hill(
        R,
        K,
        n,
    )

    M = normalized_hill(
        E * R,
        K,
        n,
    )

    return R, Q, M


# ============================================================
# Fit ONLY R and Q.
#
# M is deliberately excluded from the objective.
# ============================================================

def fit_upstream(E, R_target, Q_target):

    def objective(theta):

        logK, logn = theta

        K = np.exp(logK)
        n = np.exp(logn)

        R, Q, _ = predict(
            E,
            K,
            n,
        )

        return (
            (R - R_target) ** 2
            +
            (Q - Q_target) ** 2
        )


    starts = [
        (np.log(0.1), np.log(0.75)),
        (np.log(0.5), np.log(1.0)),
        (np.log(1.0), np.log(2.0)),
        (np.log(2.0), np.log(3.0)),
        (np.log(5.0), np.log(1.0)),
    ]


    best = None


    for start in starts:

        res = minimize(
            objective,
            x0=np.array(start),
            method="L-BFGS-B",
            bounds=[
                (
                    np.log(0.02),
                    np.log(10.0),
                ),
                (
                    np.log(0.25),
                    np.log(8.0),
                ),
            ],
        )


        if (
            best is None
            or
            res.fun < best.fun
        ):
            best = res


    K = np.exp(best.x[0])
    n = np.exp(best.x[1])

    R, Q, M = predict(
        E,
        K,
        n,
    )

    return {
        "K": K,
        "n": n,
        "Pred_R": R,
        "Pred_Q": Q,
        "Pred_M_HELDOUT": M,
        "Upstream_SSE": best.fun,
    }


# ============================================================
# Observed-data fit
# ============================================================

fit = fit_upstream(
    obs["E"],
    obs["R"],
    obs["Q"],
)


linear_M = (
    obs["E"] ** 2
)

semiempirical_M = (
    obs["E"]
    *
    obs["R"]
)


# ============================================================
# Bootstrap
#
# Resample replicates within each of the four biological groups.
# Refit R/Q each time, then predict held-out M.
# ============================================================

rng = np.random.default_rng(
    20261022
)

N_BOOT = 5000


group_names = [
    "CAG_Control",
    "CAG_TKO",
    "MAP2_Control",
    "MAP2_TKO",
]


groups = {
    g: df[
        df["Group"] == g
    ].copy()
    for g in group_names
}


rows = []


for b in range(N_BOOT):

    pieces = []

    for name, sub in groups.items():

        idx = rng.integers(
            0,
            len(sub),
            size=len(sub),
        )

        temp = sub.iloc[idx].copy()

        pieces.append(temp)


    boot_df = pd.concat(
        pieces,
        ignore_index=True,
    )


    ratios = get_ratios(
        boot_df
    )


    if ratios is None:
        continue


    # Keep only biologically sensible positive backbone ratios.
    if (
        ratios["E"] <= 0
        or
        ratios["R"] <= 0
        or
        ratios["Q"] <= 0
        or
        ratios["M"] <= 0
    ):
        continue


    model = fit_upstream(
        ratios["E"],
        ratios["R"],
        ratios["Q"],
    )


    pred_nonlinear = (
        model[
            "Pred_M_HELDOUT"
        ]
    )

    pred_linear = (
        ratios["E"] ** 2
    )

    pred_semi = (
        ratios["E"]
        *
        ratios["R"]
    )


    rows.append({

        "E":
            ratios["E"],

        "R":
            ratios["R"],

        "Q":
            ratios["Q"],

        "Observed_M":
            ratios["M"],

        "K":
            model["K"],

        "n":
            model["n"],

        "Pred_M_nonlinear":
            pred_nonlinear,

        "Pred_M_linear":
            pred_linear,

        "Pred_M_E_times_R":
            pred_semi,

        "Error_nonlinear":
            pred_nonlinear
            -
            ratios["M"],

        "Error_linear":
            pred_linear
            -
            ratios["M"],

        "Error_E_times_R":
            pred_semi
            -
            ratios["M"],

        "gamma_residual":
            (
                ratios["M"]
                /
                pred_nonlinear
            ),
    })


boot = pd.DataFrame(
    rows
)


boot.to_csv(
    OUT
    / "F9B_bootstrap_heldout_predictions.csv",
    index=False,
)


# ============================================================
# Summaries
# ============================================================

def summarize(values):

    x = np.asarray(
        values,
        dtype=float,
    )

    return {
        "median": np.median(x),
        "lo": np.quantile(x, 0.025),
        "hi": np.quantile(x, 0.975),
    }


pred_nonlin = summarize(
    boot["Pred_M_nonlinear"]
)

err_nonlin = summarize(
    boot["Error_nonlinear"]
)

err_linear = summarize(
    boot["Error_linear"]
)

err_semi = summarize(
    boot["Error_E_times_R"]
)

gamma = summarize(
    boot["gamma_residual"]
)


better_than_linear = np.mean(
    np.abs(
        boot["Error_nonlinear"]
    )
    <
    np.abs(
        boot["Error_linear"]
    )
)


# ============================================================
# Export compact summary
# ============================================================

summary = pd.DataFrame([

    {
        "Model":
            "NONLINEAR_UPSTREAM_FIT",

        "Observed_M":
            obs["M"],

        "Predicted_M":
            fit["Pred_M_HELDOUT"],

        "Absolute_error":
            abs(
                fit["Pred_M_HELDOUT"]
                -
                obs["M"]
            ),
    },

    {
        "Model":
            "LINEAR_SHARED",

        "Observed_M":
            obs["M"],

        "Predicted_M":
            linear_M,

        "Absolute_error":
            abs(
                linear_M
                -
                obs["M"]
            ),
    },

    {
        "Model":
            "E_TIMES_OBSERVED_R",

        "Observed_M":
            obs["M"],

        "Predicted_M":
            semiempirical_M,

        "Absolute_error":
            abs(
                semiempirical_M
                -
                obs["M"]
            ),
    },
])


summary.to_csv(
    OUT
    / "F9B_model_comparison.csv",
    index=False,
)


# ============================================================
# Console
# ============================================================

print()
print("=" * 100)
print("F9B — HELD-OUT MITOTIC PREDICTION")
print("=" * 100)

print()
print("OBSERVED RATIOS")
print("-" * 100)

for k in [
    "E",
    "R",
    "Q",
    "M",
]:
    print(
        f"{k} = {obs[k]:.4f}"
    )


print()
print("UPSTREAM-ONLY NONLINEAR FIT")
print("-" * 100)

print(
    f"K = {fit['K']:.4f}"
)

print(
    f"n = {fit['n']:.4f}"
)

print(
    f"Predicted R = {fit['Pred_R']:.4f}"
)

print(
    f"Observed  R = {obs['R']:.4f}"
)

print(
    f"Predicted Q = {fit['Pred_Q']:.4f}"
)

print(
    f"Observed  Q = {obs['Q']:.4f}"
)

print()
print(
    f"HELD-OUT predicted M = "
    f"{fit['Pred_M_HELDOUT']:.4f}"
)

print(
    f"Observed M           = "
    f"{obs['M']:.4f}"
)


print()
print("BENCHMARKS")
print("-" * 100)

print(
    summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}",
    )
)


print()
print("BOOTSTRAP")
print("-" * 100)

print(
    "Valid bootstrap fits:",
    len(boot),
)

print(
    "Nonlinear predicted M: "
    f"median={pred_nonlin['median']:.4f}, "
    f"95%=[{pred_nonlin['lo']:.4f}, "
    f"{pred_nonlin['hi']:.4f}]"
)

print(
    "Nonlinear prediction error: "
    f"median={err_nonlin['median']:+.4f}, "
    f"95%=[{err_nonlin['lo']:+.4f}, "
    f"{err_nonlin['hi']:+.4f}]"
)

print(
    "Linear prediction error: "
    f"median={err_linear['median']:+.4f}, "
    f"95%=[{err_linear['lo']:+.4f}, "
    f"{err_linear['hi']:+.4f}]"
)

print(
    "E*R prediction error: "
    f"median={err_semi['median']:+.4f}, "
    f"95%=[{err_semi['lo']:+.4f}, "
    f"{err_semi['hi']:+.4f}]"
)

print(
    "Residual competence gamma = Mobs/Mpred: "
    f"median={gamma['median']:.4f}, "
    f"95%=[{gamma['lo']:.4f}, "
    f"{gamma['hi']:.4f}]"
)

print(
    "Fraction bootstrap samples where nonlinear "
    "beats linear:",
    f"{better_than_linear:.3f}"
)


print()
print("Outputs written to:")
print(OUT)

print()
print(
    "[NEXT] If held-out prediction succeeds, freeze the "
    "Oshikawa shared backbone and move to the adult-NSC "
    "residual/fate-routing problem."
)