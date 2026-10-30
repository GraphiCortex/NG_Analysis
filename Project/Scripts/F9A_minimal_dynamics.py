from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# F9A — Minimal nonlinear dynamical model
#
# Question:
#
# Can the CAG-vs-MAP2 mitotic/replication imbalance be explained
# by a SHARED nonlinear regulatory architecture?
#
# Or does MAP2 require an additional context-specific reduction
# in mitotic competence?
#
#
# Dynamical architecture:
#
#   dR/dt = H(E)    - R
#   dQ/dt = H(R)    - Q
#   dM/dt = H(E*R)  - M
#
# with
#
#   H(x) = x^n / (K^n + x^n)
#
#
# We currently fit only steady-state RESPONSE RATIOS because
# static transcriptomics cannot identify kinetic time constants.
#
# H0:
#   same K,n in CAG and MAP2
#
# H1:
#   same architecture, but MAP2 mitotic output receives an
#   additional multiplicative competence factor gamma.
#
# ============================================================


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
    / "F9A_MinimalDynamics"
)

OUT.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 1. Load empirical sample-level state scores
# ============================================================

df = pd.read_csv(INPUT)


MODULES = {
    "E": "E2F_CELL_CYCLE",
    "R": "REPLICATION_ATR",
    "Q": "DNA_REPAIR_CHECKPOINT",
    "M": "MITOTIC_EXECUTION",
}


# ============================================================
# 2. Group-effect helper
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

    return (
        tko.mean()
        -
        ctrl.mean()
    )


effects = {}

for short, module in MODULES.items():

    effects[short] = {
        "CAG": effect(
            df,
            module,
            "CAG",
        ),

        "MAP2": effect(
            df,
            module,
            "MAP2",
        ),
    }


# ============================================================
# 3. Empirical ratios
#
# CAG is reference = 1.
# ============================================================

ratios = {}

for variable in MODULES:

    ratios[variable] = (
        effects[variable]["MAP2"]
        /
        effects[variable]["CAG"]
    )


print()
print("=" * 100)
print("F9A — MINIMAL NONLINEAR DYNAMICS")
print("=" * 100)

print()
print("EMPIRICAL TKO EFFECTS")
print("-" * 100)

for v in MODULES:

    print(
        f"{v}: "
        f"CAG={effects[v]['CAG']:+.4f}  "
        f"MAP2={effects[v]['MAP2']:+.4f}  "
        f"MAP2/CAG={ratios[v]:.4f}"
    )


# ============================================================
# 4. Bootstrap empirical uncertainty
#
# Resample biological replicates WITHIN each of the four groups.
# ============================================================

rng = np.random.default_rng(
    20261022
)

N_BOOT = 20000


groups = {
    group: df[
        df["Group"] == group
    ].copy()
    for group in [
        "CAG_Control",
        "CAG_TKO",
        "MAP2_Control",
        "MAP2_TKO",
    ]
}


boot_ratios = {
    v: []
    for v in MODULES
}


for _ in range(N_BOOT):

    sampled = {}

    for name, sub in groups.items():

        idx = rng.integers(
            0,
            len(sub),
            size=len(sub),
        )

        sampled[name] = (
            sub.iloc[idx]
        )


    for short, module in MODULES.items():

        cag_effect = (
            sampled["CAG_TKO"][module].mean()
            -
            sampled["CAG_Control"][module].mean()
        )

        map2_effect = (
            sampled["MAP2_TKO"][module].mean()
            -
            sampled["MAP2_Control"][module].mean()
        )


        # Ratios become unstable only if denominator approaches zero.
        # That is not expected for these strongly activated backbone
        # programs, but guard explicitly anyway.

        if abs(cag_effect) > 1e-8:

            boot_ratios[short].append(
                map2_effect
                /
                cag_effect
            )


bootstrap_summary = []


for variable in MODULES:

    arr = np.asarray(
        boot_ratios[variable],
        dtype=float,
    )

    arr = arr[
        np.isfinite(arr)
    ]

    bootstrap_summary.append({

        "Variable":
            variable,

        "Observed_ratio":
            ratios[variable],

        "Bootstrap_median":
            np.median(arr),

        "Bootstrap_SD":
            np.std(
                arr,
                ddof=1,
            ),

        "CI025":
            np.quantile(
                arr,
                0.025,
            ),

        "CI975":
            np.quantile(
                arr,
                0.975,
            ),
    })


boot = pd.DataFrame(
    bootstrap_summary
)

boot.to_csv(
    OUT
    / "F9A_bootstrap_response_ratios.csv",
    index=False,
)


print()
print()
print("BOOTSTRAP RESPONSE RATIOS")
print("-" * 100)

print(
    boot.to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}",
    )
)


# ============================================================
# 5. Data targets used by model
#
# E is treated as an empirically measured INPUT.
#
# We ask whether the observed E attenuation is sufficient to
# explain R, Q and M under a shared nonlinear architecture.
# ============================================================

E_RATIO = ratios["E"]

targets = np.array([
    ratios["R"],
    ratios["Q"],
    ratios["M"],
])


# Bootstrap uncertainty for standardized residuals

boot_lookup = (
    boot
    .set_index("Variable")
)

target_sd = np.array([
    boot_lookup.loc["R", "Bootstrap_SD"],
    boot_lookup.loc["Q", "Bootstrap_SD"],
    boot_lookup.loc["M", "Bootstrap_SD"],
])


# Avoid pathological zero uncertainty

target_sd = np.maximum(
    target_sd,
    1e-6,
)


# ============================================================
# 6. Normalized Hill response
#
# Normalization ensures:
#
#   normalized_hill(1) = 1
#
# so CAG is the reference state.
# ============================================================

def hill(x, K, n):

    x = np.maximum(
        x,
        1e-12,
    )

    return (
        x ** n
        /
        (
            K ** n
            +
            x ** n
        )
    )


def normalized_hill(
    x,
    K,
    n,
):

    return (
        hill(
            x,
            K,
            n,
        )
        /
        hill(
            1.0,
            K,
            n,
        )
    )


# ============================================================
# 7. Model H0
#
# Shared nonlinear network.
#
# E_MAP2 / E_CAG is empirical input.
#
# R:
#   H(E)
#
# Q:
#   H(R)
#
# M:
#   H(E * R)
#
# The SAME K and Hill coefficient n are used throughout this
# deliberately minimal first model.
# ============================================================

def shared_model(
    K,
    n,
):

    rR = normalized_hill(
        E_RATIO,
        K,
        n,
    )

    rQ = normalized_hill(
        rR,
        K,
        n,
    )

    rM = normalized_hill(
        E_RATIO * rR,
        K,
        n,
    )

    return np.array([
        rR,
        rQ,
        rM,
    ])


# ============================================================
# 8. Grid search H0
#
# K and n are intentionally scanned broadly rather than
# optimized to a single arbitrary initialization.
# ============================================================

K_grid = np.logspace(
    np.log10(0.05),
    np.log10(5.0),
    500,
)

n_grid = np.linspace(
    0.5,
    6.0,
    350,
)


h0_rows = []

best_h0 = None


for K in K_grid:

    for n in n_grid:

        pred = shared_model(
            K,
            n,
        )

        standardized = (
            pred
            -
            targets
        ) / target_sd

        chi2 = float(
            np.sum(
                standardized ** 2
            )
        )


        row = {
            "K":
                K,

            "n":
                n,

            "Pred_R":
                pred[0],

            "Pred_Q":
                pred[1],

            "Pred_M":
                pred[2],

            "Chi2":
                chi2,
        }


        if (
            best_h0 is None
            or
            chi2
            <
            best_h0["Chi2"]
        ):
            best_h0 = row


# ============================================================
# 9. Model H1
#
# Same nonlinear network, but MAP2 receives an additional
# mitotic-competence multiplier gamma:
#
#   M_MAP2 = gamma * M_shared
#
# gamma < 1:
#   MAP2-specific loss of mitotic competence.
#
# gamma ~ 1:
#   no extra context-specific competence term is needed.
# ============================================================

best_h1 = None


for K in K_grid:

    for n in n_grid:

        base = shared_model(
            K,
            n,
        )

        base_M = base[2]


        if base_M <= 0:
            continue


        # Best gamma for M given K,n.
        gamma = (
            targets[2]
            /
            base_M
        )

        gamma = np.clip(
            gamma,
            0.2,
            1.5,
        )


        pred = base.copy()

        pred[2] = (
            gamma
            *
            base_M
        )


        standardized = (
            pred
            -
            targets
        ) / target_sd

        chi2 = float(
            np.sum(
                standardized ** 2
            )
        )


        row = {
            "K":
                K,

            "n":
                n,

            "gamma_MAP2":
                gamma,

            "Pred_R":
                pred[0],

            "Pred_Q":
                pred[1],

            "Pred_M":
                pred[2],

            "Chi2":
                chi2,
        }


        if (
            best_h1 is None
            or
            chi2
            <
            best_h1["Chi2"]
        ):
            best_h1 = row


# ============================================================
# 10. Parameter-free linear benchmark
#
# No nonlinear amplification:
#
#   R ratio = E ratio
#   Q ratio = R ratio
#   M ratio = E * R
# ============================================================

linear_pred = np.array([
    E_RATIO,
    E_RATIO,
    E_RATIO ** 2,
])

linear_chi2 = float(
    np.sum(
        (
            (
                linear_pred
                -
                targets
            )
            /
            target_sd
        ) ** 2
    )
)


# ============================================================
# 11. Check H0 against bootstrap intervals
# ============================================================

intervals = {
    variable: (
        boot_lookup.loc[
            variable,
            "CI025"
        ],
        boot_lookup.loc[
            variable,
            "CI975"
        ],
    )
    for variable in [
        "R",
        "Q",
        "M",
    ]
}


h0_inside = {}

for variable, prediction in zip(
    [
        "R",
        "Q",
        "M",
    ],
    [
        best_h0["Pred_R"],
        best_h0["Pred_Q"],
        best_h0["Pred_M"],
    ],
):

    lo, hi = intervals[
        variable
    ]

    h0_inside[
        variable
    ] = (
        lo
        <=
        prediction
        <=
        hi
    )


# ============================================================
# 12. Export model comparison
# ============================================================

comparison = pd.DataFrame([

    {
        "Model":
            "LINEAR_SHARED",

        "K":
            np.nan,

        "n":
            np.nan,

        "gamma_MAP2":
            1.0,

        "Pred_R":
            linear_pred[0],

        "Pred_Q":
            linear_pred[1],

        "Pred_M":
            linear_pred[2],

        "Chi2":
            linear_chi2,
    },

    {
        "Model":
            "H0_SHARED_NONLINEAR",

        "K":
            best_h0["K"],

        "n":
            best_h0["n"],

        "gamma_MAP2":
            1.0,

        "Pred_R":
            best_h0["Pred_R"],

        "Pred_Q":
            best_h0["Pred_Q"],

        "Pred_M":
            best_h0["Pred_M"],

        "Chi2":
            best_h0["Chi2"],
    },

    {
        "Model":
            "H1_CONTEXT_MITOTIC_COMPETENCE",

        "K":
            best_h1["K"],

        "n":
            best_h1["n"],

        "gamma_MAP2":
            best_h1["gamma_MAP2"],

        "Pred_R":
            best_h1["Pred_R"],

        "Pred_Q":
            best_h1["Pred_Q"],

        "Pred_M":
            best_h1["Pred_M"],

        "Chi2":
            best_h1["Chi2"],
    },
])


comparison.to_csv(
    OUT
    / "F9A_model_comparison.csv",
    index=False,
)


# ============================================================
# 13. Console output
# ============================================================

print()
print()
print("=" * 100)
print("MODEL COMPARISON")
print("=" * 100)


print()
print("Observed MAP2/CAG ratios:")
print(
    f"    E = {E_RATIO:.4f}"
)

print(
    f"    R = {targets[0]:.4f}"
)

print(
    f"    Q = {targets[1]:.4f}"
)

print(
    f"    M = {targets[2]:.4f}"
)


print()
print(
    comparison.to_string(
        index=False,
        float_format=lambda x: f"{x:.5f}",
    )
)


print()
print("H0 predictions inside empirical 95% bootstrap intervals:")

for variable in [
    "R",
    "Q",
    "M",
]:

    lo, hi = intervals[
        variable
    ]

    pred = best_h0[
        f"Pred_{variable}"
    ]

    print(
        f"    {variable}: "
        f"prediction={pred:.4f}  "
        f"CI=[{lo:.4f}, {hi:.4f}]  "
        f"inside={h0_inside[variable]}"
    )


print()
print("Interpretation guide:")
print(
    "  If H0 reproduces R/Q/M within empirical uncertainty, "
    "the MAP2 imbalance does NOT require an explicit "
    "context-specific mitotic defect."
)

print(
    "  In that case, modest upstream attenuation plus shared "
    "nonlinear dynamics can generate the large mitotic deficit."
)

print(
    "  If H0 fails but H1 fits and gamma_MAP2 << 1, an explicit "
    "context-dependent mitotic-competence term is supported."
)


print()
print("Outputs written to:")
print(OUT)