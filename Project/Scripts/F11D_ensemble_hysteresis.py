from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp


# ============================================================
# F11D — Ensemble hysteresis and switching thresholds
#
# Goals
# -----
# 1. Re-run slow forward/reverse Rb-loss ramps across EVERY
#    F11B parameter set classified as bistable.
#
# 2. Estimate:
#
#       rho_up   = low-D -> high-D switching threshold
#       rho_down = high-D -> low-D switching threshold
#
# 3. Quantify:
#
#       hysteresis width = rho_up - rho_down
#
# 4. Determine whether history dependence is a robust property
#    of the admissible bistable ensemble or merely a feature
#    of the selected representative model.
#
# IMPORTANT
# ---------
# These thresholds are dimensionless model quantities.
# They are not yet experimentally calibrated doses/times.
# ============================================================


ROOT = Path(".")

F11B = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F11B_AttractorEnsemble"
)

OUT = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F11D_EnsembleHysteresis"
)

FIG = OUT / "figures"

OUT.mkdir(
    parents=True,
    exist_ok=True,
)

FIG.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# 1. Load ensemble
# ============================================================

ensemble = pd.read_csv(
    F11B
    / "F11B_parameter_ensemble.csv"
)


def as_bool(x):

    if isinstance(x, bool):
        return x

    return (
        str(x)
        .strip()
        .lower()
        ==
        "true"
    )


ensemble["Bistable"] = (
    ensemble["Bistable"]
    .apply(as_bool)
)


bistable = (
    ensemble[
        ensemble["Bistable"]
    ]
    .copy()
    .reset_index(drop=True)
)


if bistable.empty:

    raise RuntimeError(
        "No bistable parameter sets were found."
    )


print()
print("=" * 105)
print("F11D — ENSEMBLE HYSTERESIS")
print("=" * 105)

print()
print(
    "Bistable models to test:",
    len(bistable),
)


# ============================================================
# 2. Parameter names
# ============================================================

PARAMETER_COLUMNS = [

    "bU",
    "bD",
    "bT",

    "aR",

    "aUU",
    "aDU",
    "aUD",
    "aDD",
    "aTT",

    "KUU",
    "KDU",
    "KUD",
    "KDD",
    "KTT",
    "KDT",

    "nUU",
    "nDU",
    "nUD",
    "nDD",
    "nTT",
    "nDT",
]


# ============================================================
# 3. Regulatory functions
# ============================================================

def hill(x, K, n):

    x = np.maximum(
        x,
        0.0,
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


def repression(x, K, n):

    x = np.maximum(
        x,
        0.0,
    )

    return (
        K ** n
        /
        (
            K ** n
            +
            x ** n
        )
    )


# ============================================================
# 4. Model
# ============================================================

def rhs(t, state, p, rho):

    U, D, T = state


    dU = (
        p["bU"]
        +
        p["aR"] * rho
        +
        p["aUU"]
        * hill(
            U,
            p["KUU"],
            p["nUU"],
        )
        +
        p["aDU"]
        * hill(
            D,
            p["KDU"],
            p["nDU"],
        )
        -
        U
    )


    dD = (
        p["bD"]
        +
        p["aUD"]
        * hill(
            U,
            p["KUD"],
            p["nUD"],
        )
        +
        p["aDD"]
        * hill(
            D,
            p["KDD"],
            p["nDD"],
        )
        -
        D
    )


    dT = (
        p["bT"]
        +
        p["aTT"]
        * hill(
            T,
            p["KTT"],
            p["nTT"],
        )
        * repression(
            D,
            p["KDT"],
            p["nDT"],
        )
        -
        T
    )


    return np.array([
        dU,
        dD,
        dT,
    ])


# ============================================================
# 5. Integrate to equilibrium
# ============================================================

def evolve(
    state,
    p,
    rho,
    duration=70,
):

    sol = solve_ivp(

        lambda t, x:
            rhs(
                t,
                x,
                p,
                rho,
            ),

        (
            0,
            duration,
        ),

        state,

        rtol=1e-7,
        atol=1e-9,

        max_step=0.3,
    )


    return sol.y[
        :,
        -1
    ]


# ============================================================
# 6. rho scan
# ============================================================

rho_forward = np.linspace(
    0.0,
    1.30,
    201,
)

rho_reverse = (
    rho_forward[::-1]
)


# ============================================================
# 7. Threshold helper
#
# We identify the largest discrete state jump during each
# quasistatic continuation.
# ============================================================

def switch_threshold(
    rho,
    values,
):

    rho = np.asarray(
        rho,
        dtype=float,
    )

    values = np.asarray(
        values,
        dtype=float,
    )


    changes = np.abs(
        np.diff(
            values
        )
    )


    if len(changes) == 0:

        return (
            np.nan,
            np.nan,
        )


    i = int(
        np.argmax(
            changes
        )
    )


    threshold = (
        rho[i]
        +
        rho[i + 1]
    ) / 2


    jump = changes[i]


    return (
        threshold,
        jump,
    )


# ============================================================
# 8. Run every bistable model
# ============================================================

summary_rows = []

scan_rows = []


for _, row in bistable.iterrows():

    parameter_id = int(
        row[
            "Parameter_set"
        ]
    )


    p = {
        key: float(
            row[key]
        )
        for key in PARAMETER_COLUMNS
    }


    # --------------------------------------------------------
    # Forward ramp
    #
    # Start from previously identified control sink.
    # --------------------------------------------------------

    state = np.array([

        row[
            "Control_U"
        ],

        row[
            "Control_D"
        ],

        row[
            "Control_T"
        ],
    ], dtype=float)


    forward_states = []


    for rho in rho_forward:

        state = evolve(
            state,
            p,
            rho,
        )

        forward_states.append(
            state.copy()
        )


    forward_states = np.asarray(
        forward_states
    )


    # --------------------------------------------------------
    # Reverse ramp
    #
    # Start from high-drive endpoint of forward continuation.
    # --------------------------------------------------------

    state = (
        forward_states[
            -1
        ]
        .copy()
    )


    reverse_states = []


    for rho in rho_reverse:

        state = evolve(
            state,
            p,
            rho,
        )

        reverse_states.append(
            state.copy()
        )


    reverse_states = np.asarray(
        reverse_states
    )


    # --------------------------------------------------------
    # Switching thresholds
    # --------------------------------------------------------

    rho_up, forward_jump = (
        switch_threshold(
            rho_forward,
            forward_states[
                :,
                1
            ],
        )
    )


    rho_down, reverse_jump = (
        switch_threshold(
            rho_reverse,
            reverse_states[
                :,
                1
            ],
        )
    )


    hysteresis_width = (
        rho_up
        -
        rho_down
    )


    # --------------------------------------------------------
    # Compare forward/reverse states at identical rho
    # --------------------------------------------------------

    reverse_ascending = (
        reverse_states[
            ::-1
        ]
    )


    delta_D = (
        reverse_ascending[
            :,
            1
        ]
        -
        forward_states[
            :,
            1
        ]
    )


    delta_T = (
        reverse_ascending[
            :,
            2
        ]
        -
        forward_states[
            :,
            2
        ]
    )


    max_D_gap = np.max(
        np.abs(
            delta_D
        )
    )


    max_T_gap = np.max(
        np.abs(
            delta_T
        )
    )


    strong_hysteresis = (
        hysteresis_width
        >
        0.03
        and
        max_D_gap
        >
        0.20
        and
        forward_jump
        >
        0.20
        and
        reverse_jump
        >
        0.20
    )


    summary_rows.append({

        "Parameter_set":
            parameter_id,

        "rho_up":
            rho_up,

        "rho_down":
            rho_down,

        "Hysteresis_width":
            hysteresis_width,

        "Forward_D_jump":
            forward_jump,

        "Reverse_D_jump":
            reverse_jump,

        "Max_D_gap":
            max_D_gap,

        "Max_T_gap":
            max_T_gap,

        "Strong_hysteresis":
            strong_hysteresis,
    })


    # --------------------------------------------------------
    # Save scan
    # --------------------------------------------------------

    for i, rho in enumerate(
        rho_forward
    ):

        scan_rows.append({

            "Parameter_set":
                parameter_id,

            "rho":
                rho,

            "D_forward":
                forward_states[
                    i,
                    1
                ],

            "T_forward":
                forward_states[
                    i,
                    2
                ],

            "D_reverse":
                reverse_ascending[
                    i,
                    1
                ],

            "T_reverse":
                reverse_ascending[
                    i,
                    2
                ],
        })


summary = pd.DataFrame(
    summary_rows
)

scans = pd.DataFrame(
    scan_rows
)


summary.to_csv(
    OUT
    / "F11D_hysteresis_summary.csv",
    index=False,
)

scans.to_csv(
    OUT
    / "F11D_all_hysteresis_scans.csv",
    index=False,
)


# ============================================================
# 9. Ensemble statistics
# ============================================================

strong = summary[
    summary[
        "Strong_hysteresis"
    ]
]


fraction_strong = (
    len(strong)
    /
    len(summary)
)


print()
print("ENSEMBLE HYSTERESIS RESULTS")
print("-" * 105)

print(
    "Bistable models tested:",
    len(summary),
)

print(
    "Strongly hysteretic models:",
    len(strong),
)

print(
    "Fraction of bistable models with strong hysteresis:",
    f"{fraction_strong:.3f}",
)


print()
print("ALL BISTABLE MODELS")

print(
    summary.sort_values(
        "Hysteresis_width",
        ascending=False,
    )
    .to_string(
        index=False,
        float_format=lambda x: f"{x:.4f}",
    )
)


# ============================================================
# 10. Robust threshold summaries
# ============================================================

valid = summary[
    (
        summary[
            "rho_up"
        ]
        .notna()
    )
    &
    (
        summary[
            "rho_down"
        ]
        .notna()
    )
]


print()
print("SWITCHING-THRESHOLD DISTRIBUTION")
print("-" * 105)


for col in [
    "rho_up",
    "rho_down",
    "Hysteresis_width",
    "Max_D_gap",
]:

    vals = valid[
        col
    ].values


    print(
        f"{col:20s} "
        f"median={np.median(vals):.4f}  "
        f"IQR=["
        f"{np.quantile(vals, 0.25):.4f}, "
        f"{np.quantile(vals, 0.75):.4f}]"
    )


# ============================================================
# 11. Threshold scatter
# ============================================================

plt.figure(
    figsize=(7, 7)
)

plt.scatter(
    valid[
        "rho_down"
    ],
    valid[
        "rho_up"
    ],
    s=55,
)

maximum = max(
    valid[
        [
            "rho_down",
            "rho_up",
        ]
    ].max()
)

plt.plot(
    [
        0,
        maximum,
    ],
    [
        0,
        maximum,
    ],
    linestyle="--",
)

plt.xlabel(
    "Reverse switching threshold rho_down"
)

plt.ylabel(
    "Forward switching threshold rho_up"
)

plt.title(
    "Hysteresis thresholds across bistable models"
)

plt.tight_layout()

plt.savefig(
    FIG
    / "01_switching_thresholds.png",
    dpi=250,
)

plt.close()


# ============================================================
# 12. Hysteresis-width distribution
# ============================================================

plt.figure(
    figsize=(8, 6)
)

plt.hist(
    valid[
        "Hysteresis_width"
    ],
    bins=min(
        12,
        len(valid),
    ),
)

plt.xlabel(
    "Hysteresis width"
)

plt.ylabel(
    "Number of models"
)

plt.title(
    "Distribution of history dependence"
)

plt.tight_layout()

plt.savefig(
    FIG
    / "02_hysteresis_width_distribution.png",
    dpi=250,
)

plt.close()


# ============================================================
# 13. Ensemble-average hysteresis curves
# ============================================================

grouped = (
    scans.groupby(
        "rho"
    )
    .agg({

        "D_forward":
            [
                "median",
                lambda x:
                    np.quantile(
                        x,
                        0.25
                    ),
                lambda x:
                    np.quantile(
                        x,
                        0.75
                    ),
            ],

        "D_reverse":
            [
                "median",
                lambda x:
                    np.quantile(
                        x,
                        0.25
                    ),
                lambda x:
                    np.quantile(
                        x,
                        0.75
                    ),
            ],
    })
)


grouped.columns = [

    "D_forward_median",
    "D_forward_q25",
    "D_forward_q75",

    "D_reverse_median",
    "D_reverse_q25",
    "D_reverse_q75",
]


grouped = (
    grouped.reset_index()
)


grouped.to_csv(
    OUT
    / "F11D_ensemble_hysteresis_curve.csv",
    index=False,
)


plt.figure(
    figsize=(8, 6)
)


plt.plot(
    grouped["rho"],
    grouped[
        "D_forward_median"
    ],
    label="Increasing drive",
)

plt.fill_between(
    grouped["rho"],
    grouped[
        "D_forward_q25"
    ],
    grouped[
        "D_forward_q75"
    ],
    alpha=0.2,
)


plt.plot(
    grouped["rho"],
    grouped[
        "D_reverse_median"
    ],
    label="Decreasing drive",
)

plt.fill_between(
    grouped["rho"],
    grouped[
        "D_reverse_q25"
    ],
    grouped[
        "D_reverse_q75"
    ],
    alpha=0.2,
)


plt.xlabel(
    "Rb-loss drive rho"
)

plt.ylabel(
    "Death-state activity D"
)

plt.title(
    "Ensemble history dependence"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    FIG
    / "03_ensemble_hysteresis.png",
    dpi=250,
)

plt.close()


# ============================================================
# 14. Final interpretation
# ============================================================

print()
print("INTERPRETATION")
print("-" * 105)


if fraction_strong >= 0.5:

    print(
        "Hysteresis is common across the admissible "
        "bistable ensemble."
    )

    print(
        "History dependence is therefore a robust dynamical "
        "consequence of a substantial region of parameter space."
    )

elif fraction_strong >= 0.2:

    print(
        "Hysteresis occurs in a meaningful subset of the "
        "admissible bistable ensemble."
    )

    print(
        "It is structurally plausible but parameter-dependent."
    )

else:

    print(
        "Strong hysteresis is rare across the bistable ensemble."
    )

    print(
        "The representative model should therefore not be "
        "treated as a robust hysteretic prediction."
    )


print()
print("Figures written to:")
print(FIG)

print()
print(
    "[NEXT] Use the threshold distribution to screen transient "
    "and sustained perturbations of U->D, D->D, D->U, U->U "
    "and rho for their ability to switch the high-D attractor "
    "back into the low-D basin."
)