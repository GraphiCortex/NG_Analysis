from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp
from scipy.optimize import root


# ============================================================
# F11C — Basin geometry and hysteresis
#
# Uses the representative F11B model.
#
# Goals:
#
# 1. Map initial conditions into their final attractor.
# 2. Visualize the separatrix between low-D and high-D sinks.
# 3. Slowly increase and decrease Rb-loss drive rho.
# 4. Test for hysteresis / history dependence.
#
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
    / "F11C_BasinsHysteresis"
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
# 1. Load representative parameters
# ============================================================

p = (
    pd.read_csv(
        F11B
        / "F11B_representative_parameters.csv"
    )
    .iloc[0]
    .to_dict()
)


# ============================================================
# 2. Find representative rho automatically
# ============================================================

ensemble = pd.read_csv(
    F11B
    / "F11B_parameter_ensemble.csv"
)


parameter_columns = [
    k
    for k in p
    if k in ensemble.columns
]


distance = np.zeros(
    len(ensemble)
)


for col in parameter_columns:

    distance += (
        ensemble[col]
        -
        p[col]
    ) ** 2


rep_row = ensemble.iloc[
    np.argmin(distance)
]


rho_phase = (
    rep_row[
        "Best_bistable_rho"
    ]
)


if pd.isna(rho_phase):

    rho_phase = 0.5


rho_phase = float(
    rho_phase
)


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
# 4. ODE system
# ============================================================

def rhs(t, state, rho):

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
# 5. Find fixed points
# ============================================================

def find_fixed_points(rho):

    guesses = np.linspace(
        0.05,
        2.8,
        6,
    )


    fixed = []


    for U0 in guesses:

        for D0 in guesses:

            for T0 in guesses:

                sol = root(
                    lambda x:
                        rhs(
                            0,
                            x,
                            rho,
                        ),
                    [
                        U0,
                        D0,
                        T0,
                    ],
                )


                if not sol.success:
                    continue


                x = sol.x


                if (
                    np.any(
                        x < -1e-6
                    )
                    or
                    np.any(
                        x > 4
                    )
                ):
                    continue


                if np.linalg.norm(
                    rhs(
                        0,
                        x,
                        rho,
                    )
                ) > 1e-6:

                    continue


                x = np.maximum(
                    x,
                    0,
                )


                if any(
                    np.linalg.norm(
                        x - y
                    )
                    < 1e-3
                    for y in fixed
                ):
                    continue


                fixed.append(
                    x
                )


    fixed.sort(
        key=lambda x:
            x[1]
    )


    return fixed


# ============================================================
# 6. Integrate to equilibrium
# ============================================================

def evolve(
    x0,
    rho,
    duration=100,
):

    sol = solve_ivp(

        lambda t, x:
            rhs(
                t,
                x,
                rho,
            ),

        (
            0,
            duration,
        ),

        x0,

        rtol=1e-8,
        atol=1e-10,

        max_step=0.25,
    )


    return (
        sol.y[:, -1],
        sol,
    )


# ============================================================
# 7. Determine sinks at bistable rho
#
# Rather than reimplementing the Jacobian, identify which fixed
# points attract nearby trajectories.
# ============================================================

fps = find_fixed_points(
    rho_phase
)


stable_points = []


for fp in fps:

    perturbations = [
        np.array([
            1e-3,
            0,
            0,
        ]),
        np.array([
            0,
            1e-3,
            0,
        ]),
        np.array([
            0,
            0,
            1e-3,
        ]),
    ]


    stable = True


    for eps in perturbations:

        final, _ = evolve(
            fp + eps,
            rho_phase,
            duration=50,
        )


        if (
            np.linalg.norm(
                final - fp
            )
            >
            1e-2
        ):

            stable = False
            break


    if stable:

        stable_points.append(
            fp
        )


stable_points = sorted(
    stable_points,
    key=lambda x:
        x[1]
)


if len(stable_points) < 2:

    raise RuntimeError(
        "Representative rho does not contain at least "
        "two stable attractors."
    )


low_sink = stable_points[0]

high_sink = stable_points[-1]


print()
print("=" * 100)
print("F11C — BASINS AND HYSTERESIS")
print("=" * 100)

print()
print(
    f"Representative rho = "
    f"{rho_phase:.4f}"
)

print(
    "Low-D sink:",
    np.round(
        low_sink,
        4,
    )
)

print(
    "High-D sink:",
    np.round(
        high_sink,
        4,
    )
)


# ============================================================
# 8. Basin map
#
# Scan initial U and D.
#
# Start all trajectories from the same intermediate T.
# ============================================================

U_grid = np.linspace(
    0,
    2.5,
    61,
)

D_grid = np.linspace(
    0,
    2.5,
    61,
)


T_initial = 0.60


basin_rows = []


for U0 in U_grid:

    for D0 in D_grid:

        initial = np.array([
            U0,
            D0,
            T_initial,
        ])


        final, _ = evolve(
            initial,
            rho_phase,
            duration=80,
        )


        d_low = np.linalg.norm(
            final - low_sink
        )

        d_high = np.linalg.norm(
            final - high_sink
        )


        basin = (
            0
            if d_low < d_high
            else 1
        )


        basin_rows.append({

            "Initial_U":
                U0,

            "Initial_D":
                D0,

            "Initial_T":
                T_initial,

            "Final_U":
                final[0],

            "Final_D":
                final[1],

            "Final_T":
                final[2],

            "Basin":
                basin,
        })


basins = pd.DataFrame(
    basin_rows
)


basins.to_csv(
    OUT
    / "F11C_basin_map.csv",
    index=False,
)


# ============================================================
# 9. Basin figure
# ============================================================

plt.figure(
    figsize=(8, 7)
)


plt.scatter(
    basins["Initial_U"],
    basins["Initial_D"],
    c=basins["Basin"],
    s=12,
)


plt.scatter(
    low_sink[0],
    low_sink[1],
    marker="o",
    s=160,
    edgecolors="black",
)

plt.scatter(
    high_sink[0],
    high_sink[1],
    marker="o",
    s=160,
    edgecolors="black",
)


plt.xlabel(
    "Initial U"
)

plt.ylabel(
    "Initial D"
)

plt.title(
    f"Basins of attraction at rho={rho_phase:.2f}"
)

plt.tight_layout()

plt.savefig(
    FIG
    / "01_basin_map.png",
    dpi=250,
)

plt.close()


# ============================================================
# 10. Basin fractions
# ============================================================

basin_fraction = (
    basins[
        "Basin"
    ]
    .value_counts(
        normalize=True
    )
    .sort_index()
)


print()
print("BASIN FRACTIONS")
print("-" * 100)

print(
    "Low-D basin:",
    f"{basin_fraction.get(0, 0):.3f}",
)

print(
    "High-D basin:",
    f"{basin_fraction.get(1, 0):.3f}",
)


# ============================================================
# 11. Slow forward rho ramp
#
# Each final state becomes the next initial state.
# This approximates quasistatic continuation.
# ============================================================

rho_forward = np.linspace(
    0,
    1.30,
    161,
)


forward_rows = []


state = np.array([
    0.05,
    0.05,
    0.8,
])


for rho in rho_forward:

    state, _ = evolve(
        state,
        rho,
        duration=60,
    )


    forward_rows.append({

        "Direction":
            "forward",

        "rho":
            rho,

        "U":
            state[0],

        "D":
            state[1],

        "T":
            state[2],
    })


# ============================================================
# 12. Reverse rho ramp
# ============================================================

rho_reverse = rho_forward[
    ::-1
]


reverse_rows = []


# Start from the final high-drive state.

state = np.array(
    forward_rows[-1][
        ["U", "D", "T"]
    ]
    if False
    else [
        forward_rows[-1]["U"],
        forward_rows[-1]["D"],
        forward_rows[-1]["T"],
    ]
)


for rho in rho_reverse:

    state, _ = evolve(
        state,
        rho,
        duration=60,
    )


    reverse_rows.append({

        "Direction":
            "reverse",

        "rho":
            rho,

        "U":
            state[0],

        "D":
            state[1],

        "T":
            state[2],
    })


hysteresis = pd.DataFrame(
    forward_rows
    +
    reverse_rows
)


hysteresis.to_csv(
    OUT
    / "F11C_hysteresis_scan.csv",
    index=False,
)


# ============================================================
# 13. Hysteresis plot — D
# ============================================================

forward = hysteresis[
    hysteresis[
        "Direction"
    ]
    ==
    "forward"
]

reverse = hysteresis[
    hysteresis[
        "Direction"
    ]
    ==
    "reverse"
]


plt.figure(
    figsize=(8, 6)
)

plt.plot(
    forward["rho"],
    forward["D"],
    label="Increasing rho",
)

plt.plot(
    reverse["rho"],
    reverse["D"],
    label="Decreasing rho",
)

plt.xlabel(
    "Rb-loss drive rho"
)

plt.ylabel(
    "Death-state activity D"
)

plt.title(
    "History dependence of the death state"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    FIG
    / "02_death_hysteresis.png",
    dpi=250,
)

plt.close()


# ============================================================
# 14. Hysteresis plot — trophic state
# ============================================================

plt.figure(
    figsize=(8, 6)
)

plt.plot(
    forward["rho"],
    forward["T"],
    label="Increasing rho",
)

plt.plot(
    reverse["rho"],
    reverse["T"],
    label="Decreasing rho",
)

plt.xlabel(
    "Rb-loss drive rho"
)

plt.ylabel(
    "Trophic-support activity T"
)

plt.title(
    "History dependence of trophic support"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    FIG
    / "03_trophic_hysteresis.png",
    dpi=250,
)

plt.close()


# ============================================================
# 15. Quantify forward/reverse separation
# ============================================================

reverse_sorted = (
    reverse.sort_values(
        "rho"
    )
)


comparison = forward.merge(
    reverse_sorted,
    on="rho",
    suffixes=(
        "_forward",
        "_reverse",
    ),
)


comparison[
    "Delta_D"
] = (
    comparison[
        "D_reverse"
    ]
    -
    comparison[
        "D_forward"
    ]
)


comparison[
    "Delta_T"
] = (
    comparison[
        "T_reverse"
    ]
    -
    comparison[
        "T_forward"
    ]
)


comparison.to_csv(
    OUT
    / "F11C_hysteresis_difference.csv",
    index=False,
)


max_D_gap = (
    comparison[
        "Delta_D"
    ]
    .abs()
    .max()
)


max_T_gap = (
    comparison[
        "Delta_T"
    ]
    .abs()
    .max()
)


hysteretic = (
    max_D_gap > 0.20
)


print()
print("HYSTERESIS")
print("-" * 100)

print(
    f"Maximum forward/reverse D separation: "
    f"{max_D_gap:.4f}"
)

print(
    f"Maximum forward/reverse T separation: "
    f"{max_T_gap:.4f}"
)

print(
    "Strong history dependence detected:",
    hysteretic,
)


print()
print("Figures written to:")
print(FIG)

print()
print(
    "[NEXT] If hysteresis is present, determine the switching "
    "thresholds and test which molecular perturbations move "
    "the system across the separatrix."
)