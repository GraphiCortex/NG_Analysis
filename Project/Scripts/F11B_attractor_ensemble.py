from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.optimize import root, brentq
from scipy.integrate import solve_ivp


# ============================================================
# F11B — Sign-constrained adult-fate dynamical ensemble
#
# States
# ------
# U = conserved Rb/E2F-cycle/replication backbone
# D = adult death/stress program
# T = trophic/quiescence-support state
#
# External control
# ----------------
# rho = Rb-loss drive
#
# Empirically supported signs used here:
#
#   U -> U   positive
#   U -> D   positive
#   D -> U   positive
#   D -> D   positive
#   D -| T   negative
#   T -> T   positive
#
# We DO NOT assume bistability.
#
# Instead:
#   1. sample many admissible parameter sets
#   2. solve fixed points
#   3. classify Jacobian stability
#   4. ask how often multiple sinks occur
#   5. visualize one representative admissible system
#
# This establishes dynamical PLAUSIBILITY, not biological proof
# of kinetic parameters.
# ============================================================


ROOT = Path(".")

OUT = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F11B_AttractorEnsemble"
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


RNG = np.random.default_rng(
    20261022
)


# ============================================================
# 1. Regulatory functions
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


def hill_derivative(x, K, n):

    x = max(
        float(x),
        1e-12,
    )

    return (
        n
        *
        K ** n
        *
        x ** (n - 1)
        /
        (
            K ** n
            +
            x ** n
        ) ** 2
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


def repression_derivative(x, K, n):

    x = max(
        float(x),
        1e-12,
    )

    return (
        -n
        *
        K ** n
        *
        x ** (n - 1)
        /
        (
            K ** n
            +
            x ** n
        ) ** 2
    )


# ============================================================
# 2. ODE system
#
# Degradation rates are scaled to one.
# Thus time is measured in arbitrary decay-time units.
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
# 3. Jacobian
# ============================================================

def jacobian(state, p, rho):

    U, D, T = state


    J = np.zeros(
        (3, 3)
    )


    J[0, 0] = (
        p["aUU"]
        *
        hill_derivative(
            U,
            p["KUU"],
            p["nUU"],
        )
        -
        1
    )


    J[0, 1] = (
        p["aDU"]
        *
        hill_derivative(
            D,
            p["KDU"],
            p["nDU"],
        )
    )


    J[1, 0] = (
        p["aUD"]
        *
        hill_derivative(
            U,
            p["KUD"],
            p["nUD"],
        )
    )


    J[1, 1] = (
        p["aDD"]
        *
        hill_derivative(
            D,
            p["KDD"],
            p["nDD"],
        )
        -
        1
    )


    J[2, 1] = (
        p["aTT"]
        *
        hill(
            T,
            p["KTT"],
            p["nTT"],
        )
        *
        repression_derivative(
            D,
            p["KDT"],
            p["nDT"],
        )
    )


    J[2, 2] = (
        p["aTT"]
        *
        hill_derivative(
            T,
            p["KTT"],
            p["nTT"],
        )
        *
        repression(
            D,
            p["KDT"],
            p["nDT"],
        )
        -
        1
    )


    return J


# ============================================================
# 4. Fixed-point finder
# ============================================================

def find_fixed_points(
    p,
    rho,
    max_state=4.0,
):

    guesses = np.linspace(
        0.05,
        3.0,
        5,
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
                            p,
                            rho,
                        ),
                    x0=np.array([
                        U0,
                        D0,
                        T0,
                    ]),
                )


                if not sol.success:
                    continue


                x = sol.x


                if not np.all(
                    np.isfinite(x)
                ):
                    continue


                if np.any(
                    x < -1e-6
                ):
                    continue


                if np.any(
                    x > max_state
                ):
                    continue


                if np.linalg.norm(
                    rhs(
                        0,
                        x,
                        p,
                        rho,
                    )
                ) > 1e-6:
                    continue


                x = np.maximum(
                    x,
                    0
                )


                duplicate = any(
                    np.linalg.norm(
                        x - y["state"]
                    )
                    < 1e-3
                    for y in fixed
                )


                if duplicate:
                    continue


                J = jacobian(
                    x,
                    p,
                    rho,
                )

                eig = np.linalg.eigvals(
                    J
                )


                max_real = np.max(
                    np.real(eig)
                )


                if max_real < -1e-5:

                    stability = (
                        "sink"
                    )

                elif np.all(
                    np.real(eig)
                    > 1e-5
                ):

                    stability = (
                        "source"
                    )

                else:

                    stability = (
                        "saddle_or_nonhyperbolic"
                    )


                fixed.append({

                    "state":
                        x,

                    "eig":
                        eig,

                    "max_real":
                        max_real,

                    "stability":
                        stability,
                })


    fixed.sort(
        key=lambda z:
            (
                z["state"][0],
                z["state"][1],
                z["state"][2],
            )
    )


    return fixed


# ============================================================
# 5. Parameter sampler
#
# These ranges are broad exploratory priors.
#
# They are NOT measured biochemical constants.
# ============================================================

def sample_parameters():

    return {

        "bU":
            RNG.uniform(
                0.01,
                0.15,
            ),

        "bD":
            RNG.uniform(
                0.01,
                0.15,
            ),

        "bT":
            RNG.uniform(
                0.05,
                0.30,
            ),

        "aR":
            RNG.uniform(
                0.30,
                1.20,
            ),


        "aUU":
            RNG.uniform(
                0.30,
                1.60,
            ),

        "aDU":
            RNG.uniform(
                0.30,
                1.60,
            ),

        "aUD":
            RNG.uniform(
                0.30,
                1.60,
            ),

        "aDD":
            RNG.uniform(
                0.30,
                1.60,
            ),

        "aTT":
            RNG.uniform(
                0.70,
                2.00,
            ),


        "KUU":
            RNG.uniform(
                0.20,
                1.20,
            ),

        "KDU":
            RNG.uniform(
                0.20,
                1.20,
            ),

        "KUD":
            RNG.uniform(
                0.20,
                1.20,
            ),

        "KDD":
            RNG.uniform(
                0.20,
                1.20,
            ),

        "KTT":
            RNG.uniform(
                0.20,
                1.20,
            ),

        "KDT":
            RNG.uniform(
                0.20,
                1.20,
            ),


        "nUU":
            RNG.uniform(
                1.2,
                4.0,
            ),

        "nDU":
            RNG.uniform(
                1.2,
                4.0,
            ),

        "nUD":
            RNG.uniform(
                1.2,
                4.0,
            ),

        "nDD":
            RNG.uniform(
                1.2,
                4.0,
            ),

        "nTT":
            RNG.uniform(
                1.2,
                4.0,
            ),

        "nDT":
            RNG.uniform(
                1.2,
                4.0,
            ),
    }


# ============================================================
# 6. Qualitative biological consistency
#
# We require only:
#
#   Control-like rho=0:
#       low U
#       low D
#       relatively high T
#
#   Rb-loss rho=1:
#       higher U
#       higher D
#       lower T
#
# Thresholds are deliberately loose.
# ============================================================

def evaluate_parameter_set(p):

    control_fp = find_fixed_points(
        p,
        rho=0.0,
    )

    control_sinks = [
        x
        for x in control_fp
        if x["stability"]
        ==
        "sink"
    ]


    if len(control_sinks) != 1:

        return None


    ctrl = (
        control_sinks[0]["state"]
    )

    U0, D0, T0 = ctrl


    if not (
        U0 < 0.6
        and
        D0 < 0.6
        and
        T0 > 0.25
    ):

        return None


    loss_fp = find_fixed_points(
        p,
        rho=1.0,
    )

    loss_sinks = [
        x
        for x in loss_fp
        if x["stability"]
        ==
        "sink"
    ]


    if len(loss_sinks) == 0:

        return None


    candidate_loss = []


    for fp in loss_sinks:

        U1, D1, T1 = (
            fp["state"]
        )


        if (
            U1
            >
            U0 + 0.35
            and
            D1
            >
            D0 + 0.35
            and
            T1
            <
            T0 - 0.10
        ):

            candidate_loss.append(
                fp
            )


    if not candidate_loss:

        return None


    # --------------------------------------------------------
    # Scan rho for multistability
    # --------------------------------------------------------

    rho_scan = np.linspace(
        0,
        1,
        11,
    )


    max_sinks = 0

    best_rho = None

    max_D_separation = 0


    for rho in rho_scan:

        fps = find_fixed_points(
            p,
            rho,
        )

        sinks = [
            x
            for x in fps
            if x["stability"]
            ==
            "sink"
        ]


        max_sinks = max(
            max_sinks,
            len(sinks),
        )


        if len(sinks) >= 2:

            Dvals = [
                x["state"][1]
                for x in sinks
            ]

            separation = (
                max(Dvals)
                -
                min(Dvals)
            )


            if (
                separation
                >
                max_D_separation
            ):

                max_D_separation = (
                    separation
                )

                best_rho = rho


    bistable = (
        max_sinks >= 2
        and
        max_D_separation >= 0.30
    )


    # Score for choosing a clear representative model.

    best_loss = max(
        candidate_loss,
        key=lambda x:
            x["state"][1],
    )


    U1, D1, T1 = (
        best_loss["state"]
    )


    score = (
        (U1 - U0)
        +
        (D1 - D0)
        +
        (T0 - T1)
        +
        2.0
        *
        max_D_separation
    )


    return {

        "control":
            ctrl,

        "loss":
            best_loss["state"],

        "max_sinks":
            max_sinks,

        "bistable":
            bistable,

        "best_rho":
            best_rho,

        "max_D_separation":
            max_D_separation,

        "score":
            score,
    }


# ============================================================
# 7. Parameter ensemble
# ============================================================

N_PARAMETER_SETS = 1200


accepted_rows = []

accepted_models = []


print()
print("=" * 100)
print("F11B — SIGN-CONSTRAINED ATTRACTOR ENSEMBLE")
print("=" * 100)

print()
print(
    "Sampling",
    N_PARAMETER_SETS,
    "kinetic parameter sets..."
)


for i in range(
    N_PARAMETER_SETS
):

    p = sample_parameters()

    result = (
        evaluate_parameter_set(
            p
        )
    )


    if result is None:
        continue


    row = {
        "Parameter_set":
            i,

        **p,

        "Control_U":
            result[
                "control"
            ][0],

        "Control_D":
            result[
                "control"
            ][1],

        "Control_T":
            result[
                "control"
            ][2],

        "Loss_U":
            result[
                "loss"
            ][0],

        "Loss_D":
            result[
                "loss"
            ][1],

        "Loss_T":
            result[
                "loss"
            ][2],

        "Max_sinks":
            result[
                "max_sinks"
            ],

        "Bistable":
            result[
                "bistable"
            ],

        "Best_bistable_rho":
            result[
                "best_rho"
            ],

        "Max_D_separation":
            result[
                "max_D_separation"
            ],

        "Score":
            result[
                "score"
            ],
    }


    accepted_rows.append(
        row
    )

    accepted_models.append(
        (
            p,
            result,
            row,
        )
    )


ensemble = pd.DataFrame(
    accepted_rows
)


if ensemble.empty:

    raise RuntimeError(
        "No sampled parameter set satisfied the loose "
        "control/Rb-loss qualitative constraints. "
        "Broaden parameter ranges before interpreting."
    )


ensemble.to_csv(
    OUT
    / "F11B_parameter_ensemble.csv",
    index=False,
)


# ============================================================
# 8. Ensemble summary
# ============================================================

n_accepted = len(
    ensemble
)

n_bistable = int(
    ensemble[
        "Bistable"
    ].sum()
)

bistable_fraction = (
    n_bistable
    /
    n_accepted
)


print()
print("ENSEMBLE RESULTS")
print("-" * 100)

print(
    "Accepted qualitative models:",
    n_accepted,
)

print(
    "Models exhibiting >=2 distinct sinks:",
    n_bistable,
)

print(
    "Bistable fraction among accepted models:",
    f"{bistable_fraction:.3f}",
)


print(
    "Median maximum sink count:",
    f"{ensemble['Max_sinks'].median():.2f}",
)


# ============================================================
# 9. Choose representative model
#
# Prefer a bistable model if any exist.
# Otherwise use the highest-scoring admissible monostable model.
# ============================================================

if n_bistable > 0:

    candidates = ensemble[
        ensemble[
            "Bistable"
        ]
    ]

else:

    candidates = ensemble


representative_row = (
    candidates.sort_values(
        "Score",
        ascending=False,
    )
    .iloc[0]
)


parameter_index = int(
    representative_row[
        "Parameter_set"
    ]
)


representative = None


for p, result, row in accepted_models:

    if int(
        row[
            "Parameter_set"
        ]
    ) == parameter_index:

        representative = (
            p,
            result,
        )

        break


p_rep, res_rep = (
    representative
)


pd.DataFrame(
    [
        p_rep
    ]
).to_csv(
    OUT
    / "F11B_representative_parameters.csv",
    index=False,
)


print()
print("REPRESENTATIVE MODEL")
print("-" * 100)

print(
    "Parameter set:",
    parameter_index,
)

print(
    "Bistable:",
    res_rep[
        "bistable"
    ],
)

print(
    "Best bistable rho:",
    res_rep[
        "best_rho"
    ],
)

print(
    "Control state:",
    np.round(
        res_rep[
            "control"
        ],
        4,
    ),
)

print(
    "Rb-loss state:",
    np.round(
        res_rep[
            "loss"
        ],
        4,
    ),
)


# ============================================================
# 10. Dense bifurcation scan
# ============================================================

rho_dense = np.linspace(
    0,
    1.30,
    161,
)


branch_rows = []


for rho in rho_dense:

    fps = find_fixed_points(
        p_rep,
        rho,
    )


    for j, fp in enumerate(
        fps
    ):

        U, D, T = (
            fp[
                "state"
            ]
        )


        branch_rows.append({

            "rho":
                rho,

            "Fixed_point":
                j,

            "U":
                U,

            "D":
                D,

            "T":
                T,

            "Stability":
                fp[
                    "stability"
                ],

            "Max_real_eigenvalue":
                fp[
                    "max_real"
                ],
        })


branches = pd.DataFrame(
    branch_rows
)


branches.to_csv(
    OUT
    / "F11B_representative_fixed_points.csv",
    index=False,
)


# ============================================================
# 11. Bifurcation-style plots
# ============================================================

stable = branches[
    branches[
        "Stability"
    ]
    ==
    "sink"
]

unstable = branches[
    branches[
        "Stability"
    ]
    !=
    "sink"
]


plt.figure(
    figsize=(8, 6)
)

plt.scatter(
    stable["rho"],
    stable["D"],
    s=12,
    label="Stable sink",
)

plt.scatter(
    unstable["rho"],
    unstable["D"],
    s=12,
    marker="x",
    label="Unstable / saddle",
)

plt.xlabel(
    "Rb-loss drive rho"
)

plt.ylabel(
    "Death-state activity D"
)

plt.title(
    "Candidate death-state bifurcation structure"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    FIG
    / "01_death_bifurcation.png",
    dpi=220,
)

plt.close()


plt.figure(
    figsize=(8, 6)
)

plt.scatter(
    stable["rho"],
    stable["T"],
    s=12,
    label="Stable sink",
)

plt.scatter(
    unstable["rho"],
    unstable["T"],
    s=12,
    marker="x",
    label="Unstable / saddle",
)

plt.xlabel(
    "Rb-loss drive rho"
)

plt.ylabel(
    "Trophic-support activity T"
)

plt.title(
    "Trophic-support attractor structure"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    FIG
    / "02_trophic_bifurcation.png",
    dpi=220,
)

plt.close()


# ============================================================
# 12. Choose rho for phase portrait
# ============================================================

if (
    res_rep[
        "bistable"
    ]
    and
    res_rep[
        "best_rho"
    ]
    is not None
):

    rho_phase = float(
        res_rep[
            "best_rho"
        ]
    )

else:

    rho_phase = 0.5


fps_phase = find_fixed_points(
    p_rep,
    rho_phase,
)


# ============================================================
# 13. U-D phase portrait
#
# T is not required for the U-D vector field because the
# conservative signed topology contains no T -> U or T -> D
# edge.
# ============================================================

uvals = np.linspace(
    0,
    2.8,
    31,
)

dvals = np.linspace(
    0,
    2.8,
    31,
)


UU, DD = np.meshgrid(
    uvals,
    dvals,
)


DU = np.zeros_like(
    UU
)

DDOT = np.zeros_like(
    DD
)


for i in range(
    UU.shape[0]
):

    for j in range(
        UU.shape[1]
    ):

        state = np.array([
            UU[i, j],
            DD[i, j],
            0.5,
        ])

        vector = rhs(
            0,
            state,
            p_rep,
            rho_phase,
        )

        DU[i, j] = (
            vector[0]
        )

        DDOT[i, j] = (
            vector[1]
        )


speed = np.sqrt(
    DU ** 2
    +
    DDOT ** 2
)

speed[
    speed == 0
] = 1


plt.figure(
    figsize=(8, 7)
)


plt.streamplot(
    UU,
    DD,
    DU,
    DDOT,
    density=1.0,
)


# U-nullcline

plt.contour(
    UU,
    DD,
    DU,
    levels=[0],
    linewidths=2,
)


# D-nullcline

plt.contour(
    UU,
    DD,
    DDOT,
    levels=[0],
    linewidths=2,
    linestyles="--",
)


for fp in fps_phase:

    U, D, T = (
        fp[
            "state"
        ]
    )


    if fp[
        "stability"
    ] == "sink":

        marker = "o"

    else:

        marker = "x"


    plt.scatter(
        U,
        D,
        s=90,
        marker=marker,
    )


plt.xlabel(
    "U — conserved Rb/E2F stress state"
)

plt.ylabel(
    "D — adult death state"
)

plt.title(
    f"U-D phase portrait at rho={rho_phase:.2f}"
)

plt.tight_layout()

plt.savefig(
    FIG
    / "03_UD_phase_portrait.png",
    dpi=240,
)

plt.close()


# ============================================================
# 14. Trajectories from different initial conditions
# ============================================================

initials = []

for U0 in [
    0.05,
    0.6,
    1.2,
    2.0,
]:

    for D0 in [
        0.05,
        0.6,
        1.2,
        2.0,
    ]:

        initials.append([
            U0,
            D0,
            0.8,
        ])


trajectory_rows = []


plt.figure(
    figsize=(8, 7)
)


for k, x0 in enumerate(
    initials
):

    sol = solve_ivp(

        fun=lambda t, x:
            rhs(
                t,
                x,
                p_rep,
                rho_phase,
            ),

        t_span=(
            0,
            80,
        ),

        y0=x0,

        max_step=0.2,

        rtol=1e-7,

        atol=1e-9,
    )


    plt.plot(
        sol.y[0],
        sol.y[1],
        linewidth=1,
    )


    for t, U, D, T in zip(
        sol.t,
        sol.y[0],
        sol.y[1],
        sol.y[2],
    ):

        trajectory_rows.append({

            "Trajectory":
                k,

            "time":
                t,

            "U":
                U,

            "D":
                D,

            "T":
                T,
        })


for fp in fps_phase:

    U, D, T = (
        fp[
            "state"
        ]
    )


    marker = (
        "o"
        if fp[
            "stability"
        ]
        ==
        "sink"
        else
        "x"
    )


    plt.scatter(
        U,
        D,
        s=100,
        marker=marker,
    )


plt.xlabel(
    "U — conserved Rb/E2F stress state"
)

plt.ylabel(
    "D — adult death state"
)

plt.title(
    f"Trajectories and attractors at rho={rho_phase:.2f}"
)

plt.tight_layout()

plt.savefig(
    FIG
    / "04_UD_trajectories.png",
    dpi=240,
)

plt.close()


pd.DataFrame(
    trajectory_rows
).to_csv(
    OUT
    / "F11B_representative_trajectories.csv",
    index=False,
)


# ============================================================
# 15. Final fixed-point summary
# ============================================================

print()
print("FIXED POINTS AT PHASE-PORTRAIT rho")
print("-" * 100)

print(
    f"rho = {rho_phase:.3f}"
)

print()


for i, fp in enumerate(
    fps_phase
):

    U, D, T = (
        fp[
            "state"
        ]
    )

    eig = (
        fp[
            "eig"
        ]
    )


    print(
        f"FP{i+1}: "
        f"U={U:.4f}, "
        f"D={D:.4f}, "
        f"T={T:.4f}   "
        f"{fp['stability']}"
    )

    print(
        "      eigenvalues =",
        np.round(
            eig,
            5,
        )
    )


print()
print("Figures written to:")
print(FIG)

print()
print(
    "IMPORTANT: These phase portraits demonstrate dynamics "
    "permitted by the empirically signed topology. They are "
    "not yet evidence that the fitted kinetic parameters are "
    "those used by real cells."
)

print()
print(
    "[NEXT] Inspect the bifurcation and phase-portrait figures. "
    "Then test whether attractor structure is robust across the "
    "parameter ensemble and derive experimentally testable "
    "switching predictions."
)