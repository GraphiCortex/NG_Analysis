from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp


# ============================================================
# F11E — Attractor rescue / separatrix-crossing screen
#
# Question:
#
# Starting from the committed high-D attractor, which
# perturbations can:
#
#   1. suppress the high-D state while treatment is present?
#   2. push the system across the separatrix so that it remains
#      in the low-D basin AFTER treatment is removed?
#
# Tested interventions:
#
#   U -> D inhibition
#   D -> D inhibition
#   D -> U inhibition
#   U -> U inhibition
#   reduction of Rb-loss drive rho
#
#   combinations of feedback-arm inhibition
#
#   trophic boosting = mechanistic negative control
#
# IMPORTANT:
# Parameter strengths and times are dimensionless.
# ============================================================


ROOT = Path(".")

F11B = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F11B_AttractorEnsemble"
)

F11D = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F11D_EnsembleHysteresis"
)

OUT = (
    ROOT
    / "Data"
    / "Processed"
    / "Dynamics"
    / "F11E_RescueScreen"
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
# 1. Load models
# ============================================================

ensemble = pd.read_csv(
    F11B
    / "F11B_parameter_ensemble.csv"
)

hyst = pd.read_csv(
    F11D
    / "F11D_hysteresis_summary.csv"
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


models = (
    ensemble[
        ensemble["Bistable"]
    ]
    .merge(
        hyst,
        on="Parameter_set",
        how="inner",
    )
    .copy()
)


print()
print("=" * 105)
print("F11E — ATTRACTOR RESCUE SCREEN")
print("=" * 105)

print()
print(
    "Bistable/hysteretic models:",
    len(models),
)


# ============================================================
# 2. Model parameter columns
# ============================================================

PARAMS = [

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
# 4. ODE
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
# 5. Integrator
# ============================================================

def evolve(
    state,
    p,
    rho,
    duration,
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

        max_step=0.25,
    )


    return sol.y[
        :,
        -1
    ]


# ============================================================
# 6. Intervention definitions
#
# strength = fractional intervention:
#
# 0.0 = none
# 1.0 = complete inhibition for inhibitory interventions
#
# Trophic boost uses:
#
# bT -> bT * (1 + 4*strength)
#
# ============================================================

INTERVENTIONS = [

    "INHIBIT_U_TO_D",

    "INHIBIT_D_SELF",

    "INHIBIT_D_TO_U",

    "INHIBIT_U_SELF",

    "REDUCE_RB_DRIVE",

    "INHIBIT_U_TO_D_PLUS_D_SELF",

    "INHIBIT_U_TO_D_PLUS_D_TO_U",

    "INHIBIT_FULL_UD_LOOP",

    "BOOST_TROPHIC",
]


def perturb(
    p,
    rho,
    intervention,
    strength,
):

    q = dict(
        p
    )

    r = float(
        rho
    )


    inhibition = (
        1.0
        -
        strength
    )


    if intervention == "INHIBIT_U_TO_D":

        q["aUD"] *= inhibition


    elif intervention == "INHIBIT_D_SELF":

        q["aDD"] *= inhibition


    elif intervention == "INHIBIT_D_TO_U":

        q["aDU"] *= inhibition


    elif intervention == "INHIBIT_U_SELF":

        q["aUU"] *= inhibition


    elif intervention == "REDUCE_RB_DRIVE":

        r *= inhibition


    elif intervention == "INHIBIT_U_TO_D_PLUS_D_SELF":

        q["aUD"] *= inhibition

        q["aDD"] *= inhibition


    elif intervention == "INHIBIT_U_TO_D_PLUS_D_TO_U":

        q["aUD"] *= inhibition

        q["aDU"] *= inhibition


    elif intervention == "INHIBIT_FULL_UD_LOOP":

        q["aUD"] *= inhibition

        q["aDU"] *= inhibition

        q["aDD"] *= inhibition

        q["aUU"] *= inhibition


    elif intervention == "BOOST_TROPHIC":

        q["bT"] *= (
            1.0
            +
            4.0
            *
            strength
        )


    else:

        raise ValueError(
            intervention
        )


    return (
        q,
        r,
    )


# ============================================================
# 7. Screen design
# ============================================================

STRENGTHS = np.linspace(
    0.10,
    1.00,
    10,
)


PULSE_DURATIONS = [
    2,
    5,
    10,
    20,
    40,
]


HELD_DURATION = 100

RECOVERY_DURATION = 120


# ============================================================
# 8. Helpers
# ============================================================

def distance_UD(
    x,
    ref,
):

    return np.linalg.norm(
        np.asarray(
            x[:2]
        )
        -
        np.asarray(
            ref[:2]
        )
    )


def is_low_basin(
    state,
    low_ref,
    high_ref,
):

    return (
        distance_UD(
            state,
            low_ref,
        )
        <
        distance_UD(
            state,
            high_ref,
        )
    )


# ============================================================
# 9. Recover low/high states for each model
#
# We test at the midpoint of that model's hysteresis window:
#
# rho_test = (rho_up + rho_down)/2
#
# Thus both attractors should coexist.
# ============================================================

model_refs = {}

valid_model_ids = []


for _, row in models.iterrows():

    model_id = int(
        row[
            "Parameter_set"
        ]
    )


    p = {
        key:
            float(
                row[key]
            )
        for key in PARAMS
    }


    rho_test = 0.5 * (
        float(
            row[
                "rho_up"
            ]
        )
        +
        float(
            row[
                "rho_down"
            ]
        )
    )


    # Low branch:
    # begin from the known control state.

    low_start = np.array([

        row[
            "Control_U"
        ],

        row[
            "Control_D"
        ],

        row[
            "Control_T"
        ],
    ])


    low_ref = evolve(
        low_start,
        p,
        rho_test,
        150,
    )


    # High branch:
    # begin from the known rho=1 Rb-loss state.

    high_start = np.array([

        row[
            "Loss_U"
        ],

        row[
            "Loss_D"
        ],

        row[
            "Loss_T"
        ],
    ])


    high_ref = evolve(
        high_start,
        p,
        rho_test,
        150,
    )


    if (
        abs(
            high_ref[1]
            -
            low_ref[1]
        )
        <
        0.20
    ):

        print(
            "[WARNING] "
            f"Model {model_id} did not retain "
            "well-separated low/high states at "
            "the hysteresis midpoint."
        )

        continue


    model_refs[
        model_id
    ] = {

        "p":
            p,

        "rho":
            rho_test,

        "low":
            low_ref,

        "high":
            high_ref,
    }


    valid_model_ids.append(
        model_id
    )


print()
print(
    "Models with validated low/high references:",
    len(
        valid_model_ids
    ),
)


# ============================================================
# 10. Sustained intervention screen
#
# A sustained rescue means the perturbed system moves into the
# low-U/low-D side of the original bistable landscape while
# the intervention is continuously present.
# ============================================================

sustained_rows = []


for model_id in valid_model_ids:

    ref = model_refs[
        model_id
    ]

    p = ref["p"]

    rho = ref["rho"]

    low = ref["low"]

    high = ref["high"]


    midpoint_U = 0.5 * (
        low[0]
        +
        high[0]
    )

    midpoint_D = 0.5 * (
        low[1]
        +
        high[1]
    )


    for intervention in INTERVENTIONS:

        for strength in STRENGTHS:

            p2, rho2 = perturb(

                p,
                rho,
                intervention,
                strength,
            )


            final = evolve(

                high,
                p2,
                rho2,

                HELD_DURATION,
            )


            rescued = (
                final[0]
                <
                midpoint_U
                and
                final[1]
                <
                midpoint_D
            )


            sustained_rows.append({

                "Parameter_set":
                    model_id,

                "rho_test":
                    rho,

                "Intervention":
                    intervention,

                "Strength":
                    strength,

                "Final_U":
                    final[0],

                "Final_D":
                    final[1],

                "Final_T":
                    final[2],

                "Sustained_rescue":
                    rescued,
            })


sustained = pd.DataFrame(
    sustained_rows
)


sustained.to_csv(
    OUT
    / "F11E_sustained_screen.csv",
    index=False,
)


# ============================================================
# 11. Transient pulse screen
#
# Pulse intervention -> restore original model -> recovery.
#
# Persistent rescue means the recovered state lies in the
# original low-D basin.
# ============================================================

pulse_rows = []


for model_id in valid_model_ids:

    ref = model_refs[
        model_id
    ]

    p = ref["p"]

    rho = ref["rho"]

    low = ref["low"]

    high = ref["high"]


    for intervention in INTERVENTIONS:

        for duration in PULSE_DURATIONS:

            for strength in STRENGTHS:

                p2, rho2 = perturb(

                    p,
                    rho,
                    intervention,
                    strength,
                )


                after_pulse = evolve(

                    high,
                    p2,
                    rho2,

                    duration,
                )


                recovered = evolve(

                    after_pulse,
                    p,
                    rho,

                    RECOVERY_DURATION,
                )


                switched = is_low_basin(

                    recovered,
                    low,
                    high,
                )


                pulse_rows.append({

                    "Parameter_set":
                        model_id,

                    "rho_test":
                        rho,

                    "Intervention":
                        intervention,

                    "Pulse_duration":
                        duration,

                    "Strength":
                        strength,

                    "Postpulse_U":
                        after_pulse[0],

                    "Postpulse_D":
                        after_pulse[1],

                    "Postpulse_T":
                        after_pulse[2],

                    "Recovered_U":
                        recovered[0],

                    "Recovered_D":
                        recovered[1],

                    "Recovered_T":
                        recovered[2],

                    "Persistent_switch":
                        switched,
                })


pulse = pd.DataFrame(
    pulse_rows
)


pulse.to_csv(
    OUT
    / "F11E_transient_screen.csv",
    index=False,
)


# ============================================================
# 12. Sustained summary
# ============================================================

sustained_summary = (

    sustained.groupby(
        [
            "Intervention",
            "Strength",
        ]
    )[
        "Sustained_rescue"
    ]
    .mean()
    .reset_index(
        name="Rescue_fraction"
    )
)


sustained_summary.to_csv(
    OUT
    / "F11E_sustained_summary.csv",
    index=False,
)


# Minimum successful sustained strength per model

min_sustained = (

    sustained[
        sustained[
            "Sustained_rescue"
        ]
    ]

    .groupby(
        [
            "Parameter_set",
            "Intervention",
        ]
    )[
        "Strength"
    ]

    .min()

    .reset_index(
        name="Minimum_strength"
    )
)


# ============================================================
# 13. Transient summary
# ============================================================

pulse_summary = (

    pulse.groupby(
        [
            "Intervention",
            "Pulse_duration",
            "Strength",
        ]
    )[
        "Persistent_switch"
    ]

    .mean()

    .reset_index(
        name="Switch_fraction"
    )
)


pulse_summary.to_csv(
    OUT
    / "F11E_transient_summary.csv",
    index=False,
)


min_pulse = (

    pulse[
        pulse[
            "Persistent_switch"
        ]
    ]

    .groupby(
        [
            "Parameter_set",
            "Intervention",
            "Pulse_duration",
        ]
    )[
        "Strength"
    ]

    .min()

    .reset_index(
        name="Minimum_strength"
    )
)


min_pulse.to_csv(
    OUT
    / "F11E_minimum_pulse_strength.csv",
    index=False,
)


# ============================================================
# 14. Intervention-level summary
# ============================================================

summary_rows = []


for intervention in INTERVENTIONS:

    sub_sustained = sustained[
        sustained[
            "Intervention"
        ]
        ==
        intervention
    ]


    max_strength_sustained = sub_sustained[
        np.isclose(
            sub_sustained[
                "Strength"
            ],
            STRENGTHS.max(),
        )
    ]


    sustained_fraction_max = (
        max_strength_sustained[
            "Sustained_rescue"
        ]
        .mean()
    )


    mins = min_sustained[
        min_sustained[
            "Intervention"
        ]
        ==
        intervention
    ]


    median_min_strength = (
        mins[
            "Minimum_strength"
        ]
        .median()
        if len(mins)
        else np.nan
    )


    pulse_sub = pulse[
        pulse[
            "Intervention"
        ]
        ==
        intervention
    ]


    max_persistent_fraction = (

        pulse_sub.groupby(
            [
                "Pulse_duration",
                "Strength",
            ]
        )[
            "Persistent_switch"
        ]
        .mean()
        .max()
    )


    summary_rows.append({

        "Intervention":
            intervention,

        "Sustained_rescue_fraction_at_max_strength":
            sustained_fraction_max,

        "Models_with_any_sustained_rescue":
            mins[
                "Parameter_set"
            ]
            .nunique(),

        "Median_minimum_sustained_strength":
            median_min_strength,

        "Best_persistent_switch_fraction":
            max_persistent_fraction,
    })


intervention_summary = pd.DataFrame(
    summary_rows
)


intervention_summary = (
    intervention_summary.sort_values(
        [
            "Best_persistent_switch_fraction",
            "Sustained_rescue_fraction_at_max_strength",
        ],
        ascending=False,
    )
)


intervention_summary.to_csv(
    OUT
    / "F11E_intervention_summary.csv",
    index=False,
)


# ============================================================
# 15. Console
# ============================================================

print()
print("INTERVENTION SUMMARY")
print("-" * 105)


print(
    intervention_summary.to_string(
        index=False,
        float_format=lambda x: f"{x:.3f}",
    )
)


# ============================================================
# 16. Best transient intervention by duration
# ============================================================

print()
print("PERSISTENT SWITCHING BY PULSE DURATION")
print("-" * 105)


best_rows = []


for intervention in INTERVENTIONS:

    for duration in PULSE_DURATIONS:

        sub = pulse_summary[
            (
                pulse_summary[
                    "Intervention"
                ]
                ==
                intervention
            )
            &
            (
                pulse_summary[
                    "Pulse_duration"
                ]
                ==
                duration
            )
        ]


        if sub.empty:
            continue


        best = sub.loc[
            sub[
                "Switch_fraction"
            ]
            .idxmax()
        ]


        best_rows.append({

            "Intervention":
                intervention,

            "Pulse_duration":
                duration,

            "Best_strength":
                best[
                    "Strength"
                ],

            "Switch_fraction":
                best[
                    "Switch_fraction"
                ],
        })


best_pulse = pd.DataFrame(
    best_rows
)


print(
    best_pulse.to_string(
        index=False,
        float_format=lambda x: f"{x:.3f}",
    )
)


best_pulse.to_csv(
    OUT
    / "F11E_best_pulse_by_duration.csv",
    index=False,
)


# ============================================================
# 17. Plot sustained dose-response
# ============================================================

plt.figure(
    figsize=(10, 7)
)


for intervention in INTERVENTIONS:

    x = sustained_summary[
        sustained_summary[
            "Intervention"
        ]
        ==
        intervention
    ]


    plt.plot(

        x[
            "Strength"
        ],

        x[
            "Rescue_fraction"
        ],

        marker="o",

        label=intervention,
    )


plt.xlabel(
    "Intervention strength"
)

plt.ylabel(
    "Fraction of models rescued"
)

plt.title(
    "Sustained rescue across the hysteretic ensemble"
)

plt.legend(
    fontsize=7,
)

plt.tight_layout()

plt.savefig(
    FIG
    / "01_sustained_rescue_dose_response.png",
    dpi=250,
)

plt.close()


# ============================================================
# 18. Plot best pulse performance
# ============================================================

plot_data = (

    best_pulse.pivot(

        index="Pulse_duration",

        columns="Intervention",

        values="Switch_fraction",
    )
)


plt.figure(
    figsize=(10, 7)
)


for intervention in plot_data.columns:

    plt.plot(

        plot_data.index,

        plot_data[
            intervention
        ],

        marker="o",

        label=intervention,
    )


plt.xlabel(
    "Pulse duration"
)

plt.ylabel(
    "Best persistent-switch fraction"
)

plt.title(
    "Transient perturbations that permanently cross the separatrix"
)

plt.legend(
    fontsize=7,
)

plt.tight_layout()

plt.savefig(
    FIG
    / "02_transient_rescue_duration.png",
    dpi=250,
)

plt.close()


# ============================================================
# 19. Mechanistic negative-control check
# ============================================================

trophic = intervention_summary[
    intervention_summary[
        "Intervention"
    ]
    ==
    "BOOST_TROPHIC"
]


print()
print("TROPHIC NEGATIVE CONTROL")
print("-" * 105)


if not trophic.empty:

    row = trophic.iloc[0]

    print(
        "Sustained rescue fraction:",
        f"{row['Sustained_rescue_fraction_at_max_strength']:.3f}",
    )

    print(
        "Best persistent switching fraction:",
        f"{row['Best_persistent_switch_fraction']:.3f}",
    )


print()
print("Outputs written to:")
print(OUT)

print()
print(
    "[NEXT] Translate the most robust separatrix-crossing "
    "interventions back from state-level edges to molecular "
    "targets in the F5/F10 network."
)