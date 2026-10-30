#!/usr/bin/env Rscript

# ============================================================
# F11C_minimal_dynamics_predictions.R
# ============================================================
#
# PURPOSE
# -------
#
# Convert the empirically constrained F11B state architecture
# into the minimal dynamical model:
#
#     dE/dt = kE * (UE - E)
#     dS/dt = kS * (US - S)
#
# where:
#
#   E = E2F / proliferative program
#   S = immune / stress program
#
# and UE, US are perturbation- and context-dependent target
# states taken directly from the null-centered F11B geometry.
#
#
# IMPORTANT
# ---------
#
# No kinetic parameter is estimated from the data.
#
# Time is NORMALIZED / ARBITRARY.
#
# kE and kS are explored over a broad sensitivity range.
#
# The model does NOT assume:
#
#   - E -> S coupling
#   - S -> E coupling
#   - positive feedback
#   - bistability
#   - hysteresis
#   - molecular master regulators
#
#
# COUNTERFACTUAL INTERVENTIONS
# ----------------------------
#
# Adult Rb loss + selective E normalization:
#
#   E target -> 0
#   S target remains at the adult Rb-loss value
#
# Adult Rb loss + selective S normalization:
#
#   S target -> 0
#   E target remains at the adult Rb-loss value
#
# These are PROGRAM-LEVEL interventions.
# They are model predictions, not observed experiments.
#
#
# PERTURBATION REMOVAL
# --------------------
#
# The minimal model predicts return to baseline after removal.
#
# Persistence after removal would falsify the memoryless model
# and motivate feedback / memory terms.
#
# ============================================================


# ============================================================
# PACKAGES
# ============================================================

required <- c(
    "readr"
)


missing <- required[
    !vapply(
        required,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]


if (length(missing) > 0) {

    stop(
        paste(
            "Missing package(s):",
            paste(
                missing,
                collapse = ", "
            )
        )
    )
}


# ============================================================
# PATHS
# ============================================================

INPUT_FILE <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F11B_StateSpaceCoupling",
    "F11B_state_space_coordinates.csv"
)


OUTPUT_DIR <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F11C_MinimalDynamics"
)


dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# PARAMETERS
# ============================================================

DT <- 0.02

DEFAULT_KE <- 1

DEFAULT_KS <- 1


RATE_GRID <- c(
    0.25,
    0.5,
    1,
    2,
    4
)


OBSERVED_CONTEXTS <- c(
    "Adult_RbTKO",
    "CAG_RbLoss",
    "MAP2_RbLoss",
    "E2F13_DKO"
)


# ============================================================
# LOAD F11B COORDINATES
# ============================================================

message(
    "Loading null-centered F11B state coordinates..."
)


coordinates <- readr::read_csv(
    INPUT_FILE,
    show_col_types = FALSE
)


required_columns <- c(
    "Definition",
    "Context",
    "E_score",
    "S_score"
)


missing_columns <- setdiff(
    required_columns,
    names(
        coordinates
    )
)


if (length(missing_columns) > 0) {

    stop(
        paste(
            "F11B coordinate file missing:",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    )
}


# ============================================================
# TARGET LOOKUP
# ============================================================

get_target <- function(
    definition,
    context
) {

    row <- coordinates[
        coordinates$Definition ==
            definition &
        coordinates$Context ==
            context,
        ,
        drop = FALSE
    ]


    if (nrow(row) != 1) {

        stop(
            paste(
                "Expected one coordinate for",
                definition,
                context
            )
        )
    }


    data.frame(

        E_target =
            as.numeric(
                row$E_score[
                    1
                ]
            ),

        S_target =
            as.numeric(
                row$S_score[
                    1
                ]
            ),

        stringsAsFactors = FALSE
    )
}


# ============================================================
# EXACT SOLUTION FOR ONE CONSTANT-INPUT SEGMENT
# ============================================================
#
# dX/dt = k * (target - X)
#
# has exact solution:
#
# X(t) = target + (X0 - target) * exp(-k * t)
#
# No numerical ODE solver is necessary.
# ============================================================

simulate_segment <- function(
    start_time,
    end_time,
    start_E,
    start_S,
    target_E,
    target_S,
    kE,
    kS,
    dt,
    phase
) {

    time <- seq(
        start_time,
        end_time,
        by = dt
    )


    if (
        length(time) == 0 ||
        tail(
            time,
            1
        ) <
            end_time
    ) {

        time <- c(
            time,
            end_time
        )
    }


    elapsed <- (
        time -
        start_time
    )


    E <- (
        target_E +
        (
            start_E -
            target_E
        ) *
        exp(
            -kE *
            elapsed
        )
    )


    S <- (
        target_S +
        (
            start_S -
            target_S
        ) *
        exp(
            -kS *
            elapsed
        )
    )


    data.frame(

        Time =
            time,

        E =
            E,

        S =
            S,

        E_target =
            target_E,

        S_target =
            target_S,

        Phase =
            phase,

        stringsAsFactors = FALSE
    )
}


# ============================================================
# PIECEWISE SIMULATION
# ============================================================

simulate_schedule <- function(
    schedule,
    kE,
    kS,
    dt
) {

    E0 <- 0

    S0 <- 0


    output <- data.frame()


    for (i in seq_len(
        nrow(
            schedule
        )
    )) {

        segment <- schedule[
            i,
            ,
            drop = FALSE
        ]


        simulation <- simulate_segment(

            start_time =
                segment$Start[
                    1
                ],

            end_time =
                segment$End[
                    1
                ],

            start_E =
                E0,

            start_S =
                S0,

            target_E =
                segment$E_target[
                    1
                ],

            target_S =
                segment$S_target[
                    1
                ],

            kE =
                kE,

            kS =
                kS,

            dt =
                dt,

            phase =
                segment$Phase[
                    1
                ]
        )


        if (
            i > 1 &&
            nrow(
                simulation
            ) > 0
        ) {

            simulation <- simulation[
                -1,
                ,
                drop = FALSE
            ]
        }


        output <- rbind(
            output,
            simulation
        )


        E0 <- tail(
            simulation$E,
            1
        )


        S0 <- tail(
            simulation$S,
            1
        )
    }


    rownames(
        output
    ) <- NULL


    output
}


# ============================================================
# BUILD OBSERVED PERTURBATION SCHEDULE
# ============================================================

make_onset_schedule <- function(
    E_target,
    S_target
) {

    data.frame(

        Start =
            0,

        End =
            8,

        E_target =
            E_target,

        S_target =
            S_target,

        Phase =
            "Perturbation_on",

        stringsAsFactors = FALSE
    )
}


# ============================================================
# ADULT COUNTERFACTUAL SCHEDULES
# ============================================================

make_adult_pulse_schedule <- function(
    E_target,
    S_target
) {

    data.frame(

        Start =
            c(
                0,
                4
            ),

        End =
            c(
                4,
                12
            ),

        E_target =
            c(
                E_target,
                0
            ),

        S_target =
            c(
                S_target,
                0
            ),

        Phase =
            c(
                "Adult_Rb_loss_on",
                "Perturbation_removed"
            ),

        stringsAsFactors = FALSE
    )
}


make_E_normalization_schedule <- function(
    E_target,
    S_target
) {

    data.frame(

        Start =
            c(
                0,
                4
            ),

        End =
            c(
                4,
                12
            ),

        E_target =
            c(
                E_target,
                0
            ),

        S_target =
            c(
                S_target,
                S_target
            ),

        Phase =
            c(
                "Adult_Rb_loss_on",
                "Selective_E_normalization"
            ),

        stringsAsFactors = FALSE
    )
}


make_S_normalization_schedule <- function(
    E_target,
    S_target
) {

    data.frame(

        Start =
            c(
                0,
                4
            ),

        End =
            c(
                4,
                12
            ),

        E_target =
            c(
                E_target,
                E_target
            ),

        S_target =
            c(
                S_target,
                0
            ),

        Phase =
            c(
                "Adult_Rb_loss_on",
                "Selective_S_normalization"
            ),

        stringsAsFactors = FALSE
    )
}


# ============================================================
# OBSERVED-STATE RELAXATION SIMULATIONS
# ============================================================

message(
    "Simulating empirical perturbation targets..."
)


observed_simulations <- data.frame()


for (definition in unique(
    coordinates$Definition
)) {

    for (context in OBSERVED_CONTEXTS) {

        target <- get_target(
            definition,
            context
        )


        schedule <- make_onset_schedule(

            target$E_target[
                1
            ],

            target$S_target[
                1
            ]
        )


        simulation <- simulate_schedule(

            schedule =
                schedule,

            kE =
                DEFAULT_KE,

            kS =
                DEFAULT_KS,

            dt =
                DT
        )


        simulation$Definition <-
            definition


        simulation$Scenario <-
            context


        simulation$kE <-
            DEFAULT_KE


        simulation$kS <-
            DEFAULT_KS


        observed_simulations <- rbind(
            observed_simulations,
            simulation
        )
    }
}


# ============================================================
# ADULT COUNTERFACTUAL SIMULATIONS
# ============================================================

message(
    "Simulating adult counterfactual interventions..."
)


counterfactual_simulations <- data.frame()


for (definition in unique(
    coordinates$Definition
)) {

    target <- get_target(
        definition,
        "Adult_RbTKO"
    )


    adult_E <- target$E_target[
        1
    ]


    adult_S <- target$S_target[
        1
    ]


    scenario_names <- c(
        "Adult_pulse_then_release",
        "Adult_selective_E_normalization",
        "Adult_selective_S_normalization"
    )


    schedules <- list(

        make_adult_pulse_schedule(
            adult_E,
            adult_S
        ),

        make_E_normalization_schedule(
            adult_E,
            adult_S
        ),

        make_S_normalization_schedule(
            adult_E,
            adult_S
        )
    )


    for (i in seq_along(
        scenario_names
    )) {

        schedule <- getElement(
            schedules,
            i
        )


        simulation <- simulate_schedule(

            schedule =
                schedule,

            kE =
                DEFAULT_KE,

            kS =
                DEFAULT_KS,

            dt =
                DT
        )


        simulation$Definition <-
            definition


        simulation$Scenario <-
            scenario_names[
                i
            ]


        simulation$kE <-
            DEFAULT_KE


        simulation$kS <-
            DEFAULT_KS


        counterfactual_simulations <- rbind(
            counterfactual_simulations,
            simulation
        )
    }
}


# ============================================================
# RATE SENSITIVITY
# ============================================================
#
# kE and kS are NOT fitted.
#
# We vary them broadly to test whether qualitative predictions
# depend on assumed timescale.
# ============================================================

message(
    "Running relaxation-rate sensitivity..."
)


rate_sensitivity <- data.frame()


diagnostic_times <- c(
    2,
    4,
    6,
    8,
    12
)


for (definition in unique(
    coordinates$Definition
)) {

    target <- get_target(
        definition,
        "Adult_RbTKO"
    )


    adult_E <- target$E_target[
        1
    ]


    adult_S <- target$S_target[
        1
    ]


    scenario_names <- c(
        "Adult_pulse_then_release",
        "Adult_selective_E_normalization",
        "Adult_selective_S_normalization"
    )


    schedules <- list(

        make_adult_pulse_schedule(
            adult_E,
            adult_S
        ),

        make_E_normalization_schedule(
            adult_E,
            adult_S
        ),

        make_S_normalization_schedule(
            adult_E,
            adult_S
        )
    )


    for (kE in RATE_GRID) {

        for (kS in RATE_GRID) {

            for (scenario_index in seq_along(
                scenario_names
            )) {

                schedule <- getElement(
                    schedules,
                    scenario_index
                )


                simulation <- simulate_schedule(

                    schedule =
                        schedule,

                    kE =
                        kE,

                    kS =
                        kS,

                    dt =
                        DT
                )


                for (diagnostic_time in diagnostic_times) {

                    nearest <- which.min(
                        abs(
                            simulation$Time -
                            diagnostic_time
                        )
                    )


                    one_row <- data.frame(

                        Definition =
                            definition,

                        Scenario =
                            scenario_names[
                                scenario_index
                            ],

                        kE =
                            kE,

                        kS =
                            kS,

                        Diagnostic_time =
                            diagnostic_time,

                        E =
                            simulation$E[
                                nearest
                            ],

                        S =
                            simulation$S[
                                nearest
                            ],

                        stringsAsFactors = FALSE
                    )


                    rate_sensitivity <- rbind(
                        rate_sensitivity,
                        one_row
                    )
                }
            }
        }
    }
}


# ============================================================
# CHARACTERISTIC RELAXATION TIMES
# ============================================================
#
# For dX/dt = k(target-X):
#
# fraction remaining after time t = exp(-k*t)
#
# t90 = log(10)/k
# t95 = -log(0.05)/k
#
# Again: normalized units only.
# ============================================================

rate_table <- data.frame(

    k =
        RATE_GRID,

    Time_to_90_percent =
        log(
            10
        ) /
        RATE_GRID,

    Time_to_95_percent =
        -log(
            0.05
        ) /
        RATE_GRID,

    stringsAsFactors = FALSE
)


# ============================================================
# FALSIFIABLE PREDICTIONS
# ============================================================

prediction_table <- data.frame(

    Prediction_ID =
        c(
            "P1",
            "P2",
            "P3",
            "P4",
            "P5",
            "P6"
        ),

    Prediction =
        c(
            paste(
                "Adult Rb loss should drive both E and S",
                "toward positive excess program states."
            ),

            paste(
                "Embryonic cortical Rb loss should retain",
                "a positive E state without acquiring the",
                "adult positive S state."
            ),

            paste(
                "A selective intervention that normalizes E",
                "during adult Rb loss need not normalize S."
            ),

            paste(
                "A selective intervention that normalizes S",
                "during adult Rb loss need not normalize E."
            ),

            paste(
                "Removal of the initiating perturbation should",
                "return E and S toward baseline in the",
                "memoryless minimal model."
            ),

            paste(
                "Persistent E or S after perturbation removal",
                "would require additional feedback, memory,",
                "or state-retention mechanisms."
            )
        ),

    Status =
        c(
            "Anchored to observed endpoint",
            "Anchored to observed cross-context endpoints",
            "Model-generated counterfactual prediction",
            "Model-generated counterfactual prediction",
            "Model-generated falsifiable prediction",
            "Explicit falsification criterion"
        ),

    Falsified_if =
        c(
            "Future measurements fail to reproduce positive adult E and S.",
            "CAG/MAP2 Rb loss reproducibly acquires an adult-like positive S state.",
            "Selective E normalization necessarily abolishes S.",
            "Selective S normalization necessarily abolishes E.",
            "Substantial program activity persists after perturbation removal.",
            "Not applicable: this statement defines the extension criterion."
        ),

    stringsAsFactors = FALSE
)


# ============================================================
# MODEL-ASSUMPTION TABLE
# ============================================================

assumption_table <- data.frame(

    Component =
        c(
            "E target states",
            "S target states",
            "Baseline",
            "E-S coupling",
            "Relaxation rates",
            "Feedback",
            "Perturbation memory"
        ),

    Source =
        c(
            "Observed / inferred from F11B null-centered program scores",
            "Observed / inferred from F11B null-centered program scores",
            "Matched-null centering",
            "Set to zero in the minimal model",
            "Not identified; sensitivity parameter only",
            "Absent from minimal model",
            "Absent from minimal model"
        ),

    Interpretation =
        c(
            "Empirically constrained endpoint",
            "Empirically constrained endpoint",
            "Program excess equals zero",
            "Model assumption justified by F11B non-identifiability",
            "Controls timing only, not endpoint",
            "Would require additional evidence",
            "Persistence would falsify the minimal model"
        ),

    stringsAsFactors = FALSE
)


# ============================================================
# PLOT: EMPIRICAL STATE-SPACE TARGETS
# ============================================================

plot_state_targets <- function(
    definition,
    output_file
) {

    table <- coordinates[
        coordinates$Definition ==
            definition,
        ,
        drop = FALSE
    ]


    x_range <- range(
        c(
            0,
            table$E_score
        ),
        na.rm = TRUE
    )


    y_range <- range(
        c(
            0,
            table$S_score
        ),
        na.rm = TRUE
    )


    x_margin <- 0.15 *
        diff(
            x_range
        )


    y_margin <- 0.15 *
        diff(
            y_range
        )


    if (
        !is.finite(
            x_margin
        ) ||
        x_margin == 0
    ) {

        x_margin <- 0.1
    }


    if (
        !is.finite(
            y_margin
        ) ||
        y_margin == 0
    ) {

        y_margin <- 0.1
    }


    png(
        output_file,
        width = 1800,
        height = 1400,
        res = 220
    )


    par(
        mar = c(
            5.5,
            5.5,
            4,
            2
        ),
        las = 1
    )


    plot(

        table$E_score,

        table$S_score,

        type = "n",

        xlim =
            c(
                x_range[
                    1
                ] -
                x_margin,

                x_range[
                    2
                ] +
                x_margin
            ),

        ylim =
            c(
                y_range[
                    1
                ] -
                y_margin,

                y_range[
                    2
                ] +
                y_margin
            ),

        xlab =
            "Null-centered E score",

        ylab =
            "Null-centered S score",

        main =
            paste(
                "F11C empirical targets:",
                definition
            )
    )


    abline(
        h = 0,
        lty = 3
    )


    abline(
        v = 0,
        lty = 3
    )


    for (i in seq_len(
        nrow(
            table
        )
    )) {

        arrows(

            0,
            0,

            table$E_score[
                i
            ],

            table$S_score[
                i
            ],

            length = 0.08
        )
    }


    points(
        table$E_score,
        table$S_score,
        pch = 19,
        cex = 1.4
    )


    labels <- table$Context


    labels[
        labels ==
            "Adult_RbTKO"
    ] <- "Adult Rb-TKO"


    labels[
        labels ==
            "CAG_RbLoss"
    ] <- "CAG Rb loss"


    labels[
        labels ==
            "MAP2_RbLoss"
    ] <- "MAP2 Rb loss"


    labels[
        labels ==
            "E2F13_DKO"
    ] <- "E2F1/3 DKO"


    text(
        table$E_score,
        table$S_score,
        labels = labels,
        pos = 4,
        offset = 0.6
    )


    dev.off()
}


# ============================================================
# PLOT: ADULT COUNTERFACTUAL TIME COURSES
# ============================================================

plot_counterfactuals <- function(
    definition,
    output_file
) {

    table <- counterfactual_simulations[
        counterfactual_simulations$Definition ==
            definition,
        ,
        drop = FALSE
    ]


    scenarios <- unique(
        table$Scenario
    )


    png(
        output_file,
        width = 1800,
        height = 1600,
        res = 220
    )


    par(
        mfrow = c(
            2,
            1
        ),
        mar = c(
            4.5,
            5,
            3,
            2
        )
    )


    E_range <- range(
        table$E,
        na.rm = TRUE
    )


    plot(

        NA,

        xlim =
            range(
                table$Time,
                na.rm = TRUE
            ),

        ylim =
            E_range,

        xlab =
            "Normalized time",

        ylab =
            "E",

        main =
            paste(
                definition,
                "- E dynamics"
            )
    )


    abline(
        h = 0,
        lty = 3
    )


    abline(
        v = 4,
        lty = 3
    )


    line_types <- seq_along(
        scenarios
    )


    for (i in seq_along(
        scenarios
    )) {

        one <- table[
            table$Scenario ==
                scenarios[
                    i
                ],
            ,
            drop = FALSE
        ]


        lines(
            one$Time,
            one$E,
            lty =
                line_types[
                    i
                ],
            lwd = 2
        )
    }


    legend(
        "bottomright",
        legend =
            scenarios,
        lty =
            line_types,
        lwd = 2,
        cex = 0.8,
        bty = "n"
    )


    S_range <- range(
        table$S,
        na.rm = TRUE
    )


    plot(

        NA,

        xlim =
            range(
                table$Time,
                na.rm = TRUE
            ),

        ylim =
            S_range,

        xlab =
            "Normalized time",

        ylab =
            "S",

        main =
            paste(
                definition,
                "- S dynamics"
            )
    )


    abline(
        h = 0,
        lty = 3
    )


    abline(
        v = 4,
        lty = 3
    )


    for (i in seq_along(
        scenarios
    )) {

        one <- table[
            table$Scenario ==
                scenarios[
                    i
                ],
            ,
            drop = FALSE
        ]


        lines(
            one$Time,
            one$S,
            lty =
                line_types[
                    i
                ],
            lwd = 2
        )
    }


    legend(
        "bottomright",
        legend =
            scenarios,
        lty =
            line_types,
        lwd = 2,
        cex = 0.8,
        bty = "n"
    )


    dev.off()
}


# ============================================================
# GENERATE FIGURES
# ============================================================

plot_state_targets(

    "Primary",

    file.path(
        OUTPUT_DIR,
        "F11C_empirical_targets_primary.png"
    )
)


plot_state_targets(

    "Sensitivity",

    file.path(
        OUTPUT_DIR,
        "F11C_empirical_targets_sensitivity.png"
    )
)


plot_counterfactuals(

    "Primary",

    file.path(
        OUTPUT_DIR,
        "F11C_counterfactual_dynamics_primary.png"
    )
)


plot_counterfactuals(

    "Sensitivity",

    file.path(
        OUTPUT_DIR,
        "F11C_counterfactual_dynamics_sensitivity.png"
    )
)


# ============================================================
# SAVE OUTPUTS
# ============================================================

readr::write_csv(

    observed_simulations,

    file.path(
        OUTPUT_DIR,
        "F11C_observed_target_relaxation.csv"
    )
)


readr::write_csv(

    counterfactual_simulations,

    file.path(
        OUTPUT_DIR,
        "F11C_counterfactual_simulations.csv"
    )
)


readr::write_csv(

    rate_sensitivity,

    file.path(
        OUTPUT_DIR,
        "F11C_rate_sensitivity.csv"
    )
)


readr::write_csv(

    rate_table,

    file.path(
        OUTPUT_DIR,
        "F11C_normalized_relaxation_times.csv"
    )
)


readr::write_csv(

    prediction_table,

    file.path(
        OUTPUT_DIR,
        "F11C_falsifiable_predictions.csv"
    )
)


readr::write_csv(

    assumption_table,

    file.path(
        OUTPUT_DIR,
        "F11C_model_assumptions.csv"
    )
)


# ============================================================
# CONSOLE SUMMARY
# ============================================================

message("")

message(
    "============================================================"
)

message(
    "F11C — MINIMAL DYNAMICS AND FALSIFIABLE PREDICTIONS"
)

message(
    "============================================================"
)


message("")

message(
    "EMPIRICAL NULL-CENTERED TARGETS"
)

message(
    "-------------------------------"
)


print(

    coordinates[
        ,
        c(
            "Definition",
            "Context",
            "E_score",
            "S_score"
        ),
        drop = FALSE
    ],

    row.names = FALSE
)


message("")

message(
    "NORMALIZED RATE SENSITIVITY"
)

message(
    "---------------------------"
)


print(
    rate_table,
    row.names = FALSE
)


message("")

message(
    "FALSIFIABLE PREDICTIONS"
)

message(
    "-----------------------"
)


print(
    prediction_table,
    row.names = FALSE
)


message("")

message(
    "MODEL ASSUMPTIONS"
)

message(
    "-----------------"
)


print(
    assumption_table,
    row.names = FALSE
)


message("")

message(
    "Outputs written to:"
)


message(
    OUTPUT_DIR
)


message("")

message(
    "INTERPRETATION GUARDRAILS:"
)


message(
    paste(
        "1. Trajectories are model-generated.",
        "Only their target states are empirically constrained."
    )
)


message(
    paste(
        "2. Normalized time has no direct mapping to",
        "minutes, hours, or days."
    )
)


message(
    paste(
        "3. kE and kS are sensitivity parameters,",
        "not estimated biological rate constants."
    )
)


message(
    paste(
        "4. Selective E and S normalization scenarios",
        "are counterfactual program-level interventions."
    )
)


message(
    paste(
        "5. Return to baseline after perturbation removal",
        "is a prediction of the memoryless model."
    )
)


message(
    paste(
        "6. Persistence, hysteresis, or obligatory E-S coupling",
        "would falsify the minimal architecture and justify",
        "a more complex model."
    )
)