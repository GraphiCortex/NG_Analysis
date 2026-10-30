#!/usr/bin/env Rscript

# ============================================================
# F11D_experimental_falsification_matrix.R
# ============================================================
#
# PURPOSE
# -------
#
# Close the F11-v2 modeling section by translating the
# empirically constrained architecture into an explicit
# experimental discrimination / falsification framework.
#
#
# MODELS
# ------
#
# M0:
#   Independent context-routed branches
#
#       perturbation -> E
#       adult context + perturbation -> S
#
#
# M1:
#   E -> S coupling
#
#
# M2:
#   S -> E coupling
#
#
# M3:
#   Bidirectional E <-> S coupling
#
#
# M4:
#   State-memory / feedback architecture
#
#
# IMPORTANT
# ---------
#
# These are competing MODEL ARCHITECTURES.
#
# The prediction matrix below is not new biological evidence.
#
# It formalizes:
#
#   - what each model predicts;
#   - which experiments discriminate them;
#   - what observations would falsify the current minimal model.
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

F11B_COORDINATES <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F11B_StateSpaceCoupling",
    "F11B_state_space_coordinates.csv"
)


F11B_CONSTRAINTS <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F11B_StateSpaceCoupling",
    "F11B_structural_constraints.csv"
)


F11C_PREDICTIONS <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F11C_MinimalDynamics",
    "F11C_falsifiable_predictions.csv"
)


OUTPUT_DIR <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F11D_ExperimentalFalsification"
)


dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# LOAD PRIOR RESULTS
# ============================================================

message(
    "Loading F11B/F11C outputs..."
)


coordinates <- readr::read_csv(
    F11B_COORDINATES,
    show_col_types = FALSE
)


constraints <- readr::read_csv(
    F11B_CONSTRAINTS,
    show_col_types = FALSE
)


f11c_predictions <- readr::read_csv(
    F11C_PREDICTIONS,
    show_col_types = FALSE
)


# ============================================================
# MODEL DEFINITIONS
# ============================================================

model_definitions <- data.frame(

    Model =
        c(
            "M0",
            "M1",
            "M2",
            "M3",
            "M4"
        ),

    Architecture =
        c(
            "Independent context-routed branches",
            "E drives S",
            "S drives E",
            "Bidirectional E-S coupling",
            "State memory / feedback"
        ),

    Core_equation =
        c(
            "E and S have independent perturbation/context inputs",
            "S depends obligatorily on E",
            "E depends obligatorily on S",
            "E and S mutually influence one another",
            "At least one state can persist after initiating input is removed"
        ),

    Current_status =
        c(
            "Minimal model currently sufficient",
            "Not required by endpoint data",
            "Not required by endpoint data",
            "Not required by endpoint data",
            "Not supported by endpoint data"
        ),

    stringsAsFactors = FALSE
)


# ============================================================
# EXPERIMENT DEFINITIONS
# ============================================================

experiments <- data.frame(

    Experiment_ID =
        c(
            "X1",
            "X2",
            "X3",
            "X4",
            "X5",
            "X6"
        ),

    Experiment =
        c(
            paste(
                "Adult Rb loss plus E2F-program suppression"
            ),

            paste(
                "Adult Rb loss plus selective S-program suppression"
            ),

            paste(
                "Reversible removal of the initiating perturbation"
            ),

            paste(
                "Dense time-course after perturbation onset"
            ),

            paste(
                "Matched Rb loss across adult and embryonic",
                "neural contexts"
            ),

            paste(
                "Adult E2F1/3 loss versus adult Rb loss",
                "with simultaneous E and S measurement"
            )
        ),

    Primary_question =
        c(
            paste(
                "Does S require E in the adult Rb-loss state?"
            ),

            paste(
                "Does E require the adult immune/stress state?"
            ),

            paste(
                "Do E or S persist after the initiating input disappears?"
            ),

            paste(
                "Is there obligatory temporal ordering between E and S?"
            ),

            paste(
                "Is S determined by context despite similar Rb-loss E?"
            ),

            paste(
                "Can opposite E states coexist with similar positive S?"
            )
        ),

    Feasibility =
        c(
            "HIGH",
            "MEDIUM",
            "MEDIUM",
            "HIGH",
            "HIGH",
            "HIGH"
        ),

    Primary_readout =
        c(
            "E_CORE and S_IMMUNE program scores",
            "E_CORE and S_IMMUNE program scores",
            "E_CORE and S_IMMUNE over time after withdrawal",
            "Repeated E_CORE and S_IMMUNE measurements",
            "E_CORE and S_IMMUNE across contexts",
            "E_CORE and S_IMMUNE in both adult perturbations"
        ),

    stringsAsFactors = FALSE
)


# ============================================================
# MODEL PREDICTION MATRIX
# ============================================================
#
# Codes:
#
#   E_DOWN_S_PERSISTS
#   S_DOWN_E_PERSISTS
#   BOTH_DOWN
#   RETURN_BASELINE
#   PERSISTENCE_POSSIBLE
#   E_PRECEDES_S
#   S_PRECEDES_E
#   NO_REQUIRED_ORDER
#   CONTEXT_DEPENDENT_S
#   E_AND_S_COVARY
#   OPPOSITE_E_SIMILAR_S_ALLOWED
#
# These are qualitative model predictions.
# ============================================================

prediction_matrix <- data.frame(

    Experiment_ID =
        rep(
            experiments$Experiment_ID,
            each = 5
        ),

    Model =
        rep(
            model_definitions$Model,
            times = nrow(
                experiments
            )
        ),

    Prediction =
        c(

            # ------------------------------------------------
            # X1
            # Adult Rb loss + E suppression
            # ------------------------------------------------

            "E_DOWN_S_PERSISTS",
            "BOTH_DOWN",
            "E_DOWN_S_PERSISTS",
            "BOTH_DOWN",
            "MODEL_DEPENDENT",


            # ------------------------------------------------
            # X2
            # Adult Rb loss + S suppression
            # ------------------------------------------------

            "S_DOWN_E_PERSISTS",
            "S_DOWN_E_PERSISTS",
            "BOTH_DOWN",
            "BOTH_DOWN",
            "MODEL_DEPENDENT",


            # ------------------------------------------------
            # X3
            # Perturbation removal
            # ------------------------------------------------

            "RETURN_BASELINE",
            "RETURN_BASELINE",
            "RETURN_BASELINE",
            "RETURN_BASELINE",
            "PERSISTENCE_POSSIBLE",


            # ------------------------------------------------
            # X4
            # Time-course
            # ------------------------------------------------

            "NO_REQUIRED_ORDER",
            "E_PRECEDES_S",
            "S_PRECEDES_E",
            "COUPLED_OR_COMPLEX_ORDER",
            "HISTORY_DEPENDENT_DYNAMICS",


            # ------------------------------------------------
            # X5
            # Matched Rb loss adult vs embryonic
            # ------------------------------------------------

            "CONTEXT_DEPENDENT_S",
            "CONTEXT_MODULATED_COUPLING_REQUIRED",
            "STRONG_S_TO_E_PREDICTS_E_CONTEXT_DIFFERENCE",
            "CONTEXT_MODULATED_COUPLING_REQUIRED",
            "CONTEXT_AND_HISTORY_DEPENDENT",


            # ------------------------------------------------
            # X6
            # Adult E2F DKO versus adult Rb loss
            # ------------------------------------------------

            "OPPOSITE_E_SIMILAR_S_ALLOWED",
            "E_AND_S_COVARY",
            "S_CAN_REMAIN_SIMILAR_WITH_DIFFERENT_E",
            "E_AND_S_EXPECTED_TO_COVARY",
            "MODEL_DEPENDENT"
        ),

    stringsAsFactors = FALSE
)


# ============================================================
# DISCRIMINATION POWER
# ============================================================
#
# Count the number of distinct qualitative predictions made
# across models for each experiment.
#
# This is a heuristic ranking only.
# ============================================================

discrimination <- data.frame()


for (experiment_id in experiments$Experiment_ID) {

    rows <- prediction_matrix[
        prediction_matrix$Experiment_ID ==
            experiment_id,
        ,
        drop = FALSE
    ]


    unique_predictions <- unique(
        rows$Prediction
    )


    experiment_row <- experiments[
        experiments$Experiment_ID ==
            experiment_id,
        ,
        drop = FALSE
    ]


    one_row <- data.frame(

        Experiment_ID =
            experiment_id,

        Experiment =
            experiment_row$Experiment[
                1
            ],

        Feasibility =
            experiment_row$Feasibility[
                1
            ],

        Distinct_model_predictions =
            length(
                unique_predictions
            ),

        stringsAsFactors = FALSE
    )


    discrimination <- rbind(
        discrimination,
        one_row
    )
}


# ============================================================
# PRIORITY SCORE
# ============================================================
#
# Feasibility:
#
#   HIGH   = 3
#   MEDIUM = 2
#   LOW    = 1
#
# Priority is deliberately simple:
#
#   distinct predictions * feasibility weight
#
# This is NOT a biological statistic.
# ============================================================

feasibility_weight <- rep(
    NA_real_,
    nrow(
        discrimination
    )
)


feasibility_weight[
    discrimination$Feasibility ==
        "HIGH"
] <- 3


feasibility_weight[
    discrimination$Feasibility ==
        "MEDIUM"
] <- 2


feasibility_weight[
    discrimination$Feasibility ==
        "LOW"
] <- 1


discrimination$Feasibility_weight <-
    feasibility_weight


discrimination$Priority_score <- (
    discrimination$Distinct_model_predictions *
    discrimination$Feasibility_weight
)


discrimination <- discrimination[
    order(
        -discrimination$Priority_score,
        -discrimination$Distinct_model_predictions
    ),
    ,
    drop = FALSE
]


rownames(
    discrimination
) <- NULL


# ============================================================
# CURRENT EMPIRICAL CONSTRAINTS
# ============================================================

primary <- coordinates[
    coordinates$Definition ==
        "Primary",
    ,
    drop = FALSE
]


get_coordinate <- function(
    context,
    variable
) {

    row <- primary[
        primary$Context ==
            context,
        ,
        drop = FALSE
    ]


    if (nrow(row) != 1) {

        return(
            NA_real_
        )
    }


    as.numeric(
        row[
            ,
            variable,
            drop = TRUE
        ][1]
    )
}


adult_E <- get_coordinate(
    "Adult_RbTKO",
    "E_score"
)


adult_S <- get_coordinate(
    "Adult_RbTKO",
    "S_score"
)


cag_E <- get_coordinate(
    "CAG_RbLoss",
    "E_score"
)


cag_S <- get_coordinate(
    "CAG_RbLoss",
    "S_score"
)


map2_E <- get_coordinate(
    "MAP2_RbLoss",
    "E_score"
)


map2_S <- get_coordinate(
    "MAP2_RbLoss",
    "S_score"
)


e2f_E <- get_coordinate(
    "E2F13_DKO",
    "E_score"
)


e2f_S <- get_coordinate(
    "E2F13_DKO",
    "S_score"
)


current_evidence <- data.frame(

    Observation =
        c(
            "Adult Rb loss produces positive E",
            "Adult Rb loss produces positive S",
            "CAG Rb loss produces positive E",
            "CAG Rb loss does not produce positive excess S",
            "MAP2 Rb loss produces positive E",
            "MAP2 Rb loss does not produce positive excess S",
            "Adult E2F1/3 loss produces negative E",
            "Adult E2F1/3 loss retains positive S",
            "Adult Rb loss and E2F1/3 loss have opposite E but similar S",
            "F11B does not identify obligatory E-S coupling"
        ),

    Value =
        c(
            adult_E,
            adult_S,
            cag_E,
            cag_S,
            map2_E,
            map2_S,
            e2f_E,
            e2f_S,
            abs(
                adult_S -
                e2f_S
            ),
            NA_real_
        ),

    Interpretation =
        c(
            "Observed/inferred endpoint",
            "Observed/inferred endpoint",
            "Observed/inferred endpoint",
            "Observed/inferred endpoint",
            "Observed/inferred endpoint",
            "Observed/inferred endpoint",
            "Observed/inferred endpoint",
            "Observed/inferred endpoint",
            "Strong evidence against simple monotonic E-only explanation of S",
            "Independent-branch model remains sufficient"
        ),

    stringsAsFactors = FALSE
)


# ============================================================
# FALSIFICATION RULES FOR M0
# ============================================================

m0_falsification <- data.frame(

    Test =
        c(
            "Adult Rb loss plus E suppression",
            "Adult Rb loss plus S suppression",
            "Perturbation withdrawal",
            "Time-course ordering",
            "Cross-context replication"
        ),

    M0_prediction =
        c(
            "E decreases while S can persist",
            "S decreases while E can persist",
            "Both states relax toward baseline",
            "No obligatory E-before-S or S-before-E ordering",
            "E remains broadly conserved while S remains context dependent"
        ),

    M0_falsified_if =
        c(
            "S reproducibly collapses whenever E is selectively suppressed",
            "E reproducibly collapses whenever S is selectively suppressed",
            "Either state remains persistently displaced after complete removal of input",
            "One program reproducibly precedes and is required for activation of the other",
            "Matched contexts show S is fully predicted by E without independent context effect"
        ),

    stringsAsFactors = FALSE
)


# ============================================================
# MODEL-COMPARISON SUMMARY
# ============================================================

model_summary <- data.frame(

    Model =
        c(
            "M0",
            "M1",
            "M2",
            "M3",
            "M4"
        ),

    Current_empirical_status =
        c(
            "SUFFICIENT",
            "NOT_REQUIRED",
            "NOT_REQUIRED",
            "NOT_REQUIRED",
            "NOT_SUPPORTED"
        ),

    Why =
        c(
            paste(
                "Explains conserved E and context-routed S",
                "without invoking unsupported coupling."
            ),

            paste(
                "F11B found no robust common monotonic",
                "E-to-S relationship."
            ),

            paste(
                "Large S differences coexist with relatively",
                "stable E across Rb-loss contexts."
            ),

            paste(
                "No endpoint evidence requires mutual coupling."
            ),

            paste(
                "No persistence or hysteresis experiment",
                "has demonstrated state memory."
            )
        ),

    stringsAsFactors = FALSE
)


# ============================================================
# SIMPLE DISCRIMINATION FIGURE
# ============================================================

plot_priority <- function(
    table,
    output_file
) {

    ordered <- table[
        order(
            table$Priority_score
        ),
        ,
        drop = FALSE
    ]


    labels <- paste(
        ordered$Experiment_ID,
        ordered$Experiment,
        sep = ": "
    )


    png(
        output_file,
        width = 1800,
        height = 1300,
        res = 220
    )


    par(
        mar = c(
            5,
            14,
            4,
            2
        ),
        las = 1
    )


    barplot(

        ordered$Priority_score,

        names.arg =
            labels,

        horiz =
            TRUE,

        xlab =
            "Heuristic discrimination priority",

        main =
            "F11D experimental model-discrimination priority",

        cex.names =
            0.8
    )


    dev.off()
}


plot_priority(

    discrimination,

    file.path(
        OUTPUT_DIR,
        "F11D_experiment_priority.png"
    )
)


# ============================================================
# SAVE OUTPUTS
# ============================================================

readr::write_csv(

    model_definitions,

    file.path(
        OUTPUT_DIR,
        "F11D_model_definitions.csv"
    )
)


readr::write_csv(

    experiments,

    file.path(
        OUTPUT_DIR,
        "F11D_experiment_definitions.csv"
    )
)


readr::write_csv(

    prediction_matrix,

    file.path(
        OUTPUT_DIR,
        "F11D_model_prediction_matrix.csv"
    )
)


readr::write_csv(

    discrimination,

    file.path(
        OUTPUT_DIR,
        "F11D_experiment_priority.csv"
    )
)


readr::write_csv(

    current_evidence,

    file.path(
        OUTPUT_DIR,
        "F11D_current_empirical_constraints.csv"
    )
)


readr::write_csv(

    m0_falsification,

    file.path(
        OUTPUT_DIR,
        "F11D_M0_falsification_rules.csv"
    )
)


readr::write_csv(

    model_summary,

    file.path(
        OUTPUT_DIR,
        "F11D_model_comparison_summary.csv"
    )
)


# ============================================================
# CONSOLE REPORT
# ============================================================

message("")

message(
    "============================================================"
)

message(
    "F11D — EXPERIMENTAL FALSIFICATION MATRIX"
)

message(
    "============================================================"
)


message("")

message(
    "CURRENT MODEL STATUS"
)

message(
    "--------------------"
)


print(
    model_summary,
    row.names = FALSE
)


message("")

message(
    "CURRENT EMPIRICAL CONSTRAINTS"
)

message(
    "-----------------------------"
)


print(
    current_evidence,
    row.names = FALSE
)


message("")

message(
    "EXPERIMENT PRIORITY"
)

message(
    "-------------------"
)


print(
    discrimination,
    row.names = FALSE
)


message("")

message(
    "M0 FALSIFICATION RULES"
)

message(
    "----------------------"
)


print(
    m0_falsification,
    row.names = FALSE
)


message("")

message(
    "MODEL PREDICTION MATRIX"
)

message(
    "-----------------------"
)


print(
    prediction_matrix,
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
        "1. The prediction matrix encodes model assumptions;",
        "it is not new empirical evidence."
    )
)


message(
    paste(
        "2. Priority scores are heuristic experimental-design",
        "scores, not statistical quantities."
    )
)


message(
    paste(
        "3. M0 is currently sufficient, not proven true."
    )
)


message(
    paste(
        "4. Experiments X1-X4 are explicitly designed",
        "to falsify or extend M0."
    )
)


message(
    paste(
        "5. A more complex model should only be adopted",
        "if an experiment demonstrates obligatory coupling,",
        "temporal dependence, or state persistence."
    )
)