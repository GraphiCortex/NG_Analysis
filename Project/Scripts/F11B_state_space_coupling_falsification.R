#!/usr/bin/env Rscript

# ============================================================
# F11B_state_space_coupling_falsification.R
# ============================================================
#
# Uses F11A empirical state variables to:
#
#   1) plot the four perturbations in (E,S) state space;
#
#   2) quantify the geometry of the three Rb-loss contexts;
#
#   3) test two deliberately restrictive explanations:
#
#        R1: E alone explains cross-context S variation
#        R2: S alone explains cross-context E variation
#
#   4) evaluate their stability with:
#
#        - leave-one-context-out prediction
#        - gene bootstrap
#        - pairwise required gains
#
#
# IMPORTANT
# ---------
#
# These slopes are endpoint geometry.
#
# They are NOT:
#
#   - biochemical coupling constants
#   - kinetic rate constants
#   - causal effect estimates
#
# We have perturbational endpoints, not time-series data.
#
# Therefore this analysis cannot establish:
#
#   - bistability
#   - hysteresis
#   - oscillation
#   - time constants
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

INPUT_DIR <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F11A_EmpiricalStateVariables"
)


OUTPUT_DIR <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F11B_StateSpaceCoupling"
)


dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


SCORE_FILE <- file.path(
    INPUT_DIR,
    "F11A_state_variable_scores.csv"
)


DIFF_FILE <- file.path(
    INPUT_DIR,
    "F11A_context_differences.csv"
)


MEMBERSHIP_FILE <- file.path(
    INPUT_DIR,
    "F11A_program_membership.csv"
)


STATE_FILE <- file.path(
    INPUT_DIR,
    "F11A_gene_state_evidence_matrix.csv"
)


# ============================================================
# PARAMETERS
# ============================================================

N_BOOT <- 5000

RANDOM_SEED <- 20261029


RB_CONTEXTS <- c(
    "Adult_RbTKO",
    "CAG_RbLoss",
    "MAP2_RbLoss"
)


ALL_CONTEXTS <- c(
    "Adult_RbTKO",
    "CAG_RbLoss",
    "MAP2_RbLoss",
    "E2F13_DKO"
)


set.seed(
    RANDOM_SEED
)


# ============================================================
# HELPER: R-SQUARED
# ============================================================

safe_r2 <- function(
    x,
    y
) {

    ok <- (
        is.finite(x) &
        is.finite(y)
    )


    x <- x[
        ok
    ]


    y <- y[
        ok
    ]


    if (
        length(x) < 2 ||
        sd(x) < 1e-12 ||
        sd(y) < 1e-12
    ) {

        return(
            NA_real_
        )
    }


    cor(
        x,
        y
    ) ^ 2
}


# ============================================================
# HELPER: SIMPLE LINEAR GEOMETRY
# ============================================================

linear_fit <- function(
    x,
    y
) {

    ok <- (
        is.finite(x) &
        is.finite(y)
    )


    x <- x[
        ok
    ]


    y <- y[
        ok
    ]


    if (
        length(x) < 2 ||
        var(x) < 1e-12
    ) {

        return(
            data.frame(
                Slope = NA_real_,
                Intercept = NA_real_,
                R2 = NA_real_,
                RMSE = NA_real_,
                stringsAsFactors = FALSE
            )
        )
    }


    slope <- (
        cov(
            x,
            y
        ) /
        var(
            x
        )
    )


    intercept <- (
        mean(
            y
        ) -
        slope *
        mean(
            x
        )
    )


    fitted <- (
        intercept +
        slope *
        x
    )


    rmse <- sqrt(
        mean(
            (
                y -
                fitted
            ) ^ 2
        )
    )


    data.frame(
        Slope = slope,
        Intercept = intercept,
        R2 = safe_r2(
            x,
            y
        ),
        RMSE = rmse,
        stringsAsFactors = FALSE
    )
}


# ============================================================
# HELPER: LEAVE-ONE-CONTEXT-OUT
# ============================================================

loo_fit <- function(
    x,
    y,
    labels
) {

    n <- length(
        x
    )


    prediction <- rep(
        NA_real_,
        n
    )


    error <- rep(
        NA_real_,
        n
    )


    training_slope <- rep(
        NA_real_,
        n
    )


    for (i in seq_len(
        n
    )) {

        train <- seq_len(
            n
        )


        train <- train[
            train != i
        ]


        if (length(train) != 2) {

            next
        }


        dx <- (
            x[
                train[2]
            ] -
            x[
                train[1]
            ]
        )


        if (
            !is.finite(dx) ||
            abs(dx) < 1e-12
        ) {

            next
        }


        slope <- (
            y[
                train[2]
            ] -
            y[
                train[1]
            ]
        ) / dx


        intercept <- (
            y[
                train[1]
            ] -
            slope *
            x[
                train[1]
            ]
        )


        prediction[
            i
        ] <- (
            intercept +
            slope *
            x[
                i
            ]
        )


        error[
            i
        ] <- (
            prediction[
                i
            ] -
            y[
                i
            ]
        )


        training_slope[
            i
        ] <- slope
    }


    finite_error <- error[
        is.finite(
            error
        )
    ]


    if (length(finite_error) > 0) {

        rmse <- sqrt(
            mean(
                finite_error ^ 2
            )
        )


        mae <- mean(
            abs(
                finite_error
            )
        )


        max_error <- max(
            abs(
                finite_error
            )
        )

    } else {

        rmse <- NA_real_

        mae <- NA_real_

        max_error <- NA_real_
    }


    outcome_range <- diff(
        range(
            y,
            na.rm = TRUE
        )
    )


    normalized_rmse <- if (
        is.finite(
            outcome_range
        ) &&
        outcome_range > 1e-12
    ) {

        rmse /
            outcome_range

    } else {

        NA_real_
    }


    list(

        details =
            data.frame(

                Held_out =
                    labels,

                Prediction =
                    prediction,

                Observed =
                    y,

                Error =
                    error,

                Training_slope =
                    training_slope,

                stringsAsFactors = FALSE
            ),


        summary =
            data.frame(

                LOO_RMSE =
                    rmse,

                LOO_MAE =
                    mae,

                LOO_MaxAbsError =
                    max_error,

                Outcome_range =
                    outcome_range,

                LOO_Normalized_RMSE =
                    normalized_rmse,

                stringsAsFactors = FALSE
            )
    )
}


# ============================================================
# HELPER: PAIRWISE REQUIRED GAIN
# ============================================================

pairwise_gains <- function(
    x,
    y,
    labels,
    direction_name
) {

    pairs <- combn(
        seq_along(
            labels
        ),
        2
    )


    out <- data.frame()


    for (i in seq_len(
        ncol(
            pairs
        )
    )) {

        a <- pairs[
            1,
            i
        ]


        b <- pairs[
            2,
            i
        ]


        dx <- (
            x[
                a
            ] -
            x[
                b
            ]
        )


        dy <- (
            y[
                a
            ] -
            y[
                b
            ]
        )


        gain <- if (
            is.finite(
                dx
            ) &&
            abs(dx) > 1e-12
        ) {

            dy /
                dx

        } else {

            NA_real_
        }


        one_row <- data.frame(

            Direction =
                direction_name,

            Context_1 =
                labels[
                    a
                ],

            Context_2 =
                labels[
                    b
                ],

            Delta_input =
                dx,

            Delta_output =
                dy,

            Pairwise_gain =
                gain,

            stringsAsFactors = FALSE
        )


        out <- rbind(
            out,
            one_row
        )
    }


    out
}


# ============================================================
# HELPER: SAFE QUANTILE
# ============================================================

q_safe <- function(
    x,
    p
) {

    x <- x[
        is.finite(
            x
        )
    ]


    if (length(x) == 0) {

        return(
            NA_real_
        )
    }


    as.numeric(
        quantile(
            x,
            p,
            na.rm = TRUE
        )
    )
}


# ============================================================
# LOAD F11A OUTPUTS
# ============================================================

message(
    "Loading F11A outputs..."
)


scores <- readr::read_csv(
    SCORE_FILE,
    show_col_types = FALSE
)


differences <- readr::read_csv(
    DIFF_FILE,
    show_col_types = FALSE
)


membership <- readr::read_csv(
    MEMBERSHIP_FILE,
    show_col_types = FALSE
)


state_data <- readr::read_csv(
    STATE_FILE,
    show_col_types = FALSE
)


required_state <- c(
    "Gene_key",
    ALL_CONTEXTS
)


missing_state <- setdiff(
    required_state,
    names(
        state_data
    )
)


if (length(missing_state) > 0) {

    stop(
        paste(
            "State matrix missing:",
            paste(
                missing_state,
                collapse = ", "
            )
        )
    )
}


state_matrix <- as.matrix(
    state_data[
        ,
        ALL_CONTEXTS,
        drop = FALSE
    ]
)


storage.mode(
    state_matrix
) <- "double"


rownames(
    state_matrix
) <- as.character(
    state_data$Gene_key
)


# ============================================================
# SCORE LOOKUP
# ============================================================

get_score_row <- function(
    program,
    context
) {

    out <- scores[
        scores$Program ==
            program &
        scores$Context ==
            context,
        ,
        drop = FALSE
    ]


    if (nrow(out) != 1) {

        stop(
            paste(
                "Expected one score row for",
                program,
                context
            )
        )
    }


    out
}


# ============================================================
# BUILD STATE-SPACE COORDINATES
# ============================================================

make_coordinates <- function(
    E_program,
    S_program,
    definition
) {

    out <- data.frame()


    for (context in ALL_CONTEXTS) {

        e <- get_score_row(
            E_program,
            context
        )


        s <- get_score_row(
            S_program,
            context
        )


        one_row <- data.frame(

    Definition =
        definition,

    Context =
        context,

    E_program =
        E_program,

    E_raw_score =
        e$Mean_score,

    E_null_mean =
        e$Null_mean,

    E_score =
        e$Mean_score -
        e$Null_mean,

    E_CI_low =
        e$Bootstrap_CI_low -
        e$Null_mean,

    E_CI_high =
        e$Bootstrap_CI_high -
        e$Null_mean,

    S_program =
        S_program,

    S_raw_score =
        s$Mean_score,

    S_null_mean =
        s$Null_mean,

    S_score =
        s$Mean_score -
        s$Null_mean,

    S_CI_low =
        s$Bootstrap_CI_low -
        s$Null_mean,

    S_CI_high =
        s$Bootstrap_CI_high -
        s$Null_mean,

    stringsAsFactors = FALSE
)

        out <- rbind(
            out,
            one_row
        )
    }


    out
}


primary_coordinates <- make_coordinates(
    "E_CORE",
    "S_IMMUNE",
    "Primary"
)


sensitivity_coordinates <- make_coordinates(
    "E_BROAD",
    "S_INNATE",
    "Sensitivity"
)


coordinates <- rbind(
    primary_coordinates,
    sensitivity_coordinates
)


# ============================================================
# RESTRICTED ENDPOINT MODELS
# ============================================================

analyze_definition <- function(
    table,
    definition
) {

    rb <- table[
        table$Context %in%
            RB_CONTEXTS,
        ,
        drop = FALSE
    ]


    rb <- rb[
        match(
            RB_CONTEXTS,
            rb$Context
        ),
        ,
        drop = FALSE
    ]


    E <- rb$E_score

    S <- rb$S_score

    labels <- rb$Context


    # --------------------------------------------------------
    # R1:
    #
    # E alone explains S.
    # --------------------------------------------------------

    fit_es <- linear_fit(
        E,
        S
    )


    loo_es <- loo_fit(
        E,
        S,
        labels
    )


    # --------------------------------------------------------
    # R2:
    #
    # S alone explains E.
    # --------------------------------------------------------

    fit_se <- linear_fit(
        S,
        E
    )


    loo_se <- loo_fit(
        S,
        E,
        labels
    )


    summary <- rbind(

        data.frame(

            Definition =
                definition,

            Restricted_model =
                "E_ONLY_EXPLAINS_S",

            Input =
                "E",

            Outcome =
                "S",

            Slope =
                fit_es$Slope,

            Intercept =
                fit_es$Intercept,

            R2 =
                fit_es$R2,

            InSample_RMSE =
                fit_es$RMSE,

            LOO_RMSE =
                loo_es$summary$LOO_RMSE,

            LOO_MAE =
                loo_es$summary$LOO_MAE,

            LOO_MaxAbsError =
                loo_es$summary$LOO_MaxAbsError,

            Outcome_range =
                loo_es$summary$Outcome_range,

            LOO_Normalized_RMSE =
                loo_es$summary$LOO_Normalized_RMSE,

            stringsAsFactors = FALSE
        ),


        data.frame(

            Definition =
                definition,

            Restricted_model =
                "S_ONLY_EXPLAINS_E",

            Input =
                "S",

            Outcome =
                "E",

            Slope =
                fit_se$Slope,

            Intercept =
                fit_se$Intercept,

            R2 =
                fit_se$R2,

            InSample_RMSE =
                fit_se$RMSE,

            LOO_RMSE =
                loo_se$summary$LOO_RMSE,

            LOO_MAE =
                loo_se$summary$LOO_MAE,

            LOO_MaxAbsError =
                loo_se$summary$LOO_MaxAbsError,

            Outcome_range =
                loo_se$summary$Outcome_range,

            LOO_Normalized_RMSE =
                loo_se$summary$LOO_Normalized_RMSE,

            stringsAsFactors = FALSE
        )
    )


    loo_es$details$Definition <-
        definition


    loo_es$details$Restricted_model <-
        "E_ONLY_EXPLAINS_S"


    loo_se$details$Definition <-
        definition


    loo_se$details$Restricted_model <-
        "S_ONLY_EXPLAINS_E"


    details <- rbind(
        loo_es$details,
        loo_se$details
    )


    gain_es <- pairwise_gains(
        E,
        S,
        labels,
        "E_TO_S"
    )


    gain_se <- pairwise_gains(
        S,
        E,
        labels,
        "S_TO_E"
    )


    gain_es$Definition <-
        definition


    gain_se$Definition <-
        definition


    list(

        summary =
            summary,

        details =
            details,

        gains =
            rbind(
                gain_es,
                gain_se
            )
    )
}


primary_geometry <- analyze_definition(
    primary_coordinates,
    "Primary"
)


sensitivity_geometry <- analyze_definition(
    sensitivity_coordinates,
    "Sensitivity"
)


geometry_summary <- rbind(
    primary_geometry$summary,
    sensitivity_geometry$summary
)


loo_details <- rbind(
    primary_geometry$details,
    sensitivity_geometry$details
)


gain_table <- rbind(
    primary_geometry$gains,
    sensitivity_geometry$gains
)


# ============================================================
# GENE-BOOTSTRAP STABILITY
# ============================================================
#
# Genes are resampled within each program.
#
# This quantifies program-level gene stability.
#
# It does NOT quantify biological-replicate uncertainty.
# ============================================================

bootstrap_definition <- function(
    E_program,
    S_program,
    definition
) {

    if (!E_program %in%
        names(
            membership
        )) {

        stop(
            paste(
                "Missing membership column:",
                E_program
            )
        )
    }


    if (!S_program %in%
        names(
            membership
        )) {

        stop(
            paste(
                "Missing membership column:",
                S_program
            )
        )
    }


    E_member <- as.logical(
        membership[
            ,
            E_program,
            drop = TRUE
        ]
    )


    S_member <- as.logical(
        membership[
            ,
            S_program,
            drop = TRUE
        ]
    )


    E_genes <- as.character(
        membership$Gene_key[
            !is.na(
                E_member
            ) &
            E_member
        ]
    )


    S_genes <- as.character(
        membership$Gene_key[
            !is.na(
                S_member
            ) &
            S_member
        ]
    )


    E_genes <- intersect(
        E_genes,
        rownames(
            state_matrix
        )
    )


    S_genes <- intersect(
        S_genes,
        rownames(
            state_matrix
        )
    )


    if (
        length(E_genes) < 5 ||
        length(S_genes) < 5
    ) {

        stop(
            paste(
                "Too few program genes for",
                definition
            )
        )
    }


    # --------------------------------------------------------
    # MATCH Rb-loss contexts to the state matrix
    # --------------------------------------------------------

    rb_cols <- match(
        RB_CONTEXTS,
        colnames(
            state_matrix
        )
    )


    if (any(
        is.na(
            rb_cols
        )
    )) {

        stop(
            "One or more RB_CONTEXTS are absent from state_matrix."
        )
    }


    # --------------------------------------------------------
    # NULL-CENTERING VALUES
    #
    # For each program and context:
    #
    #   excess score = observed program score - matched-null mean
    #
    # This removes context-specific background offsets before
    # evaluating cross-context E-S geometry.
    # --------------------------------------------------------

    E_null <- rep(
        NA_real_,
        length(
            RB_CONTEXTS
        )
    )


    S_null <- rep(
        NA_real_,
        length(
            RB_CONTEXTS
        )
    )


    for (j in seq_along(
        RB_CONTEXTS
    )) {

        context <- RB_CONTEXTS[
            j
        ]


        E_row <- get_score_row(
            E_program,
            context
        )


        S_row <- get_score_row(
            S_program,
            context
        )


        E_null[
            j
        ] <- E_row$Null_mean


        S_null[
            j
        ] <- S_row$Null_mean
    }


    # --------------------------------------------------------
    # OUTPUT TABLE
    # --------------------------------------------------------

    out <- data.frame(

        Iteration =
            seq_len(
                N_BOOT
            ),

        Definition =
            definition,

        E_to_S_slope =
            NA_real_,

        E_to_S_R2 =
            NA_real_,

        E_to_S_LOO_Normalized_RMSE =
            NA_real_,

        S_to_E_slope =
            NA_real_,

        S_to_E_R2 =
            NA_real_,

        S_to_E_LOO_Normalized_RMSE =
            NA_real_,

        stringsAsFactors = FALSE
    )


    # --------------------------------------------------------
    # GENE BOOTSTRAP
    # --------------------------------------------------------

    for (i in seq_len(
        N_BOOT
    )) {

        sampled_E <- sample(

            E_genes,

            size =
                length(
                    E_genes
                ),

            replace =
                TRUE
        )


        sampled_S <- sample(

            S_genes,

            size =
                length(
                    S_genes
                ),

            replace =
                TRUE
        )


        # ----------------------------------------------------
        # Raw bootstrap program means
        # ----------------------------------------------------

        E_score <- colMeans(

            state_matrix[
                sampled_E,
                rb_cols,
                drop = FALSE
            ],

            na.rm = TRUE
        )


        S_score <- colMeans(

            state_matrix[
                sampled_S,
                rb_cols,
                drop = FALSE
            ],

            na.rm = TRUE
        )


        # ----------------------------------------------------
        # NULL-CENTER THE BOOTSTRAP SCORES
        #
        # These are now excess program scores relative to the
        # expression-matched expectation from F11A.
        # ----------------------------------------------------

        E_score <- (
            E_score -
            E_null
        )


        S_score <- (
            S_score -
            S_null
        )


        # ----------------------------------------------------
        # RESTRICTED MODEL 1:
        #
        # E alone explains S
        # ----------------------------------------------------

        fit_es <- linear_fit(
            E_score,
            S_score
        )


        loo_es <- loo_fit(
            E_score,
            S_score,
            RB_CONTEXTS
        )


        # ----------------------------------------------------
        # RESTRICTED MODEL 2:
        #
        # S alone explains E
        # ----------------------------------------------------

        fit_se <- linear_fit(
            S_score,
            E_score
        )


        loo_se <- loo_fit(
            S_score,
            E_score,
            RB_CONTEXTS
        )


        # ----------------------------------------------------
        # STORE RESULTS
        # ----------------------------------------------------

        out$E_to_S_slope[
            i
        ] <- fit_es$Slope


        out$E_to_S_R2[
            i
        ] <- fit_es$R2


        out$E_to_S_LOO_Normalized_RMSE[
            i
        ] <- loo_es$summary$LOO_Normalized_RMSE


        out$S_to_E_slope[
            i
        ] <- fit_se$Slope


        out$S_to_E_R2[
            i
        ] <- fit_se$R2


        out$S_to_E_LOO_Normalized_RMSE[
            i
        ] <- loo_se$summary$LOO_Normalized_RMSE
    }


    out
}

message(
    "Bootstrapping primary state-space geometry..."
)


boot_primary <- bootstrap_definition(
    "E_CORE",
    "S_IMMUNE",
    "Primary"
)


message(
    "Bootstrapping sensitivity state-space geometry..."
)


boot_sensitivity <- bootstrap_definition(
    "E_BROAD",
    "S_INNATE",
    "Sensitivity"
)


bootstrap_all <- rbind(
    boot_primary,
    boot_sensitivity
)


# ============================================================
# SUMMARIZE BOOTSTRAP
# ============================================================

summarize_boot <- function(
    df,
    definition,
    direction,
    slope_col,
    r2_col,
    loo_col
) {

    slope <- df[
        ,
        slope_col,
        drop = TRUE
    ]


    r2 <- df[
        ,
        r2_col,
        drop = TRUE
    ]


    loo <- df[
        ,
        loo_col,
        drop = TRUE
    ]


    finite_slope <- slope[
        is.finite(
            slope
        )
    ]


    finite_r2 <- r2[
        is.finite(
            r2
        )
    ]


    finite_loo <- loo[
        is.finite(
            loo
        )
    ]


    data.frame(

        Definition =
            definition,

        Direction =
            direction,

        Finite_fraction =
            mean(
                is.finite(
                    slope
                )
            ),

        Slope_median =
            median(
                finite_slope,
                na.rm = TRUE
            ),

        Slope_CI_low =
            q_safe(
                finite_slope,
                0.025
            ),

        Slope_CI_high =
            q_safe(
                finite_slope,
                0.975
            ),

        Fraction_slope_positive =
            mean(
                finite_slope > 0,
                na.rm = TRUE
            ),

        R2_median =
            median(
                finite_r2,
                na.rm = TRUE
            ),

        LOO_Normalized_RMSE_median =
            median(
                finite_loo,
                na.rm = TRUE
            ),

        LOO_Normalized_RMSE_CI_low =
            q_safe(
                finite_loo,
                0.025
            ),

        LOO_Normalized_RMSE_CI_high =
            q_safe(
                finite_loo,
                0.975
            ),

        stringsAsFactors = FALSE
    )
}


bootstrap_summary <- rbind(

    summarize_boot(
        boot_primary,
        "Primary",
        "E_TO_S",
        "E_to_S_slope",
        "E_to_S_R2",
        "E_to_S_LOO_Normalized_RMSE"
    ),


    summarize_boot(
        boot_primary,
        "Primary",
        "S_TO_E",
        "S_to_E_slope",
        "S_to_E_R2",
        "S_to_E_LOO_Normalized_RMSE"
    ),


    summarize_boot(
        boot_sensitivity,
        "Sensitivity",
        "E_TO_S",
        "E_to_S_slope",
        "E_to_S_R2",
        "E_to_S_LOO_Normalized_RMSE"
    ),


    summarize_boot(
        boot_sensitivity,
        "Sensitivity",
        "S_TO_E",
        "S_to_E_slope",
        "S_to_E_R2",
        "S_to_E_LOO_Normalized_RMSE"
    )
)


# ============================================================
# STRUCTURAL CONSTRAINTS FROM F11A
# ============================================================

get_difference <- function(
    program,
    contrast
) {

    out <- differences[
        differences$Program ==
            program &
        differences$Contrast ==
            contrast,
        ,
        drop = FALSE
    ]


    if (nrow(out) != 1) {

        stop(
            paste(
                "Expected one difference row for",
                program,
                contrast
            )
        )
    }


    out
}


E_AC <- get_difference(
    "E_CORE",
    "Adult_minus_CAG"
)


E_AM <- get_difference(
    "E_CORE",
    "Adult_minus_MAP2"
)


S_AC <- get_difference(
    "S_IMMUNE",
    "Adult_minus_CAG"
)


S_AM <- get_difference(
    "S_IMMUNE",
    "Adult_minus_MAP2"
)


E_AC_separated <- (
    E_AC$Bootstrap_CI_low > 0 ||
    E_AC$Bootstrap_CI_high < 0
)


E_AM_separated <- (
    E_AM$Bootstrap_CI_low > 0 ||
    E_AM$Bootstrap_CI_high < 0
)


S_AC_positive <- (
    S_AC$Bootstrap_CI_low > 0
)


S_AM_positive <- (
    S_AM$Bootstrap_CI_low > 0
)


structural_constraints <- data.frame(

    Statement =
        c(
            "Core E differs between Adult and CAG Rb loss",
            "Core E differs between Adult and MAP2 Rb loss",
            "Immune S is higher in Adult than CAG Rb loss",
            "Immune S is higher in Adult than MAP2 Rb loss",
            "Context-routed S input is sufficient without E-S coupling",
            "E-only explanation of S is identifiable from endpoints",
            "S-only explanation of E is required by endpoints",
            "Bistability or hysteresis is identified by endpoints"
        ),


    Supported =
        c(
            E_AC_separated,
            E_AM_separated,
            S_AC_positive,
            S_AM_positive,
            TRUE,
            FALSE,
            FALSE,
            FALSE
        ),


    Interpretation =
        c(
            "F11A gene-bootstrap CI",
            "F11A gene-bootstrap CI",
            "Adult-context S separation",
            "Adult-context S separation",
            "Minimal independent-branch architecture can reproduce endpoint routing",
            "Only three Rb-loss endpoints; direct context input remains an alternative",
            "Core E remains similar despite substantial S variation",
            "No switching/time-series/hysteresis experiment is available"
        ),


    stringsAsFactors = FALSE
)


# ============================================================
# STATE-SPACE FIGURES
# ============================================================

pretty_context <- function(x) {

    out <- x


    out[
        x == "Adult_RbTKO"
    ] <- "Adult Rb-TKO"


    out[
        x == "CAG_RbLoss"
    ] <- "CAG Rb loss"


    out[
        x == "MAP2_RbLoss"
    ] <- "MAP2 Rb loss"


    out[
        x == "E2F13_DKO"
    ] <- "E2F1/3 DKO"


    out
}


plot_state_space <- function(
    table,
    output_file,
    xlab_text,
    ylab_text,
    title_text
) {

    x_values <- c(
        table$E_CI_low,
        table$E_CI_high
    )


    y_values <- c(
        table$S_CI_low,
        table$S_CI_high
    )


    x_range <- range(
        x_values,
        na.rm = TRUE
    )


    y_range <- range(
        y_values,
        na.rm = TRUE
    )


    x_margin <- (
        0.12 *
        diff(
            x_range
        )
    )


    y_margin <- (
        0.12 *
        diff(
            y_range
        )
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
            xlab_text,

        ylab =
            ylab_text,

        main =
            title_text
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

        segments(
            table$E_CI_low[
                i
            ],
            table$S_score[
                i
            ],
            table$E_CI_high[
                i
            ],
            table$S_score[
                i
            ]
        )


        segments(
            table$E_score[
                i
            ],
            table$S_CI_low[
                i
            ],
            table$E_score[
                i
            ],
            table$S_CI_high[
                i
            ]
        )
    }


    symbols <- c(
        19,
        17,
        15,
        18
    )


    symbols <- symbols[
        match(
            table$Context,
            ALL_CONTEXTS
        )
    ]


    points(
        table$E_score,
        table$S_score,
        pch = symbols,
        cex = 1.5
    )


    text(
        table$E_score,
        table$S_score,
        labels = pretty_context(
            table$Context
        ),
        pos = 4,
        offset = 0.7,
        cex = 0.9
    )


    dev.off()
}


plot_state_space(

    primary_coordinates,

    file.path(
        OUTPUT_DIR,
        "F11B_state_space_primary.png"
    ),

    "E_CORE score",

    "S_IMMUNE score",

    "F11B primary empirical state space"
)


plot_state_space(

    sensitivity_coordinates,

    file.path(
        OUTPUT_DIR,
        "F11B_state_space_sensitivity.png"
    ),

    "E_BROAD score",

    "S_INNATE score",

    "F11B sensitivity empirical state space"
)


# ============================================================
# SAVE OUTPUTS
# ============================================================

readr::write_csv(

    coordinates,

    file.path(
        OUTPUT_DIR,
        "F11B_state_space_coordinates.csv"
    )
)


readr::write_csv(

    geometry_summary,

    file.path(
        OUTPUT_DIR,
        "F11B_restricted_coupling_models.csv"
    )
)


readr::write_csv(

    loo_details,

    file.path(
        OUTPUT_DIR,
        "F11B_leave_one_context_out.csv"
    )
)


readr::write_csv(

    gain_table,

    file.path(
        OUTPUT_DIR,
        "F11B_pairwise_required_gains.csv"
    )
)


readr::write_csv(

    bootstrap_all,

    file.path(
        OUTPUT_DIR,
        "F11B_gene_bootstrap_coupling.csv"
    )
)


readr::write_csv(

    bootstrap_summary,

    file.path(
        OUTPUT_DIR,
        "F11B_gene_bootstrap_summary.csv"
    )
)


readr::write_csv(

    structural_constraints,

    file.path(
        OUTPUT_DIR,
        "F11B_structural_constraints.csv"
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
    "F11B — STATE-SPACE GEOMETRY AND COUPLING FALSIFICATION"
)

message(
    "============================================================"
)


message("")

message(
    "STATE-SPACE COORDINATES"
)

message(
    "-----------------------"
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
    "RESTRICTED COUPLING MODELS"
)

message(
    "--------------------------"
)


print(
    geometry_summary,
    row.names = FALSE
)


message("")

message(
    "LEAVE-ONE-CONTEXT-OUT"
)

message(
    "---------------------"
)


print(
    loo_details,
    row.names = FALSE
)


message("")

message(
    "PAIRWISE REQUIRED GAINS"
)

message(
    "-----------------------"
)


print(
    gain_table,
    row.names = FALSE
)


message("")

message(
    "GENE-BOOTSTRAP COUPLING STABILITY"
)

message(
    "---------------------------------"
)


print(
    bootstrap_summary,
    row.names = FALSE
)


message("")

message(
    "STRUCTURAL CONSTRAINTS"
)

message(
    "----------------------"
)


print(
    structural_constraints,
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
    "1. Slopes are endpoint geometry, not biochemical rate constants."
)


message(
    "2. With only three Rb-loss contexts, in-sample R2 is weak evidence."
)


message(
    "3. Leave-one-context-out stability is more informative than R2 here."
)


message(
    "4. Gene bootstrap measures program-level gene stability, not replicate uncertainty."
)


message(
    "5. Direct context input to S remains an alternative to E->S coupling."
)


message(
    "6. Endpoints cannot establish bistability, hysteresis, or kinetic parameters."
)