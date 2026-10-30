#!/usr/bin/env Rscript

# ============================================================
# F10K_REL_conditional_test.R
# ============================================================
#
# PURPOSE
# -------
# F10J identified REL_ALT as the only GTRD TF-target set
# significantly enriched in Up_A2 after multiple testing.
#
# F10K asks whether that enrichment survives after removing
# genes already annotated to the broad immune/stress axis.
#
# This addresses the key alternative explanation:
#
#   Up_A2 is immune enriched
#   +
#   REL targets are immune enriched
#   --------------------------------
#   apparent REL enrichment
#
# PRIMARY TEST
# ------------
# Among NON-IMMUNE adult-Up genes:
#
#       Up_A2 vs Up_A1
#
# is REL_ALT target membership still enriched?
#
# SECONDARY TESTS
# ---------------
# REL, RELB, and NFKB1 are evaluated in parallel.
#
# IMPORTANT
# ---------
# Persistence of enrichment would support regulatory structure
# beyond simple immune annotation overlap.
#
# It would still NOT prove causal TF activity.
# ============================================================


suppressPackageStartupMessages({
    library(msigdbr)
    library(dplyr)
    library(readr)
    library(stringr)
})


# ============================================================
# PATHS
# ============================================================

F10G_FILE <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F10G_EmpiricalResponseArchetypes",
    "F10G_cluster_assignments.csv"
)

F10I_FILE <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F10I_UpA2FunctionalDecomposition",
    "F10I_gene_functional_membership.csv"
)

OUTPUT_DIR <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F10K_RELConditionalTest"
)

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# PARAMETERS
# ============================================================

REL_FAMILY <- c(
    "REL_ALT",
    "REL",
    "RELB",
    "NFKB1"
)


# ============================================================
# HELPERS
# ============================================================

safe_fraction <- function(numerator, denominator) {

    if (denominator == 0) {
        return(NA_real_)
    }

    numerator / denominator
}


run_fisher <- function(
    data,
    target_keys,
    label,
    condition_label
) {

    is_target <- data$Gene_key %in% target_keys

    is_a2 <- data$Archetype_key == "Up_A2"
    is_a1 <- data$Archetype_key == "Up_A1"

    a2_yes <- sum(
        is_a2 & is_target
    )

    a2_no <- sum(
        is_a2 & !is_target
    )

    a1_yes <- sum(
        is_a1 & is_target
    )

    a1_no <- sum(
        is_a1 & !is_target
    )

    contingency <- matrix(
        c(
            a2_yes,
            a2_no,
            a1_yes,
            a1_no
        ),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
            Archetype = c(
                "Up_A2",
                "Up_A1"
            ),
            Target = c(
                "Yes",
                "No"
            )
        )
    )

    if (any(rowSums(contingency) == 0)) {
        stop(
            paste0(
                "Cannot run Fisher test for ",
                label,
                " under ",
                condition_label,
                ": one archetype has zero genes."
            )
        )
    }

    test <- fisher.test(
        contingency
    )

    odds_ratio <- unname(
        test$estimate
    )

    ci_low <- unname(
        test$conf.int[1]
    )

    ci_high <- unname(
        test$conf.int[2]
    )

    tibble(
        TF_symbol = label,
        Condition = condition_label,

        Up_A2_targets = a2_yes,
        Up_A2_total = a2_yes + a2_no,

        Up_A1_targets = a1_yes,
        Up_A1_total = a1_yes + a1_no,

        Up_A2_fraction = safe_fraction(
            a2_yes,
            a2_yes + a2_no
        ),

        Up_A1_fraction = safe_fraction(
            a1_yes,
            a1_yes + a1_no
        ),

        Odds_ratio = odds_ratio,

        OR_CI95_low = ci_low,

        OR_CI95_high = ci_high,

        P_value = test$p.value
    )
}


# ============================================================
# LOAD F10G
# ============================================================

message(
    "Loading F10G..."
)

f10g <- read_csv(
    F10G_FILE,
    show_col_types = FALSE
) %>%
    filter(
        RbTKO_Direction == "Up"
    ) %>%
    mutate(
        Gene_key = toupper(
            trimws(
                Gene
            )
        )
    )


if (!setequal(
    unique(
        f10g$Archetype_key
    ),
    c(
        "Up_A1",
        "Up_A2"
    )
)) {

    stop(
        "F10K requires exactly Up_A1 and Up_A2 among adult-Up genes."
    )
}


message(
    "Adult-Up genes: ",
    nrow(
        f10g
    )
)

message(
    "  Up_A1: ",
    sum(
        f10g$Archetype_key == "Up_A1"
    )
)

message(
    "  Up_A2: ",
    sum(
        f10g$Archetype_key == "Up_A2"
    )
)


# ============================================================
# LOAD F10I FUNCTIONAL MEMBERSHIP
# ============================================================

message("")
message(
    "Loading F10I immune annotations..."
)

f10i <- read_csv(
    F10I_FILE,
    show_col_types = FALSE
)

required_f10i <- c(
    "Gene",
    "Immune_any",
    "Death_any"
)

missing_f10i <- setdiff(
    required_f10i,
    colnames(
        f10i
    )
)

if (length(
    missing_f10i
) > 0) {

    stop(
        paste0(
            "F10I membership file missing columns: ",
            paste(
                missing_f10i,
                collapse = ", "
            )
        )
    )
}


f10i <- f10i %>%
    mutate(
        Gene_key = toupper(
            trimws(
                Gene
            )
        )
    ) %>%
    select(
        Gene_key,
        Immune_any,
        Death_any
    ) %>%
    distinct(
        Gene_key,
        .keep_all = TRUE
    )


analysis_df <- f10g %>%
    left_join(
        f10i,
        by = "Gene_key"
    )


if (any(
    is.na(
        analysis_df$Immune_any
    )
)) {

    missing_annotation <- analysis_df %>%
        filter(
            is.na(
                Immune_any
            )
        )

    stop(
        paste0(
            nrow(
                missing_annotation
            ),
            " adult-Up genes are missing F10I immune annotations."
        )
    )
}


message(
    "Immune-annotated adult-Up genes: ",
    sum(
        analysis_df$Immune_any
    )
)

message(
    "Nonimmune adult-Up genes: ",
    sum(
        !analysis_df$Immune_any
    )
)


# ============================================================
# LOAD NATIVE MOUSE GTRD
# ============================================================

message("")
message(
    "Loading native mouse MSigDB M3/GTRD..."
)

gtrd_raw <- msigdbr(
    db_species = "MM",
    species = "Mus musculus",
    collection = "M3",
    subcollection = "GTRD"
)


gtrd <- gtrd_raw %>%
    transmute(
        Gene_set = gs_name,

        TF_symbol = str_remove(
            gs_name,
            "_TARGET_GENES$"
        ),

        Target_gene = gene_symbol,

        Target_key = toupper(
            trimws(
                gene_symbol
            )
        )
    ) %>%
    distinct(
        TF_symbol,
        Target_key
    )


# ============================================================
# CHECK REQUESTED REL SETS
# ============================================================

available_rel <- intersect(
    REL_FAMILY,
    unique(
        gtrd$TF_symbol
    )
)

missing_rel <- setdiff(
    REL_FAMILY,
    available_rel
)


message(
    "REL-family sets available: ",
    paste(
        available_rel,
        collapse = ", "
    )
)


if (length(
    missing_rel
) > 0) {

    message(
        "REL-family sets absent from GTRD: ",
        paste(
            missing_rel,
            collapse = ", "
        )
    )
}


if (length(
    available_rel
) == 0) {

    stop(
        "None of the requested REL-family GTRD sets are available."
    )
}


# ============================================================
# ANALYSIS STRATA
# ============================================================

all_adult_up <- analysis_df

nonimmune <- analysis_df %>%
    filter(
        !Immune_any
    )

nondeath <- analysis_df %>%
    filter(
        !Death_any
    )

nonimmune_nondeath <- analysis_df %>%
    filter(
        !Immune_any,
        !Death_any
    )


message("")
message(
    "Analysis strata:"
)

message(
    "  All adult-Up: ",
    nrow(
        all_adult_up
    )
)

message(
    "  Nonimmune: ",
    nrow(
        nonimmune
    )
)

message(
    "  Nondeath: ",
    nrow(
        nondeath
    )
)

message(
    "  Nonimmune + nondeath: ",
    nrow(
        nonimmune_nondeath
    )
)


# ============================================================
# RUN REL-FAMILY TESTS
# ============================================================

results <- list()


for (tf in available_rel) {

    target_keys <- gtrd %>%
        filter(
            TF_symbol == tf
        ) %>%
        pull(
            Target_key
        ) %>%
        unique()


    results[[
        paste0(
            tf,
            "__ALL"
        )
    ]] <- run_fisher(
        all_adult_up,
        target_keys,
        tf,
        "ALL_ADULT_UP"
    )


    results[[
        paste0(
            tf,
            "__NONIMMUNE"
        )
    ]] <- run_fisher(
        nonimmune,
        target_keys,
        tf,
        "AFTER_REMOVING_IMMUNE"
    )


    results[[
        paste0(
            tf,
            "__NONDEATH"
        )
    ]] <- run_fisher(
        nondeath,
        target_keys,
        tf,
        "AFTER_REMOVING_DEATH"
    )


    results[[
        paste0(
            tf,
            "__NONIMMUNE_NONDEATH"
        )
    ]] <- run_fisher(
        nonimmune_nondeath,
        target_keys,
        tf,
        "AFTER_REMOVING_IMMUNE_AND_DEATH"
    )
}


result_table <- bind_rows(
    results
)


# ============================================================
# MULTIPLE-TEST CORRECTION
#
# Two FDRs are reported:
#
#   FDR_global:
#       correction across every TF x condition test
#
#   FDR_within_condition:
#       correction across TFs inside each analysis stratum
#
# The primary pre-specified test is:
#
#   REL_ALT / AFTER_REMOVING_IMMUNE
#
# Its raw Fisher p-value should therefore also be inspected
# directly rather than relying only on exploratory family-wide
# correction.
# ============================================================

result_table <- result_table %>%
    mutate(
        FDR_global = p.adjust(
            P_value,
            method = "BH"
        )
    ) %>%
    group_by(
        Condition
    ) %>%
    mutate(
        FDR_within_condition = p.adjust(
            P_value,
            method = "BH"
        )
    ) %>%
    ungroup() %>%
    arrange(
        Condition,
        FDR_within_condition,
        P_value
    )


# ============================================================
# TARGET OVERLAP DECOMPOSITION FOR REL_ALT
# ============================================================

rel_alt_targets <- gtrd %>%
    filter(
        TF_symbol == "REL_ALT"
    ) %>%
    pull(
        Target_key
    ) %>%
    unique()


rel_alt_gene_audit <- analysis_df %>%
    mutate(
        REL_ALT_target = Gene_key %in% rel_alt_targets
    ) %>%
    filter(
        REL_ALT_target
    ) %>%
    select(
        Gene,
        Archetype_key,
        Immune_any,
        Death_any,
        RbTKO_LFC,
        E2F13_LFC,
        CAG_LFC,
        MAP2_LFC,
        Adult_residual_z
    ) %>%
    arrange(
        Archetype_key,
        desc(
            Adult_residual_z
        )
    )


rel_alt_overlap_summary <- rel_alt_gene_audit %>%
    count(
        Archetype_key,
        Immune_any,
        Death_any,
        name = "N_genes"
    ) %>%
    group_by(
        Archetype_key
    ) %>%
    mutate(
        Fraction_within_REL_ALT_targets =
            N_genes / sum(
                N_genes
            )
    ) %>%
    ungroup()


# ============================================================
# SAVE
# ============================================================

write_csv(
    result_table,
    file.path(
        OUTPUT_DIR,
        "F10K_REL_family_conditional_tests.csv"
    )
)

write_csv(
    rel_alt_gene_audit,
    file.path(
        OUTPUT_DIR,
        "F10K_REL_ALT_gene_audit.csv"
    )
)

write_csv(
    rel_alt_overlap_summary,
    file.path(
        OUTPUT_DIR,
        "F10K_REL_ALT_annotation_overlap.csv"
    )
)


# ============================================================
# CONSOLE OUTPUT
# ============================================================

message("")
message(
    "============================================================"
)

message(
    "F10K — REL CONDITIONAL TEST"
)

message(
    "============================================================"
)

message("")

print(
    as.data.frame(
        result_table %>%
            select(
                TF_symbol,
                Condition,
                Up_A2_targets,
                Up_A2_total,
                Up_A1_targets,
                Up_A1_total,
                Up_A2_fraction,
                Up_A1_fraction,
                Odds_ratio,
                OR_CI95_low,
                OR_CI95_high,
                P_value,
                FDR_within_condition,
                FDR_global
            )
    ),
    row.names = FALSE
)


message("")
message(
    "REL_ALT TARGET ANNOTATION COMPOSITION"
)

message(
    "-------------------------------------"
)

print(
    as.data.frame(
        rel_alt_overlap_summary
    ),
    row.names = FALSE
)


message("")
message(
    "PRIMARY TEST"
)

message(
    "------------"
)

primary_result <- result_table %>%
    filter(
        TF_symbol == "REL_ALT",
        Condition == "AFTER_REMOVING_IMMUNE"
    )


print(
    as.data.frame(
        primary_result
    ),
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
    "INTERPRETATION GUARDRAIL:"
)

message(
    paste(
        "Persistence of REL_ALT enrichment after immune genes are",
        "removed would support regulatory structure beyond simple",
        "immune-annotation overlap, but would still not demonstrate",
        "causal c-Rel activity in adult neural stem cells."
    )
)