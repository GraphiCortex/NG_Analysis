#!/usr/bin/env Rscript

# ============================================================
# F10J_regulatory_driver_screen.R
# ============================================================
#
# PURPOSE
# -------
# Identify transcription-factor target programs that distinguish
# the empirically derived adult-Up response archetypes:
#
#   Up_A1 = conserved E2F/cell-cycle response
#   Up_A2 = adult-context immune/stimulus response
#
# PRIMARY QUESTION
# ----------------
# Among genes already upregulated after adult Rb-family loss,
# which GTRD TF-target sets are preferentially represented in
# Up_A2 versus Up_A1?
#
# IMPORTANT INTERPRETATION
# ------------------------
# GTRD target-set enrichment indicates enrichment of genes with
# promoter binding evidence for a TF.
#
# It does NOT by itself establish:
#   - TF activation,
#   - edge direction,
#   - causal regulation,
#   - direct regulation in adult neural stem cells.
#
# Therefore TF-target enrichment and TF expression behavior are
# kept as separate evidence layers.
#
# ============================================================


# ============================================================
# PACKAGES
# ============================================================

required_packages <- c(
    "msigdbr",
    "dplyr",
    "readr",
    "stringr"
)

missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (length(missing_packages) > 0) {
    stop(
        paste0(
            "Missing R package(s): ",
            paste(
                missing_packages,
                collapse = ", "
            ),
            "\nInstall them before running F10J."
        )
    )
}

suppressPackageStartupMessages({
    library(msigdbr)
    library(dplyr)
    library(readr)
    library(stringr)
})


# ============================================================
# PATHS
# ============================================================

INPUT_FILE <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F10G_EmpiricalResponseArchetypes",
    "F10G_cluster_assignments.csv"
)

OUTPUT_DIR <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F10J_RegulatoryDriverScreen"
)

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# PARAMETERS
# ============================================================

MIN_TARGETS_IN_UNIVERSE <- 10

MAX_TARGETS_IN_UNIVERSE <- 500

FDR_THRESHOLD <- 0.05


# ============================================================
# HELPERS
# ============================================================

safe_fraction <- function(numerator, denominator) {
    if (denominator == 0) {
        return(NA_real_)
    }

    numerator / denominator
}


safe_collapse <- function(x) {
    x <- sort(
        unique(
            na.omit(
                x
            )
        )
    )

    if (length(x) == 0) {
        return("")
    }

    paste(
        x,
        collapse = ";"
    )
}


# ============================================================
# LOAD F10G
# ============================================================

message(
    "Loading F10G cluster assignments..."
)

df <- read_csv(
    INPUT_FILE,
    show_col_types = FALSE
)

required_columns <- c(
    "Gene",
    "RbTKO_Direction",
    "Archetype_key",
    "RbTKO_LFC",
    "E2F13_LFC",
    "CAG_LFC",
    "MAP2_LFC",
    "Adult_residual_z",
    "Adult_E2F_reciprocal_strong",
    "Rb_loss_sign_conserved"
)

missing_columns <- setdiff(
    required_columns,
    colnames(
        df
    )
)

if (length(missing_columns) > 0) {
    stop(
        paste0(
            "F10G table is missing required columns: ",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    )
}


# ============================================================
# NORMALIZE GENE KEYS
# ============================================================

df <- df %>%
    mutate(
        Gene_key = toupper(
            trimws(
                Gene
            )
        )
    )


# ============================================================
# PRIMARY UNIVERSE = ALL ADULT-UP COMPLETE-CONTEXT GENES
# ============================================================

adult_up <- df %>%
    filter(
        RbTKO_Direction == "Up"
    )

if (!setequal(
    unique(
        adult_up$Archetype_key
    ),
    c(
        "Up_A1",
        "Up_A2"
    )
)) {
    stop(
        paste0(
            "F10J expects exactly Up_A1 and Up_A2 among adult-Up genes. ",
            "Observed: ",
            paste(
                sort(
                    unique(
                        adult_up$Archetype_key
                    )
                ),
                collapse = ", "
            )
        )
    )
}

n_a1 <- sum(
    adult_up$Archetype_key == "Up_A1"
)

n_a2 <- sum(
    adult_up$Archetype_key == "Up_A2"
)

message(
    "Adult-Up analysis universe: ",
    nrow(
        adult_up
    ),
    " genes"
)

message(
    "  Up_A1: ",
    n_a1
)

message(
    "  Up_A2: ",
    n_a2
)


# ============================================================
# CHECK CURRENT MOUSE MSigDB GTRD COLLECTION
# ============================================================

message("")
message(
    "Checking native mouse MSigDB collections..."
)

collections <- msigdbr_collections(
    db_species = "MM"
)

# Native Mouse MSigDB uses:
#
#   M3 = regulatory target gene sets
#   GTRD = transcription-factor target subcollection
#
# Human MSigDB instead uses C3 / TFT:GTRD.
# We explicitly use the native mouse M3/GTRD collection here.

gtrd_collection_rows <- collections %>%
    filter(
        gs_collection == "M3",
        gs_subcollection == "GTRD"
    )

if (nrow(gtrd_collection_rows) == 0) {

    message("")
    message(
        "Available mouse regulatory collections containing 'GTRD':"
    )

    diagnostic_rows <- collections %>%
        filter(
            grepl(
                "GTRD",
                gs_subcollection,
                ignore.case = TRUE
            )
        )

    print(
        as.data.frame(
            diagnostic_rows
        ),
        row.names = FALSE
    )

    stop(
        paste(
            "Native mouse MSigDB M3/GTRD was not found.",
            "See the diagnostic collection table above."
        )
    )
}

message(
    "Native mouse GTRD collection found: M3 / GTRD"
)

if ("num_genesets" %in% colnames(
    gtrd_collection_rows
)) {

    message(
        "Reported GTRD gene sets: ",
        paste(
            unique(
                gtrd_collection_rows$num_genesets
            ),
            collapse = "; "
        )
    )
}

if ("db_version" %in% colnames(
    gtrd_collection_rows
)) {

    message(
        "MSigDB version: ",
        paste(
            unique(
                gtrd_collection_rows$db_version
            ),
            collapse = "; "
        )
    )
}


# ============================================================
# LOAD NATIVE-MOUSE GTRD TF-TARGET SETS
# ============================================================

message("")
message(
    "Loading native mouse M3/GTRD gene sets..."
)

gtrd_raw <- msigdbr(
    db_species = "MM",
    species = "Mus musculus",
    collection = "M3",
    subcollection = "GTRD"
)

if (nrow(gtrd_raw) == 0) {
    stop(
        "No native-mouse M3/GTRD gene-set records were returned."
    )
}

required_gtrd_columns <- c(
    "gene_symbol",
    "gs_name"
)

missing_gtrd_columns <- setdiff(
    required_gtrd_columns,
    colnames(
        gtrd_raw
    )
)

if (length(missing_gtrd_columns) > 0) {
    stop(
        paste0(
            "Unexpected msigdbr schema. Missing: ",
            paste(
                missing_gtrd_columns,
                collapse = ", "
            )
        )
    )
}

message(
    "Raw GTRD memberships: ",
    nrow(
        gtrd_raw
    )
)

message(
    "GTRD gene sets: ",
    n_distinct(
        gtrd_raw$gs_name
    )
)

if (nrow(gtrd_raw) == 0) {
    stop(
        "No GTRD gene-set records were returned."
    )
}

required_gtrd_columns <- c(
    "gene_symbol",
    "gs_name"
)

missing_gtrd_columns <- setdiff(
    required_gtrd_columns,
    colnames(
        gtrd_raw
    )
)

if (length(missing_gtrd_columns) > 0) {
    stop(
        paste0(
            "Unexpected msigdbr schema. Missing: ",
            paste(
                missing_gtrd_columns,
                collapse = ", "
            )
        )
    )
}

message(
    "Raw GTRD memberships: ",
    nrow(
        gtrd_raw
    )
)

message(
    "GTRD gene sets: ",
    n_distinct(
        gtrd_raw$gs_name
    )
)


# ============================================================
# PREPARE GTRD TABLE
# ============================================================

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
        ),

        Source_species = if (
            "gs_source_species" %in% colnames(
                gtrd_raw
            )
        ) {
            gs_source_species
        } else {
            NA_character_
        },

        Description = if (
            "gs_description" %in% colnames(
                gtrd_raw
            )
        ) {
            gs_description
        } else {
            NA_character_
        },

        Exact_source = if (
            "gs_exact_source" %in% colnames(
                gtrd_raw
            )
        ) {
            gs_exact_source
        } else {
            NA_character_
        },

        DB_version = if (
            "db_version" %in% colnames(
                gtrd_raw
            )
        ) {
            db_version
        } else {
            NA_character_
        }
    ) %>%
    distinct(
        Gene_set,
        TF_symbol,
        Target_key,
        .keep_all = TRUE
    )


# ============================================================
# LIMIT TARGET SETS TO THE ACTUAL ADULT-UP UNIVERSE
# ============================================================

universe_keys <- unique(
    adult_up$Gene_key
)

gtrd_universe <- gtrd %>%
    filter(
        Target_key %in% universe_keys
    )

set_sizes <- gtrd_universe %>%
    group_by(
        Gene_set,
        TF_symbol
    ) %>%
    summarise(
        N_targets_in_universe = n_distinct(
            Target_key
        ),

        Source_species = safe_collapse(
            Source_species
        ),

        DB_version = safe_collapse(
            DB_version
        ),

        .groups = "drop"
    )

eligible_sets <- set_sizes %>%
    filter(
        N_targets_in_universe >= MIN_TARGETS_IN_UNIVERSE,
        N_targets_in_universe <= MAX_TARGETS_IN_UNIVERSE
    )

message("")
message(
    "GTRD sets with >= ",
    MIN_TARGETS_IN_UNIVERSE,
    " and <= ",
    MAX_TARGETS_IN_UNIVERSE,
    " targets in the adult-Up universe: ",
    nrow(
        eligible_sets
    )
)

if (nrow(
    eligible_sets
) == 0) {
    stop(
        "No GTRD TF-target sets passed the universe-size filter."
    )
}


# ============================================================
# FISHER TEST EACH TF TARGET SET:
#
#                  TARGET     NOT TARGET
#   Up_A2             a            b
#   Up_A1             c            d
#
# OR > 1  -> targets preferentially represented in Up_A2
# OR < 1  -> targets preferentially represented in Up_A1
#
# This is a direct archetype comparison, conditional on already
# being an adult Rb-TKO-upregulated gene.
# ============================================================

message("")
message(
    "Running TF-target enrichment tests..."
)

test_one_set <- function(gene_set_name, tf_symbol) {

    targets <- gtrd_universe %>%
        filter(
            Gene_set == gene_set_name
        ) %>%
        pull(
            Target_key
        ) %>%
        unique()

    in_target <- adult_up$Gene_key %in% targets

    is_a2 <- adult_up$Archetype_key == "Up_A2"
    is_a1 <- adult_up$Archetype_key == "Up_A1"

    a2_target <- sum(
        is_a2 & in_target
    )

    a2_not <- sum(
        is_a2 & !in_target
    )

    a1_target <- sum(
        is_a1 & in_target
    )

    a1_not <- sum(
        is_a1 & !in_target
    )

    contingency <- matrix(
        c(
            a2_target,
            a2_not,
            a1_target,
            a1_not
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

    fisher_result <- fisher.test(
        contingency,
        alternative = "two.sided"
    )

    odds_ratio <- if (
        length(
            fisher_result$estimate
        ) == 1
    ) {
        unname(
            fisher_result$estimate
        )
    } else {
        NA_real_
    }

    ci_low <- if (
        length(
            fisher_result$conf.int
        ) >= 2
    ) {
        unname(
            fisher_result$conf.int[1]
        )
    } else {
        NA_real_
    }

    ci_high <- if (
        length(
            fisher_result$conf.int
        ) >= 2
    ) {
        unname(
            fisher_result$conf.int[2]
        )
    } else {
        NA_real_
    }

    a2_overlap_genes <- adult_up %>%
        filter(
            Archetype_key == "Up_A2",
            Gene_key %in% targets
        ) %>%
        pull(
            Gene
        )

    a1_overlap_genes <- adult_up %>%
        filter(
            Archetype_key == "Up_A1",
            Gene_key %in% targets
        ) %>%
        pull(
            Gene
        )

    tibble(
        Gene_set = gene_set_name,
        TF_symbol = tf_symbol,

        Targets_in_adult_up =
            length(
                targets
            ),

        Up_A2_targets =
            a2_target,

        Up_A2_total =
            n_a2,

        Up_A1_targets =
            a1_target,

        Up_A1_total =
            n_a1,

        Up_A2_fraction =
            safe_fraction(
                a2_target,
                n_a2
            ),

        Up_A1_fraction =
            safe_fraction(
                a1_target,
                n_a1
            ),

        Odds_ratio =
            odds_ratio,

        OR_CI95_low =
            ci_low,

        OR_CI95_high =
            ci_high,

        P_value =
            fisher_result$p.value,

        Up_A2_overlap_genes =
            safe_collapse(
                a2_overlap_genes
            ),

        Up_A1_overlap_genes =
            safe_collapse(
                a1_overlap_genes
            )
    )
}


test_results <- bind_rows(
    lapply(
        seq_len(
            nrow(
                eligible_sets
            )
        ),
        function(i) {
            test_one_set(
                eligible_sets$Gene_set[[i]],
                eligible_sets$TF_symbol[[i]]
            )
        }
    )
)


# ============================================================
# MULTIPLE-TEST CORRECTION + DIRECTION
# ============================================================

test_results <- test_results %>%
    mutate(
        FDR = p.adjust(
            P_value,
            method = "BH"
        ),

        Enrichment_direction = case_when(
            Odds_ratio > 1 ~
                "Up_A2",

            Odds_ratio < 1 ~
                "Up_A1",

            TRUE ~
                "Neutral"
        ),

        log2_OR = case_when(
            is.finite(
                Odds_ratio
            ) &
            Odds_ratio > 0 ~
                log2(
                    Odds_ratio
                ),

            is.infinite(
                Odds_ratio
            ) &
            Odds_ratio > 0 ~
                Inf,

            TRUE ~
                NA_real_
        )
    ) %>%
    left_join(
        eligible_sets,
        by = c(
            "Gene_set",
            "TF_symbol"
        )
    ) %>%
    arrange(
        FDR,
        P_value
    )


# ============================================================
# TF EXPRESSION / PERTURBATION AUDIT
#
# IMPORTANT:
# This is deliberately kept separate from target enrichment.
#
# Absence of TF differential expression does NOT refute TF
# activity because TFs can be regulated post-transcriptionally.
# Conversely, TF differential expression does NOT prove that
# its GTRD targets are causally regulated by it here.
# ============================================================

tf_expression_lookup <- df %>%
    transmute(
        TF_key = Gene_key,

        TF_gene_in_F10G =
            Gene,

        TF_RbTKO_direction =
            RbTKO_Direction,

        TF_archetype =
            Archetype_key,

        TF_RbTKO_LFC =
            RbTKO_LFC,

        TF_E2F13_LFC =
            E2F13_LFC,

        TF_CAG_LFC =
            CAG_LFC,

        TF_MAP2_LFC =
            MAP2_LFC,

        TF_Adult_residual_z =
            Adult_residual_z,

        TF_Adult_E2F_reciprocal_strong =
            Adult_E2F_reciprocal_strong,

        TF_Rb_loss_sign_conserved =
            Rb_loss_sign_conserved
    ) %>%
    distinct(
        TF_key,
        .keep_all = TRUE
    )


results_annotated <- test_results %>%
    mutate(
        TF_key = toupper(
            trimws(
                TF_symbol
            )
        )
    ) %>%
    left_join(
        tf_expression_lookup,
        by = "TF_key"
    ) %>%
    mutate(
        TF_expression_status = case_when(
            is.na(
                TF_gene_in_F10G
            ) ~
                "TF_not_in_F10G_universe",

            TF_archetype == "Up_A2" ~
                "TF_itself_Up_A2",

            TF_archetype == "Up_A1" ~
                "TF_itself_Up_A1",

            TF_RbTKO_direction == "Down" ~
                "TF_adult_down",

            TRUE ~
                "TF_present_other"
        )
    ) %>%
    arrange(
        FDR,
        P_value
    )


# ============================================================
# SIGNIFICANT SETS
# ============================================================

significant <- results_annotated %>%
    filter(
        FDR < FDR_THRESHOLD
    )

significant_a2 <- significant %>%
    filter(
        Enrichment_direction == "Up_A2"
    ) %>%
    arrange(
        FDR,
        desc(
            Odds_ratio
        )
    )

significant_a1 <- significant %>%
    filter(
        Enrichment_direction == "Up_A1"
    ) %>%
    arrange(
        FDR,
        Odds_ratio
    )


# ============================================================
# EMPIRICAL PRIORITY TABLE
#
# This is NOT a causal-driver score.
#
# Priority means only:
#   1. target set is significantly enriched in Up_A2;
#   2. enrichment is reasonably sized;
#   3. TF's own perturbation behavior is displayed for inspection.
#
# No bonus is assigned merely because a TF itself is DE.
# ============================================================

priority_table <- significant_a2 %>%
    mutate(
        Evidence_class = case_when(
            FDR < 0.001 &
            Odds_ratio >= 2 ~
                "Strong_target_enrichment",

            FDR < 0.05 &
            Odds_ratio >= 1.5 ~
                "Moderate_target_enrichment",

            TRUE ~
                "Significant_target_enrichment"
        )
    ) %>%
    select(
        TF_symbol,
        Gene_set,
        Evidence_class,
        Targets_in_adult_up,
        Up_A2_targets,
        Up_A1_targets,
        Up_A2_fraction,
        Up_A1_fraction,
        Odds_ratio,
        OR_CI95_low,
        OR_CI95_high,
        P_value,
        FDR,
        TF_expression_status,
        TF_RbTKO_direction,
        TF_archetype,
        TF_RbTKO_LFC,
        TF_E2F13_LFC,
        TF_CAG_LFC,
        TF_MAP2_LFC,
        TF_Adult_residual_z,
        TF_Adult_E2F_reciprocal_strong,
        TF_Rb_loss_sign_conserved,
        Source_species,
        DB_version,
        Up_A2_overlap_genes
    )


# ============================================================
# SAVE OUTPUTS
# ============================================================

write_csv(
    results_annotated,
    file.path(
        OUTPUT_DIR,
        "F10J_all_GTRD_TF_tests.csv"
    )
)

write_csv(
    significant_a2,
    file.path(
        OUTPUT_DIR,
        "F10J_significant_UpA2_TF_targets.csv"
    )
)

write_csv(
    significant_a1,
    file.path(
        OUTPUT_DIR,
        "F10J_significant_UpA1_TF_targets.csv"
    )
)

write_csv(
    priority_table,
    file.path(
        OUTPUT_DIR,
        "F10J_UpA2_regulatory_priority.csv"
    )
)

write_csv(
    gtrd_universe,
    file.path(
        OUTPUT_DIR,
        "F10J_GTRD_membership_in_adult_up_universe.csv"
    )
)


# ============================================================
# SUMMARY
# ============================================================

db_versions <- safe_collapse(
    gtrd$DB_version
)

summary_lines <- c(
    "F10J — REGULATORY DRIVER SCREEN",
    "================================",
    "",
    paste0(
        "Adult-Up universe: ",
        nrow(
            adult_up
        )
    ),
    paste0(
        "Up_A1 genes: ",
        n_a1
    ),
    paste0(
        "Up_A2 genes: ",
        n_a2
    ),
    paste0(
        "MSigDB version(s): ",
        db_versions
    ),
    paste0(
        "Eligible GTRD TF-target sets tested: ",
        nrow(
            test_results
        )
    ),
    paste0(
        "Significant Up_A2-enriched TF target sets (FDR < ",
        FDR_THRESHOLD,
        "): ",
        nrow(
            significant_a2
        )
    ),
    paste0(
        "Significant Up_A1-enriched TF target sets (FDR < ",
        FDR_THRESHOLD,
        "): ",
        nrow(
            significant_a1
        )
    ),
    "",
    "Interpretation guardrail:",
    paste(
        "TF-target enrichment is binding-target enrichment.",
        "It is not evidence by itself that the TF is activated,",
        "that it causally drives Up_A2, or that GTRD edges are",
        "functional in adult neural stem cells."
    )
)

writeLines(
    summary_lines,
    con = file.path(
        OUTPUT_DIR,
        "F10J_summary.txt"
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
    "F10J — REGULATORY DRIVER SCREEN"
)

message(
    "============================================================"
)

message("")

message(
    "MSigDB version(s): ",
    db_versions
)

message(
    "Eligible GTRD TF-target sets tested: ",
    nrow(
        test_results
    )
)

message(
    "Significant Up_A2-enriched sets: ",
    nrow(
        significant_a2
    )
)

message(
    "Significant Up_A1-enriched sets: ",
    nrow(
        significant_a1
    )
)

message("")


# ------------------------------------------------------------
# TOP Up_A2 TARGET PROGRAMS
# ------------------------------------------------------------

message(
    "TOP Up_A2-ENRICHED TF TARGET PROGRAMS"
)

message(
    "-------------------------------------"
)

if (nrow(
    significant_a2
) == 0) {

    message(
        "No GTRD TF-target sets significant at FDR < ",
        FDR_THRESHOLD,
        "."
    )

} else {

    print(
        as.data.frame(
            significant_a2 %>%
                select(
                    TF_symbol,
                    Targets_in_adult_up,
                    Up_A2_targets,
                    Up_A1_targets,
                    Odds_ratio,
                    OR_CI95_low,
                    OR_CI95_high,
                    P_value,
                    FDR,
                    TF_expression_status,
                    TF_RbTKO_LFC,
                    TF_E2F13_LFC,
                    TF_CAG_LFC,
                    TF_MAP2_LFC
                ) %>%
                slice_head(
                    n = 30
                )
        ),
        row.names = FALSE
    )
}


message("")


# ------------------------------------------------------------
# TOP Up_A1 TARGET PROGRAMS
# ------------------------------------------------------------

message(
    "TOP Up_A1-ENRICHED TF TARGET PROGRAMS"
)

message(
    "-------------------------------------"
)

if (nrow(
    significant_a1
) == 0) {

    message(
        "No GTRD TF-target sets significant at FDR < ",
        FDR_THRESHOLD,
        "."
    )

} else {

    print(
        as.data.frame(
            significant_a1 %>%
                select(
                    TF_symbol,
                    Targets_in_adult_up,
                    Up_A2_targets,
                    Up_A1_targets,
                    Odds_ratio,
                    OR_CI95_low,
                    OR_CI95_high,
                    P_value,
                    FDR,
                    TF_expression_status,
                    TF_RbTKO_LFC,
                    TF_E2F13_LFC,
                    TF_CAG_LFC,
                    TF_MAP2_LFC
                ) %>%
                slice_head(
                    n = 20
                )
        ),
        row.names = FALSE
    )
}


message("")
message(
    "Outputs written to:"
)

message(
    OUTPUT_DIR
)

message("")

message(
    "IMPORTANT:"
)

message(
    paste(
        "Do not call the top TF a driver from this analysis alone.",
        "F10J identifies regulatory hypotheses whose target sets",
        "distinguish Up_A2 from Up_A1."
    )
)