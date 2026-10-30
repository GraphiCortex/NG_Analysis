#!/usr/bin/env Rscript

# ============================================================
# F10I_upA2_functional_decomposition.R
# ============================================================
#
# PURPOSE
# -------
# Decompose the data-derived adult-Up archetype Up_A2 into broad,
# pre-specified functional annotation axes:
#
#   - immune / innate immune / inflammatory
#   - programmed cell death / apoptosis / extrinsic apoptosis
#   - neurogenesis / neuronal differentiation / glial differentiation
#
# PRIMARY COMPARISON
# ------------------
# Up_A2 versus Up_A1, restricted to genes already classified as
# adult Rb-TKO Up in F10G.
#
# CENTRAL TEST
# ------------
# Does death annotation remain enriched in Up_A2 after removing
# all immune-annotated genes?
#
# Conversely, does immune annotation remain enriched in Up_A2
# after removing all death-annotated genes?
#
# IMPORTANT
# ---------
# This is an annotation-level decomposition. It does NOT establish
# causal regulatory relationships or prove that annotated death
# pathways are functionally active.
# ============================================================


suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Mm.eg.db)
    library(dplyr)
    library(readr)
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
    "F10I_UpA2FunctionalDecomposition"
)

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# PRE-SPECIFIED GO ANCHORS
# ============================================================
#
# GOALL contains inherited/ancestor annotations. Using broad,
# pre-specified anchors lets us ask whether Up_A2 is enriched for
# these functional axes without selecting individual GO terms
# after seeing the F10H enrichment output.
# ============================================================

GO_ANCHORS <- tibble(
    Functional_axis = c(
        "IMMUNE",
        "INNATE_IMMUNE",
        "INFLAMMATORY",
        "PROGRAMMED_DEATH",
        "APOPTOSIS",
        "EXTRINSIC_APOPTOSIS",
        "NEUROGENESIS",
        "NEURON_DIFFERENTIATION",
        "GLIAL_DIFFERENTIATION"
    ),
    GO_ID = c(
        "GO:0006955",  # immune response
        "GO:0045087",  # innate immune response
        "GO:0006954",  # inflammatory response
        "GO:0012501",  # programmed cell death
        "GO:0006915",  # apoptotic process
        "GO:0097191",  # extrinsic apoptotic signaling pathway
        "GO:0022008",  # neurogenesis
        "GO:0030182",  # neuron differentiation
        "GO:0010001"   # glial cell differentiation
    )
)


# ============================================================
# PRE-SPECIFIED CANDIDATES
# ============================================================

CANDIDATES <- c(
    "MYC",
    "DDIT3",
    "MAP3K5",
    "FAS",
    "TNFRSF10B",
    "PMAIP1",
    "BBC3",
    "ZBP1",
    "PYCARD",
    "MLKL",
    "EIF2AK2",
    "NR4A1",
    "NR4A2",
    "NR4A3",
    "CDKN1A",
    "FGF2",
    "KDR",
    "CXCL12",
    "BDNF"
)


# ============================================================
# HELPERS
# ============================================================

stop_if_missing_columns <- function(data, required, label) {
    missing <- setdiff(required, colnames(data))

    if (length(missing) > 0) {
        stop(
            paste0(
                label,
                " is missing required columns: ",
                paste(missing, collapse = ", ")
            )
        )
    }
}


safe_fraction <- function(numerator, denominator) {
    if (denominator == 0) {
        return(NA_real_)
    }

    numerator / denominator
}


fisher_for_flag <- function(data, flag_column, label) {
    if (!flag_column %in% colnames(data)) {
        stop(
            paste0(
                "Unknown functional flag: ",
                flag_column
            )
        )
    }

    if (!all(data$Archetype_key %in% c("Up_A1", "Up_A2"))) {
        stop(
            "fisher_for_flag() received rows outside Up_A1/Up_A2."
        )
    }

    flag <- as.logical(
        data[[flag_column]]
    )

    flag[is.na(flag)] <- FALSE

    is_a2 <- data$Archetype_key == "Up_A2"
    is_a1 <- data$Archetype_key == "Up_A1"

    a2_yes <- sum(
        is_a2 & flag
    )

    a2_no <- sum(
        is_a2 & !flag
    )

    a1_yes <- sum(
        is_a1 & flag
    )

    a1_no <- sum(
        is_a1 & !flag
    )

    table_2x2 <- matrix(
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
            Annotation = c(
                "Yes",
                "No"
            )
        )
    )

    if (any(rowSums(table_2x2) == 0)) {
        stop(
            paste0(
                "Cannot run Fisher test for ",
                label,
                ": one archetype has zero genes after filtering."
            )
        )
    }

    test <- fisher.test(
        table_2x2
    )

    odds_ratio <- if (length(test$estimate) == 1) {
        unname(
            test$estimate
        )
    } else {
        NA_real_
    }

    ci_lower <- if (length(test$conf.int) >= 2) {
        unname(
            test$conf.int[1]
        )
    } else {
        NA_real_
    }

    ci_upper <- if (length(test$conf.int) >= 2) {
        unname(
            test$conf.int[2]
        )
    } else {
        NA_real_
    }

    tibble(
        Functional_axis = label,
        Up_A2_positive = a2_yes,
        Up_A2_total = a2_yes + a2_no,
        Up_A1_positive = a1_yes,
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
        OR_CI95_low = ci_lower,
        OR_CI95_high = ci_upper,
        P_value = test$p.value
    )
}


# ============================================================
# LOAD F10G ASSIGNMENTS
# ============================================================

message(
    "Loading F10G assignments..."
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
    "Adult_residual_z"
)

stop_if_missing_columns(
    df,
    required_columns,
    "F10G cluster assignment table"
)


# ============================================================
# RESTRICT TO ADULT-UP GENES
# ============================================================

adult_up <- df %>%
    filter(
        RbTKO_Direction == "Up"
    )

message(
    "Adult-Up genes: ",
    nrow(adult_up)
)

message(
    "Up_A1 genes: ",
    sum(
        adult_up$Archetype_key == "Up_A1"
    )
)

message(
    "Up_A2 genes: ",
    sum(
        adult_up$Archetype_key == "Up_A2"
    )
)

observed_up_archetypes <- sort(
    unique(
        adult_up$Archetype_key
    )
)

expected_up_archetypes <- c(
    "Up_A1",
    "Up_A2"
)

if (!setequal(
    observed_up_archetypes,
    expected_up_archetypes
)) {
    stop(
        paste0(
            "Expected adult-Up archetypes {Up_A1, Up_A2}, but found: ",
            paste(
                observed_up_archetypes,
                collapse = ", "
            ),
            ". F10I assumes the F10G adult-Up k=2 solution."
        )
    )
}


# ============================================================
# CASE-INSENSITIVE OFFICIAL MOUSE SYMBOL MAPPING
# ============================================================

message(
    "Building official mouse-symbol map..."
)

all_symbols <- AnnotationDbi::keys(
    org.Mm.eg.db,
    keytype = "SYMBOL"
)

symbol_map <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = all_symbols,
    keytype = "SYMBOL",
    columns = c(
        "SYMBOL",
        "ENTREZID"
    )
)

symbol_map <- symbol_map %>%
    filter(
        !is.na(SYMBOL),
        !is.na(ENTREZID)
    ) %>%
    mutate(
        Gene_match_key = toupper(
            trimws(
                SYMBOL
            )
        )
    ) %>%
    distinct(
        Gene_match_key,
        SYMBOL,
        ENTREZID
    )

adult_up <- adult_up %>%
    mutate(
        Gene_match_key = toupper(
            trimws(
                Gene
            )
        )
    )

requested_keys <- unique(
    adult_up$Gene_match_key
)

collision_table <- symbol_map %>%
    filter(
        Gene_match_key %in% requested_keys
    ) %>%
    group_by(
        Gene_match_key
    ) %>%
    summarise(
        N_symbols = n_distinct(
            SYMBOL
        ),
        N_entrez = n_distinct(
            ENTREZID
        ),
        Symbols = paste(
            sort(
                unique(
                    SYMBOL
                )
            ),
            collapse = ";"
        ),
        Entrez_IDs = paste(
            sort(
                unique(
                    ENTREZID
                )
            ),
            collapse = ";"
        ),
        .groups = "drop"
    ) %>%
    filter(
        N_symbols > 1 |
        N_entrez > 1
    )

write_csv(
    collision_table,
    file.path(
        OUTPUT_DIR,
        "F10I_symbol_case_collisions.csv"
    )
)

message(
    "Ambiguous uppercase keys among adult-Up genes: ",
    nrow(
        collision_table
    )
)

ambiguous_keys <- collision_table$Gene_match_key

symbol_map_clean <- symbol_map %>%
    filter(
        !Gene_match_key %in% ambiguous_keys
    ) %>%
    distinct(
        Gene_match_key,
        .keep_all = TRUE
    )

adult_up <- adult_up %>%
    left_join(
        symbol_map_clean %>%
            select(
                Gene_match_key,
                Official_symbol = SYMBOL,
                ENTREZID
            ),
        by = "Gene_match_key"
    )

mapped_n <- sum(
    !is.na(
        adult_up$Official_symbol
    )
)

mapping_fraction <- mapped_n / nrow(
    adult_up
)

message(
    "Mapped adult-Up genes: ",
    mapped_n,
    "/",
    nrow(adult_up),
    " (",
    round(
        100 * mapping_fraction,
        1
    ),
    "%)"
)

mapping_audit <- adult_up %>%
    transmute(
        Gene,
        Gene_match_key,
        Official_symbol,
        ENTREZID,
        Archetype_key,
        Mapping_status = case_when(
            Gene_match_key %in% ambiguous_keys ~
                "Ambiguous_uppercase_key",
            !is.na(Official_symbol) ~
                "Mapped",
            TRUE ~
                "Unmapped"
        )
    )

write_csv(
    mapping_audit,
    file.path(
        OUTPUT_DIR,
        "F10I_gene_mapping_audit.csv"
    )
)

if (mapping_fraction < 0.80) {
    stop(
        paste0(
            "Only ",
            round(
                100 * mapping_fraction,
                1
            ),
            "% of adult-Up genes mapped. ",
            "Do not interpret F10I until mapping is repaired."
        )
    )
}


# ============================================================
# RETRIEVE GOALL BIOLOGICAL-PROCESS ANNOTATIONS
# ============================================================

official_symbols <- unique(
    na.omit(
        adult_up$Official_symbol
    )
)

message(
    "Retrieving GOALL annotations for ",
    length(
        official_symbols
    ),
    " official symbols..."
)

go_annotations <- AnnotationDbi::select(
    org.Mm.eg.db,
    keys = official_symbols,
    keytype = "SYMBOL",
    columns = c(
        "GOALL",
        "ONTOLOGYALL"
    )
)

go_annotations <- go_annotations %>%
    filter(
        !is.na(GOALL),
        ONTOLOGYALL == "BP"
    ) %>%
    distinct(
        SYMBOL,
        GOALL
    )

message(
    "Unique SYMBOL-GOALL BP annotations retained: ",
    nrow(
        go_annotations
    )
)


# ============================================================
# BUILD GENE-LEVEL FUNCTIONAL MEMBERSHIP
# ============================================================

membership <- adult_up %>%
    select(
        Gene,
        Official_symbol,
        ENTREZID,
        Archetype_key,
        RbTKO_LFC,
        E2F13_LFC,
        CAG_LFC,
        MAP2_LFC,
        Adult_residual_z
    )

for (i in seq_len(
    nrow(
        GO_ANCHORS
    )
)) {
    axis_name <- GO_ANCHORS$Functional_axis[[i]]

    go_id <- GO_ANCHORS$GO_ID[[i]]

    positive_symbols <- unique(
        go_annotations$SYMBOL[
            go_annotations$GOALL == go_id
        ]
    )

    membership[[axis_name]] <- (
        !is.na(
            membership$Official_symbol
        ) &
        membership$Official_symbol %in% positive_symbols
    )
}

membership <- membership %>%
    mutate(
        Immune_any = (
            IMMUNE |
            INNATE_IMMUNE |
            INFLAMMATORY
        ),
        Death_any = (
            PROGRAMMED_DEATH |
            APOPTOSIS |
            EXTRINSIC_APOPTOSIS
        ),
        Neural_any = (
            NEUROGENESIS |
            NEURON_DIFFERENTIATION |
            GLIAL_DIFFERENTIATION
        ),
        Immune_death_compartment = case_when(
            Immune_any & Death_any ~
                "Immune_and_death",
            Immune_any & !Death_any ~
                "Immune_only",
            !Immune_any & Death_any ~
                "Death_only",
            TRUE ~
                "Neither"
        )
    )

write_csv(
    membership,
    file.path(
        OUTPUT_DIR,
        "F10I_gene_functional_membership.csv"
    )
)


# ============================================================
# PRIMARY ANCHOR ENRICHMENT TESTS
# ============================================================

anchor_tests <- bind_rows(
    lapply(
        GO_ANCHORS$Functional_axis,
        function(axis_name) {
            fisher_for_flag(
                membership,
                axis_name,
                axis_name
            )
        }
    ),

    fisher_for_flag(
        membership,
        "Immune_any",
        "IMMUNE_ANY"
    ),

    fisher_for_flag(
        membership,
        "Death_any",
        "DEATH_ANY"
    ),

    fisher_for_flag(
        membership,
        "Neural_any",
        "NEURAL_ANY"
    )
) %>%
    mutate(
        FDR = p.adjust(
            P_value,
            method = "BH"
        )
    ) %>%
    arrange(
        FDR,
        P_value
    )

write_csv(
    anchor_tests,
    file.path(
        OUTPUT_DIR,
        "F10I_anchor_enrichment_tests.csv"
    )
)


# ============================================================
# IMMUNE / DEATH OVERLAP WITHIN Up_A2
# ============================================================

up_a2 <- membership %>%
    filter(
        Archetype_key == "Up_A2"
    )

compartment_counts <- up_a2 %>%
    count(
        Immune_death_compartment,
        name = "N_genes"
    ) %>%
    mutate(
        Fraction_Up_A2 = N_genes / sum(
            N_genes
        )
    ) %>%
    arrange(
        desc(
            N_genes
        )
    )

write_csv(
    compartment_counts,
    file.path(
        OUTPUT_DIR,
        "F10I_UpA2_immune_death_overlap.csv"
    )
)

immune_genes <- up_a2 %>%
    filter(
        Immune_any
    ) %>%
    pull(
        Gene
    ) %>%
    unique()

death_genes <- up_a2 %>%
    filter(
        Death_any
    ) %>%
    pull(
        Gene
    ) %>%
    unique()

overlap_genes <- intersect(
    immune_genes,
    death_genes
)

union_genes <- union(
    immune_genes,
    death_genes
)

death_only_genes <- setdiff(
    death_genes,
    immune_genes
)

immune_only_genes <- setdiff(
    immune_genes,
    death_genes
)

jaccard <- if (length(
    union_genes
) > 0) {
    length(
        overlap_genes
    ) / length(
        union_genes
    )
} else {
    NA_real_
}

overlap_summary <- tibble(
    Up_A2_total = nrow(
        up_a2
    ),
    Immune_genes = length(
        immune_genes
    ),
    Death_genes = length(
        death_genes
    ),
    Immune_and_death = length(
        overlap_genes
    ),
    Death_only = length(
        death_only_genes
    ),
    Immune_only = length(
        immune_only_genes
    ),
    Immune_or_death = length(
        union_genes
    ),
    Jaccard = jaccard
)

write_csv(
    overlap_summary,
    file.path(
        OUTPUT_DIR,
        "F10I_overlap_summary.csv"
    )
)

write_csv(
    tibble(
        Gene = sort(
            immune_genes
        )
    ),
    file.path(
        OUTPUT_DIR,
        "F10I_UpA2_immune_genes.csv"
    )
)

write_csv(
    tibble(
        Gene = sort(
            death_genes
        )
    ),
    file.path(
        OUTPUT_DIR,
        "F10I_UpA2_death_genes.csv"
    )
)

write_csv(
    tibble(
        Gene = sort(
            overlap_genes
        )
    ),
    file.path(
        OUTPUT_DIR,
        "F10I_UpA2_immune_and_death_genes.csv"
    )
)

write_csv(
    tibble(
        Gene = sort(
            death_only_genes
        )
    ),
    file.path(
        OUTPUT_DIR,
        "F10I_UpA2_death_only_genes.csv"
    )
)

write_csv(
    tibble(
        Gene = sort(
            immune_only_genes
        )
    ),
    file.path(
        OUTPUT_DIR,
        "F10I_UpA2_immune_only_genes.csv"
    )
)


# ============================================================
# CRITICAL CONDITIONAL / INDEPENDENCE TESTS
# ============================================================
#
# Test 1:
# Remove every immune-annotated gene, then ask whether the
# remaining death-annotated genes are still enriched in Up_A2.
#
# Test 2:
# Remove every death-annotated gene, then ask whether the
# remaining immune-annotated genes are still enriched in Up_A2.
#
# These tests establish annotation-level separability, not causal
# independence between biological pathways.
# ============================================================

nonimmune <- membership %>%
    filter(
        !Immune_any
    )

nondeath <- membership %>%
    filter(
        !Death_any
    )

message(
    "Nonimmune adult-Up genes retained: ",
    nrow(
        nonimmune
    )
)

message(
    "Nondeath adult-Up genes retained: ",
    nrow(
        nondeath
    )
)

death_after_removing_immune <- fisher_for_flag(
    nonimmune,
    "Death_any",
    "DEATH_AFTER_REMOVING_IMMUNE"
)

immune_after_removing_death <- fisher_for_flag(
    nondeath,
    "Immune_any",
    "IMMUNE_AFTER_REMOVING_DEATH"
)

independence_tests <- bind_rows(
    death_after_removing_immune,
    immune_after_removing_death
) %>%
    mutate(
        FDR = p.adjust(
            P_value,
            method = "BH"
        )
    ) %>%
    arrange(
        FDR,
        P_value
    )

write_csv(
    independence_tests,
    file.path(
        OUTPUT_DIR,
        "F10I_independence_tests.csv"
    )
)


# ============================================================
# RESPONSE SUMMARY BY Up_A2 FUNCTIONAL COMPARTMENT
# ============================================================

response_summary <- up_a2 %>%
    group_by(
        Immune_death_compartment
    ) %>%
    summarise(
        N = n(),
        Median_RbTKO_LFC = median(
            RbTKO_LFC,
            na.rm = TRUE
        ),
        Median_E2F13_LFC = median(
            E2F13_LFC,
            na.rm = TRUE
        ),
        Median_CAG_LFC = median(
            CAG_LFC,
            na.rm = TRUE
        ),
        Median_MAP2_LFC = median(
            MAP2_LFC,
            na.rm = TRUE
        ),
        Median_Adult_residual_z = median(
            Adult_residual_z,
            na.rm = TRUE
        ),
        .groups = "drop"
    ) %>%
    arrange(
        desc(
            N
        )
    )

write_csv(
    response_summary,
    file.path(
        OUTPUT_DIR,
        "F10I_compartment_response_summary.csv"
    )
)


# ============================================================
# CANDIDATE FUNCTIONAL AUDIT
# ============================================================
#
# Candidates outside adult-Up are retained in the output as
# "Not_in_adult_up_analysis" rather than silently disappearing.
# ============================================================

candidate_audit <- tibble(
    Gene = CANDIDATES
) %>%
    left_join(
        membership %>%
            select(
                Gene,
                Archetype_key,
                Immune_any,
                Death_any,
                Neural_any,
                Immune_death_compartment,
                PROGRAMMED_DEATH,
                APOPTOSIS,
                EXTRINSIC_APOPTOSIS,
                IMMUNE,
                INNATE_IMMUNE,
                INFLAMMATORY,
                RbTKO_LFC,
                E2F13_LFC,
                CAG_LFC,
                MAP2_LFC,
                Adult_residual_z
            ),
        by = "Gene"
    ) %>%
    mutate(
        Analysis_status = if_else(
            is.na(
                Archetype_key
            ),
            "Not_in_adult_up_analysis",
            "Included"
        )
    )

write_csv(
    candidate_audit,
    file.path(
        OUTPUT_DIR,
        "F10I_candidate_functional_audit.csv"
    )
)


# ============================================================
# HUMAN-READABLE SUMMARY
# ============================================================

summary_lines <- c(
    "F10I — Up_A2 FUNCTIONAL DECOMPOSITION",
    "======================================",
    "",
    paste0(
        "Adult-Up genes: ",
        nrow(
            membership
        )
    ),
    paste0(
        "Mapped adult-Up genes: ",
        mapped_n,
        "/",
        nrow(
            adult_up
        ),
        " (",
        round(
            100 * mapping_fraction,
            1
        ),
        "%)"
    ),
    paste0(
        "Up_A2 immune genes: ",
        length(
            immune_genes
        )
    ),
    paste0(
        "Up_A2 death genes: ",
        length(
            death_genes
        )
    ),
    paste0(
        "Up_A2 immune-and-death genes: ",
        length(
            overlap_genes
        )
    ),
    paste0(
        "Up_A2 death-only genes: ",
        length(
            death_only_genes
        )
    ),
    paste0(
        "Up_A2 immune-only genes: ",
        length(
            immune_only_genes
        )
    ),
    paste0(
        "Immune/death Jaccard overlap: ",
        signif(
            jaccard,
            4
        )
    ),
    "",
    "Critical conditional tests:",
    paste0(
        "  DEATH_AFTER_REMOVING_IMMUNE: OR=",
        signif(
            death_after_removing_immune$Odds_ratio,
            4
        ),
        ", p=",
        signif(
            death_after_removing_immune$P_value,
            4
        )
    ),
    paste0(
        "  IMMUNE_AFTER_REMOVING_DEATH: OR=",
        signif(
            immune_after_removing_death$Odds_ratio,
            4
        ),
        ", p=",
        signif(
            immune_after_removing_death$P_value,
            4
        )
    ),
    "",
    "Interpretation guardrail:",
    paste(
        "These tests assess annotation-level separability between",
        "immune and death GO categories. They do not establish",
        "causal independence, pathway activation, or cell death."
    )
)

writeLines(
    summary_lines,
    con = file.path(
        OUTPUT_DIR,
        "F10I_summary.txt"
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
    "F10I — Up_A2 FUNCTIONAL DECOMPOSITION"
)
message(
    "============================================================"
)
message("")

message(
    "ANCHOR ENRICHMENT: Up_A2 versus Up_A1"
)
message(
    "--------------------------------------"
)

print(
    as.data.frame(
        anchor_tests %>%
            select(
                Functional_axis,
                Up_A2_positive,
                Up_A2_total,
                Up_A1_positive,
                Up_A1_total,
                Up_A2_fraction,
                Up_A1_fraction,
                Odds_ratio,
                OR_CI95_low,
                OR_CI95_high,
                P_value,
                FDR
            )
    ),
    row.names = FALSE
)

message("")
message(
    "Up_A2 IMMUNE / DEATH OVERLAP"
)
message(
    "----------------------------"
)

print(
    as.data.frame(
        compartment_counts
    ),
    row.names = FALSE
)

message("")

print(
    as.data.frame(
        overlap_summary
    ),
    row.names = FALSE
)

message("")
message(
    "CRITICAL CONDITIONAL TESTS"
)
message(
    "--------------------------"
)

print(
    as.data.frame(
        independence_tests
    ),
    row.names = FALSE
)

message("")
message(
    "Up_A2 RESPONSE BY FUNCTIONAL COMPARTMENT"
)
message(
    "----------------------------------------"
)

print(
    as.data.frame(
        response_summary
    ),
    row.names = FALSE
)

message("")
message(
    "CANDIDATE FUNCTIONAL AUDIT"
)
message(
    "--------------------------"
)

print(
    as.data.frame(
        candidate_audit
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
    "INTERPRETATION:"
)

message(
    paste(
        "The key quantity is DEATH_AFTER_REMOVING_IMMUNE.",
        "If its odds ratio remains > 1 with convincing statistical",
        "support, death annotation is enriched in Up_A2 even after",
        "immune-annotated genes are removed."
    )
)

message(
    paste(
        "This would support annotation-level separability, not",
        "causal independence or proof of active cell death."
    )
)