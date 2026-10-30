#!/usr/bin/env Rscript

# ============================================================
# F10L_adult_signed_TF_activity.R
# ============================================================
#
# PURPOSE
# -------
# Infer signed transcription-factor activities after adult
# Rb-family loss using the FULL GSE190765 processed RNA-seq
# dataset rather than the significant-gene subset.
#
# TWO INDEPENDENT EVIDENCE LAYERS
# -------------------------------
#
# 1. CONTRAST-LEVEL ACTIVITY
#
#    Use the DESeq2 Wald statistic ("stat") for every sufficiently
#    expressed gene as the molecular readout.
#
#    Positive TF activity:
#      activated targets tend to move up and/or repressed targets
#      tend to move down after Rb-TKO.
#
#    Negative TF activity:
#      the converse pattern.
#
#
# 2. REPLICATE-LEVEL ACTIVITY
#
#    Use log2(normalized_count + 1) for the six biological samples,
#    center each gene across samples, and infer TF activity in each
#    biological replicate.
#
#    We then compare:
#
#      RB_TKO_1/2/3
#          versus
#      RB_THC_1/2/3
#
#
# IMPORTANT
# ---------
# CollecTRI is prior regulatory knowledge.
#
# Inferred TF activity is NOT direct protein-activity measurement.
# It is an inference from concordant signed target behavior.
#
# ============================================================


# ============================================================
# PACKAGES
# ============================================================

required_packages <- c(
    "decoupleR",
    "dplyr",
    "readr",
    "tibble"
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
            "Missing package(s): ",
            paste(
                missing_packages,
                collapse = ", "
            )
        )
    )
}

suppressPackageStartupMessages({
    library(decoupleR)
    library(dplyr)
    library(readr)
    library(tibble)
})


# ============================================================
# PATHS
# ============================================================

INPUT_FILE <- file.path(
    "Data",
    "Raw",
    "RNAseq",
    "Fong2022",
    "GSE190765",
    "GSE190765_FDSC_FC.Batch.IHW.csv.gz"
)

OUTPUT_DIR <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F10L_AdultSignedTFActivity"
)

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# PARAMETERS
# ============================================================

PRIMARY_BASEMEAN <- 10

SENSITIVITY_BASEMEAN <- 1

MIN_REGULON_SIZE <- 5

SAMPLE_COLUMNS <- c(
    "RB_THC_1",
    "RB_THC_2",
    "RB_THC_3",
    "RB_TKO_1",
    "RB_TKO_2",
    "RB_TKO_3"
)

CANDIDATE_TFS <- c(
    "E2F1",
    "E2F2",
    "E2F3",
    "E2F4",

    "REL",
    "RELA",
    "RELB",
    "NFKB1",
    "NFKB2",

    "STAT1",
    "STAT2",
    "STAT3",
    "STAT5A",
    "STAT5B",

    "IRF1",
    "IRF2",
    "IRF3",
    "IRF7",
    "IRF8",
    "IRF9",

    "MYC",
    "TP53",
    "JUN",
    "FOS",
    "ATF4",
    "CEBPD",
    "NR4A1",
    "NR4A2",
    "NR4A3"
)


# ============================================================
# HELPERS
# ============================================================

safe_fdr <- function(p) {
    p.adjust(
        p,
        method = "BH"
    )
}


safe_mean <- function(x) {
    if (all(is.na(x))) {
        return(NA_real_)
    }

    mean(
        x,
        na.rm = TRUE
    )
}


safe_sd <- function(x) {
    x <- x[
        is.finite(
            x
        )
    ]

    if (length(x) < 2) {
        return(NA_real_)
    }

    sd(x)
}


safe_effect_size <- function(
    tko,
    thc
) {
    tko <- tko[
        is.finite(
            tko
        )
    ]

    thc <- thc[
        is.finite(
            thc
        )
    ]

    n1 <- length(tko)
    n0 <- length(thc)

    if (n1 < 2 || n0 < 2) {
        return(
            c(
                Cohen_d = NA_real_,
                Hedges_g = NA_real_
            )
        )
    }

    pooled_var <- (
        ((n1 - 1) * var(tko)) +
        ((n0 - 1) * var(thc))
    ) / (
        n1 + n0 - 2
    )

    if (
        !is.finite(
            pooled_var
        ) ||
        pooled_var <= 0
    ) {
        return(
            c(
                Cohen_d = NA_real_,
                Hedges_g = NA_real_
            )
        )
    }

    pooled_sd <- sqrt(
        pooled_var
    )

    d <- (
        mean(tko) -
        mean(thc)
    ) / pooled_sd

    correction <- 1 - (
        3 / (
            4 * (n1 + n0) - 9
        )
    )

    g <- correction * d

    c(
        Cohen_d = d,
        Hedges_g = g
    )
}


exact_group_permutation_p <- function(
    scores,
    groups
) {
    scores <- as.numeric(
        scores
    )

    groups <- as.character(
        groups
    )

    keep <- (
        is.finite(scores) &
        groups %in% c(
            "THC",
            "TKO"
        )
    )

    scores <- scores[keep]
    groups <- groups[keep]

    if (
        sum(groups == "THC") != 3 ||
        sum(groups == "TKO") != 3
    ) {
        return(NA_real_)
    }

    observed <- (
        mean(
            scores[
                groups == "TKO"
            ]
        ) -
        mean(
            scores[
                groups == "THC"
            ]
        )
    )

    all_tko_assignments <- combn(
        seq_along(
            scores
        ),
        3,
        simplify = FALSE
    )

    null_differences <- vapply(
        all_tko_assignments,
        function(tko_idx) {
            thc_idx <- setdiff(
                seq_along(
                    scores
                ),
                tko_idx
            )

            mean(
                scores[tko_idx]
            ) -
            mean(
                scores[thc_idx]
            )
        },
        numeric(1)
    )

    mean(
        abs(
            null_differences
        ) >= (
            abs(
                observed
            ) - 1e-12
        )
    )
}


prepare_ulm_output <- function(
    x,
    score_name,
    p_name
) {
    required <- c(
        "source",
        "condition",
        "score"
    )

    missing <- setdiff(
        required,
        colnames(x)
    )

    if (length(missing) > 0) {
        stop(
            paste0(
                "Unexpected decoupleR ULM output. Missing: ",
                paste(
                    missing,
                    collapse = ", "
                )
            )
        )
    }

    if (!"p_value" %in% colnames(x)) {
        x$p_value <- NA_real_
    }

    x %>%
        transmute(
            TF = source,
            Condition = condition,
            !!score_name := score,
            !!p_name := p_value
        )
}


# ============================================================
# INPUT
# ============================================================

message(
    "Loading complete GSE190765 processed table..."
)

raw <- read_csv(
    INPUT_FILE,
    show_col_types = FALSE
)

message(
    "Raw rows: ",
    nrow(raw)
)

message(
    "Raw columns: ",
    ncol(raw)
)


required_columns <- c(
    "Row.names",
    "name",
    "baseMean",
    "log2FoldChange",
    "stat",
    "pvalue",
    "padj",
    SAMPLE_COLUMNS
)

missing_columns <- setdiff(
    required_columns,
    colnames(raw)
)

if (length(missing_columns) > 0) {
    stop(
        paste0(
            "Input table missing required columns: ",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    )
}


# ============================================================
# CLEAN GENE SYMBOLS
# ============================================================
#
# Multiple Ensembl rows can occasionally map to the same symbol.
#
# For a duplicated symbol we retain the row with the largest
# baseMean. This rule is independent of differential direction
# and significance, avoiding selection by logFC/p-value.
# ============================================================

clean <- raw %>%
    mutate(
        Gene = trimws(
            as.character(
                name
            )
        ),

        Gene_key = toupper(
            Gene
        ),

        baseMean_numeric = as.numeric(
            baseMean
        )
    ) %>%
    filter(
        !is.na(
            Gene
        ),
        Gene != "",
        Gene != "NA",
        is.finite(
            baseMean_numeric
        )
    ) %>%
    arrange(
        Gene_key,
        desc(
            baseMean_numeric
        )
    ) %>%
    group_by(
        Gene_key
    ) %>%
    slice_head(
        n = 1
    ) %>%
    ungroup()


message(
    "Unique gene symbols retained: ",
    nrow(clean)
)

message(
    "Duplicate/invalid rows removed: ",
    nrow(raw) - nrow(clean)
)

# ============================================================
# COLLECTRI — LOCAL STATIC MOUSE RESOURCE
# ============================================================
#
# decoupleR::get_collectri() currently fails in some OmnipathR
# versions because organism lookup depends on external Ensembl
# metadata.
#
# We therefore use the mouse-specific static OmniPath resource
# directly:
#
#   interactions_collectri_10090.tsv.gz
#
# 10090 = Mus musculus.
#
# This is the same static resource used by OmniPathR's own
# fallback mechanism.
# ============================================================

COLLECTRI_FILE <- file.path(
    "Data",
    "Raw",
    "Regulatory",
    "CollecTRI",
    "interactions_collectri_10090.tsv.gz"
)


if (!file.exists(
    COLLECTRI_FILE
)) {
    stop(
        paste0(
            "Mouse CollecTRI static file not found:\n",
            COLLECTRI_FILE
        )
    )
}


message("")
message(
    "Loading local mouse CollecTRI static resource..."
)


collectri_raw <- read_tsv(
    COLLECTRI_FILE,
    show_col_types = FALSE,
    progress = FALSE
)


message(
    "Raw CollecTRI rows: ",
    nrow(
        collectri_raw
    )
)


message(
    "Raw CollecTRI columns: ",
    ncol(
        collectri_raw
    )
)


# ------------------------------------------------------------
# Validate expected OmniPath static-table schema
# ------------------------------------------------------------

required_network_columns <- c(
    "source_genesymbol",
    "target_genesymbol",
    "is_stimulation",
    "is_inhibition"
)


missing_network_columns <- setdiff(
    required_network_columns,
    colnames(
        collectri_raw
    )
)


if (length(
    missing_network_columns
) > 0) {

    stop(
        paste0(
            "Unexpected static CollecTRI schema. Missing: ",
            paste(
                missing_network_columns,
                collapse = ", "
            ),
            "\nAvailable columns:\n",
            paste(
                colnames(
                    collectri_raw
                ),
                collapse = ", "
            )
        )
    )
}


# ============================================================
# REPRODUCE decoupleR::get_collectri() FORMATTING
# ============================================================
#
# decoupleR:
#
#   - keeps TF complexes intact when split_complexes = FALSE
#   - collapses JUN/FOS complexes to AP1
#   - collapses REL/NFKB complexes to NFKB
#   - converts stimulation/inhibition into signed mor values
#
# We reproduce that logic here while keeping only the fields
# needed for ULM.
# ============================================================


collectri_work <- collectri_raw %>%
    mutate(
        source_genesymbol = as.character(
            source_genesymbol
        ),

        target_genesymbol = as.character(
            target_genesymbol
        ),

        source_id = if (
            "source" %in%
                colnames(
                    collectri_raw
                )
        ) {
            as.character(
                source
            )
        } else {
            source_genesymbol
        }
    )


# ------------------------------------------------------------
# Non-complex interactions
# ------------------------------------------------------------

collectri_interactions <- collectri_work %>%
    filter(
        !grepl(
            "COMPLEX",
            source_id,
            fixed = TRUE
        )
    )


# ------------------------------------------------------------
# Complex TF interactions
# ------------------------------------------------------------

collectri_complex <- collectri_work %>%
    filter(
        grepl(
            "COMPLEX",
            source_id,
            fixed = TRUE
        )
    ) %>%
    mutate(
        source_genesymbol = case_when(

            grepl(
                "JUN",
                source_genesymbol,
                ignore.case = TRUE
            ) |
            grepl(
                "FOS",
                source_genesymbol,
                ignore.case = TRUE
            ) ~
                "AP1",

            grepl(
                "REL",
                source_genesymbol,
                ignore.case = TRUE
            ) |
            grepl(
                "NFKB",
                source_genesymbol,
                ignore.case = TRUE
            ) ~
                "NFKB",

            TRUE ~
                source_genesymbol
        )
    )


collectri_formatted <- bind_rows(
    collectri_interactions,
    collectri_complex
)


# ------------------------------------------------------------
# Convert to decoupleR network schema
# ------------------------------------------------------------

network <- collectri_formatted %>%
    transmute(
        source = toupper(
            trimws(
                source_genesymbol
            )
        ),

        target = toupper(
            trimws(
                target_genesymbol
            )
        ),

        is_stimulation = as.numeric(
            is_stimulation
        ),

        is_inhibition = as.numeric(
            is_inhibition
        )
    ) %>%
    filter(
        !is.na(
            source
        ),
        source != "",
        !is.na(
            target
        ),
        target != "",
        is_stimulation %in% c(
            0,
            1
        )
    ) %>%
    mutate(
        mor = if_else(
            is_stimulation == 1,
            1,
            -1
        )
    ) %>%
    select(
        source,
        target,
        mor
    ) %>%
    distinct(
        source,
        target,
        .keep_all = TRUE
    )


message(
    "Formatted CollecTRI edges retained: ",
    nrow(
        network
    )
)


message(
    "Formatted CollecTRI regulators: ",
    n_distinct(
        network$source
    )
)


if (nrow(
    network
) == 0) {
    stop(
        "Formatted CollecTRI network contains zero edges."
    )
}
primary <- clean %>%
    filter(
        baseMean_numeric >= PRIMARY_BASEMEAN,
        is.finite(
            stat
        )
    )


sensitivity <- clean %>%
    filter(
        baseMean_numeric >= SENSITIVITY_BASEMEAN,
        is.finite(
            stat
        )
    )


message("")
message(
    "Primary genes (baseMean >= ",
    PRIMARY_BASEMEAN,
    ", finite DESeq2 stat): ",
    nrow(primary)
)

message(
    "Sensitivity genes (baseMean >= ",
    SENSITIVITY_BASEMEAN,
    ", finite DESeq2 stat): ",
    nrow(sensitivity)
)


# ============================================================
# NETWORK COVERAGE
# ============================================================

primary_gene_keys <- unique(
    primary$Gene_key
)


network_coverage <- network %>%
    mutate(
        Target_in_primary =
            target %in%
            primary_gene_keys
    ) %>%
    group_by(
        source
    ) %>%
    summarise(
        Total_edges =
            n(),

        Targets_in_primary =
            sum(
                Target_in_primary
            ),

        .groups =
            "drop"
    ) %>%
    arrange(
        desc(
            Targets_in_primary
        )
    )


write_csv(
    network_coverage,
    file.path(
        OUTPUT_DIR,
        "F10L_CollecTRI_network_coverage.csv"
    )
)


eligible_regulators <- network_coverage %>%
    filter(
        Targets_in_primary >= MIN_REGULON_SIZE
    ) %>%
    pull(
        source
    )


message(
    "Regulators with >= ",
    MIN_REGULON_SIZE,
    " measured targets: ",
    length(
        eligible_regulators
    )
)


network_primary <- network %>%
    filter(
        source %in% eligible_regulators
    )


# ============================================================
# PRIMARY CONTRAST-LEVEL ACTIVITY
#
# Input = full DESeq2 Wald statistic vector.
# ============================================================

message("")
message(
    "Inferring primary contrast-level TF activity..."
)


contrast_mat <- matrix(
    as.numeric(
        primary$stat
    ),
    ncol = 1,
    dimnames = list(
        primary$Gene_key,
        "RB_TKO_vs_RB_THC_DESeq2_stat"
    )
)


contrast_ulm_raw <- decoupleR::run_ulm(
    mat = contrast_mat,
    network = network_primary,
    .source = "source",
    .target = "target",
    .mor = "mor",
    center = FALSE,
    minsize = MIN_REGULON_SIZE
)


contrast_activity <- prepare_ulm_output(
    contrast_ulm_raw,
    "Contrast_activity",
    "Contrast_activity_p"
) %>%
    mutate(
        Contrast_activity_FDR =
            safe_fdr(
                Contrast_activity_p
            )
    ) %>%
    arrange(
        Contrast_activity_FDR,
        desc(
            abs(
                Contrast_activity
            )
        )
    )


# ============================================================
# SENSITIVITY CONTRAST
#
# Lower detectability threshold: baseMean >= 1.
# ============================================================

message(
    "Running baseMean >= 1 sensitivity contrast..."
)


sensitivity_mat <- matrix(
    as.numeric(
        sensitivity$stat
    ),
    ncol = 1,
    dimnames = list(
        sensitivity$Gene_key,
        "RB_TKO_vs_RB_THC_DESeq2_stat"
    )
)


sensitivity_network <- network %>%
    filter(
        target %in%
            sensitivity$Gene_key
    )


sensitivity_ulm_raw <- decoupleR::run_ulm(
    mat = sensitivity_mat,
    network = sensitivity_network,
    .source = "source",
    .target = "target",
    .mor = "mor",
    center = FALSE,
    minsize = MIN_REGULON_SIZE
)


sensitivity_activity <- prepare_ulm_output(
    sensitivity_ulm_raw,
    "Sensitivity_activity",
    "Sensitivity_activity_p"
) %>%
    mutate(
        Sensitivity_activity_FDR =
            safe_fdr(
                Sensitivity_activity_p
            )
    )


# ============================================================
# REPLICATE-LEVEL TF ACTIVITY
# ============================================================
#
# Use log2(normalized count + 1).
#
# ULM center=TRUE subtracts each gene's mean across the six
# samples, so activity scores represent within-experiment
# differences rather than absolute baseline expression.
# ============================================================

message("")
message(
    "Inferring replicate-level TF activity..."
)


sample_counts <- as.matrix(
    primary[
        SAMPLE_COLUMNS
    ]
)

storage.mode(
    sample_counts
) <- "double"


if (any(
    !is.finite(
        sample_counts
    )
)) {
    stop(
        "Non-finite normalized count detected in sample matrix."
    )
}


if (any(
    sample_counts < 0
)) {
    stop(
        "Negative values detected in normalized counts."
    )
}


sample_log <- log2(
    sample_counts + 1
)

rownames(
    sample_log
) <- primary$Gene_key


sample_ulm_raw <- decoupleR::run_ulm(
    mat = sample_log,
    network = network_primary,
    .source = "source",
    .target = "target",
    .mor = "mor",
    center = TRUE,
    minsize = MIN_REGULON_SIZE
)


sample_activity <- prepare_ulm_output(
    sample_ulm_raw,
    "Sample_activity",
    "Sample_activity_p"
) %>%
    mutate(
        Group = case_when(
            grepl(
                "^RB_TKO_",
                Condition
            ) ~
                "TKO",

            grepl(
                "^RB_THC_",
                Condition
            ) ~
                "THC",

            TRUE ~
                NA_character_
        )
    )


if (any(
    is.na(
        sample_activity$Group
    )
)) {
    stop(
        "Could not assign all activity samples to THC/TKO."
    )
}


# ============================================================
# REPLICATE GROUP COMPARISON
# ============================================================

message(
    "Comparing TKO versus control activity across replicates..."
)


regulators <- sort(
    unique(
        sample_activity$TF
    )
)


sample_group_results <- lapply(
    regulators,
    function(tf_name) {

        x <- sample_activity %>%
            filter(
                TF == tf_name
            )

        tko <- x %>%
            filter(
                Group == "TKO"
            ) %>%
            pull(
                Sample_activity
            )

        thc <- x %>%
            filter(
                Group == "THC"
            ) %>%
            pull(
                Sample_activity
            )

        effect <- safe_effect_size(
            tko,
            thc
        )

        tibble(
            TF =
                tf_name,

            Mean_TKO_activity =
                safe_mean(
                    tko
                ),

            Mean_THC_activity =
                safe_mean(
                    thc
                ),

            Delta_activity =
                safe_mean(
                    tko
                ) -
                safe_mean(
                    thc
                ),

            SD_TKO =
                safe_sd(
                    tko
                ),

            SD_THC =
                safe_sd(
                    thc
                ),

            Cohen_d =
                unname(
                    effect[
                        "Cohen_d"
                    ]
                ),

            Hedges_g =
                unname(
                    effect[
                        "Hedges_g"
                    ]
                ),

            Exact_permutation_p =
                exact_group_permutation_p(
                    x$Sample_activity,
                    x$Group
                ),

            Min_TKO_activity =
                min(
                    tko,
                    na.rm = TRUE
                ),

            Max_TKO_activity =
                max(
                    tko,
                    na.rm = TRUE
                ),

            Min_THC_activity =
                min(
                    thc,
                    na.rm = TRUE
                ),

            Max_THC_activity =
                max(
                    thc,
                    na.rm = TRUE
                )
        )
    }
) %>%
    bind_rows()


# ============================================================
# INTEGRATE CONTRAST + REPLICATE + SENSITIVITY EVIDENCE
# ============================================================

integrated <- contrast_activity %>%
    select(
        TF,
        Contrast_activity,
        Contrast_activity_p,
        Contrast_activity_FDR
    ) %>%
    left_join(
        sample_group_results,
        by = "TF"
    ) %>%
    left_join(
        sensitivity_activity %>%
            select(
                TF,
                Sensitivity_activity,
                Sensitivity_activity_p,
                Sensitivity_activity_FDR
            ),
        by = "TF"
    ) %>%
    left_join(
        network_coverage %>%
            rename(
                TF = source
            ),
        by = "TF"
    ) %>%
    mutate(
        Contrast_replicate_sign_concordant =
            case_when(
                is.na(
                    Contrast_activity
                ) |
                is.na(
                    Delta_activity
                ) ~
                    NA,

                sign(
                    Contrast_activity
                ) ==
                sign(
                    Delta_activity
                ) ~
                    TRUE,

                TRUE ~
                    FALSE
            ),

        Primary_sensitivity_sign_concordant =
            case_when(
                is.na(
                    Contrast_activity
                ) |
                is.na(
                    Sensitivity_activity
                ) ~
                    NA,

                sign(
                    Contrast_activity
                ) ==
                sign(
                    Sensitivity_activity
                ) ~
                    TRUE,

                TRUE ~
                    FALSE
            ),

        Activity_direction =
            case_when(
                Contrast_activity > 0 ~
                    "Activated_in_RbTKO",

                Contrast_activity < 0 ~
                    "Reduced_in_RbTKO",

                TRUE ~
                    "Neutral"
            )
    ) %>%
    arrange(
        Contrast_activity_FDR,
        desc(
            abs(
                Contrast_activity
            )
        )
    )


# ============================================================
# CANDIDATE AUDIT
# ============================================================

candidate_audit <- integrated %>%
    mutate(
        TF_key = toupper(
            TF
        )
    ) %>%
    filter(
        TF_key %in%
            CANDIDATE_TFS
    ) %>%
    arrange(
        Contrast_activity_FDR,
        desc(
            abs(
                Contrast_activity
            )
        )
    )


# ============================================================
# TOP POSITIVE / NEGATIVE ACTIVITIES
# ============================================================

top_positive <- integrated %>%
    filter(
        Contrast_activity > 0
    ) %>%
    arrange(
        Contrast_activity_FDR,
        desc(
            Contrast_activity
        )
    ) %>%
    slice_head(
        n = 30
    )


top_negative <- integrated %>%
    filter(
        Contrast_activity < 0
    ) %>%
    arrange(
        Contrast_activity_FDR,
        Contrast_activity
    ) %>%
    slice_head(
        n = 30
    )


# ============================================================
# SAVE
# ============================================================

write_csv(
    contrast_activity,
    file.path(
        OUTPUT_DIR,
        "F10L_contrast_TF_activity_primary.csv"
    )
)


write_csv(
    sensitivity_activity,
    file.path(
        OUTPUT_DIR,
        "F10L_contrast_TF_activity_sensitivity_baseMean1.csv"
    )
)


write_csv(
    sample_activity,
    file.path(
        OUTPUT_DIR,
        "F10L_sample_TF_activity.csv"
    )
)


write_csv(
    sample_group_results,
    file.path(
        OUTPUT_DIR,
        "F10L_sample_group_activity_comparison.csv"
    )
)


write_csv(
    integrated,
    file.path(
        OUTPUT_DIR,
        "F10L_integrated_TF_activity.csv"
    )
)


write_csv(
    candidate_audit,
    file.path(
        OUTPUT_DIR,
        "F10L_candidate_TF_audit.csv"
    )
)


write_csv(
    top_positive,
    file.path(
        OUTPUT_DIR,
        "F10L_top_activated_TFs.csv"
    )
)


write_csv(
    top_negative,
    file.path(
        OUTPUT_DIR,
        "F10L_top_reduced_TFs.csv"
    )
)


# ============================================================
# CONSOLE
# ============================================================

message("")
message(
    "============================================================"
)

message(
    "F10L — ADULT SIGNED TF ACTIVITY"
)

message(
    "============================================================"
)

message("")

message(
    "Primary genes: ",
    nrow(
        primary
    )
)

message(
    "Eligible regulators: ",
    length(
        eligible_regulators
    )
)

message("")


message(
    "TOP INFERRED ACTIVATED TFs"
)

message(
    "--------------------------"
)


print(
    as.data.frame(
        top_positive %>%
            select(
                TF,
                Contrast_activity,
                Contrast_activity_FDR,
                Delta_activity,
                Hedges_g,
                Exact_permutation_p,
                Contrast_replicate_sign_concordant,
                Sensitivity_activity,
                Primary_sensitivity_sign_concordant,
                Targets_in_primary
            )
    ),
    row.names = FALSE
)


message("")


message(
    "TOP INFERRED REDUCED TFs"
)

message(
    "------------------------"
)


print(
    as.data.frame(
        top_negative %>%
            select(
                TF,
                Contrast_activity,
                Contrast_activity_FDR,
                Delta_activity,
                Hedges_g,
                Exact_permutation_p,
                Contrast_replicate_sign_concordant,
                Sensitivity_activity,
                Primary_sensitivity_sign_concordant,
                Targets_in_primary
            )
    ),
    row.names = FALSE
)


message("")


message(
    "PRE-SPECIFIED CANDIDATE AUDIT"
)

message(
    "-----------------------------"
)


if (nrow(
    candidate_audit
) == 0) {

    message(
        "None of the pre-specified candidate TFs were present."
    )

} else {

    print(
        as.data.frame(
            candidate_audit %>%
                select(
                    TF,
                    Contrast_activity,
                    Contrast_activity_p,
                    Contrast_activity_FDR,
                    Delta_activity,
                    Hedges_g,
                    Exact_permutation_p,
                    Contrast_replicate_sign_concordant,
                    Sensitivity_activity,
                    Sensitivity_activity_FDR,
                    Primary_sensitivity_sign_concordant,
                    Targets_in_primary
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
    "INTERPRETATION GUARDRAILS:"
)

message(
    paste(
        "Contrast_activity is a CollecTRI/ULM inference from the",
        "full DESeq2 Wald-statistic profile; its FDR is model-level",
        "target-pattern evidence, not biological replicate FDR."
    )
)

message(
    paste(
        "The replicate-level analysis is an independent consistency",
        "check. With only 3+3 samples, the exact two-sided permutation",
        "test has very coarse resolution, so effect size and sign",
        "consistency matter more than its nominal p-value."
    )
)