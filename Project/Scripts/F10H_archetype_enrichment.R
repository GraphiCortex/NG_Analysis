#!/usr/bin/env Rscript

# ============================================================
# F10H_archetype_enrichment.R
# ============================================================
#
# Functional enrichment of the EMPIRICALLY discovered F10G
# response archetypes.
#
# PRIMARY QUESTION
# ----------------
# Once genes are grouped solely by empirical perturbational
# response profile, do those groups show biological coherence?
#
# IMPORTANT:
# The enrichment is performed AFTER clustering.
# Biological annotations did not contribute to cluster formation.
#
# PRIMARY BACKGROUND:
#   - adult-Up archetypes -> all adult-Up complete-context genes
#   - adult-Down archetypes -> all adult-Down complete-context genes
#
# SECONDARY BACKGROUND:
#   - all 1,206 complete-context genes
#
# This prevents adult-direction selection from trivially driving
# the functional enrichment.
#
# OUTPUT:
# Data/Processed/Dynamics/F10H_ArchetypeEnrichment/
#
# ============================================================


suppressPackageStartupMessages({

    library(clusterProfiler)
    library(org.Mm.eg.db)
    library(AnnotationDbi)
    library(dplyr)
    library(readr)
    library(ggplot2)
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
    "F10H_ArchetypeEnrichment"
)

FIGURE_DIR <- file.path(
    OUTPUT_DIR,
    "figures"
)

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    FIGURE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# PARAMETERS
# ============================================================

PRIMARY_PADJ_THRESHOLD <- 0.05

MIN_GS_SIZE <- 10

MAX_GS_SIZE <- 500

TOP_TERMS_TO_PRINT <- 15

TOP_TERMS_TO_PLOT <- 12


# ============================================================
# LOAD F10G ASSIGNMENTS
# ============================================================

message("Loading F10G cluster assignments...")

df <- read_csv(
    INPUT_FILE,
    show_col_types = FALSE
)


required_columns <- c(
    "Gene",
    "RbTKO_Direction",
    "Archetype_ID",
    "Archetype_key",
    "RbTKO_LFC",
    "E2F13_LFC",
    "CAG_LFC",
    "MAP2_LFC"
)

missing_columns <- setdiff(
    required_columns,
    colnames(df)
)

if (length(missing_columns) > 0) {

    stop(
        paste(
            "Missing required columns:",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    )
}


message(
    "Input genes: ",
    nrow(df)
)


# ============================================================
# MAP SYMBOLS -> ENTREZ
# ============================================================
#
# IMPORTANT:
# Gene symbols in the computational pipeline are largely stored
# in uppercase, whereas official mouse symbols in org.Mm.eg.db
# use mouse-specific capitalization (e.g. E2f1, Myc, Cdkn1a).
#
# AnnotationDbi::select() is case-sensitive, so direct lookup
# would incorrectly lose most genes.
#
# We therefore:
#
#   1. retrieve the complete official mouse SYMBOL -> ENTREZID
#      annotation table;
#   2. create an uppercase matching key on BOTH sides;
#   3. match case-insensitively;
#   4. explicitly audit any uppercase-key collisions.
#
# This changes only identifier matching. It does NOT modify the
# biological gene symbols stored in F10G.
# ============================================================

message(
    "Mapping mouse gene symbols to Entrez IDs "
)

message(
    "using case-insensitive official-symbol matching..."
)


# ------------------------------------------------------------
# Retrieve complete official mouse annotation
# ------------------------------------------------------------

all_mouse_symbols <- AnnotationDbi::keys(
    org.Mm.eg.db,
    keytype = "SYMBOL"
)


mouse_symbol_map <- AnnotationDbi::select(

    org.Mm.eg.db,

    keys = all_mouse_symbols,

    keytype = "SYMBOL",

    columns = c(
        "SYMBOL",
        "ENTREZID"
    )
)


mouse_symbol_map <- mouse_symbol_map %>%

    filter(
        !is.na(SYMBOL),
        !is.na(ENTREZID)
    ) %>%

    mutate(

        Gene_match_key =
            toupper(
                trimws(
                    SYMBOL
                )
            )
    )


# ------------------------------------------------------------
# Audit uppercase-key collisions
#
# Example problem:
# two different official mouse symbols could theoretically map
# to the same uppercase representation.
#
# We do not silently ignore this.
# ------------------------------------------------------------

collision_table <- mouse_symbol_map %>%

    distinct(
        Gene_match_key,
        SYMBOL,
        ENTREZID
    ) %>%

    group_by(
        Gene_match_key
    ) %>%

    summarise(

        N_symbols =
            n_distinct(
                SYMBOL
            ),

        N_entrez =
            n_distinct(
                ENTREZID
            ),

        Symbols =
            paste(
                sort(
                    unique(
                        SYMBOL
                    )
                ),
                collapse = ";"
            ),

        Entrez_IDs =
            paste(
                sort(
                    unique(
                        ENTREZID
                    )
                ),
                collapse = ";"
            ),

        .groups =
            "drop"
    ) %>%

    filter(
        N_symbols > 1
        |
        N_entrez > 1
    )


write_csv(

    collision_table,

    file.path(
        OUTPUT_DIR,
        "F10H_symbol_case_collisions.csv"
    )
)


message(
    "Uppercase annotation-key collisions: ",
    nrow(
        collision_table
    )
)


# ------------------------------------------------------------
# For unambiguous matching keys, retain one SYMBOL -> ENTREZ ID.
#
# If an uppercase key is genuinely ambiguous, exclude it rather
# than choosing an arbitrary annotation.
# ------------------------------------------------------------

ambiguous_keys <- collision_table$Gene_match_key


mouse_symbol_map_clean <- mouse_symbol_map %>%

    filter(
        !Gene_match_key %in%
            ambiguous_keys
    ) %>%

    distinct(
        Gene_match_key,
        .keep_all = TRUE
    )


# ------------------------------------------------------------
# Build matching key for the F10G genes
# ------------------------------------------------------------

df <- df %>%

    mutate(

        Gene_match_key =
            toupper(
                trimws(
                    Gene
                )
            )
    )


# ------------------------------------------------------------
# Join case-insensitively
# ------------------------------------------------------------

df_mapped <- df %>%

    left_join(

        mouse_symbol_map_clean %>%

            select(
                Gene_match_key,
                Official_mouse_symbol =
                    SYMBOL,
                ENTREZID
            ),

        by =
            "Gene_match_key"
    )


# ------------------------------------------------------------
# Mapping audit
# ------------------------------------------------------------

mapping_summary <- tibble(

    total_genes =
        nrow(
            df
        ),

    mapped_genes =
        sum(
            !is.na(
                df_mapped$ENTREZID
            )
        ),

    unmapped_genes =
        sum(
            is.na(
                df_mapped$ENTREZID
            )
        ),

    mapped_fraction =
        mean(
            !is.na(
                df_mapped$ENTREZID
            )
        ),

    ambiguous_case_keys_excluded =
        sum(
            df_mapped$Gene_match_key
            %in%
            ambiguous_keys
        )
)


write_csv(

    mapping_summary,

    file.path(
        OUTPUT_DIR,
        "F10H_mapping_summary.csv"
    )
)


# ------------------------------------------------------------
# Save gene-level mapping audit
# ------------------------------------------------------------

mapping_audit <- df_mapped %>%

    select(

        Gene,

        Gene_match_key,

        Official_mouse_symbol,

        ENTREZID,

        RbTKO_Direction,

        Archetype_key
    ) %>%

    mutate(

        Mapping_status =
            case_when(

                Gene_match_key
                %in%
                ambiguous_keys
                ~
                "Ambiguous_uppercase_key",

                !is.na(
                    ENTREZID
                )
                ~
                "Mapped",

                TRUE
                ~
                "Unmapped"
            )
    )


write_csv(

    mapping_audit,

    file.path(
        OUTPUT_DIR,
        "F10H_gene_mapping_audit.csv"
    )
)


message(
    "Mapped ",
    mapping_summary$mapped_genes,
    "/",
    mapping_summary$total_genes,
    " genes (",
    round(
        100 *
        mapping_summary$mapped_fraction,
        1
    ),
    "%)."
)

if (
    mapping_summary$mapped_fraction
    <
    0.80
) {

    warning(
        paste0(
            "Only ",
            round(
                100 *
                mapping_summary$mapped_fraction,
                1
            ),
            "% of F10G genes mapped to Entrez IDs. ",
            "Do not interpret enrichment until the mapping problem is understood."
        )
    )
}

# ============================================================
# HELPER: RUN GO-BP ENRICHMENT
# ============================================================

run_go_bp <- function(
    target_entrez,
    universe_entrez,
    archetype_key,
    background_name
) {

    target_entrez <- unique(
        na.omit(
            target_entrez
        )
    )

    universe_entrez <- unique(
        na.omit(
            universe_entrez
        )
    )


    # Ensure target lies inside its declared universe.
    target_entrez <- intersect(
        target_entrez,
        universe_entrez
    )


    if (length(target_entrez) < MIN_GS_SIZE) {

        return(
            NULL
        )
    }


    result <- enrichGO(

        gene = target_entrez,

        universe = universe_entrez,

        OrgDb = org.Mm.eg.db,

        keyType = "ENTREZID",

        ont = "BP",

        pAdjustMethod = "BH",

        pvalueCutoff = 1,

        qvalueCutoff = 1,

        minGSSize = MIN_GS_SIZE,

        maxGSSize = MAX_GS_SIZE,

        readable = TRUE
    )


    result_df <- as.data.frame(
        result
    )


    if (nrow(result_df) == 0) {

        return(
            NULL
        )
    }


    result_df <- result_df %>%

        mutate(

            Archetype_key =
                archetype_key,

            Background =
                background_name,

            Target_n =
                length(
                    target_entrez
                ),

            Universe_n =
                length(
                    universe_entrez
                )
        )


    return(
        result_df
    )
}


# ============================================================
# PRIMARY ENRICHMENT
#
# Background = same adult Rb-TKO direction
# ============================================================

message("")
message(
    "Running PRIMARY within-direction GO-BP enrichment..."
)


primary_results <- list()


archetypes <- sort(
    unique(
        df_mapped$Archetype_key
    )
)


for (archetype in archetypes) {

    group <- df_mapped %>%

        filter(
            Archetype_key ==
                archetype
        )


    adult_direction <- unique(
        group$RbTKO_Direction
    )


    if (length(adult_direction) != 1) {

        stop(
            paste(
                "Archetype",
                archetype,
                "contains multiple adult directions."
            )
        )
    }


    stratum_background <- df_mapped %>%

        filter(
            RbTKO_Direction ==
                adult_direction
        )


    result <- run_go_bp(

        target_entrez =
            group$ENTREZID,

        universe_entrez =
            stratum_background$ENTREZID,

        archetype_key =
            archetype,

        background_name =
            paste0(
                "Adult_",
                adult_direction,
                "_stratum"
            )
    )


    if (!is.null(result)) {

        primary_results[[
            archetype
        ]] <- result
    }
}


if (
    length(
        primary_results
    ) > 0
) {

    primary_table <- bind_rows(
        primary_results
    )

} else {

    primary_table <- tibble()
}


write_csv(
    primary_table,
    file.path(
        OUTPUT_DIR,
        "F10H_GO_BP_primary_within_direction.csv"
    )
)


# ============================================================
# SECONDARY ENRICHMENT
#
# Background = all 1,206 genes
# ============================================================

message(
    "Running SECONDARY whole-F10G-universe enrichment..."
)


global_universe <- unique(
    na.omit(
        df_mapped$ENTREZID
    )
)


secondary_results <- list()


for (archetype in archetypes) {

    group <- df_mapped %>%

        filter(
            Archetype_key ==
                archetype
        )


    result <- run_go_bp(

        target_entrez =
            group$ENTREZID,

        universe_entrez =
            global_universe,

        archetype_key =
            archetype,

        background_name =
            "All_complete_context_genes"
    )


    if (!is.null(result)) {

        secondary_results[[
            archetype
        ]] <- result
    }
}


if (
    length(
        secondary_results
    ) > 0
) {

    secondary_table <- bind_rows(
        secondary_results
    )

} else {

    secondary_table <- tibble()
}


write_csv(
    secondary_table,
    file.path(
        OUTPUT_DIR,
        "F10H_GO_BP_secondary_global.csv"
    )
)


# ============================================================
# PRIMARY SIGNIFICANT TERMS
# ============================================================

if (
    nrow(
        primary_table
    ) > 0
) {

    primary_significant <- primary_table %>%

        filter(
            p.adjust <
                PRIMARY_PADJ_THRESHOLD
        ) %>%

        arrange(
            Archetype_key,
            p.adjust
        )

} else {

    primary_significant <- tibble()
}


write_csv(
    primary_significant,
    file.path(
        OUTPUT_DIR,
        "F10H_GO_BP_primary_significant_FDR05.csv"
    )
)


# ============================================================
# SUMMARY PER ARCHETYPE
# ============================================================

summary_rows <- list()


for (archetype in archetypes) {

    genes <- df_mapped %>%

        filter(
            Archetype_key ==
                archetype
        )


    if (
        nrow(
            primary_table
        ) > 0
    ) {

        terms <- primary_table %>%

            filter(
                Archetype_key ==
                    archetype
            ) %>%

            arrange(
                p.adjust
            )

    } else {

        terms <- tibble()
    }


    significant_terms <- terms %>%

        filter(
            p.adjust <
                PRIMARY_PADJ_THRESHOLD
        )


    top_descriptions <- if (
        nrow(
            significant_terms
        ) > 0
    ) {

        paste(
            head(
                significant_terms$Description,
                10
            ),
            collapse = "; "
        )

    } else {

        ""
    }


    summary_rows[[
        archetype
    ]] <- tibble(

        Archetype_key =
            archetype,

        Adult_direction =
            first(
                genes$RbTKO_Direction
            ),

        N_genes =
            nrow(
                genes
            ),

        N_mapped =
            sum(
                !is.na(
                    genes$ENTREZID
                )
            ),

        N_significant_GO_BP =
            nrow(
                significant_terms
            ),

        Best_GO_term =
            if (
                nrow(
                    terms
                ) > 0
            ) {
                first(
                    terms$Description
                )
            } else {
                NA_character_
            },

        Best_GO_padj =
            if (
                nrow(
                    terms
                ) > 0
            ) {
                first(
                    terms$p.adjust
                )
            } else {
                NA_real_
            },

        Top_significant_terms =
            top_descriptions
    )
}


archetype_summary <- bind_rows(
    summary_rows
)


write_csv(
    archetype_summary,
    file.path(
        OUTPUT_DIR,
        "F10H_archetype_enrichment_summary.csv"
    )
)


# ============================================================
# PLOT TOP TERMS FOR EACH ARCHETYPE
# ============================================================

if (
    nrow(
        primary_table
    ) > 0
) {

    for (archetype in archetypes) {

        plot_df <- primary_table %>%

            filter(
                Archetype_key ==
                    archetype
            ) %>%

            arrange(
                p.adjust
            ) %>%

            slice_head(
                n =
                    TOP_TERMS_TO_PLOT
            )


        if (
            nrow(
                plot_df
            ) == 0
        ) {

            next
        }


        plot_df <- plot_df %>%

            mutate(

                minus_log10_FDR =
                    -log10(
                        pmax(
                            p.adjust,
                            1e-300
                        )
                    ),

                Description =
                    factor(
                        Description,
                        levels =
                            rev(
                                Description
                            )
                    )
            )


        p <- ggplot(
            plot_df,
            aes(
                x =
                    minus_log10_FDR,

                y =
                    Description
            )
        ) +

            geom_point(
                aes(
                    size =
                        Count
                )
            ) +

            labs(
                title =
                    paste0(
                        archetype,
                        ": GO biological process"
                    ),

                subtitle =
                    "Primary background: same adult Rb-TKO direction",

                x =
                    expression(
                        -log[10](
                            "BH FDR"
                        )
                    ),

                y =
                    NULL,

                size =
                    "Genes"
            ) +

            theme_bw(
                base_size = 11
            ) +

            theme(
                plot.title.position =
                    "plot"
            )


        ggsave(

            filename =
                file.path(
                    FIGURE_DIR,
                    paste0(
                        "F10H_",
                        archetype,
                        "_GO_BP.png"
                    )
                ),

            plot =
                p,

            width =
                9,

            height =
                6,

            dpi =
                300
        )
    }
}


# ============================================================
# CONSOLE OUTPUT
# ============================================================

message("")
message(
    "============================================================"
)

message(
    "F10H PRIMARY WITHIN-DIRECTION ENRICHMENT"
)

message(
    "============================================================"
)

message("")


for (archetype in archetypes) {

    message(
        archetype
    )

    message(
        paste(
            rep(
                "-",
                nchar(
                    archetype
                )
            ),
            collapse = ""
        )
    )


    terms <- primary_table %>%

        filter(
            Archetype_key ==
                archetype
        ) %>%

        arrange(
            p.adjust
        )


    if (
        nrow(
            terms
        ) == 0
    ) {

        message(
            "No GO-BP terms returned."
        )

        message("")

        next
    }


    console_terms <- terms %>%

        select(
            ID,
            Description,
            GeneRatio,
            BgRatio,
            p.adjust,
            Count
        ) %>%

        slice_head(
            n =
                TOP_TERMS_TO_PRINT
        )


    print(
        as.data.frame(
            console_terms
        ),
        row.names = FALSE
    )


    n_significant <- sum(
        terms$p.adjust
        <
        PRIMARY_PADJ_THRESHOLD
    )


    message(
        "Significant GO-BP terms at BH FDR < 0.05: ",
        n_significant
    )

    message("")
}


message(
    "============================================================"
)

message(
    "ARCHETYPE SUMMARY"
)

message(
    "============================================================"
)

print(
    as.data.frame(
        archetype_summary
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
    "IMPORTANT:"
)

message(
    "Do not biologically name an archetype from a single redundant "
)

message(
    "GO term. Interpret coherent term families, effect direction, "
)

message(
    "and the empirical perturbation profile together."
)