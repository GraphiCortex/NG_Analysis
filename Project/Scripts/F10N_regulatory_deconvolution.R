#!/usr/bin/env Rscript

# ============================================================
# F10N_regulatory_deconvolution.R
# ============================================================
#
# Adult Rb-TKO regulatory deconvolution using the full DESeq2
# Wald-statistic profile and frozen mouse CollecTRI.
#
# Tests:
#
#   FULL
#
#   NO_IMMUNE_GO
#       remove GO:0006955 (immune response) plus GOALL descendants
#
#   NO_CELL_CYCLE_GO
#       remove GO:0007049 (cell cycle) plus GOALL descendants
#
#   NO_E2F_FAMILY_TARGETS
#       remove the measured union of CollecTRI targets belonging to
#       E2F1 / E2F2 / E2F3 / E2F4 / TFDP1
#
# Every condition is analyzed using:
#
#   ULM = univariate regulator fitting
#   MLM = multivariate joint regulator fitting
#
# Main biological questions:
#
#   1. Does NF-kB / IRF / STAT-like activity survive removal
#      of immune-response genes?
#
#   2. Does MYC activity survive removal of cell-cycle genes?
#
#   3. Does MYC activity survive removal of the E2F-family
#      target universe itself?
#
# IMPORTANT:
# MLM reduces shared-target confounding by fitting regulators
# jointly. It does not establish causal TF activity.
#
# ============================================================


# ============================================================
# PACKAGES
# ============================================================

required <- c(
    "decoupleR",
    "readr",
    "AnnotationDbi",
    "org.Mm.eg.db"
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


suppressPackageStartupMessages({

    library(AnnotationDbi)

    library(org.Mm.eg.db)
})


# ============================================================
# PATHS
# ============================================================

ADULT_FILE <- file.path(
    "Data",
    "Raw",
    "RNAseq",
    "Fong2022",
    "GSE190765",
    "GSE190765_FDSC_FC.Batch.IHW.csv.gz"
)


COLLECTRI_FILE <- file.path(
    "Data",
    "Raw",
    "Regulatory",
    "CollecTRI",
    "interactions_collectri_10090.tsv.gz"
)


OUTPUT_DIR <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F10N_RegulatoryDeconvolution"
)


dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# PARAMETERS
# ============================================================

BASEMEAN_THRESHOLD <- 10

MIN_REGULON_SIZE <- 5


IMMUNE_GO <- "GO:0006955"

CELL_CYCLE_GO <- "GO:0007049"


E2F_FAMILY <- c(
    "E2F1",
    "E2F2",
    "E2F3",
    "E2F4",
    "TFDP1"
)


CANDIDATE_TFS <- c(

    "E2F1",
    "E2F2",
    "E2F3",
    "E2F4",
    "TFDP1",
    "FOXM1",

    "MYC",

    "REL",
    "RELA",
    "RELB",
    "NFKB",
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
    "IRF4",
    "IRF5",
    "IRF7",
    "IRF8",
    "IRF9",

    "JUN",
    "FOS",
    "AP1",
    "CEBPD",
    "CEBPB",
    "DDIT3",
    "ATF4",

    "TP53",

    "NR4A1",
    "NR4A2",
    "NR4A3"
)


# ============================================================
# BASIC HELPERS
# ============================================================

safe_fdr <- function(p) {

    p.adjust(
        as.numeric(p),
        method = "BH"
    )
}


sign_same <- function(
    x,
    y
) {

    answer <- rep(
        NA,
        length(x)
    )


    valid <- (
        is.finite(x) &
        is.finite(y)
    )


    answer[
        valid
    ] <- (
        sign(
            x[
                valid
            ]
        ) ==
        sign(
            y[
                valid
            ]
        )
    )


    answer
}


safe_ratio <- function(
    numerator,
    denominator
) {

    answer <- rep(
        NA_real_,
        length(numerator)
    )


    valid <- (
        is.finite(numerator) &
        is.finite(denominator) &
        abs(
            denominator
        ) > 1e-12
    )


    answer[
        valid
    ] <- (
        numerator[
            valid
        ] /
        denominator[
            valid
        ]
    )


    answer
}


# ============================================================
# CLEAN ADULT RNA-SEQ
# ============================================================

clean_rnaseq <- function(
    df,
    threshold
) {

    needed <- c(
        "name",
        "baseMean",
        "stat"
    )


    absent <- setdiff(
        needed,
        names(df)
    )


    if (length(absent) > 0) {

        stop(
            paste(
                "RNA-seq table missing:",
                paste(
                    absent,
                    collapse = ", "
                )
            )
        )
    }


    gene <- trimws(
        as.character(
            df[
                ,
                "name",
                drop = TRUE
            ]
        )
    )


    gene_key <- toupper(
        gene
    )


    base_mean <- suppressWarnings(
        as.numeric(
            df[
                ,
                "baseMean",
                drop = TRUE
            ]
        )
    )


    stat_value <- suppressWarnings(
        as.numeric(
            df[
                ,
                "stat",
                drop = TRUE
            ]
        )
    )


    keep <- (
        !is.na(gene) &
        gene != "" &
        gene != "NA" &
        is.finite(base_mean)
    )


    out <- data.frame(

        Gene =
            gene[
                keep
            ],

        Gene_key =
            gene_key[
                keep
            ],

        baseMean =
            base_mean[
                keep
            ],

        stat =
            stat_value[
                keep
            ],

        stringsAsFactors = FALSE
    )


    # --------------------------------------------------------
    # Reproduce F10L/F10M deduplication:
    #
    # highest-baseMean row per case-insensitive symbol.
    #
    # This selection is independent of differential direction
    # and statistical significance.
    # --------------------------------------------------------

    ordering <- order(
        out$Gene_key,
        -out$baseMean
    )


    out <- out[
        ordering,
        ,
        drop = FALSE
    ]


    out <- out[
        !duplicated(
            out$Gene_key
        ),
        ,
        drop = FALSE
    ]


    out <- out[
        is.finite(
            out$stat
        ) &
        out$baseMean >=
            threshold,
        ,
        drop = FALSE
    ]


    rownames(
        out
    ) <- NULL


    out
}


# ============================================================
# LOAD AND FORMAT FROZEN MOUSE COLLECTRI
# ============================================================

load_collectri <- function(
    path
) {

    if (!file.exists(
        path
    )) {

        stop(
            paste(
                "Missing CollecTRI file:",
                path
            )
        )
    }


    raw <- readr::read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE
    )


    needed <- c(
        "source_genesymbol",
        "target_genesymbol",
        "is_stimulation",
        "is_inhibition"
    )


    absent <- setdiff(
        needed,
        names(raw)
    )


    if (length(absent) > 0) {

        stop(
            paste(
                "Unexpected CollecTRI schema; missing:",
                paste(
                    absent,
                    collapse = ", "
                )
            )
        )
    }


    source_symbol <- as.character(
        raw[
            ,
            "source_genesymbol",
            drop = TRUE
        ]
    )


    target_symbol <- as.character(
        raw[
            ,
            "target_genesymbol",
            drop = TRUE
        ]
    )


    if ("source" %in%
        names(raw)) {

        source_id <- as.character(
            raw[
                ,
                "source",
                drop = TRUE
            ]
        )

    } else {

        source_id <- source_symbol
    }


    # --------------------------------------------------------
    # Collapse known TF complexes in the same conceptual way
    # used previously.
    # --------------------------------------------------------

    is_complex <- grepl(
        "COMPLEX",
        source_id,
        fixed = TRUE
    )


    is_ap1 <- (
        is_complex &
        (
            grepl(
                "JUN",
                source_symbol,
                ignore.case = TRUE
            ) |
            grepl(
                "FOS",
                source_symbol,
                ignore.case = TRUE
            )
        )
    )


    is_nfkb <- (
        is_complex &
        (
            grepl(
                "REL",
                source_symbol,
                ignore.case = TRUE
            ) |
            grepl(
                "NFKB",
                source_symbol,
                ignore.case = TRUE
            )
        )
    )


    source_symbol[
        is_ap1
    ] <- "AP1"


    source_symbol[
        is_nfkb
    ] <- "NFKB"


    stimulation <- suppressWarnings(
        as.numeric(
            raw[
                ,
                "is_stimulation",
                drop = TRUE
            ]
        )
    )


    inhibition <- suppressWarnings(
        as.numeric(
            raw[
                ,
                "is_inhibition",
                drop = TRUE
            ]
        )
    )


    sign_known <- (
        is.finite(
            stimulation
        ) &
        is.finite(
            inhibition
        )
    )


    activation_only <- (
        sign_known &
        stimulation == 1 &
        inhibition == 0
    )


    inhibition_only <- (
        sign_known &
        stimulation == 0 &
        inhibition == 1
    )


    ambiguous <- (
        sign_known &
        !(
            activation_only |
            inhibition_only
        )
    )


    unknown <- !sign_known


    message(
        "CollecTRI sign audit: activation=",
        sum(
            activation_only
        ),
        ", inhibition=",
        sum(
            inhibition_only
        ),
        ", ambiguous/unsigned=",
        sum(
            ambiguous
        ),
        ", unknown=",
        sum(
            unknown
        )
    )


    source_key <- toupper(
        trimws(
            source_symbol
        )
    )


    target_key <- toupper(
        trimws(
            target_symbol
        )
    )


    valid <- (
        !is.na(
            source_key
        ) &
        source_key != "" &
        !is.na(
            target_key
        ) &
        target_key != "" &
        (
            activation_only |
            inhibition_only
        )
    )


    mor <- rep(
        NA_real_,
        length(
            source_key
        )
    )


    mor[
        activation_only
    ] <- 1


    mor[
        inhibition_only
    ] <- -1


    network <- data.frame(

        source =
            source_key[
                valid
            ],

        target =
            target_key[
                valid
            ],

        mor =
            mor[
                valid
            ],

        stringsAsFactors = FALSE
    )


    network <- unique(
        network
    )


    # --------------------------------------------------------
    # Remove TF-target pairs that appear with conflicting signs.
    # --------------------------------------------------------

    pair_key <- paste(
        network$source,
        network$target,
        sep = "|||"
    )


    groups <- split(
        seq_len(
            nrow(
                network
            )
        ),
        pair_key
    )


    non_conflicting <- vapply(

        groups,

        function(index) {

            length(
                unique(
                    network$mor[
                        index
                    ]
                )
            ) == 1
        },

        logical(1)
    )


    valid_pairs <- names(
        non_conflicting[
            non_conflicting
        ]
    )


    network <- network[
        pair_key %in%
            valid_pairs,
        ,
        drop = FALSE
    ]


    pair_key <- paste(
        network$source,
        network$target,
        sep = "|||"
    )


    network <- network[
        !duplicated(
            pair_key
        ),
        ,
        drop = FALSE
    ]


    rownames(
        network
    ) <- NULL


    if (nrow(
        network
    ) == 0) {

        stop(
            "No usable signed CollecTRI edges remain."
        )
    }


    network
}


# ============================================================
# BUILD GOALL CONDITIONING SETS
# ============================================================

build_go_sets <- function(
    gene_keys
) {

    available_columns <- AnnotationDbi::columns(
        org.Mm.eg.db
    )


    needed <- c(
        "GOALL",
        "ONTOLOGYALL"
    )


    if (!all(
        needed %in%
            available_columns
    )) {

        stop(
            paste(
                "org.Mm.eg.db lacks:",
                paste(
                    setdiff(
                        needed,
                        available_columns
                    ),
                    collapse = ", "
                )
            )
        )
    }


    official_symbols <- AnnotationDbi::keys(
        org.Mm.eg.db,
        keytype = "SYMBOL"
    )


    symbol_map <- data.frame(

        SYMBOL =
            official_symbols,

        Gene_key =
            toupper(
                official_symbols
            ),

        stringsAsFactors = FALSE
    )


    # --------------------------------------------------------
    # Avoid ambiguous case-insensitive symbol mappings.
    # --------------------------------------------------------

    counts <- table(
        symbol_map$Gene_key
    )


    unique_upper <- names(
        counts[
            counts == 1
        ]
    )


    symbol_map <- symbol_map[
        symbol_map$Gene_key %in%
            unique_upper,
        ,
        drop = FALSE
    ]


    mapped <- symbol_map$SYMBOL[
        match(
            gene_keys,
            symbol_map$Gene_key
        )
    ]


    mapped <- unique(
        mapped[
            !is.na(
                mapped
            )
        ]
    )


    message(
        "GO mapping: ",
        length(
            mapped
        ),
        " / ",
        length(
            unique(
                gene_keys
            )
        ),
        " adult genes mapped"
    )


    go_table <- suppressMessages(

        AnnotationDbi::select(

            org.Mm.eg.db,

            keys =
                mapped,

            keytype =
                "SYMBOL",

            columns =
                c(
                    "GOALL",
                    "ONTOLOGYALL"
                )
        )
    )


    immune_symbol <- unique(
        as.character(
            go_table$SYMBOL[
                go_table$GOALL ==
                    IMMUNE_GO &
                go_table$ONTOLOGYALL ==
                    "BP"
            ]
        )
    )


    cycle_symbol <- unique(
        as.character(
            go_table$SYMBOL[
                go_table$GOALL ==
                    CELL_CYCLE_GO &
                go_table$ONTOLOGYALL ==
                    "BP"
            ]
        )
    )


    list(

        immune =
            unique(
                toupper(
                    immune_symbol
                )
            ),

        cell_cycle =
            unique(
                toupper(
                    cycle_symbol
                )
            )
    )
}


# ============================================================
# RUN ONE ACTIVITY ANALYSIS
# ============================================================

run_activity <- function(
    df,
    network,
    condition_name,
    method_name
) {

    genes <- as.character(
        df[
            ,
            "Gene_key",
            drop = TRUE
        ]
    )


    evidence <- as.numeric(
        df[
            ,
            "stat",
            drop = TRUE
        ]
    )


    keep <- (
        !is.na(
            genes
        ) &
        genes != "" &
        is.finite(
            evidence
        )
    )


    genes <- genes[
        keep
    ]


    evidence <- evidence[
        keep
    ]


    first_occurrence <- !duplicated(
        genes
    )


    genes <- genes[
        first_occurrence
    ]


    evidence <- evidence[
        first_occurrence
    ]


    context_network <- network[
        network$target %in%
            genes,
        ,
        drop = FALSE
    ]


    regulon_sizes <- table(
        context_network$source
    )


    eligible <- names(
        regulon_sizes[
            regulon_sizes >=
                MIN_REGULON_SIZE
        ]
    )


    context_network <- context_network[
        context_network$source %in%
            eligible,
        ,
        drop = FALSE
    ]


    if (nrow(
        context_network
    ) == 0) {

        stop(
            paste(
                "No eligible regulons for",
                condition_name,
                method_name
            )
        )
    }


    matrix_input <- matrix(
        evidence,
        ncol = 1
    )


    rownames(
        matrix_input
    ) <- genes


    colnames(
        matrix_input
    ) <- condition_name


    if (method_name == "ULM") {

        fit <- decoupleR::run_ulm(

            mat =
                matrix_input,

            network =
                context_network,

            .source =
                "source",

            .target =
                "target",

            .mor =
                "mor",

            center =
                FALSE,

            minsize =
                MIN_REGULON_SIZE
        )

    } else if (
        method_name == "MLM"
    ) {

        fit <- decoupleR::run_mlm(

            mat =
                matrix_input,

            network =
                context_network,

            .source =
                "source",

            .target =
                "target",

            .mor =
                "mor",

            center =
                FALSE,

            minsize =
                MIN_REGULON_SIZE
        )

    } else {

        stop(
            paste(
                "Unknown method:",
                method_name
            )
        )
    }


    if (!"p_value" %in%
        names(fit)) {

        fit$p_value <- NA_real_
    }


    tf <- as.character(
        fit[
            ,
            "source",
            drop = TRUE
        ]
    )


    activity <- as.numeric(
        fit[
            ,
            "score",
            drop = TRUE
        ]
    )


    p_value <- as.numeric(
        fit[
            ,
            "p_value",
            drop = TRUE
        ]
    )


    result <- data.frame(

        TF =
            tf,

        Condition =
            condition_name,

        Method =
            method_name,

        Activity =
            activity,

        P_value =
            p_value,

        FDR =
            safe_fdr(
                p_value
            ),

        Targets_measured =
            as.numeric(
                regulon_sizes[
                    tf
                ]
            ),

        stringsAsFactors = FALSE
    )


    result <- result[
        order(
            result$FDR,
            -abs(
                result$Activity
            ),
            na.last = TRUE
        ),
        ,
        drop = FALSE
    ]


    rownames(
        result
    ) <- NULL


    result
}


# ============================================================
# EXTRACT ONE RESULT INTO A WIDE TABLE
# ============================================================

extract_condition <- function(
    long_table,
    condition_name,
    method_name,
    prefix
) {

    out <- long_table[

        long_table$Condition ==
            condition_name &
        long_table$Method ==
            method_name,

        c(
            "TF",
            "Activity",
            "P_value",
            "FDR",
            "Targets_measured"
        ),

        drop = FALSE
    ]


    names(
        out
    ) <- c(
        "TF",
        paste0(
            prefix,
            "_Activity"
        ),
        paste0(
            prefix,
            "_P"
        ),
        paste0(
            prefix,
            "_FDR"
        ),
        paste0(
            prefix,
            "_Targets"
        )
    )


    out
}


# ============================================================
# CANDIDATE REGULON OVERLAP
# ============================================================

make_overlap_table <- function(
    network,
    universe
) {

    network <- network[
        network$target %in%
            universe,
        ,
        drop = FALSE
    ]


    regulators <- intersect(
        CANDIDATE_TFS,
        unique(
            network$source
        )
    )


    if (length(
        regulators
    ) < 2) {

        return(
            data.frame()
        )
    }


    pairs <- combn(
        regulators,
        2
    )


    out <- data.frame(

        TF1 =
            character(0),

        TF2 =
            character(0),

        TF1_targets =
            integer(0),

        TF2_targets =
            integer(0),

        Shared_targets =
            integer(0),

        Jaccard =
            numeric(0),

        stringsAsFactors = FALSE
    )


    for (i in seq_len(
        ncol(
            pairs
        )
    )) {

        tf1 <- pairs[
            1,
            i
        ]


        tf2 <- pairs[
            2,
            i
        ]


        targets_1 <- unique(
            network$target[
                network$source ==
                    tf1
            ]
        )


        targets_2 <- unique(
            network$target[
                network$source ==
                    tf2
            ]
        )


        shared <- length(
            intersect(
                targets_1,
                targets_2
            )
        )


        union_n <- length(
            union(
                targets_1,
                targets_2
            )
        )


        jaccard <- if (
            union_n > 0
        ) {

            shared /
                union_n

        } else {

            NA_real_
        }


        one_row <- data.frame(

            TF1 =
                tf1,

            TF2 =
                tf2,

            TF1_targets =
                length(
                    targets_1
                ),

            TF2_targets =
                length(
                    targets_2
                ),

            Shared_targets =
                shared,

            Jaccard =
                jaccard,

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
# LOAD ADULT PROFILE
# ============================================================

message(
    "Loading adult full-profile RNA-seq..."
)


adult_raw <- readr::read_csv(
    ADULT_FILE,
    show_col_types = FALSE
)


adult <- clean_rnaseq(
    adult_raw,
    BASEMEAN_THRESHOLD
)


message(
    "Adult primary genes: ",
    nrow(
        adult
    )
)


# ============================================================
# LOAD COLLECTRI
# ============================================================

message("")

message(
    "Loading frozen mouse CollecTRI..."
)


network <- load_collectri(
    COLLECTRI_FILE
)


message(
    "Usable signed CollecTRI edges: ",
    nrow(
        network
    )
)


message(
    "Usable regulators: ",
    length(
        unique(
            network$source
        )
    )
)


# ============================================================
# BUILD CONDITIONING SETS
# ============================================================

message("")

message(
    "Building GOALL conditioning sets..."
)


program_sets <- build_go_sets(
    adult$Gene_key
)


immune_keys <- program_sets$immune


cell_cycle_keys <- program_sets$cell_cycle


e2f_target_keys <- unique(
    network$target[
        network$source %in%
            E2F_FAMILY &
        network$target %in%
            adult$Gene_key
    ]
)


membership <- data.frame(

    Gene =
        adult$Gene,

    Gene_key =
        adult$Gene_key,

    Immune_GO =
        adult$Gene_key %in%
            immune_keys,

    Cell_cycle_GO =
        adult$Gene_key %in%
            cell_cycle_keys,

    E2F_family_target =
        adult$Gene_key %in%
            e2f_target_keys,

    stringsAsFactors = FALSE
)


message(
    "Immune GO genes in adult universe: ",
    sum(
        membership$Immune_GO
    )
)


message(
    "Cell-cycle GO genes in adult universe: ",
    sum(
        membership$Cell_cycle_GO
    )
)


message(
    "Measured E2F-family target union: ",
    sum(
        membership$E2F_family_target
    )
)


# ============================================================
# CONDITIONED EXPRESSION UNIVERSES
# ============================================================

full_df <- adult


no_immune_df <- adult[
    !membership$Immune_GO,
    ,
    drop = FALSE
]


no_cell_cycle_df <- adult[
    !membership$Cell_cycle_GO,
    ,
    drop = FALSE
]


no_e2f_targets_df <- adult[
    !membership$E2F_family_target,
    ,
    drop = FALSE
]


removal_summary <- data.frame(

    Condition =
        c(
            "FULL",
            "NO_IMMUNE_GO",
            "NO_CELL_CYCLE_GO",
            "NO_E2F_FAMILY_TARGETS"
        ),

    Genes_retained =
        c(
            nrow(
                full_df
            ),
            nrow(
                no_immune_df
            ),
            nrow(
                no_cell_cycle_df
            ),
            nrow(
                no_e2f_targets_df
            )
        ),

    Genes_removed =
        c(
            0,

            nrow(
                full_df
            ) -
            nrow(
                no_immune_df
            ),

            nrow(
                full_df
            ) -
            nrow(
                no_cell_cycle_df
            ),

            nrow(
                full_df
            ) -
            nrow(
                no_e2f_targets_df
            )
        ),

    stringsAsFactors = FALSE
)


# ============================================================
# RUN ULM + MLM
# ============================================================

message("")

message(
    "Running FULL / ULM..."
)


full_ulm <- run_activity(
    full_df,
    network,
    "FULL",
    "ULM"
)


message(
    "Running FULL / MLM..."
)


full_mlm <- run_activity(
    full_df,
    network,
    "FULL",
    "MLM"
)


message(
    "Running NO_IMMUNE_GO / ULM..."
)


no_immune_ulm <- run_activity(
    no_immune_df,
    network,
    "NO_IMMUNE_GO",
    "ULM"
)


message(
    "Running NO_IMMUNE_GO / MLM..."
)


no_immune_mlm <- run_activity(
    no_immune_df,
    network,
    "NO_IMMUNE_GO",
    "MLM"
)


message(
    "Running NO_CELL_CYCLE_GO / ULM..."
)


no_cell_ulm <- run_activity(
    no_cell_cycle_df,
    network,
    "NO_CELL_CYCLE_GO",
    "ULM"
)


message(
    "Running NO_CELL_CYCLE_GO / MLM..."
)


no_cell_mlm <- run_activity(
    no_cell_cycle_df,
    network,
    "NO_CELL_CYCLE_GO",
    "MLM"
)


message(
    "Running NO_E2F_FAMILY_TARGETS / ULM..."
)


no_e2f_ulm <- run_activity(
    no_e2f_targets_df,
    network,
    "NO_E2F_FAMILY_TARGETS",
    "ULM"
)


message(
    "Running NO_E2F_FAMILY_TARGETS / MLM..."
)


no_e2f_mlm <- run_activity(
    no_e2f_targets_df,
    network,
    "NO_E2F_FAMILY_TARGETS",
    "MLM"
)


results_long <- rbind(

    full_ulm,
    full_mlm,

    no_immune_ulm,
    no_immune_mlm,

    no_cell_ulm,
    no_cell_mlm,

    no_e2f_ulm,
    no_e2f_mlm
)


# ============================================================
# BUILD WIDE DECONVOLUTION TABLE
# ============================================================

ulm_full <- extract_condition(
    results_long,
    "FULL",
    "ULM",
    "ULM_FULL"
)


ulm_no_immune <- extract_condition(
    results_long,
    "NO_IMMUNE_GO",
    "ULM",
    "ULM_NO_IMMUNE"
)


ulm_no_cell <- extract_condition(
    results_long,
    "NO_CELL_CYCLE_GO",
    "ULM",
    "ULM_NO_CELL_CYCLE"
)


ulm_no_e2f <- extract_condition(
    results_long,
    "NO_E2F_FAMILY_TARGETS",
    "ULM",
    "ULM_NO_E2F_TARGETS"
)


mlm_full <- extract_condition(
    results_long,
    "FULL",
    "MLM",
    "MLM_FULL"
)


mlm_no_immune <- extract_condition(
    results_long,
    "NO_IMMUNE_GO",
    "MLM",
    "MLM_NO_IMMUNE"
)


mlm_no_cell <- extract_condition(
    results_long,
    "NO_CELL_CYCLE_GO",
    "MLM",
    "MLM_NO_CELL_CYCLE"
)


mlm_no_e2f <- extract_condition(
    results_long,
    "NO_E2F_FAMILY_TARGETS",
    "MLM",
    "MLM_NO_E2F_TARGETS"
)


wide <- merge(
    ulm_full,
    ulm_no_immune,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


wide <- merge(
    wide,
    ulm_no_cell,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


wide <- merge(
    wide,
    ulm_no_e2f,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


wide <- merge(
    wide,
    mlm_full,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


wide <- merge(
    wide,
    mlm_no_immune,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


wide <- merge(
    wide,
    mlm_no_cell,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


wide <- merge(
    wide,
    mlm_no_e2f,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


# ============================================================
# RETENTION METRICS
# ============================================================

wide$ULM_immune_retention <- safe_ratio(
    wide$ULM_NO_IMMUNE_Activity,
    wide$ULM_FULL_Activity
)


wide$MLM_immune_retention <- safe_ratio(
    wide$MLM_NO_IMMUNE_Activity,
    wide$MLM_FULL_Activity
)


wide$ULM_cell_cycle_retention <- safe_ratio(
    wide$ULM_NO_CELL_CYCLE_Activity,
    wide$ULM_FULL_Activity
)


wide$MLM_cell_cycle_retention <- safe_ratio(
    wide$MLM_NO_CELL_CYCLE_Activity,
    wide$MLM_FULL_Activity
)


wide$ULM_no_e2f_target_retention <- safe_ratio(
    wide$ULM_NO_E2F_TARGETS_Activity,
    wide$ULM_FULL_Activity
)


wide$MLM_no_e2f_target_retention <- safe_ratio(
    wide$MLM_NO_E2F_TARGETS_Activity,
    wide$MLM_FULL_Activity
)


wide$ULM_MLM_full_sign_concordant <- sign_same(
    wide$ULM_FULL_Activity,
    wide$MLM_FULL_Activity
)


wide$MLM_immune_sign_retained <- sign_same(
    wide$MLM_FULL_Activity,
    wide$MLM_NO_IMMUNE_Activity
)


wide$MLM_cell_cycle_sign_retained <- sign_same(
    wide$MLM_FULL_Activity,
    wide$MLM_NO_CELL_CYCLE_Activity
)


wide$MLM_no_e2f_target_sign_retained <- sign_same(
    wide$MLM_FULL_Activity,
    wide$MLM_NO_E2F_TARGETS_Activity
)


# ============================================================
# CANDIDATE AUDIT
# ============================================================

candidate_audit <- wide[
    toupper(
        wide$TF
    ) %in%
        CANDIDATE_TFS,
    ,
    drop = FALSE
]


candidate_audit <- candidate_audit[
    order(
        candidate_audit$MLM_FULL_FDR,
        -abs(
            candidate_audit$MLM_FULL_Activity
        ),
        na.last = TRUE
    ),
    ,
    drop = FALSE
]


rownames(
    candidate_audit
) <- NULL


# ============================================================
# REGULON OVERLAP AUDIT
# ============================================================

overlap_table <- make_overlap_table(
    network,
    adult$Gene_key
)


if (nrow(
    overlap_table
) > 0) {

    overlap_table <- overlap_table[
        order(
            -overlap_table$Jaccard,
            -overlap_table$Shared_targets,
            na.last = TRUE
        ),
        ,
        drop = FALSE
    ]


    rownames(
        overlap_table
    ) <- NULL
}


# ============================================================
# SAVE
# ============================================================

readr::write_csv(
    results_long,
    file.path(
        OUTPUT_DIR,
        "F10N_activity_all_conditions_long.csv"
    )
)


readr::write_csv(
    wide,
    file.path(
        OUTPUT_DIR,
        "F10N_activity_deconvolution_wide.csv"
    )
)


readr::write_csv(
    candidate_audit,
    file.path(
        OUTPUT_DIR,
        "F10N_candidate_TF_audit.csv"
    )
)


readr::write_csv(
    overlap_table,
    file.path(
        OUTPUT_DIR,
        "F10N_candidate_regulon_overlap.csv"
    )
)


readr::write_csv(
    removal_summary,
    file.path(
        OUTPUT_DIR,
        "F10N_conditioning_gene_counts.csv"
    )
)


readr::write_csv(
    membership,
    file.path(
        OUTPUT_DIR,
        "F10N_gene_conditioning_membership.csv"
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
    "F10N — REGULATORY DECONVOLUTION"
)

message(
    "============================================================"
)


message("")

message(
    "CONDITIONING SUMMARY"
)

message(
    "--------------------"
)


print(
    removal_summary,
    row.names = FALSE
)


console_columns <- c(

    "TF",

    "ULM_FULL_Activity",
    "ULM_FULL_FDR",

    "MLM_FULL_Activity",
    "MLM_FULL_FDR",

    "MLM_NO_IMMUNE_Activity",
    "MLM_NO_IMMUNE_FDR",
    "MLM_immune_retention",

    "MLM_NO_CELL_CYCLE_Activity",
    "MLM_NO_CELL_CYCLE_FDR",
    "MLM_cell_cycle_retention",

    "MLM_NO_E2F_TARGETS_Activity",
    "MLM_NO_E2F_TARGETS_FDR",
    "MLM_no_e2f_target_retention",

    "ULM_MLM_full_sign_concordant"
)


message("")

message(
    "PRE-SPECIFIED CANDIDATE TFs"
)

message(
    "---------------------------"
)


print(
    candidate_audit[
        ,
        console_columns,
        drop = FALSE
    ],
    row.names = FALSE
)


if (nrow(
    overlap_table
) > 0) {

    message("")

    message(
        "TOP CANDIDATE REGULON OVERLAPS"
    )

    message(
        "------------------------------"
    )


    n_show <- min(
        30,
        nrow(
            overlap_table
        )
    )


    print(
        overlap_table[
            seq_len(
                n_show
            ),
            ,
            drop = FALSE
        ],
        row.names = FALSE
    )
}


message("")

message(
    "FOCAL INTERPRETATION TARGETS"
)

message(
    "----------------------------"
)


focal_names <- c(

    "MYC",

    "E2F1",
    "E2F2",
    "E2F3",
    "E2F4",
    "TFDP1",

    "REL",
    "RELA",
    "NFKB1",

    "STAT2",
    "STAT3",

    "IRF1",
    "IRF2",
    "IRF9",

    "JUN",
    "CEBPD"
)


focal <- candidate_audit[
    candidate_audit$TF %in%
        focal_names,
    console_columns,
    drop = FALSE
]


print(
    focal,
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
        "1. MLM jointly fits regulators and therefore reduces,",
        "but does not eliminate, shared-target confounding."
    )
)


message(
    paste(
        "2. NO_IMMUNE_GO asks whether inflammatory TF activity",
        "survives removal of GO:0006955 immune-response genes."
    )
)


message(
    paste(
        "3. NO_CELL_CYCLE_GO and NO_E2F_FAMILY_TARGETS test",
        "whether MYC activity is separable from the",
        "E2F/proliferative program."
    )
)


message(
    paste(
        "4. Persistence supports regulatory separability;",
        "collapse supports shared-program explanations.",
        "Neither establishes causality."
    )
)