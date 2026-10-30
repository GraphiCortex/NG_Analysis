#!/usr/bin/env Rscript

# ============================================================
# F11A_empirical_state_variables.R
# ============================================================
#
# PURPOSE
# -------
#
# Operationalize the two state variables that will enter the
# F11-v2 dynamical model:
#
#   E = E2F / proliferative program
#   S = immune / stress program
#
#
# PRIMARY PROGRAM DEFINITIONS
# ---------------------------
#
# E_CORE:
#
#   GO:0007049 cell-cycle genes
#       INTERSECT
#   signed CollecTRI targets of:
#
#       E2F1
#       E2F2
#       E2F3
#       E2F4
#       TFDP1
#
#
# S_IMMUNE:
#
#   GO:0006955 immune-response genes
#
#
# SENSITIVITY DEFINITIONS
# -----------------------
#
# E_BROAD:
#
#   all GO:0007049 cell-cycle genes
#
#
# S_INNATE:
#
#   GO:0045087 innate-immune-response genes
#
#
# EVIDENCE
# --------
#
# Adult Rb-TKO:
#   DESeq2 Wald statistic
#
# E2F1/3-DKO:
#   DESeq2 Wald statistic
#
# CAG / MAP2:
#   signed normal score reconstructed from limma P-value
#
#
# All four contexts are restricted to a common gene universe,
# then independently transformed using the same sign-preserving
# rank-Gaussian transform.
#
#
# NULL MODEL
# ----------
#
# Each program is compared with random gene sets matched on:
#
#   adult RNA-seq baseline-expression quintile
#   E2F-DKO baseline-expression quintile
#
# The matched null therefore controls gene-set size and major
# RNA-seq detectability structure.
#
#
# IMPORTANT
# ---------
#
# This script does NOT infer dynamics.
#
# It constructs and tests empirical state variables that may
# subsequently enter a dynamical model.
#
# ============================================================


# ============================================================
# PACKAGES
# ============================================================

required <- c(
    "readr",
    "readxl",
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


E2F_FILE <- file.path(
    "Data",
    "Raw",
    "microarray",
    "Fong2022",
    "Fong 2022 E2F1.3 microarray data.xlsx"
)


E2F_SHEET <- "GSE190764_FDSD_FC.IHW"


OSH_FILE <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F8A_Oshikawa",
    "oshikawa_gene_response_matrix.csv"
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
    "F11A_EmpiricalStateVariables"
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

N_NULL <- 5000

N_BOOT <- 3000

RANDOM_SEED <- 20261028


CELL_CYCLE_GO <- "GO:0007049"

IMMUNE_GO <- "GO:0006955"

INNATE_GO <- "GO:0045087"


E2F_FAMILY <- c(
    "E2F1",
    "E2F2",
    "E2F3",
    "E2F4",
    "TFDP1"
)


set.seed(
    RANDOM_SEED
)


# ============================================================
# HELPERS
# ============================================================

signed_z_from_p <- function(
    logfc,
    pvalue
) {

    logfc <- as.numeric(
        logfc
    )

    pvalue <- as.numeric(
        pvalue
    )


    result <- rep(
        NA_real_,
        length(logfc)
    )


    valid <- (
        is.finite(logfc) &
        is.finite(pvalue)
    )


    if (!any(valid)) {

        return(
            result
        )
    }


    p_clean <- pvalue[
        valid
    ]


    p_clean <- pmax(
        p_clean,
        1e-300
    )


    p_clean <- pmin(
        p_clean,
        1
    )


    z_abs <- qnorm(
        p_clean / 2,
        lower.tail = FALSE
    )


    result[
        valid
    ] <- (
        sign(
            logfc[
                valid
            ]
        ) *
        z_abs
    )


    zero_fc <- (
        valid &
        logfc == 0
    )


    result[
        zero_fc
    ] <- 0


    result
}


signed_rank_gaussian <- function(x) {

    x <- as.numeric(
        x
    )


    result <- rep(
        NA_real_,
        length(x)
    )


    valid <- is.finite(
        x
    )


    if (sum(valid) < 2) {

        return(
            result
        )
    }


    magnitude <- abs(
        x[
            valid
        ]
    )


    ranks <- rank(
        magnitude,
        ties.method = "average"
    )


    percentile <- (
        ranks - 0.5
    ) /
        length(
            ranks
        )


    transformed <- qnorm(
        0.5 +
            0.5 *
            percentile
    )


    result[
        valid
    ] <- (
        sign(
            x[
                valid
            ]
        ) *
        transformed
    )


    result[
        valid &
        x == 0
    ] <- 0


    result
}


make_quintile <- function(x) {

    x <- as.numeric(
        x
    )


    ranks <- rank(
        x,
        ties.method = "average",
        na.last = "keep"
    )


    percentile <- ranks /
        sum(
            is.finite(
                ranks
            )
        )


    result <- ceiling(
        percentile * 5
    )


    result[
        result < 1
    ] <- 1


    result[
        result > 5
    ] <- 5


    result
}


# ============================================================
# CLEAN RNA-SEQ TABLE
# ============================================================

clean_rnaseq <- function(
    df,
    threshold
) {

    required <- c(
        "name",
        "baseMean",
        "stat"
    )


    absent <- setdiff(
        required,
        names(df)
    )


    if (length(absent) > 0) {

        stop(
            paste(
                "RNA-seq input missing:",
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


    stat <- suppressWarnings(
        as.numeric(
            df[
                ,
                "stat",
                drop = TRUE
            ]
        )
    )


    valid <- (
        !is.na(gene) &
        gene != "" &
        gene != "NA" &
        is.finite(base_mean)
    )


    out <- data.frame(

        Gene =
            gene[
                valid
            ],

        Gene_key =
            gene_key[
                valid
            ],

        baseMean =
            base_mean[
                valid
            ],

        Evidence =
            stat[
                valid
            ],

        stringsAsFactors = FALSE
    )


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
            out$Evidence
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
# LOAD SIGNED COLLECTRI
# ============================================================

load_collectri <- function(path) {

    raw <- readr::read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE
    )


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


    activation_only <- (
        is.finite(
            stimulation
        ) &
        is.finite(
            inhibition
        ) &
        stimulation == 1 &
        inhibition == 0
    )


    inhibition_only <- (
        is.finite(
            stimulation
        ) &
        is.finite(
            inhibition
        ) &
        stimulation == 0 &
        inhibition == 1
    )


    valid <- (
        activation_only |
        inhibition_only
    )


    network <- data.frame(

        source =
            toupper(
                trimws(
                    source_symbol[
                        valid
                    ]
                )
            ),

        target =
            toupper(
                trimws(
                    target_symbol[
                        valid
                    ]
                )
            ),

        mor =
            ifelse(
                activation_only[
                    valid
                ],
                1,
                -1
            ),

        stringsAsFactors = FALSE
    )


    network <- unique(
        network
    )


    pair <- paste(
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
        pair
    )


    keep_pair <- vapply(

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
        keep_pair[
            keep_pair
        ]
    )


    network <- network[
        pair %in%
            valid_pairs,
        ,
        drop = FALSE
    ]


    pair <- paste(
        network$source,
        network$target,
        sep = "|||"
    )


    network <- network[
        !duplicated(
            pair
        ),
        ,
        drop = FALSE
    ]


    rownames(
        network
    ) <- NULL


    network
}


# ============================================================
# GO SETS
# ============================================================

build_go_sets <- function(
    universe
) {

    official <- AnnotationDbi::keys(
        org.Mm.eg.db,
        keytype = "SYMBOL"
    )


    map <- data.frame(

        SYMBOL =
            official,

        Gene_key =
            toupper(
                official
            ),

        stringsAsFactors = FALSE
    )


    counts <- table(
        map$Gene_key
    )


    unique_keys <- names(
        counts[
            counts == 1
        ]
    )


    map <- map[
        map$Gene_key %in%
            unique_keys,
        ,
        drop = FALSE
    ]


    symbol <- map$SYMBOL[
        match(
            universe,
            map$Gene_key
        )
    ]


    symbol <- unique(
        symbol[
            !is.na(
                symbol
            )
        ]
    )


    annotation <- suppressMessages(

        AnnotationDbi::select(

            org.Mm.eg.db,

            keys =
                symbol,

            keytype =
                "SYMBOL",

            columns =
                c(
                    "GOALL",
                    "ONTOLOGYALL"
                )
        )
    )


    get_set <- function(go_id) {

        genes <- annotation$SYMBOL[
            annotation$GOALL ==
                go_id &
            annotation$ONTOLOGYALL ==
                "BP"
        ]


        unique(
            toupper(
                as.character(
                    genes
                )
            )
        )
    }


    list(

        cell_cycle =
            get_set(
                CELL_CYCLE_GO
            ),

        immune =
            get_set(
                IMMUNE_GO
            ),

        innate =
            get_set(
                INNATE_GO
            )
    )
}


# ============================================================
# PROGRAM SCORE
# ============================================================

score_program <- function(
    program_name,
    genes,
    state_matrix,
    strata,
    n_null,
    n_boot
) {

    genes <- intersect(
        genes,
        rownames(
            state_matrix
        )
    )


    if (length(genes) < 5) {

        stop(
            paste(
                program_name,
                "contains fewer than 5 common genes."
            )
        )
    }


    program_matrix <- state_matrix[
        genes,
        ,
        drop = FALSE
    ]


    observed <- colMeans(
        program_matrix,
        na.rm = TRUE
    )


    median_score <- apply(
        program_matrix,
        2,
        median,
        na.rm = TRUE
    )


    positive_fraction <- colMeans(
        program_matrix > 0,
        na.rm = TRUE
    )


    # --------------------------------------------------------
    # BOOTSTRAP CI
    # --------------------------------------------------------

    bootstrap <- matrix(
        NA_real_,
        nrow = n_boot,
        ncol = ncol(
            program_matrix
        )
    )


    for (i in seq_len(
        n_boot
    )) {

        sample_index <- sample(
            seq_len(
                nrow(
                    program_matrix
                )
            ),
            size =
                nrow(
                    program_matrix
                ),
            replace = TRUE
        )


        bootstrap[
            i,
        ] <- colMeans(
            program_matrix[
                sample_index,
                ,
                drop = FALSE
            ],
            na.rm = TRUE
        )
    }


    ci_low <- apply(
        bootstrap,
        2,
        quantile,
        probs = 0.025,
        na.rm = TRUE
    )


    ci_high <- apply(
        bootstrap,
        2,
        quantile,
        probs = 0.975,
        na.rm = TRUE
    )


    # --------------------------------------------------------
    # EXPRESSION-MATCHED RANDOM-SET NULL
    # --------------------------------------------------------

    program_strata <- strata[
        genes
    ]


    stratum_counts <- table(
        program_strata
    )


    universe <- rownames(
        state_matrix
    )


    null_scores <- matrix(
        NA_real_,
        nrow = n_null,
        ncol = ncol(
            state_matrix
        )
    )


    for (i in seq_len(
        n_null
    )) {

        sampled <- character(0)


        for (stratum in names(
            stratum_counts
        )) {

            required_n <- as.integer(
                stratum_counts[
                    stratum
                ]
            )


            pool <- universe[
                strata[
                    universe
                ] ==
                    stratum &
                !universe %in%
                    genes
            ]


            if (length(pool) == 0) {

                pool <- universe[
                    strata[
                        universe
                    ] ==
                        stratum
                ]
            }


            replace_flag <- (
                length(pool) <
                    required_n
            )


            sampled <- c(
                sampled,
                sample(
                    pool,
                    size =
                        required_n,
                    replace =
                        replace_flag
                )
            )
        }


        null_scores[
            i,
        ] <- colMeans(
            state_matrix[
                sampled,
                ,
                drop = FALSE
            ],
            na.rm = TRUE
        )
    }


    null_mean <- colMeans(
        null_scores,
        na.rm = TRUE
    )


    null_sd <- apply(
        null_scores,
        2,
        sd,
        na.rm = TRUE
    )


    null_z <- (
        observed -
        null_mean
    ) /
        null_sd


    p_two_sided <- rep(
        NA_real_,
        length(
            observed
        )
    )


    p_positive <- rep(
        NA_real_,
        length(
            observed
        )
    )


    p_negative <- rep(
        NA_real_,
        length(
            observed
        )
    )


    for (j in seq_along(
        observed
    )) {

        p_two_sided[
            j
        ] <- (
            1 +
            sum(
                abs(
                    null_scores[
                        ,
                        j
                    ] -
                    null_mean[
                        j
                    ]
                ) >=
                    abs(
                        observed[
                            j
                        ] -
                        null_mean[
                            j
                        ]
                    )
            )
        ) / (
            n_null +
            1
        )


        p_positive[
            j
        ] <- (
            1 +
            sum(
                null_scores[
                    ,
                    j
                ] >=
                    observed[
                        j
                    ]
            )
        ) / (
            n_null +
            1
        )


        p_negative[
            j
        ] <- (
            1 +
            sum(
                null_scores[
                    ,
                    j
                ] <=
                    observed[
                        j
                    ]
            )
        ) / (
            n_null +
            1
        )
    }


    summary <- data.frame(

        Program =
            program_name,

        Context =
            colnames(
                state_matrix
            ),

        N_genes =
            length(
                genes
            ),

        Mean_score =
            as.numeric(
                observed
            ),

        Median_score =
            as.numeric(
                median_score
            ),

        Fraction_positive =
            as.numeric(
                positive_fraction
            ),

        Bootstrap_CI_low =
            as.numeric(
                ci_low
            ),

        Bootstrap_CI_high =
            as.numeric(
                ci_high
            ),

        Null_mean =
            as.numeric(
                null_mean
            ),

        Null_SD =
            as.numeric(
                null_sd
            ),

        Null_Z =
            as.numeric(
                null_z
            ),

        Null_P_two_sided =
            as.numeric(
                p_two_sided
            ),

        Null_P_positive =
            as.numeric(
                p_positive
            ),

        Null_P_negative =
            as.numeric(
                p_negative
            ),

        stringsAsFactors = FALSE
    )


    list(

        summary =
            summary,

        genes =
            genes,

        null =
            null_scores,

        bootstrap =
            bootstrap
    )
}


# ============================================================
# LOAD ADULT RNA-SEQ
# ============================================================

message(
    "Loading adult Rb-TKO..."
)


adult_raw <- readr::read_csv(
    ADULT_FILE,
    show_col_types = FALSE
)


adult <- clean_rnaseq(
    adult_raw,
    BASEMEAN_THRESHOLD
)


# ============================================================
# LOAD E2F1/3-DKO
# ============================================================

message(
    "Loading E2F1/3-DKO..."
)


e2f_raw <- readxl::read_excel(
    E2F_FILE,
    sheet = E2F_SHEET
)


e2f <- clean_rnaseq(
    e2f_raw,
    BASEMEAN_THRESHOLD
)


# ============================================================
# LOAD CAG / MAP2
# ============================================================

message(
    "Loading CAG/MAP2..."
)


osh <- readr::read_csv(
    OSH_FILE,
    show_col_types = FALSE
)


gene <- toupper(
    trimws(
        as.character(
            osh[
                ,
                "Gene",
                drop = TRUE
            ]
        )
    )
)


cag_evidence <- signed_z_from_p(

    osh[
        ,
        "logFC_CAG",
        drop = TRUE
    ],

    osh[
        ,
        "P.Value_CAG",
        drop = TRUE
    ]
)


map2_evidence <- signed_z_from_p(

    osh[
        ,
        "logFC_MAP2",
        drop = TRUE
    ],

    osh[
        ,
        "P.Value_MAP2",
        drop = TRUE
    ]
)


cag <- data.frame(

    Gene_key =
        gene,

    Evidence =
        cag_evidence,

    stringsAsFactors = FALSE
)


map2 <- data.frame(

    Gene_key =
        gene,

    Evidence =
        map2_evidence,

    stringsAsFactors = FALSE
)


cag <- cag[
    !is.na(
        cag$Gene_key
    ) &
    cag$Gene_key != "" &
    is.finite(
        cag$Evidence
    ),
    ,
    drop = FALSE
]


map2 <- map2[
    !is.na(
        map2$Gene_key
    ) &
    map2$Gene_key != "" &
    is.finite(
        map2$Evidence
    ),
    ,
    drop = FALSE
]


cag <- cag[
    !duplicated(
        cag$Gene_key
    ),
    ,
    drop = FALSE
]


map2 <- map2[
    !duplicated(
        map2$Gene_key
    ),
    ,
    drop = FALSE
]


# ============================================================
# COMMON UNIVERSE
# ============================================================

common <- Reduce(

    intersect,

    list(

        adult$Gene_key,

        e2f$Gene_key,

        cag$Gene_key,

        map2$Gene_key
    )
)


common <- sort(
    unique(
        common
    )
)


message(
    "Common four-context universe: ",
    length(
        common
    )
)


# ============================================================
# ALIGN EVIDENCE
# ============================================================

adult_index <- match(
    common,
    adult$Gene_key
)


e2f_index <- match(
    common,
    e2f$Gene_key
)


cag_index <- match(
    common,
    cag$Gene_key
)


map2_index <- match(
    common,
    map2$Gene_key
)


adult_vector <- adult$Evidence[
    adult_index
]


e2f_vector <- e2f$Evidence[
    e2f_index
]


cag_vector <- cag$Evidence[
    cag_index
]


map2_vector <- map2$Evidence[
    map2_index
]


state_matrix <- cbind(

    Adult_RbTKO =
        signed_rank_gaussian(
            adult_vector
        ),

    E2F13_DKO =
        signed_rank_gaussian(
            e2f_vector
        ),

    CAG_RbLoss =
        signed_rank_gaussian(
            cag_vector
        ),

    MAP2_RbLoss =
        signed_rank_gaussian(
            map2_vector
        )
)


rownames(
    state_matrix
) <- common


# ============================================================
# DETECTABILITY STRATA
# ============================================================

adult_base <- adult$baseMean[
    adult_index
]


e2f_base <- e2f$baseMean[
    e2f_index
]


adult_quintile <- make_quintile(
    adult_base
)


e2f_quintile <- make_quintile(
    e2f_base
)


strata <- paste0(
    "A",
    adult_quintile,
    "_E",
    e2f_quintile
)


names(
    strata
) <- common


# ============================================================
# BUILD PRIOR-DEFINED PROGRAMS
# ============================================================

message(
    "Building GO program definitions..."
)


go_sets <- build_go_sets(
    common
)


message(
    "Loading corrected signed CollecTRI..."
)


network <- load_collectri(
    COLLECTRI_FILE
)


e2f_targets <- unique(
    network$target[
        network$source %in%
            E2F_FAMILY
    ]
)


E_CORE <- Reduce(

    intersect,

    list(
        common,
        go_sets$cell_cycle,
        e2f_targets
    )
)


E_BROAD <- intersect(
    common,
    go_sets$cell_cycle
)


S_IMMUNE <- intersect(
    common,
    go_sets$immune
)


S_INNATE <- intersect(
    common,
    go_sets$innate
)


message(
    "E_CORE genes: ",
    length(
        E_CORE
    )
)


message(
    "E_BROAD genes: ",
    length(
        E_BROAD
    )
)


message(
    "S_IMMUNE genes: ",
    length(
        S_IMMUNE
    )
)


message(
    "S_INNATE genes: ",
    length(
        S_INNATE
    )
)


# ============================================================
# SCORE PROGRAMS
# ============================================================

message("")

message(
    "Scoring E_CORE..."
)


score_E_core <- score_program(

    "E_CORE",

    E_CORE,

    state_matrix,

    strata,

    N_NULL,

    N_BOOT
)


message(
    "Scoring E_BROAD..."
)


score_E_broad <- score_program(

    "E_BROAD",

    E_BROAD,

    state_matrix,

    strata,

    N_NULL,

    N_BOOT
)


message(
    "Scoring S_IMMUNE..."
)


score_S_immune <- score_program(

    "S_IMMUNE",

    S_IMMUNE,

    state_matrix,

    strata,

    N_NULL,

    N_BOOT
)


message(
    "Scoring S_INNATE..."
)


score_S_innate <- score_program(

    "S_INNATE",

    S_INNATE,

    state_matrix,

    strata,

    N_NULL,

    N_BOOT
)


summary_table <- rbind(

    score_E_core$summary,

    score_E_broad$summary,

    score_S_immune$summary,

    score_S_innate$summary
)


# ============================================================
# CONTEXT-DIFFERENCE TABLE
# ============================================================
#
# Because every state score uses the same genes and the same
# common four-context universe, paired context differences can
# be calculated gene-by-gene.
# ============================================================

program_list <- list(

    E_CORE =
        E_CORE,

    E_BROAD =
        E_BROAD,

    S_IMMUNE =
        S_IMMUNE,

    S_INNATE =
        S_INNATE
)


contrast_definitions <- list(

    Adult_minus_CAG =
        c(
            "Adult_RbTKO",
            "CAG_RbLoss"
        ),

    Adult_minus_MAP2 =
        c(
            "Adult_RbTKO",
            "MAP2_RbLoss"
        ),

    Adult_minus_E2FDKO =
        c(
            "Adult_RbTKO",
            "E2F13_DKO"
        )
)


contrast_table <- data.frame()


for (program_name in names(
    program_list
)) {

    genes <- getElement(
        program_list,
        program_name
    )


    genes <- intersect(
        genes,
        rownames(
            state_matrix
        )
    )


    for (contrast_name in names(
        contrast_definitions
    )) {

        contexts <- getElement(
            contrast_definitions,
            contrast_name
        )


        context_1 <- contexts[
            1
        ]


        context_2 <- contexts[
            2
        ]


        x <- state_matrix[
            genes,
            context_1,
            drop = TRUE
        ]


        y <- state_matrix[
            genes,
            context_2,
            drop = TRUE
        ]


        difference <- (
            x -
            y
        )


        boot_mean <- rep(
            NA_real_,
            N_BOOT
        )


        for (i in seq_len(
            N_BOOT
        )) {

            index <- sample(
                seq_along(
                    difference
                ),
                size =
                    length(
                        difference
                    ),
                replace = TRUE
            )


            boot_mean[
                i
            ] <- mean(
                difference[
                    index
                ],
                na.rm = TRUE
            )
        }


        one_row <- data.frame(

            Program =
                program_name,

            Contrast =
                contrast_name,

            N_genes =
                length(
                    difference
                ),

            Mean_difference =
                mean(
                    difference,
                    na.rm = TRUE
                ),

            Median_difference =
                median(
                    difference,
                    na.rm = TRUE
                ),

            Bootstrap_CI_low =
                as.numeric(
                    quantile(
                        boot_mean,
                        0.025,
                        na.rm = TRUE
                    )
                ),

            Bootstrap_CI_high =
                as.numeric(
                    quantile(
                        boot_mean,
                        0.975,
                        na.rm = TRUE
                    )
                ),

            Fraction_context1_greater =
                mean(
                    difference > 0,
                    na.rm = TRUE
                ),

            stringsAsFactors = FALSE
        )


        contrast_table <- rbind(
            contrast_table,
            one_row
        )
    }
}
# ============================================================
# MODEL-CONSTRAINT CHECKS
# ============================================================

get_score <- function(
    program,
    context
) {

    value <- summary_table$Mean_score[
        summary_table$Program ==
            program &
        summary_table$Context ==
            context
    ]


    if (length(
        value
    ) != 1) {

        return(
            NA_real_
        )
    }


    value
}


E_adult <- get_score(
    "E_CORE",
    "Adult_RbTKO"
)


E_e2f <- get_score(
    "E_CORE",
    "E2F13_DKO"
)


E_cag <- get_score(
    "E_CORE",
    "CAG_RbLoss"
)


E_map2 <- get_score(
    "E_CORE",
    "MAP2_RbLoss"
)


S_adult <- get_score(
    "S_IMMUNE",
    "Adult_RbTKO"
)


S_e2f <- get_score(
    "S_IMMUNE",
    "E2F13_DKO"
)


S_cag <- get_score(
    "S_IMMUNE",
    "CAG_RbLoss"
)


S_map2 <- get_score(
    "S_IMMUNE",
    "MAP2_RbLoss"
)


constraint_table <- data.frame(

    Constraint =
        c(
            "E positive after adult Rb loss",
            "E negative after E2F1/3 loss",
            "E positive after CAG Rb loss",
            "E positive after MAP2 Rb loss",
            "S positive after adult Rb loss",
            "S adult greater than CAG",
            "S adult greater than MAP2"
        ),

    Passed =
        c(
            E_adult > 0,
            E_e2f < 0,
            E_cag > 0,
            E_map2 > 0,
            S_adult > 0,
            S_adult > S_cag,
            S_adult > S_map2
        ),

    stringsAsFactors = FALSE
)


# ============================================================
# PROGRAM MEMBERSHIP
# ============================================================

membership <- data.frame(

    Gene_key =
        common,

    E_CORE =
        common %in%
            E_CORE,

    E_BROAD =
        common %in%
            E_BROAD,

    S_IMMUNE =
        common %in%
            S_IMMUNE,

    S_INNATE =
        common %in%
            S_INNATE,

    Adult_baseMean =
        adult_base,

    E2F_baseMean =
        e2f_base,

    Detectability_stratum =
        strata[
            common
        ],

    stringsAsFactors = FALSE
)


# ============================================================
# SAVE OUTPUTS
# ============================================================

readr::write_csv(
    summary_table,
    file.path(
        OUTPUT_DIR,
        "F11A_state_variable_scores.csv"
    )
)


readr::write_csv(
    contrast_table,
    file.path(
        OUTPUT_DIR,
        "F11A_context_differences.csv"
    )
)


readr::write_csv(
    constraint_table,
    file.path(
        OUTPUT_DIR,
        "F11A_model_constraints.csv"
    )
)


readr::write_csv(
    membership,
    file.path(
        OUTPUT_DIR,
        "F11A_program_membership.csv"
    )
)


state_matrix_output <- data.frame(

    Gene_key =
        rownames(
            state_matrix
        ),

    state_matrix,

    stringsAsFactors = FALSE
)


readr::write_csv(
    state_matrix_output,
    file.path(
        OUTPUT_DIR,
        "F11A_gene_state_evidence_matrix.csv"
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
    "F11A — EMPIRICAL STATE VARIABLES"
)

message(
    "============================================================"
)


message("")

message(
    "PROGRAM DEFINITIONS"
)

message(
    "-------------------"
)


message(
    "E_CORE: ",
    length(
        E_CORE
    ),
    " genes"
)


message(
    "E_BROAD: ",
    length(
        E_BROAD
    ),
    " genes"
)


message(
    "S_IMMUNE: ",
    length(
        S_IMMUNE
    ),
    " genes"
)


message(
    "S_INNATE: ",
    length(
        S_INNATE
    ),
    " genes"
)


message("")

message(
    "STATE SCORES"
)

message(
    "------------"
)


print(
    summary_table,
    row.names = FALSE
)


message("")

message(
    "PAIRED CONTEXT DIFFERENCES"
)

message(
    "--------------------------"
)


print(
    contrast_table,
    row.names = FALSE
)


message("")

message(
    "F11-v2 MODEL CONSTRAINTS"
)

message(
    "-----------------------"
)


print(
    constraint_table,
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
        "1. E and S are empirical program scores, not",
        "direct biochemical concentrations."
    )
)


message(
    paste(
        "2. Program definitions are prior-anchored rather",
        "than selected from adult differential expression."
    )
)


message(
    paste(
        "3. Expression-matched nulls test whether the",
        "program score is unusual relative to similarly",
        "detectable genes."
    )
)


message(
    paste(
        "4. These data constrain steady perturbational",
        "states. They do not identify kinetic parameters,",
        "feedback strengths, or time constants."
    )
)


message(
    paste(
        "5. Dynamics will only be introduced after these",
        "state constraints are established."
    )
)