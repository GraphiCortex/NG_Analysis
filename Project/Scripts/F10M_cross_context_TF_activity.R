#!/usr/bin/env Rscript

# ============================================================
# F10M_cross_context_TF_activity.R
# ============================================================
#
# PURPOSE
# -------
# Signed TF-activity comparison across:
#
#   1. Adult Rb-family TKO
#   2. Adult E2F1/3 DKO
#   3. Embryonic cortical CAG Rb loss
#   4. Embryonic cortical MAP2 Rb loss
#
#
# PRIMARY EVIDENCE
# ----------------
#
# Adult Rb-TKO:
#   DESeq2 Wald statistic
#
# E2F1/3-DKO:
#   DESeq2 Wald statistic
#
# CAG / MAP2:
#   signed normal score reconstructed from limma P-value:
#
#     sign(logFC) * Phi^{-1}(1 - p/2)
#
#
# CROSS-PLATFORM SENSITIVITY
# --------------------------
#
# Each context is also analyzed after a sign-preserving
# rank-Gaussian transform.
#
# Cross-platform interpretation therefore emphasizes:
#
#   - activity direction
#   - within-context percentile
#   - FDR
#   - native/rank-normalized agreement
#
# Raw ULM magnitudes are NOT interpreted as directly
# comparable biological quantities across platforms.
#
#
# IMPORTANT
# ---------
#
# CollecTRI provides prior signed TF -> target relationships.
# ULM activity is an inference from target behavior.
#
# It does NOT directly measure:
#
#   - TF protein abundance
#   - phosphorylation
#   - nuclear localization
#   - DNA occupancy
#   - causal necessity
#
# ============================================================


# ============================================================
# PACKAGES
# ============================================================

required <- c(
    "decoupleR",
    "readr",
    "readxl"
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

F10L_FILE <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F10L_AdultSignedTFActivity",
    "F10L_integrated_TF_activity.csv"
)

OUTPUT_DIR <- file.path(
    "Data",
    "Processed",
    "Dynamics",
    "F10M_CrossContextTFActivity"
)

dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# PARAMETERS
# ============================================================

BASEMEAN_PRIMARY <- 10

BASEMEAN_SENSITIVITY <- 1

MIN_REGULON_SIZE <- 5


CANDIDATE_TFS <- c(
    "E2F1",
    "E2F2",
    "E2F3",
    "E2F4",
    "TFDP1",
    "FOXM1",

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
    "IRF7",
    "IRF8",
    "IRF9",

    "MYC",

    "CEBPD",
    "JUN",
    "FOS",
    "DDIT3",
    "ATF4",

    "TP53",

    "NR4A1",
    "NR4A2",
    "NR4A3"
)


# ============================================================
# HELPER: FDR
# ============================================================

safe_fdr <- function(p) {

    p.adjust(
        as.numeric(p),
        method = "BH"
    )
}


# ============================================================
# HELPER: SIGNED Z FROM TWO-SIDED P-VALUE
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


# ============================================================
# HELPER: SIGN-PRESERVING RANK-GAUSSIAN TRANSFORM
# ============================================================

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


    n <- length(
        ranks
    )


    percentile <- (
        ranks - 0.5
    ) / n


    gaussian_magnitude <- qnorm(
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
        gaussian_magnitude
    )


    zero_values <- (
        valid &
        x == 0
    )


    result[
        zero_values
    ] <- 0


    result
}


# ============================================================
# HELPER: SIGNED ACTIVITY PERCENTILE
# ============================================================

percentile_signed <- function(x) {

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


    if (!any(valid)) {

        return(
            result
        )
    }


    result[
        valid
    ] <- (
        rank(
            x[
                valid
            ],
            ties.method = "average"
        ) /
        sum(valid)
    )


    result
}


# ============================================================
# HELPER: ABSOLUTE ACTIVITY PERCENTILE
# ============================================================

percentile_absolute <- function(x) {

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


    if (!any(valid)) {

        return(
            result
        )
    }


    result[
        valid
    ] <- (
        rank(
            abs(
                x[
                    valid
                ]
            ),
            ties.method = "average"
        ) /
        sum(valid)
    )


    result
}


# ============================================================
# HELPER: SIGN CONCORDANCE
# ============================================================

sign_concordant <- function(
    x,
    y
) {

    result <- rep(
        NA,
        length(x)
    )


    valid <- (
        is.finite(x) &
        is.finite(y)
    )


    result[
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


    result
}


# ============================================================
# HELPER: SAFE MAXIMUM OF TWO VALUES
# ============================================================

safe_max2 <- function(
    x,
    y
) {

    result <- rep(
        NA_real_,
        length(x)
    )


    x_valid <- is.finite(
        x
    )

    y_valid <- is.finite(
        y
    )


    both <- (
        x_valid &
        y_valid
    )


    x_only <- (
        x_valid &
        !y_valid
    )


    y_only <- (
        !x_valid &
        y_valid
    )


    result[
        both
    ] <- pmax(
        x[
            both
        ],
        y[
            both
        ]
    )


    result[
        x_only
    ] <- x[
        x_only
    ]


    result[
        y_only
    ] <- y[
        y_only
    ]


    result
}


# ============================================================
# HELPER: SAFE MEAN OF TWO VALUES
# ============================================================

safe_mean2 <- function(
    x,
    y
) {

    result <- rep(
        NA_real_,
        length(x)
    )


    x_valid <- is.finite(
        x
    )

    y_valid <- is.finite(
        y
    )


    both <- (
        x_valid &
        y_valid
    )


    x_only <- (
        x_valid &
        !y_valid
    )


    y_only <- (
        !x_valid &
        y_valid
    )


    result[
        both
    ] <- (
        x[
            both
        ] +
        y[
            both
        ]
    ) / 2


    result[
        x_only
    ] <- x[
        x_only
    ]


    result[
        y_only
    ] <- y[
        y_only
    ]


    result
}


# ============================================================
# HELPER: CLEAN RNA-SEQ TABLE
# ============================================================
#
# F10L selected one row per symbol by highest baseMean before
# applying the expression/statistic filter.
#
# We reproduce that logic here so that the adult F10M activity
# can be compared directly against F10L.
# ============================================================

clean_rnaseq <- function(
    df,
    threshold
) {

    required_columns <- c(
        "name",
        "baseMean",
        "stat"
    )


    missing_columns <- setdiff(
        required_columns,
        names(df)
    )


    if (length(
        missing_columns
    ) > 0) {

        stop(
            paste(
                "RNA-seq table missing:",
                paste(
                    missing_columns,
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


    valid_gene <- (
        !is.na(gene) &
        gene != "" &
        gene != "NA" &
        is.finite(base_mean)
    )


    out <- data.frame(
        Gene =
            gene[
                valid_gene
            ],

        Gene_key =
            gene_key[
                valid_gene
            ],

        baseMean =
            base_mean[
                valid_gene
            ],

        stat =
            stat_value[
                valid_gene
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
# HELPER: LOAD FROZEN MOUSE COLLECTRI NETWORK
# ============================================================
#
# This reproduces the network construction used in F10L.
#
# F10L used:
#
#   is_stimulation == 1 -> mor = +1
#   is_stimulation == 0 -> mor = -1
#
# Individual TFs remain individual.
#
# JUN/FOS complexes -> AP1
# REL/NFKB complexes -> NFKB
#
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
# HELPER: RUN ULM FOR ONE CONTEXT
# ============================================================

run_ulm_context <- function(
    df,
    evidence_column,
    network,
    context_name,
    evidence_type
) {

    if (!"Gene_key" %in%
        names(df)) {

        stop(
            paste(
                "Gene_key missing in",
                context_name
            )
        )
    }


    if (!evidence_column %in%
        names(df)) {

        stop(
            paste(
                "Missing evidence column",
                evidence_column,
                "in",
                context_name
            )
        )
    }


    genes <- as.character(
        df[
            ,
            "Gene_key",
            drop = TRUE
        ]
    )


    evidence <- suppressWarnings(
        as.numeric(
            df[
                ,
                evidence_column,
                drop = TRUE
            ]
        )
    )


    valid <- (
        !is.na(
            genes
        ) &
        genes != "" &
        is.finite(
            evidence
        )
    )


    genes <- genes[
        valid
    ]


    evidence <- evidence[
        valid
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


    eligible_regulators <- names(
        regulon_sizes[
            regulon_sizes >=
                MIN_REGULON_SIZE
        ]
    )


    context_network <- context_network[
        context_network$source %in%
            eligible_regulators,
        ,
        drop = FALSE
    ]


    if (nrow(
        context_network
    ) == 0) {

        stop(
            paste(
                "No eligible regulons in",
                context_name
            )
        )
    }


    evidence_matrix <- matrix(
        evidence,
        ncol = 1
    )


    rownames(
        evidence_matrix
    ) <- genes


    colnames(
        evidence_matrix
    ) <- context_name


    fit <- decoupleR::run_ulm(
        mat = evidence_matrix,
        network = context_network,
        .source = "source",
        .target = "target",
        .mor = "mor",
        center = FALSE,
        minsize = MIN_REGULON_SIZE
    )


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


    targets_measured <- as.numeric(
        regulon_sizes[
            tf
        ]
    )


    result <- data.frame(
        TF =
            tf,

        Context =
            context_name,

        Evidence_type =
            evidence_type,

        Activity =
            activity,

        P_value =
            p_value,

        FDR =
            safe_fdr(
                p_value
            ),

        Activity_percentile =
            percentile_signed(
                activity
            ),

        Abs_activity_percentile =
            percentile_absolute(
                activity
            ),

        Targets_measured =
            targets_measured,

        stringsAsFactors = FALSE
    )


    ordering <- order(
        result$FDR,
        -abs(
            result$Activity
        ),
        na.last = TRUE
    )


    result <- result[
        ordering,
        ,
        drop = FALSE
    ]


    rownames(
        result
    ) <- NULL


    result
}


# ============================================================
# HELPER: CONVERT CONTEXT TABLE TO WIDE FORMAT
# ============================================================

wide_context <- function(
    context_table,
    prefix
) {

    required_columns <- c(
        "TF",
        "Activity",
        "P_value",
        "FDR",
        "Activity_percentile",
        "Abs_activity_percentile",
        "Targets_measured"
    )


    out <- context_table[
        ,
        required_columns,
        drop = FALSE
    ]


    names(out) <- c(
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
            "_Pct"
        ),

        paste0(
            prefix,
            "_AbsPct"
        ),

        paste0(
            prefix,
            "_Targets"
        )
    )


    out
}


# ============================================================
# LOAD FROZEN COLLECTRI NETWORK
# ============================================================

message(
    "Loading frozen mouse CollecTRI network..."
)


network <- load_collectri(
    COLLECTRI_FILE
)


message(
    "CollecTRI signed edges: ",
    nrow(
        network
    )
)


message(
    "CollecTRI regulators: ",
    length(
        unique(
            network$source
        )
    )
)


# ============================================================
# LOAD ADULT Rb-TKO
# ============================================================

message("")

message(
    "Loading adult Rb-TKO full profile..."
)


adult_raw <- readr::read_csv(
    ADULT_FILE,
    show_col_types = FALSE
)


adult_primary <- clean_rnaseq(
    adult_raw,
    BASEMEAN_PRIMARY
)


adult_low <- clean_rnaseq(
    adult_raw,
    BASEMEAN_SENSITIVITY
)


adult_primary$Evidence_native <-
    adult_primary$stat


adult_primary$Evidence_rank <-
    signed_rank_gaussian(
        adult_primary$Evidence_native
    )


adult_low$Evidence_native <-
    adult_low$stat


message(
    "Adult primary genes: ",
    nrow(
        adult_primary
    )
)


message(
    "Adult sensitivity genes: ",
    nrow(
        adult_low
    )
)


# ============================================================
# LOAD E2F1/3-DKO
# ============================================================

message("")

message(
    "Loading E2F1/3-DKO full profile..."
)


e2f_raw <- readxl::read_excel(
    E2F_FILE,
    sheet = E2F_SHEET
)


e2f_primary <- clean_rnaseq(
    e2f_raw,
    BASEMEAN_PRIMARY
)


e2f_low <- clean_rnaseq(
    e2f_raw,
    BASEMEAN_SENSITIVITY
)


e2f_primary$Evidence_native <-
    e2f_primary$stat


e2f_primary$Evidence_rank <-
    signed_rank_gaussian(
        e2f_primary$Evidence_native
    )


e2f_low$Evidence_native <-
    e2f_low$stat


message(
    "E2F-DKO primary genes: ",
    nrow(
        e2f_primary
    )
)


message(
    "E2F-DKO sensitivity genes: ",
    nrow(
        e2f_low
    )
)


# ============================================================
# LOAD OSHIKAWA CAG / MAP2
# ============================================================

message("")

message(
    "Loading Oshikawa CAG/MAP2 response matrix..."
)


osh <- readr::read_csv(
    OSH_FILE,
    show_col_types = FALSE
)


required_osh <- c(
    "Gene",
    "logFC_CAG",
    "P.Value_CAG",
    "logFC_MAP2",
    "P.Value_MAP2"
)


missing_osh <- setdiff(
    required_osh,
    names(osh)
)


if (length(
    missing_osh
) > 0) {

    stop(
        paste(
            "Oshikawa table missing:",
            paste(
                missing_osh,
                collapse = ", "
            )
        )
    )
}


gene <- trimws(
    as.character(
        osh[
            ,
            "Gene",
            drop = TRUE
        ]
    )
)


gene_key <- toupper(
    gene
)


valid_gene <- (
    !is.na(
        gene_key
    ) &
    gene_key != ""
)


gene <- gene[
    valid_gene
]


gene_key <- gene_key[
    valid_gene
]


osh <- osh[
    valid_gene,
    ,
    drop = FALSE
]


if (anyDuplicated(
    gene_key
) > 0) {

    stop(
        "Oshikawa matrix contains duplicated case-insensitive gene symbols."
    )
}


cag <- data.frame(
    Gene =
        gene,

    Gene_key =
        gene_key,

    logFC =
        as.numeric(
            osh[
                ,
                "logFC_CAG",
                drop = TRUE
            ]
        ),

    P_value =
        as.numeric(
            osh[
                ,
                "P.Value_CAG",
                drop = TRUE
            ]
        ),

    stringsAsFactors = FALSE
)


cag$Evidence_native <-
    signed_z_from_p(
        cag$logFC,
        cag$P_value
    )


cag$Evidence_rank <-
    signed_rank_gaussian(
        cag$Evidence_native
    )


map2 <- data.frame(
    Gene =
        gene,

    Gene_key =
        gene_key,

    logFC =
        as.numeric(
            osh[
                ,
                "logFC_MAP2",
                drop = TRUE
            ]
        ),

    P_value =
        as.numeric(
            osh[
                ,
                "P.Value_MAP2",
                drop = TRUE
            ]
        ),

    stringsAsFactors = FALSE
)


map2$Evidence_native <-
    signed_z_from_p(
        map2$logFC,
        map2$P_value
    )


map2$Evidence_rank <-
    signed_rank_gaussian(
        map2$Evidence_native
    )


message(
    "CAG/MAP2 genes: ",
    nrow(
        osh
    )
)


# ============================================================
# PRIMARY NATIVE-EVIDENCE TF ACTIVITY
# ============================================================

message("")

message(
    "Running native-evidence TF activity..."
)


adult_native <- run_ulm_context(
    adult_primary,
    "Evidence_native",
    network,
    "Adult_RbTKO",
    "Native"
)


e2f_native <- run_ulm_context(
    e2f_primary,
    "Evidence_native",
    network,
    "E2F13_DKO",
    "Native"
)


cag_native <- run_ulm_context(
    cag,
    "Evidence_native",
    network,
    "CAG_RbLoss",
    "Native"
)


map2_native <- run_ulm_context(
    map2,
    "Evidence_native",
    network,
    "MAP2_RbLoss",
    "Native"
)


native_long <- rbind(
    adult_native,
    e2f_native,
    cag_native,
    map2_native
)


# ============================================================
# RANK-NORMALIZED SENSITIVITY
# ============================================================

message(
    "Running rank-normalized sensitivity..."
)


adult_rank <- run_ulm_context(
    adult_primary,
    "Evidence_rank",
    network,
    "Adult_RbTKO",
    "RankNormalized"
)


e2f_rank <- run_ulm_context(
    e2f_primary,
    "Evidence_rank",
    network,
    "E2F13_DKO",
    "RankNormalized"
)


cag_rank <- run_ulm_context(
    cag,
    "Evidence_rank",
    network,
    "CAG_RbLoss",
    "RankNormalized"
)


map2_rank <- run_ulm_context(
    map2,
    "Evidence_rank",
    network,
    "MAP2_RbLoss",
    "RankNormalized"
)


rank_long <- rbind(
    adult_rank,
    e2f_rank,
    cag_rank,
    map2_rank
)


# ============================================================
# RNA-SEQ DETECTABILITY SENSITIVITY
# ============================================================

message(
    "Running baseMean >= 1 RNA-seq sensitivity..."
)


adult_low_result <- run_ulm_context(
    adult_low,
    "Evidence_native",
    network,
    "Adult_RbTKO_baseMean1",
    "Native_baseMean1"
)


e2f_low_result <- run_ulm_context(
    e2f_low,
    "Evidence_native",
    network,
    "E2F13_DKO_baseMean1",
    "Native_baseMean1"
)


# ============================================================
# BUILD NATIVE WIDE TABLE
# ============================================================

adult_w <- wide_context(
    adult_native,
    "Adult"
)


e2f_w <- wide_context(
    e2f_native,
    "E2F_DKO"
)


cag_w <- wide_context(
    cag_native,
    "CAG"
)


map2_w <- wide_context(
    map2_native,
    "MAP2"
)


native_wide <- merge(
    adult_w,
    e2f_w,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


native_wide <- merge(
    native_wide,
    cag_w,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


native_wide <- merge(
    native_wide,
    map2_w,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


# ============================================================
# BUILD RANK-NORMALIZED WIDE TABLE
# ============================================================

adult_rank_w <- wide_context(
    adult_rank,
    "Adult_rank"
)


e2f_rank_w <- wide_context(
    e2f_rank,
    "E2F_DKO_rank"
)


cag_rank_w <- wide_context(
    cag_rank,
    "CAG_rank"
)


map2_rank_w <- wide_context(
    map2_rank,
    "MAP2_rank"
)


rank_wide <- merge(
    adult_rank_w,
    e2f_rank_w,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


rank_wide <- merge(
    rank_wide,
    cag_rank_w,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


rank_wide <- merge(
    rank_wide,
    map2_rank_w,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


cross_context <- merge(
    native_wide,
    rank_wide,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


# ============================================================
# DERIVED CROSS-CONTEXT METRICS
# ============================================================

cross_context$Cortex_max_Pct <-
    safe_max2(
        cross_context$CAG_Pct,
        cross_context$MAP2_Pct
    )


cross_context$Cortex_mean_Pct <-
    safe_mean2(
        cross_context$CAG_Pct,
        cross_context$MAP2_Pct
    )


cross_context$Adult_minus_max_cortex_Pct <- (
    cross_context$Adult_Pct -
    cross_context$Cortex_max_Pct
)


cross_context$Adult_minus_mean_cortex_Pct <- (
    cross_context$Adult_Pct -
    cross_context$Cortex_mean_Pct
)


# ============================================================
# ADULT POSITIVE / E2F-DKO NEGATIVE
# ============================================================

cross_context$Adult_E2F_reciprocal_activity <- rep(
    NA,
    nrow(
        cross_context
    )
)


valid_reciprocal <- (
    is.finite(
        cross_context$Adult_Activity
    ) &
    is.finite(
        cross_context$E2F_DKO_Activity
    )
)


cross_context$Adult_E2F_reciprocal_activity[
    valid_reciprocal
] <- (
    cross_context$Adult_Activity[
        valid_reciprocal
    ] > 0 &
    cross_context$E2F_DKO_Activity[
        valid_reciprocal
    ] < 0
)


# ============================================================
# CONSERVED SIGN ACROSS THREE Rb-LOSS CONTEXTS
# ============================================================

cross_context$Rb_loss_activity_sign_conserved <- rep(
    NA,
    nrow(
        cross_context
    )
)


valid_rb <- (
    is.finite(
        cross_context$Adult_Activity
    ) &
    is.finite(
        cross_context$CAG_Activity
    ) &
    is.finite(
        cross_context$MAP2_Activity
    )
)


cross_context$Rb_loss_activity_sign_conserved[
    valid_rb
] <- (
    sign(
        cross_context$Adult_Activity[
            valid_rb
        ]
    ) ==
    sign(
        cross_context$CAG_Activity[
            valid_rb
        ]
    ) &
    sign(
        cross_context$Adult_Activity[
            valid_rb
        ]
    ) ==
    sign(
        cross_context$MAP2_Activity[
            valid_rb
        ]
    )
)


# ============================================================
# NATIVE VS RANK-NORMALIZED SIGN ROBUSTNESS
# ============================================================

cross_context$Adult_native_rank_sign_concordant <-
    sign_concordant(
        cross_context$Adult_Activity,
        cross_context$Adult_rank_Activity
    )


cross_context$E2F_native_rank_sign_concordant <-
    sign_concordant(
        cross_context$E2F_DKO_Activity,
        cross_context$E2F_DKO_rank_Activity
    )


cross_context$CAG_native_rank_sign_concordant <-
    sign_concordant(
        cross_context$CAG_Activity,
        cross_context$CAG_rank_Activity
    )


cross_context$MAP2_native_rank_sign_concordant <-
    sign_concordant(
        cross_context$MAP2_Activity,
        cross_context$MAP2_rank_Activity
    )


# ============================================================
# MEAN Rb-LOSS ACTIVITY PERCENTILE
# ============================================================

cross_context$Mean_RbLoss_Pct <- rep(
    NA_real_,
    nrow(
        cross_context
    )
)


for (i in seq_len(
    nrow(
        cross_context
    )
)) {

    values <- c(
        cross_context$Adult_Pct[
            i
        ],

        cross_context$CAG_Pct[
            i
        ],

        cross_context$MAP2_Pct[
            i
        ]
    )


    values <- values[
        is.finite(
            values
        )
    ]


    if (length(
        values
    ) > 0) {

        cross_context$Mean_RbLoss_Pct[
            i
        ] <- mean(
            values
        )
    }
}


# ============================================================
# F10L REPRODUCIBILITY CHECK
# ============================================================

if (file.exists(
    F10L_FILE
)) {

    message("")

    message(
        "Checking adult F10M against frozen F10L..."
    )


    old <- readr::read_csv(
        F10L_FILE,
        show_col_types = FALSE
    )


    required_f10l <- c(
        "TF",
        "Contrast_activity"
    )


    if (all(
        required_f10l %in%
            names(old)
    )) {

        old_table <- data.frame(
            TF =
                as.character(
                    old[
                        ,
                        "TF",
                        drop = TRUE
                    ]
                ),

            F10L_Adult =
                as.numeric(
                    old[
                        ,
                        "Contrast_activity",
                        drop = TRUE
                    ]
                ),

            stringsAsFactors = FALSE
        )


        current_table <- data.frame(
            TF =
                adult_native$TF,

            F10M_Adult =
                adult_native$Activity,

            stringsAsFactors = FALSE
        )


        check_table <- merge(
            current_table,
            old_table,
            by = "TF",
            all = FALSE,
            sort = FALSE
        )


        check_table$Difference <- (
            check_table$F10M_Adult -
            check_table$F10L_Adult
        )


        correlation <- cor(
            check_table$F10M_Adult,
            check_table$F10L_Adult,
            use = "complete.obs"
        )


        maximum_difference <- max(
            abs(
                check_table$Difference
            ),
            na.rm = TRUE
        )


        message(
            "Adult F10M vs F10L Pearson r: ",
            signif(
                correlation,
                8
            )
        )


        message(
            "Maximum absolute score difference: ",
            signif(
                maximum_difference,
                8
            )
        )


        readr::write_csv(
            check_table,
            file.path(
                OUTPUT_DIR,
                "F10M_F10L_adult_reproducibility_check.csv"
            )
        )
    }
}


# ============================================================
# BASEMEAN SENSITIVITY TABLE
# ============================================================

adult_threshold_primary <- data.frame(
    TF =
        adult_native$TF,

    Adult_baseMean10 =
        adult_native$Activity,

    stringsAsFactors = FALSE
)


adult_threshold_low <- data.frame(
    TF =
        adult_low_result$TF,

    Adult_baseMean1 =
        adult_low_result$Activity,

    stringsAsFactors = FALSE
)


e2f_threshold_primary <- data.frame(
    TF =
        e2f_native$TF,

    E2F_baseMean10 =
        e2f_native$Activity,

    stringsAsFactors = FALSE
)


e2f_threshold_low <- data.frame(
    TF =
        e2f_low_result$TF,

    E2F_baseMean1 =
        e2f_low_result$Activity,

    stringsAsFactors = FALSE
)


threshold_sensitivity <- merge(
    adult_threshold_primary,
    adult_threshold_low,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


threshold_sensitivity <- merge(
    threshold_sensitivity,
    e2f_threshold_primary,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


threshold_sensitivity <- merge(
    threshold_sensitivity,
    e2f_threshold_low,
    by = "TF",
    all = TRUE,
    sort = FALSE
)


threshold_sensitivity$Adult_threshold_sign_concordant <-
    sign_concordant(
        threshold_sensitivity$Adult_baseMean10,
        threshold_sensitivity$Adult_baseMean1
    )


threshold_sensitivity$E2F_threshold_sign_concordant <-
    sign_concordant(
        threshold_sensitivity$E2F_baseMean10,
        threshold_sensitivity$E2F_baseMean1
    )


# ============================================================
# CANDIDATE AUDIT
# ============================================================

candidate_index <- (
    toupper(
        cross_context$TF
    ) %in%
    CANDIDATE_TFS
)


candidate_audit <- cross_context[
    candidate_index,
    ,
    drop = FALSE
]


candidate_audit <- candidate_audit[
    order(
        -candidate_audit$Adult_Pct,
        na.last = TRUE
    ),
    ,
    drop = FALSE
]


rownames(
    candidate_audit
) <- NULL


# ============================================================
# EXPLORATORY ADULT-SELECTIVITY RANKING
# ============================================================
#
# This is NOT a formal context-interaction test.
#
# Requirements:
#
#   Adult activity > 0
#   Adult FDR < 0.05
#   Adult native/rank sign agreement
#
# Ranking:
#
#   Adult percentile
#       -
#   max(CAG percentile, MAP2 percentile)
#
# ============================================================

adult_specific_index <- (
    is.finite(
        cross_context$Adult_Activity
    ) &
    cross_context$Adult_Activity > 0 &
    is.finite(
        cross_context$Adult_FDR
    ) &
    cross_context$Adult_FDR < 0.05 &
    !is.na(
        cross_context$Adult_native_rank_sign_concordant
    ) &
    cross_context$Adult_native_rank_sign_concordant &
    is.finite(
        cross_context$Adult_minus_max_cortex_Pct
    )
)


adult_specificity <- cross_context[
    adult_specific_index,
    ,
    drop = FALSE
]


adult_specificity <- adult_specificity[
    order(
        -adult_specificity$Adult_minus_max_cortex_Pct,
        na.last = TRUE
    ),
    ,
    drop = FALSE
]


rownames(
    adult_specificity
) <- NULL


# ============================================================
# CONSERVED Rb-LOSS ACTIVITY RANKING
# ============================================================

conserved_index <- (
    is.finite(
        cross_context$Adult_Activity
    ) &
    is.finite(
        cross_context$CAG_Activity
    ) &
    is.finite(
        cross_context$MAP2_Activity
    ) &
    cross_context$Adult_Activity > 0 &
    cross_context$CAG_Activity > 0 &
    cross_context$MAP2_Activity > 0 &
    !is.na(
        cross_context$Adult_native_rank_sign_concordant
    ) &
    !is.na(
        cross_context$CAG_native_rank_sign_concordant
    ) &
    !is.na(
        cross_context$MAP2_native_rank_sign_concordant
    ) &
    cross_context$Adult_native_rank_sign_concordant &
    cross_context$CAG_native_rank_sign_concordant &
    cross_context$MAP2_native_rank_sign_concordant
)


conserved_rb <- cross_context[
    conserved_index,
    ,
    drop = FALSE
]


conserved_rb <- conserved_rb[
    order(
        -conserved_rb$Mean_RbLoss_Pct,
        na.last = TRUE
    ),
    ,
    drop = FALSE
]


rownames(
    conserved_rb
) <- NULL


# ============================================================
# SAVE OUTPUTS
# ============================================================

readr::write_csv(
    native_long,
    file.path(
        OUTPUT_DIR,
        "F10M_TF_activity_native_long.csv"
    )
)


readr::write_csv(
    rank_long,
    file.path(
        OUTPUT_DIR,
        "F10M_TF_activity_ranknormalized_long.csv"
    )
)


readr::write_csv(
    cross_context,
    file.path(
        OUTPUT_DIR,
        "F10M_cross_context_TF_activity.csv"
    )
)


readr::write_csv(
    candidate_audit,
    file.path(
        OUTPUT_DIR,
        "F10M_candidate_TF_audit.csv"
    )
)


readr::write_csv(
    adult_specificity,
    file.path(
        OUTPUT_DIR,
        "F10M_adult_specificity_ranking.csv"
    )
)


readr::write_csv(
    conserved_rb,
    file.path(
        OUTPUT_DIR,
        "F10M_conserved_Rb_loss_activity.csv"
    )
)


readr::write_csv(
    threshold_sensitivity,
    file.path(
        OUTPUT_DIR,
        "F10M_baseMean_threshold_sensitivity.csv"
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
    "F10M — CROSS-CONTEXT SIGNED TF ACTIVITY"
)

message(
    "============================================================"
)

message("")


message(
    "Adult genes: ",
    nrow(
        adult_primary
    )
)


message(
    "E2F-DKO genes: ",
    nrow(
        e2f_primary
    )
)


message(
    "CAG/MAP2 genes: ",
    nrow(
        osh
    )
)


message(
    "Adult TFs: ",
    nrow(
        adult_native
    )
)


message(
    "E2F-DKO TFs: ",
    nrow(
        e2f_native
    )
)


message(
    "CAG TFs: ",
    nrow(
        cag_native
    )
)


message(
    "MAP2 TFs: ",
    nrow(
        map2_native
    )
)


# ============================================================
# PRINT CANDIDATE AUDIT
# ============================================================

candidate_columns <- c(
    "TF",

    "Adult_Activity",
    "Adult_FDR",
    "Adult_Pct",

    "E2F_DKO_Activity",
    "E2F_DKO_FDR",
    "E2F_DKO_Pct",

    "CAG_Activity",
    "CAG_FDR",
    "CAG_Pct",

    "MAP2_Activity",
    "MAP2_FDR",
    "MAP2_Pct",

    "Adult_minus_max_cortex_Pct",

    "Adult_E2F_reciprocal_activity",

    "Rb_loss_activity_sign_conserved",

    "Adult_native_rank_sign_concordant",

    "E2F_native_rank_sign_concordant",

    "CAG_native_rank_sign_concordant",

    "MAP2_native_rank_sign_concordant"
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
        candidate_columns,
        drop = FALSE
    ],
    row.names = FALSE
)


# ============================================================
# PRINT TOP ADULT-SELECTIVE TFs
# ============================================================

adult_columns <- c(
    "TF",

    "Adult_Activity",
    "Adult_FDR",
    "Adult_Pct",

    "CAG_Activity",
    "CAG_FDR",
    "CAG_Pct",

    "MAP2_Activity",
    "MAP2_FDR",
    "MAP2_Pct",

    "Adult_minus_max_cortex_Pct",

    "E2F_DKO_Activity",
    "E2F_DKO_FDR"
)


message("")

message(
    "TOP ADULT-SELECTIVE ACTIVITY PATTERNS"
)

message(
    "------------------------------------"
)


n_show <- min(
    30,
    nrow(
        adult_specificity
    )
)


if (n_show > 0) {

    print(
        adult_specificity[
            seq_len(
                n_show
            ),
            adult_columns,
            drop = FALSE
        ],
        row.names = FALSE
    )

} else {

    message(
        "No TFs passed the adult-selectivity screen."
    )
}


# ============================================================
# PRINT TOP CONSERVED Rb-LOSS TFs
# ============================================================

conserved_columns <- c(
    "TF",

    "Adult_Activity",
    "Adult_FDR",
    "Adult_Pct",

    "CAG_Activity",
    "CAG_FDR",
    "CAG_Pct",

    "MAP2_Activity",
    "MAP2_FDR",
    "MAP2_Pct",

    "E2F_DKO_Activity",
    "E2F_DKO_FDR",

    "Adult_E2F_reciprocal_activity"
)


message("")

message(
    "TOP CONSERVED Rb-LOSS ACTIVITY PATTERNS"
)

message(
    "---------------------------------------"
)


n_show <- min(
    30,
    nrow(
        conserved_rb
    )
)


if (n_show > 0) {

    print(
        conserved_rb[
            seq_len(
                n_show
            ),
            conserved_columns,
            drop = FALSE
        ],
        row.names = FALSE
    )

} else {

    message(
        "No TFs passed the conserved-Rb-loss screen."
    )
}


# ============================================================
# FINISH
# ============================================================

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
        "1. Compare activity direction, within-context percentile,",
        "FDR, and robustness across platforms; do not interpret",
        "raw ULM magnitude as directly comparable between platforms."
    )
)


message(
    paste(
        "2. Adult_minus_max_cortex_Pct is an exploratory",
        "context-specificity ranking, not a formal interaction test."
    )
)


message(
    paste(
        "3. Related TFs can share regulon targets, so family-level",
        "regulatory conclusions can be more defensible than",
        "assigning unique causal activity to one TF."
    )
)


message(
    paste(
        "4. ULM activity is inferred from target behavior and does",
        "not directly measure TF protein activity, localization,",
        "occupancy, phosphorylation, or causal necessity."
    )
)