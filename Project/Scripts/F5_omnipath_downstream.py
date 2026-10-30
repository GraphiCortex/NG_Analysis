#!/usr/bin/env python3

"""
Phase F5 — Corrected directed OmniPath downstream analysis
==========================================================

Why this version exists
-----------------------
The first F5 diagnostic run showed two problems:

1. Every reachable source-target pair hit the 25-path cap.
   Therefore raw retained-path count was saturated and uninformative.

2. The original E4_X2 / E5_X2 / E4_X3 / E5_X3 settings were too
   permissive after adding a ~149k-edge signaling + TF network.

This corrected version separates:

    EXACT CONSTRAINED REACHABILITY
        Does ANY route exist under a given E/X constraint?

    EXACT MINIMUM EXTERNAL SUPPORT
        What is the fewest number of non-F2 intermediates required?

    EXACT REQUIRED-NODE TESTING
        If an intermediate is removed, does ANY allowed route remain?

    REPRESENTATIVE PATH ENUMERATION
        A capped set of short, experimentally anchored routes retained
        only for inspection/visualization — never used as evidence that
        no other paths exist.


Prior-knowledge network
-----------------------
SIGNALING:
    omnipath
    pathwayextra
    kinaseextra
    ligrecextra

TRANSCRIPTIONAL:
    collectri
    dorothea A/B
    tf_target

Species:
    Mus musculus, NCBI taxonomy 10090


Experimental seed universe
--------------------------
321 F2 genes.

X counts intermediate graph genes OUTSIDE the F2 experimental universe.
Source and target anchors never count toward X.


Corrected robustness settings
-----------------------------
    E4_X0   <=4 edges, 0 external intermediates
    E4_X1   <=4 edges, <=1 external intermediate
    E4_X2   <=4 edges, <=2 external intermediates
    E5_X2   <=5 edges, <=2 external intermediates

These define increasingly permissive mechanistic support tiers:

    Tier A  reachable under E4_X0
            entirely F2-anchored intermediates

    Tier B  first becomes reachable under E4_X1

    Tier C  first becomes reachable under E4_X2

    Tier D  requires E5_X2

    UNREACHABLE
            no route under any tested condition


Source families
---------------
E2F_cycle:
    E2F1, E2F2, MYC, CDK1, CDK2

replication_stress:
    PCNA, H2AFX, CHEK1, CHEK2,
    TOPBP1, CLSPN, BRCA1, EXO1


Target families
---------------
apoptosis:
    PMAIP1, BBC3, FAS, TNFRSF10B

regulated_necrosis:
    ZBP1, PYCARD, MLKL

survival_trophic:
    BDNF, FGF2, KDR, MET, IRS2, ANGPT1


Identifier resolution
---------------------
Configured anchor names are preserved.

If an anchor is absent from OmniPath, explicitly declared aliases are
tested and recorded in the anchor audit.

For example:

    H2AFX -> H2AX

No substitution is silent.


Outputs
-------
Data/Processed/FongDownstream/F5_Mechanistic/Corrected/

    F5_corrected_summary.xlsx
    anchor_audit.csv
    setting_summary.csv
    pair_setting_summary.csv
    pair_robustness.csv
    family_summary.csv
    representative_paths.csv
    exact_required_nodes.csv
    node_robustness.csv
    representative_union_edges.csv
    F5_corrected_metadata.json

    graphs/
        E4_X0.graphml
        E4_X1.graphml
        E4_X2.graphml
        E5_X2.graphml


Raw data cache
--------------
Data/Raw/Network/FongDownstream/OmniPath/

Existing signaling cache is reused.
The new transcriptional layers are cached separately.


Run
---
From Project/:

    py .\\Scripts\\F5_omnipath_downstream.py
"""

from __future__ import annotations

from pathlib import Path
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from typing import Any
import hashlib
import heapq
import io
import itertools
import json
import math
import time

import pandas as pd
import requests
import networkx as nx

from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


# ============================================================
# Project paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

F2_MASTER_FILE = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F2_CellFate"
    / "cell_fate_master.csv"
)

F4_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F4_Modules"
)

CORE_ASSIGN_FILE = (
    F4_DIR
    / "CORE"
    / "module_assignments.csv"
)

EXT_ASSIGN_FILE = (
    F4_DIR
    / "EXTENDED"
    / "module_assignments.csv"
)

F5_PARENT_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F5_Mechanistic"
)

OUTPUT_DIR = (
    F5_PARENT_DIR
    / "Corrected"
)

GRAPH_DIR = (
    OUTPUT_DIR
    / "graphs"
)

RAW_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Raw"
    / "Network"
    / "FongDownstream"
    / "OmniPath"
)

OUTPUT_XLSX = (
    OUTPUT_DIR
    / "F5_corrected_summary.xlsx"
)

OUTPUT_METADATA = (
    OUTPUT_DIR
    / "F5_corrected_metadata.json"
)


# ============================================================
# OmniPath cache paths
# ============================================================

SIGNALING_RAW = (
    RAW_DIR
    / "omnipath_mouse_signaling.tsv"
)

COLLECTRI_RAW = (
    RAW_DIR
    / "omnipath_mouse_collectri.tsv"
)

DOROTHEA_RAW = (
    RAW_DIR
    / "omnipath_mouse_dorothea_AB.tsv"
)

TF_TARGET_RAW = (
    RAW_DIR
    / "omnipath_mouse_tf_target.tsv"
)

RAW_METADATA = (
    RAW_DIR
    / "corrected_download_metadata.json"
)


# ============================================================
# OmniPath configuration
# ============================================================

OMNIPATH_URL = (
    "https://omnipathdb.org/interactions"
)

SPECIES = 10090

REQUEST_TIMEOUT = 300

FORCE_REFRESH = False

USER_AGENT = (
    "BIOL_295_FongDownstream_F5_corrected/1.0"
)


SIGNALING_PARAMS = {
    "genesymbols": "yes",
    "datasets":
        "omnipath,pathwayextra,kinaseextra,ligrecextra",
    "organisms":
        str(SPECIES),
    "fields":
        "sources,references,curation_effort",
    "license":
        "academic",
    "format":
        "tsv",
}


COLLECTRI_PARAMS = {
    "genesymbols":
        "yes",
    "datasets":
        "collectri",
    "organisms":
        str(SPECIES),
    "fields":
        "sources,references,curation_effort",
    "license":
        "academic",
    "format":
        "tsv",
}


DOROTHEA_PARAMS = {
    "genesymbols":
        "yes",
    "datasets":
        "dorothea",
    "dorothea_levels":
        "A,B",
    "organisms":
        str(SPECIES),
    "fields":
        "sources,references,curation_effort",
    "license":
        "academic",
    "format":
        "tsv",
}


TF_TARGET_PARAMS = {
    "genesymbols":
        "yes",
    "datasets":
        "tf_target",
    "organisms":
        str(SPECIES),
    "fields":
        "sources,references,curation_effort",
    "license":
        "academic",
    "format":
        "tsv",
}


# ============================================================
# Mechanistic anchors
# ============================================================

SOURCE_FAMILIES = {

    "E2F_cycle": [
        "E2F1",
        "E2F2",
        "MYC",
        "CDK1",
        "CDK2",
    ],

    "replication_stress": [
        "PCNA",
        "H2AFX",
        "CHEK1",
        "CHEK2",
        "TOPBP1",
        "CLSPN",
        "BRCA1",
        "EXO1",
    ],
}


TARGET_FAMILIES = {

    "apoptosis": [
        "PMAIP1",
        "BBC3",
        "FAS",
        "TNFRSF10B",
    ],

    "regulated_necrosis": [
        "ZBP1",
        "PYCARD",
        "MLKL",
    ],

    "survival_trophic": [
        "BDNF",
        "FGF2",
        "KDR",
        "MET",
        "IRS2",
        "ANGPT1",
    ],
}


# Explicit aliases.
#
# These are never silently substituted.
# Every resolution is written to anchor_audit.csv.
#
# H2AFX -> H2AX is especially important because mouse nomenclature
# now uses H2ax as the official symbol.

ANCHOR_ALIASES = {

    "H2AFX": [
        "H2AX",
    ],

    "PMAIP1": [
        "NOXA",
    ],

    "BBC3": [
        "PUMA",
    ],

    "ZBP1": [
        "DAI",
        "DLM1",
    ],

    "MET": [
        "HGFR",
    ],
}


# ============================================================
# Corrected E/X settings
# ============================================================

SETTINGS = {

    "E4_X0": {
        "max_edges": 4,
        "max_external": 0,
    },

    "E4_X1": {
        "max_edges": 4,
        "max_external": 1,
    },

    "E4_X2": {
        "max_edges": 4,
        "max_external": 2,
    },

    "E5_X2": {
        "max_edges": 5,
        "max_external": 2,
    },
}


SETTING_ORDER = [
    "E4_X0",
    "E4_X1",
    "E4_X2",
    "E5_X2",
]


SUPPORT_TIER = {
    "E4_X0": "A",
    "E4_X1": "B",
    "E4_X2": "C",
    "E5_X2": "D",
}


# ============================================================
# Representative path configuration
# ============================================================

# These paths are illustrative only.
# They are NOT used to determine whether a pair is reachable.

MAX_REPRESENTATIVE_PATHS = 25

# After finding:
#
#     minimal external count X*
#     shortest route length E*
#
# retain representative routes up to E* + slack.
#
# This keeps examples local rather than dumping arbitrary long routes.

REPRESENTATIVE_EDGE_SLACK = 1

MAX_REPRESENTATIVE_SEARCH_STATES = 500_000


# ============================================================
# Helpers
# ============================================================

def normalize_gene(
    value: Any,
) -> str | None:

    if pd.isna(value):
        return None

    value = str(
        value
    ).strip()

    if not value:
        return None

    if value.upper() in {
        "NA",
        "N/A",
        "NONE",
        "NULL",
    }:
        return None

    return value.upper()


def bool_value(
    value: Any,
) -> bool:

    if isinstance(
        value,
        bool,
    ):
        return value

    if pd.isna(
        value
    ):
        return False

    return (
        str(
            value
        )
        .strip()
        .upper()
        in {
            "1",
            "TRUE",
            "YES",
            "Y",
        }
    )


def find_column(
    df: pd.DataFrame,
    candidates: list[str],
) -> str | None:

    lower_lookup = {
        str(column).lower():
            column
        for column
        in df.columns
    }

    for candidate in candidates:

        key = (
            candidate.lower()
        )

        if key in lower_lookup:

            return lower_lookup[
                key
            ]

    return None


def sha256_file(
    path: Path,
) -> str:

    digest = (
        hashlib.sha256()
    )

    with open(
        path,
        "rb",
    ) as handle:

        while True:

            block = handle.read(
                1024 * 1024
            )

            if not block:
                break

            digest.update(
                block
            )

    return digest.hexdigest()


# ============================================================
# Download / cache OmniPath
# ============================================================

def download_table(
    label: str,
    params: dict,
    output_file: Path,
) -> pd.DataFrame:

    if (
        output_file.exists()
        and
        not FORCE_REFRESH
    ):

        print(
            f"[INFO] Using cached {label}:"
        )

        print(
            f"       {output_file}"
        )

        return pd.read_csv(
            output_file,
            sep="\t",
            low_memory=False,
        )

    print(
        f"[INFO] Downloading {label} from OmniPath..."
    )

    response = requests.get(
        OMNIPATH_URL,
        params=params,
        timeout=REQUEST_TIMEOUT,
        headers={
            "User-Agent":
                USER_AGENT
        },
    )

    response.raise_for_status()

    text = (
        response.text
    )

    if (
        not text.strip()
        or
        "\t"
        not in text
    ):

        raise RuntimeError(
            f"Unexpected OmniPath response for {label}:\n"
            f"{text[:1000]}"
        )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_file.write_text(
        text,
        encoding="utf-8",
    )

    result = pd.read_csv(
        io.StringIO(
            text
        ),
        sep="\t",
        low_memory=False,
    )

    print(
        f"[INFO] Downloaded "
        f"{len(result):,} {label} rows."
    )

    time.sleep(
        1.0
    )

    return result


# ============================================================
# Clean OmniPath edges
# ============================================================

def get_effect(
    row: pd.Series,
    columns: set[str],
) -> int:

    consensus_stim = (
        bool_value(
            row.get(
                "consensus_stimulation"
            )
        )
        if "consensus_stimulation"
        in columns
        else False
    )

    consensus_inhib = (
        bool_value(
            row.get(
                "consensus_inhibition"
            )
        )
        if "consensus_inhibition"
        in columns
        else False
    )

    if (
        consensus_stim
        or
        consensus_inhib
    ):

        if (
            consensus_stim
            and
            not consensus_inhib
        ):
            return 1

        if (
            consensus_inhib
            and
            not consensus_stim
        ):
            return -1

        return 0

    stim = (
        bool_value(
            row.get(
                "is_stimulation"
            )
        )
        if "is_stimulation"
        in columns
        else False
    )

    inhib = (
        bool_value(
            row.get(
                "is_inhibition"
            )
        )
        if "is_inhibition"
        in columns
        else False
    )

    if stim and not inhib:
        return 1

    if inhib and not stim:
        return -1

    return 0


def prepare_edges(
    raw: pd.DataFrame,
    layer: str,
) -> pd.DataFrame:

    source_col = find_column(
        raw,
        [
            "source_genesymbol",
            "source_genesymbols",
            "source_genesymb",
        ],
    )

    target_col = find_column(
        raw,
        [
            "target_genesymbol",
            "target_genesymbols",
            "target_genesymb",
        ],
    )

    if (
        source_col is None
        or
        target_col is None
    ):

        raise ValueError(
            "Could not identify OmniPath gene-symbol columns.\n"
            f"Columns:\n{list(raw.columns)}"
        )

    df = raw.copy()

    directed_col = find_column(
        df,
        [
            "is_directed",
            "consensus_direction",
        ],
    )

    if directed_col is not None:

        df = df[
            df[
                directed_col
            ]
            .map(
                bool_value
            )
        ].copy()

    df[
        "source_gene"
    ] = (
        df[
            source_col
        ]
        .map(
            normalize_gene
        )
    )

    df[
        "target_gene"
    ] = (
        df[
            target_col
        ]
        .map(
            normalize_gene
        )
    )

    df = df[
        df[
            "source_gene"
        ].notna()
        &
        df[
            "target_gene"
        ].notna()
    ].copy()

    df = df[
        df[
            "source_gene"
        ]
        !=
        df[
            "target_gene"
        ]
    ].copy()

    column_set = set(
        df.columns
    )

    df[
        "effect"
    ] = df.apply(
        lambda row:
            get_effect(
                row,
                column_set,
            ),
        axis=1,
    )

    df[
        "layer"
    ] = layer

    if "sources" not in df.columns:
        df[
            "sources"
        ] = ""

    if "references" not in df.columns:
        df[
            "references"
        ] = ""

    if "curation_effort" not in df.columns:
        df[
            "curation_effort"
        ] = 0

    return df[
        [
            "source_gene",
            "target_gene",
            "effect",
            "layer",
            "sources",
            "references",
            "curation_effort",
        ]
    ].copy()


def merge_text(
    values,
) -> str:

    output = set()

    for value in values:

        if pd.isna(
            value
        ):
            continue

        for item in str(
            value
        ).split(
            ";"
        ):

            item = (
                item.strip()
            )

            if item:

                output.add(
                    item
                )

    return ";".join(
        sorted(
            output
        )
    )


def aggregate_edges(
    edges: pd.DataFrame,
) -> pd.DataFrame:

    rows = []

    for (
        source,
        target
    ), group in edges.groupby(
        [
            "source_gene",
            "target_gene",
        ],
        sort=False,
    ):

        signed = {
            int(value)
            for value
            in group[
                "effect"
            ]
            if int(value) != 0
        }

        if signed == {
            1
        }:

            effect = 1
            label = (
                "stimulation"
            )

        elif signed == {
            -1
        }:

            effect = -1
            label = (
                "inhibition"
            )

        elif len(
            signed
        ) > 1:

            effect = 0
            label = (
                "ambiguous"
            )

        else:

            effect = 0
            label = (
                "unsigned"
            )

        rows.append(
            {
                "source_gene":
                    source,

                "target_gene":
                    target,

                "effect":
                    effect,

                "effect_label":
                    label,

                "layers":
                    ";".join(
                        sorted(
                            set(
                                group[
                                    "layer"
                                ]
                                .astype(
                                    str
                                )
                            )
                        )
                    ),

                "sources":
                    merge_text(
                        group[
                            "sources"
                        ]
                    ),

                "references":
                    merge_text(
                        group[
                            "references"
                        ]
                    ),

                "max_curation_effort":
                    float(
                        pd.to_numeric(
                            group[
                                "curation_effort"
                            ],
                            errors="coerce",
                        )
                        .fillna(
                            0
                        )
                        .max()
                    ),

                "supporting_records":
                    len(
                        group
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


def build_graph(
    edges: pd.DataFrame,
) -> nx.DiGraph:

    G = nx.DiGraph()

    for _, row in edges.iterrows():

        G.add_edge(
            row[
                "source_gene"
            ],
            row[
                "target_gene"
            ],

            effect=int(
                row[
                    "effect"
                ]
            ),

            effect_label=str(
                row[
                    "effect_label"
                ]
            ),

            layers=str(
                row[
                    "layers"
                ]
            ),

            sources=str(
                row[
                    "sources"
                ]
            ),

            max_curation_effort=float(
                row[
                    "max_curation_effort"
                ]
            ),

            supporting_records=int(
                row[
                    "supporting_records"
                ]
            ),
        )

    return G


# ============================================================
# Anchor resolution
# ============================================================

def resolve_anchor(
    canonical_gene: str,
    G: nx.DiGraph,
) -> tuple[
    str | None,
    str,
]:

    canonical_gene = (
        normalize_gene(
            canonical_gene
        )
    )

    if canonical_gene in G:

        return (
            canonical_gene,
            "exact",
        )

    aliases = (
        ANCHOR_ALIASES.get(
            canonical_gene,
            [],
        )
    )

    available_aliases = [
        normalize_gene(
            alias
        )
        for alias
        in aliases
        if normalize_gene(
            alias
        )
        in G
    ]

    if len(
        available_aliases
    ) == 1:

        return (
            available_aliases[
                0
            ],
            "alias",
        )

    if len(
        available_aliases
    ) > 1:

        return (
            available_aliases[
                0
            ],
            "multiple_aliases_first_used",
        )

    return (
        None,
        "missing",
    )


# ============================================================
# Reverse-distance cache
# ============================================================

REVERSE_DISTANCE_CACHE = {}


def reverse_distances(
    G: nx.DiGraph,
    target: str,
    max_edges: int,
) -> dict[str, int]:

    key = (
        target,
        max_edges,
    )

    if key not in REVERSE_DISTANCE_CACHE:

        reverse_graph = (
            G.reverse(
                copy=False
            )
        )

        REVERSE_DISTANCE_CACHE[
            key
        ] = (
            nx.single_source_shortest_path_length(
                reverse_graph,
                target,
                cutoff=max_edges,
            )
        )

    return (
        REVERSE_DISTANCE_CACHE[
            key
        ]
    )


# ============================================================
# Exact constrained reachability
# ============================================================

def constrained_target_states(
    G: nx.DiGraph,
    source: str,
    target: str,
    seed_graph_genes: set[str],
    max_edges: int,
    max_external: int,
    blocked_node: str | None = None,
) -> set[tuple[int, int]]:

    """
    EXACT bounded reachability.

    State:
        (node, edges_used, external_intermediates_used)

    Cycles do not create false positive reachability here:
    any bounded walk reaching target can have cycles removed,
    yielding an equal-or-shorter simple path with no larger X.

    Returns all reachable target states:

        (edge_count, external_count)
    """

    if (
        source not in G
        or
        target not in G
    ):

        return set()

    if blocked_node in {
        source,
        target,
    }:

        return set()

    distances = (
        reverse_distances(
            G,
            target,
            max_edges,
        )
    )

    if source not in distances:

        return set()

    queue = deque(
        [
            (
                source,
                0,
                0,
            )
        ]
    )

    seen = {
        (
            source,
            0,
            0,
        )
    }

    target_states = set()

    while queue:

        (
            current,
            edges_used,
            external_used,
        ) = queue.popleft()

        if edges_used >= max_edges:
            continue

        for neighbor in G.successors(
            current
        ):

            if (
                blocked_node
                is not None
                and
                neighbor
                ==
                blocked_node
            ):
                continue

            new_edges = (
                edges_used
                + 1
            )

            if (
                neighbor
                !=
                target
            ):

                remaining = (
                    max_edges
                    -
                    new_edges
                )

                distance = (
                    distances.get(
                        neighbor
                    )
                )

                if (
                    distance
                    is None
                    or
                    distance
                    >
                    remaining
                ):
                    continue

            new_external = (
                external_used
            )

            if (
                neighbor
                !=
                target
                and
                neighbor
                not in seed_graph_genes
            ):

                new_external += 1

            if (
                new_external
                >
                max_external
            ):
                continue

            if neighbor == target:

                target_states.add(
                    (
                        new_edges,
                        new_external,
                    )
                )

                continue

            state = (
                neighbor,
                new_edges,
                new_external,
            )

            if state in seen:
                continue

            seen.add(
                state
            )

            queue.append(
                state
            )

    return target_states


def exact_reachable(
    G: nx.DiGraph,
    source: str,
    target: str,
    seed_graph_genes: set[str],
    max_edges: int,
    max_external: int,
    blocked_node: str | None = None,
) -> bool:

    return bool(
        constrained_target_states(
            G=G,
            source=source,
            target=target,
            seed_graph_genes=seed_graph_genes,
            max_edges=max_edges,
            max_external=max_external,
            blocked_node=blocked_node,
        )
    )


# ============================================================
# Representative path enumeration
# ============================================================

def representative_paths(
    G: nx.DiGraph,
    source: str,
    target: str,
    seed_graph_genes: set[str],
    max_edges: int,
    max_external: int,
    max_paths: int,
) -> tuple[
    list[list[str]],
    bool,
]:

    """
    Enumerate short SIMPLE representative paths.

    Priority:
        1. fewer edges
        2. deterministic lexical tie break

    The caller already supplies the EXACT minimum external count,
    so all successful representatives use the strongest available
    experimental anchoring.

    These paths are illustrative only.
    """

    distances = (
        reverse_distances(
            G,
            target,
            max_edges,
        )
    )

    if source not in distances:

        return (
            [],
            False,
        )

    counter = (
        itertools.count()
    )

    heap = []

    heapq.heappush(
        heap,
        (
            0,
            next(
                counter
            ),
            [
                source
            ],
            0,
        )
    )

    paths = []

    states_examined = 0

    state_cap_hit = False

    while heap:

        (
            edges_used,
            _,
            path,
            external_used,
        ) = heapq.heappop(
            heap
        )

        states_examined += 1

        if (
            states_examined
            >
            MAX_REPRESENTATIVE_SEARCH_STATES
        ):

            state_cap_hit = True
            break

        current = (
            path[
                -1
            ]
        )

        if edges_used >= max_edges:
            continue

        for neighbor in sorted(
            G.successors(
                current
            )
        ):

            # Representative paths themselves are simple.
            if neighbor in path:
                continue

            new_edges = (
                edges_used
                + 1
            )

            if (
                neighbor
                !=
                target
            ):

                remaining = (
                    max_edges
                    -
                    new_edges
                )

                distance = (
                    distances.get(
                        neighbor
                    )
                )

                if (
                    distance
                    is None
                    or
                    distance
                    >
                    remaining
                ):
                    continue

            new_external = (
                external_used
            )

            if (
                neighbor
                !=
                target
                and
                neighbor
                not in seed_graph_genes
            ):

                new_external += 1

            if (
                new_external
                >
                max_external
            ):
                continue

            new_path = (
                path
                +
                [
                    neighbor
                ]
            )

            if neighbor == target:

                paths.append(
                    new_path
                )

                if len(
                    paths
                ) >= max_paths:

                    return (
                        paths,
                        state_cap_hit,
                    )

                continue

            heapq.heappush(
                heap,
                (
                    new_edges,
                    next(
                        counter
                    ),
                    new_path,
                    new_external,
                )
            )

    return (
        paths,
        state_cap_hit,
    )


# ============================================================
# Path annotations
# ============================================================

def annotate_path_effect(
    G: nx.DiGraph,
    path: list[str],
) -> tuple[
    str,
    int | None,
]:

    effects = []

    for a, b in zip(
        path[:-1],
        path[1:],
    ):

        effects.append(
            int(
                G[
                    a
                ][
                    b
                ].get(
                    "effect",
                    0,
                )
            )
        )

    if any(
        value == 0
        for value
        in effects
    ):

        return (
            "unknown_or_ambiguous",
            None,
        )

    product = (
        math.prod(
            effects
        )
    )

    if product > 0:

        return (
            "net_positive",
            1,
        )

    return (
        "net_negative",
        -1,
    )


def annotate_path_layers(
    G: nx.DiGraph,
    path: list[str],
) -> str:

    layers = set()

    for a, b in zip(
        path[:-1],
        path[1:],
    ):

        for layer in str(
            G[
                a
            ][
                b
            ].get(
                "layers",
                "",
            )
        ).split(
            ";"
        ):

            layer = (
                layer.strip()
            )

            if layer:

                layers.add(
                    layer
                )

    return ";".join(
        sorted(
            layers
        )
    )


# ============================================================
# Exact required-node test
# ============================================================

def exact_required_nodes(
    G: nx.DiGraph,
    source: str,
    target: str,
    seed_graph_genes: set[str],
    max_edges: int,
    max_external: int,
    representative: list[list[str]],
) -> list[dict]:

    """
    Find exact required intermediates.

    Any truly required node must occur in EVERY representative path.

    Therefore:
        1. intersect representative intermediates;
        2. only those candidates need an expensive exact deletion test;
        3. remove candidate;
        4. rerun exact constrained reachability.

    If target becomes unreachable, the node is genuinely required
    under the full E/X constraint — not merely required among the
    sampled representative paths.
    """

    if not representative:
        return []

    intermediate_sets = [
        set(
            path[
                1:-1
            ]
        )
        for path
        in representative
    ]

    candidates = (
        set.intersection(
            *intermediate_sets
        )
        if intermediate_sets
        else set()
    )

    occurrence_counts = Counter()

    for path in representative:

        occurrence_counts.update(
            path[
                1:-1
            ]
        )

    rows = []

    for node in sorted(
        candidates
    ):

        still_reachable = (
            exact_reachable(
                G=G,
                source=source,
                target=target,
                seed_graph_genes=seed_graph_genes,
                max_edges=max_edges,
                max_external=max_external,
                blocked_node=node,
            )
        )

        rows.append(
            {
                "Node":
                    node,

                "Node_class":
                    (
                        "F2_internal"
                        if node
                        in seed_graph_genes
                        else
                        "external"
                    ),

                "Representative_occurrences":
                    occurrence_counts[
                        node
                    ],

                "Representative_path_count":
                    len(
                        representative
                    ),

                "Exact_required":
                    not still_reachable,
            }
        )

    return rows


# ============================================================
# GraphML export
# ============================================================

def write_representative_graphml(
    setting: str,
    usage: dict,
    G: nx.DiGraph,
    seed_graph_genes: set[str],
    output_file: Path,
) -> None:

    H = nx.DiGraph()

    for (
        source,
        target
    ), info in usage.items():

        attrs = G[
            source
        ][
            target
        ]

        H.add_node(
            source,
            F2_seed=(
                source
                in seed_graph_genes
            ),
        )

        H.add_node(
            target,
            F2_seed=(
                target
                in seed_graph_genes
            ),
        )

        H.add_edge(
            source,
            target,

            effect=int(
                attrs.get(
                    "effect",
                    0,
                )
            ),

            effect_label=str(
                attrs.get(
                    "effect_label",
                    "",
                )
            ),

            layers=str(
                attrs.get(
                    "layers",
                    "",
                )
            ),

            representative_path_usage=int(
                info[
                    "path_usage"
                ]
            ),

            representative_pair_usage=int(
                len(
                    info[
                        "pairs"
                    ]
                )
            ),
        )

    nx.write_graphml(
        H,
        output_file,
    )


# ============================================================
# Excel formatting
# ============================================================

def format_workbook(
    path: Path,
) -> None:

    wb = load_workbook(
        path
    )

    header_fill = PatternFill(
        fill_type="solid",
        fgColor="D9E1F2",
    )

    for ws in wb.worksheets:

        ws.freeze_panes = (
            "A2"
        )

        for cell in ws[
            1
        ]:

            cell.fill = (
                header_fill
            )

            cell.font = Font(
                bold=True
            )

            cell.alignment = Alignment(
                horizontal="center",
                vertical="center",
            )

        if ws.max_row > 1:

            ws.auto_filter.ref = (
                ws.dimensions
            )

        for cells in ws.columns:

            letter = (
                get_column_letter(
                    cells[
                        0
                    ].column
                )
            )

            max_length = 0

            for cell in cells[
                :100
            ]:

                if cell.value is None:
                    continue

                max_length = max(
                    max_length,
                    len(
                        str(
                            cell.value
                        )
                    ),
                )

            ws.column_dimensions[
                letter
            ].width = min(
                max(
                    max_length + 2,
                    12,
                ),
                55,
            )

    wb.save(
        path
    )


# ============================================================
# Validate inputs
# ============================================================

print()
print(
    "=" * 82
)
print(
    "PHASE F5 — CORRECTED DIRECTED OMNIPATH ANALYSIS"
)
print(
    "=" * 82
)

for required in [
    F2_MASTER_FILE,
    CORE_ASSIGN_FILE,
    EXT_ASSIGN_FILE,
]:

    if not required.exists():

        raise FileNotFoundError(
            f"Required input not found:\n"
            f"{required}"
        )

OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

GRAPH_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

RAW_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Load F2 universe
# ============================================================

print()
print(
    "[INFO] Loading F2 experimental universe..."
)

f2 = pd.read_csv(
    F2_MASTER_FILE
)

f2[
    "Gene"
] = (
    f2[
        "Gene"
    ]
    .map(
        normalize_gene
    )
)

f2 = f2[
    f2[
        "Gene"
    ].notna()
].copy()

seed_genes = set(
    f2[
        "Gene"
    ]
)

print(
    f"[INFO] F2 experimental genes: "
    f"{len(seed_genes):,}"
)

f2_lookup = (
    f2
    .drop_duplicates(
        "Gene"
    )
    .set_index(
        "Gene"
    )
)


# ============================================================
# F4 module annotations
# ============================================================

core_assign = pd.read_csv(
    CORE_ASSIGN_FILE
)

core_assign[
    "Gene"
] = (
    core_assign[
        "Gene"
    ]
    .map(
        normalize_gene
    )
)

core_module_lookup = dict(
    zip(
        core_assign[
            "Gene"
        ],
        core_assign[
            "Module"
        ],
    )
)


ext_assign = pd.read_csv(
    EXT_ASSIGN_FILE
)

ext_assign[
    "Gene"
] = (
    ext_assign[
        "Gene"
    ]
    .map(
        normalize_gene
    )
)

ext_module_lookup = dict(
    zip(
        ext_assign[
            "Gene"
        ],
        ext_assign[
            "Module"
        ],
    )
)


# ============================================================
# Load prior knowledge
# ============================================================

print()
print(
    "[INFO] Loading signaling layer..."
)

signaling_raw = download_table(
    "mouse signaling/activity-flow",
    SIGNALING_PARAMS,
    SIGNALING_RAW,
)

print()
print(
    "[INFO] Loading CollecTRI..."
)

collectri_raw = download_table(
    "mouse CollecTRI",
    COLLECTRI_PARAMS,
    COLLECTRI_RAW,
)

print()
print(
    "[INFO] Loading DoRothEA A/B..."
)

dorothea_raw = download_table(
    "mouse DoRothEA A/B",
    DOROTHEA_PARAMS,
    DOROTHEA_RAW,
)

print()
print(
    "[INFO] Loading other TF-target interactions..."
)

tf_target_raw = download_table(
    "mouse tf_target",
    TF_TARGET_PARAMS,
    TF_TARGET_RAW,
)


# ============================================================
# Prepare and merge edges
# ============================================================

print()
print(
    "[INFO] Cleaning directed interactions..."
)

signaling_edges = prepare_edges(
    signaling_raw,
    "SIGNALING",
)

collectri_edges = prepare_edges(
    collectri_raw,
    "COLLECTRI",
)

dorothea_edges = prepare_edges(
    dorothea_raw,
    "DOROTHEA_AB",
)

tf_target_edges = prepare_edges(
    tf_target_raw,
    "TF_TARGET",
)


print(
    f"[INFO] Directed signaling rows : "
    f"{len(signaling_edges):,}"
)

print(
    f"[INFO] Directed CollecTRI rows  : "
    f"{len(collectri_edges):,}"
)

print(
    f"[INFO] Directed DoRothEA rows   : "
    f"{len(dorothea_edges):,}"
)

print(
    f"[INFO] Directed tf_target rows  : "
    f"{len(tf_target_edges):,}"
)


combined_edges = pd.concat(
    [
        signaling_edges,
        collectri_edges,
        dorothea_edges,
        tf_target_edges,
    ],
    ignore_index=True,
)


aggregated_edges = aggregate_edges(
    combined_edges
)


G = build_graph(
    aggregated_edges
)


print(
    f"[INFO] Unique directed edges    : "
    f"{G.number_of_edges():,}"
)

print(
    f"[INFO] OmniPath graph nodes     : "
    f"{G.number_of_nodes():,}"
)


# ============================================================
# Raw metadata
# ============================================================

raw_metadata = {

    "download_time_utc":
        datetime.now(
            timezone.utc
        ).isoformat(),

    "species":
        SPECIES,

    "signaling_rows":
        len(
            signaling_raw
        ),

    "collectri_rows":
        len(
            collectri_raw
        ),

    "dorothea_rows":
        len(
            dorothea_raw
        ),

    "tf_target_rows":
        len(
            tf_target_raw
        ),

    "signaling_sha256":
        sha256_file(
            SIGNALING_RAW
        ),

    "collectri_sha256":
        sha256_file(
            COLLECTRI_RAW
        ),

    "dorothea_sha256":
        sha256_file(
            DOROTHEA_RAW
        ),

    "tf_target_sha256":
        sha256_file(
            TF_TARGET_RAW
        ),
}


with open(
    RAW_METADATA,
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        raw_metadata,
        handle,
        indent=2,
    )


# ============================================================
# Resolve all anchors
# ============================================================

anchor_rows = []

resolved_sources = defaultdict(
    dict
)

resolved_targets = defaultdict(
    dict
)

resolved_alias_to_canonical = {}


for role, families in [

    (
        "SOURCE",
        SOURCE_FAMILIES,
    ),

    (
        "TARGET",
        TARGET_FAMILIES,
    ),

]:

    for family, genes in families.items():

        for canonical in genes:

            canonical = (
                normalize_gene(
                    canonical
                )
            )

            (
                resolved,
                status,
            ) = resolve_anchor(
                canonical,
                G,
            )

            if role == "SOURCE":

                resolved_sources[
                    family
                ][
                    canonical
                ] = resolved

            else:

                resolved_targets[
                    family
                ][
                    canonical
                ] = resolved

            if (
                resolved is not None
                and
                resolved != canonical
            ):

                resolved_alias_to_canonical[
                    resolved
                ] = canonical

            annotation = (
                f2_lookup.loc[
                    canonical
                ]
                if canonical
                in f2_lookup.index
                else None
            )

            anchor_rows.append(
                {
                    "Role":
                        role,

                    "Family":
                        family,

                    "Configured_gene":
                        canonical,

                    "Resolved_graph_gene":
                        resolved,

                    "Resolution_status":
                        status,

                    "Aliases_tested":
                        ";".join(
                            ANCHOR_ALIASES.get(
                                canonical,
                                [],
                            )
                        ),

                    "In_F2_seed":
                        canonical
                        in seed_genes,

                    "In_OmniPath_graph":
                        resolved
                        is not None,

                    "Direction":
                        (
                            annotation.get(
                                "Direction"
                            )
                            if annotation
                            is not None
                            else None
                        ),

                    "log2FoldChange":
                        (
                            annotation.get(
                                "log2FoldChange"
                            )
                            if annotation
                            is not None
                            else None
                        ),

                    "F2_tier":
                        (
                            annotation.get(
                                "F2_tier"
                            )
                            if annotation
                            is not None
                            else None
                        ),

                    "CORE_module":
                        core_module_lookup.get(
                            canonical
                        ),

                    "EXTENDED_module":
                        ext_module_lookup.get(
                            canonical
                        ),
                }
            )


anchor_audit = pd.DataFrame(
    anchor_rows
)

anchor_audit.to_csv(
    OUTPUT_DIR
    / "anchor_audit.csv",
    index=False,
)


# Any resolved alias corresponding to an experimental F2 gene should
# count as experimentally internal if it occurs inside another route.

seed_graph_genes = set(
    seed_genes
)

for alias, canonical in (
    resolved_alias_to_canonical.items()
):

    if canonical in seed_genes:

        seed_graph_genes.add(
            alias
        )


print()
print(
    "[INFO] Anchor resolution:"
)

for _, row in anchor_audit.iterrows():

    resolved_display = (
        row[
            "Resolved_graph_gene"
        ]
        if pd.notna(
            row[
                "Resolved_graph_gene"
            ]
        )
        else "MISSING"
    )

    print(
        f"       "
        f"{row['Role']:<6} "
        f"{row['Configured_gene']:<10} "
        f"-> "
        f"{resolved_display:<10} "
        f"[{row['Resolution_status']}]"
    )


# ============================================================
# Exact path analysis
# ============================================================

pair_setting_rows = []

representative_path_rows = []

required_node_rows = []

representative_edge_usage = {
    setting:
        defaultdict(
            lambda: {
                "path_usage": 0,
                "pairs": set(),
            }
        )
    for setting
    in SETTING_ORDER
}


print()
print(
    "=" * 82
)
print(
    "EXACT CONSTRAINED REACHABILITY"
)
print(
    "=" * 82
)


for setting_name in SETTING_ORDER:

    max_edges = (
        SETTINGS[
            setting_name
        ][
            "max_edges"
        ]
    )

    max_external = (
        SETTINGS[
            setting_name
        ][
            "max_external"
        ]
    )

    print()
    print(
        f"[{setting_name}] "
        f"max_edges={max_edges}, "
        f"max_external={max_external}"
    )

    reachable_count = 0

    testable_count = 0

    representative_count = 0

    exact_required_count = 0

    for source_family, sources in (
        SOURCE_FAMILIES.items()
    ):

        for source_canonical in sources:

            source_canonical = (
                normalize_gene(
                    source_canonical
                )
            )

            source = (
                resolved_sources[
                    source_family
                ][
                    source_canonical
                ]
            )

            for target_family, targets in (
                TARGET_FAMILIES.items()
            ):

                for target_canonical in targets:

                    target_canonical = (
                        normalize_gene(
                            target_canonical
                        )
                    )

                    target = (
                        resolved_targets[
                            target_family
                        ][
                            target_canonical
                        ]
                    )

                    pair_id = (
                        f"{source_canonical}"
                        f"->{target_canonical}"
                    )

                    if (
                        source is None
                        or
                        target is None
                    ):

                        pair_setting_rows.append(
                            {
                                "Setting":
                                    setting_name,

                                "Max_edges":
                                    max_edges,

                                "Max_external":
                                    max_external,

                                "Source_family":
                                    source_family,

                                "Source":
                                    source_canonical,

                                "Resolved_source":
                                    source,

                                "Target_family":
                                    target_family,

                                "Target":
                                    target_canonical,

                                "Resolved_target":
                                    target,

                                "Testable":
                                    False,

                                "Reachable_exact":
                                    False,

                                "Min_external_exact":
                                    None,

                                "Shortest_edges_any_exact":
                                    None,

                                "Shortest_edges_at_min_external":
                                    None,

                                "Representative_paths":
                                    0,

                                "Representative_state_cap_hit":
                                    False,
                            }
                        )

                        continue

                    testable_count += 1

                    target_states = (
                        constrained_target_states(
                            G=G,

                            source=source,

                            target=target,

                            seed_graph_genes=seed_graph_genes,

                            max_edges=max_edges,

                            max_external=max_external,
                        )
                    )

                    reachable = bool(
                        target_states
                    )

                    if not reachable:

                        pair_setting_rows.append(
                            {
                                "Setting":
                                    setting_name,

                                "Max_edges":
                                    max_edges,

                                "Max_external":
                                    max_external,

                                "Source_family":
                                    source_family,

                                "Source":
                                    source_canonical,

                                "Resolved_source":
                                    source,

                                "Target_family":
                                    target_family,

                                "Target":
                                    target_canonical,

                                "Resolved_target":
                                    target,

                                "Testable":
                                    True,

                                "Reachable_exact":
                                    False,

                                "Min_external_exact":
                                    None,

                                "Shortest_edges_any_exact":
                                    None,

                                "Shortest_edges_at_min_external":
                                    None,

                                "Representative_paths":
                                    0,

                                "Representative_state_cap_hit":
                                    False,
                            }
                        )

                        continue

                    reachable_count += 1

                    min_external = min(
                        external
                        for (
                            edges,
                            external
                        )
                        in target_states
                    )

                    shortest_any = min(
                        edges
                        for (
                            edges,
                            external
                        )
                        in target_states
                    )

                    shortest_at_min_external = min(
                        edges
                        for (
                            edges,
                            external
                        )
                        in target_states
                        if external
                        ==
                        min_external
                    )

                    representative_max_edges = min(
                        max_edges,
                        shortest_at_min_external
                        +
                        REPRESENTATIVE_EDGE_SLACK,
                    )

                    (
                        reps,
                        representative_state_cap_hit,
                    ) = representative_paths(
                        G=G,

                        source=source,

                        target=target,

                        seed_graph_genes=seed_graph_genes,

                        max_edges=representative_max_edges,

                        max_external=min_external,

                        max_paths=MAX_REPRESENTATIVE_PATHS,
                    )

                    representative_count += len(
                        reps
                    )

                    pair_setting_rows.append(
                        {
                            "Setting":
                                setting_name,

                            "Max_edges":
                                max_edges,

                            "Max_external":
                                max_external,

                            "Source_family":
                                source_family,

                            "Source":
                                source_canonical,

                            "Resolved_source":
                                source,

                            "Target_family":
                                target_family,

                            "Target":
                                target_canonical,

                            "Resolved_target":
                                target,

                            "Testable":
                                True,

                            "Reachable_exact":
                                True,

                            "Min_external_exact":
                                min_external,

                            "Shortest_edges_any_exact":
                                shortest_any,

                            "Shortest_edges_at_min_external":
                                shortest_at_min_external,

                            "Representative_paths":
                                len(
                                    reps
                                ),

                            "Representative_state_cap_hit":
                                representative_state_cap_hit,
                        }
                    )

                    # ----------------------------------------
                    # Save representative routes
                    # ----------------------------------------

                    for index, path in enumerate(
                        reps,
                        start=1,
                    ):

                        intermediates = (
                            path[
                                1:-1
                            ]
                        )

                        f2_internal = [
                            node
                            for node
                            in intermediates
                            if node
                            in seed_graph_genes
                        ]

                        external = [
                            node
                            for node
                            in intermediates
                            if node
                            not in seed_graph_genes
                        ]

                        (
                            effect_label,
                            effect_numeric,
                        ) = annotate_path_effect(
                            G,
                            path,
                        )

                        layers = (
                            annotate_path_layers(
                                G,
                                path,
                            )
                        )

                        representative_path_rows.append(
                            {
                                "Setting":
                                    setting_name,

                                "Source_family":
                                    source_family,

                                "Source":
                                    source_canonical,

                                "Resolved_source":
                                    source,

                                "Target_family":
                                    target_family,

                                "Target":
                                    target_canonical,

                                "Resolved_target":
                                    target,

                                "Path_index":
                                    index,

                                "Edges":
                                    len(
                                        path
                                    )
                                    - 1,

                                "External_intermediate_count":
                                    len(
                                        external
                                    ),

                                "F2_internal_intermediate_count":
                                    len(
                                        f2_internal
                                    ),

                                "External_intermediates":
                                    ";".join(
                                        external
                                    ),

                                "F2_internal_intermediates":
                                    ";".join(
                                        f2_internal
                                    ),

                                "Path_layers":
                                    layers,

                                "Path_effect":
                                    effect_label,

                                "Path_effect_numeric":
                                    effect_numeric,

                                "Path":
                                    " -> ".join(
                                        path
                                    ),
                            }
                        )

                        for a, b in zip(
                            path[:-1],
                            path[1:],
                        ):

                            info = (
                                representative_edge_usage[
                                    setting_name
                                ][
                                    (
                                        a,
                                        b,
                                    )
                                ]
                            )

                            info[
                                "path_usage"
                            ] += 1

                            info[
                                "pairs"
                            ].add(
                                pair_id
                            )

                    # ----------------------------------------
                    # EXACT required-node testing
                    # ----------------------------------------

                    exact_rows = (
                        exact_required_nodes(
                            G=G,

                            source=source,

                            target=target,

                            seed_graph_genes=seed_graph_genes,

                            max_edges=max_edges,

                            max_external=max_external,

                            representative=reps,
                        )
                    )

                    for required_row in exact_rows:

                        if (
                            required_row[
                                "Exact_required"
                            ]
                        ):

                            exact_required_count += 1

                        required_node_rows.append(
                            {
                                "Setting":
                                    setting_name,

                                "Source_family":
                                    source_family,

                                "Source":
                                    source_canonical,

                                "Target_family":
                                    target_family,

                                "Target":
                                    target_canonical,

                                **required_row,
                            }
                        )

    print(
        f"[{setting_name}] Testable pairs       : "
        f"{testable_count:,}"
    )

    print(
        f"[{setting_name}] Exactly reachable     : "
        f"{reachable_count:,}"
    )

    print(
        f"[{setting_name}] Representative paths  : "
        f"{representative_count:,}"
    )

    print(
        f"[{setting_name}] Exact required hits    : "
        f"{exact_required_count:,}"
    )


# ============================================================
# DataFrames
# ============================================================

pair_setting = pd.DataFrame(
    pair_setting_rows
)

representative_paths_df = pd.DataFrame(
    representative_path_rows
)

required_nodes = pd.DataFrame(
    required_node_rows
)


pair_setting.to_csv(
    OUTPUT_DIR
    / "pair_setting_summary.csv",
    index=False,
)

representative_paths_df.to_csv(
    OUTPUT_DIR
    / "representative_paths.csv",
    index=False,
)

required_nodes.to_csv(
    OUTPUT_DIR
    / "exact_required_nodes.csv",
    index=False,
)


# ============================================================
# Pair robustness / mechanistic support tiers
# ============================================================

pair_robustness_rows = []


pair_columns = [
    "Source_family",
    "Source",
    "Target_family",
    "Target",
]


for keys, group in pair_setting.groupby(
    pair_columns,
    dropna=False,
):

    (
        source_family,
        source,
        target_family,
        target,
    ) = keys

    group_by_setting = {
        row[
            "Setting"
        ]:
            row
        for _, row
        in group.iterrows()
    }

    testable = bool(
        group[
            "Testable"
        ].any()
    )

    tier = (
        "UNTESTABLE"
    )

    first_setting = None

    if testable:

        tier = (
            "UNREACHABLE"
        )

        for setting in SETTING_ORDER:

            row = (
                group_by_setting.get(
                    setting
                )
            )

            if (
                row is not None
                and
                bool(
                    row[
                        "Reachable_exact"
                    ]
                )
            ):

                first_setting = (
                    setting
                )

                tier = (
                    SUPPORT_TIER[
                        setting
                    ]
                )

                break

    reachable_rows = group[
        group[
            "Reachable_exact"
        ]
    ]

    if reachable_rows.empty:

        min_external = None
        shortest = None

    else:

        min_external = (
            pd.to_numeric(
                reachable_rows[
                    "Min_external_exact"
                ],
                errors="coerce",
            )
            .min()
        )

        best_external_rows = (
            reachable_rows[
                pd.to_numeric(
                    reachable_rows[
                        "Min_external_exact"
                    ],
                    errors="coerce",
                )
                ==
                min_external
            ]
        )

        shortest = (
            pd.to_numeric(
                best_external_rows[
                    "Shortest_edges_at_min_external"
                ],
                errors="coerce",
            )
            .min()
        )

    pair_robustness_rows.append(
        {
            "Source_family":
                source_family,

            "Source":
                source,

            "Target_family":
                target_family,

            "Target":
                target,

            "Testable":
                testable,

            "Support_tier":
                tier,

            "Strictest_reachable_setting":
                first_setting,

            "Reachable_E4_X0":
                bool(
                    group_by_setting.get(
                        "E4_X0",
                        {}
                    ).get(
                        "Reachable_exact",
                        False,
                    )
                ),

            "Reachable_E4_X1":
                bool(
                    group_by_setting.get(
                        "E4_X1",
                        {}
                    ).get(
                        "Reachable_exact",
                        False,
                    )
                ),

            "Reachable_E4_X2":
                bool(
                    group_by_setting.get(
                        "E4_X2",
                        {}
                    ).get(
                        "Reachable_exact",
                        False,
                    )
                ),

            "Reachable_E5_X2":
                bool(
                    group_by_setting.get(
                        "E5_X2",
                        {}
                    ).get(
                        "Reachable_exact",
                        False,
                    )
                ),

            "Minimum_external_any_setting":
                min_external,

            "Shortest_edges_at_min_external":
                shortest,
        }
    )


pair_robustness = pd.DataFrame(
    pair_robustness_rows
)


tier_rank = {
    "A": 0,
    "B": 1,
    "C": 2,
    "D": 3,
    "UNREACHABLE": 4,
    "UNTESTABLE": 5,
}


pair_robustness[
    "_tier_rank"
] = (
    pair_robustness[
        "Support_tier"
    ]
    .map(
        tier_rank
    )
)


pair_robustness = (
    pair_robustness
    .sort_values(
        by=[
            "_tier_rank",
            "Minimum_external_any_setting",
            "Shortest_edges_at_min_external",
            "Source",
            "Target",
        ],
        ascending=[
            True,
            True,
            True,
            True,
            True,
        ],
        na_position="last",
    )
    .drop(
        columns=[
            "_tier_rank",
        ]
    )
    .reset_index(
        drop=True
    )
)


pair_robustness.to_csv(
    OUTPUT_DIR
    / "pair_robustness.csv",
    index=False,
)


# ============================================================
# Family summary
# ============================================================

family_rows = []


for (
    source_family,
    target_family
), group in pair_robustness.groupby(
    [
        "Source_family",
        "Target_family",
    ]
):

    family_rows.append(
        {
            "Source_family":
                source_family,

            "Target_family":
                target_family,

            "Configured_pairs":
                len(
                    group
                ),

            "Testable_pairs":
                int(
                    group[
                        "Testable"
                    ]
                    .sum()
                ),

            "Tier_A":
                int(
                    (
                        group[
                            "Support_tier"
                        ]
                        ==
                        "A"
                    )
                    .sum()
                ),

            "Tier_B":
                int(
                    (
                        group[
                            "Support_tier"
                        ]
                        ==
                        "B"
                    )
                    .sum()
                ),

            "Tier_C":
                int(
                    (
                        group[
                            "Support_tier"
                        ]
                        ==
                        "C"
                    )
                    .sum()
                ),

            "Tier_D":
                int(
                    (
                        group[
                            "Support_tier"
                        ]
                        ==
                        "D"
                    )
                    .sum()
                ),

            "Unreachable":
                int(
                    (
                        group[
                            "Support_tier"
                        ]
                        ==
                        "UNREACHABLE"
                    )
                    .sum()
                ),

            "Untestable":
                int(
                    (
                        group[
                            "Support_tier"
                        ]
                        ==
                        "UNTESTABLE"
                    )
                    .sum()
                ),
        }
    )


family_summary = pd.DataFrame(
    family_rows
)

family_summary.to_csv(
    OUTPUT_DIR
    / "family_summary.csv",
    index=False,
)


# ============================================================
# Exact node robustness
# ============================================================

node_rows = []


if not required_nodes.empty:

    exact_only = required_nodes[
        required_nodes[
            "Exact_required"
        ]
    ].copy()

    if not exact_only.empty:

        exact_only[
            "_pair"
        ] = (
            exact_only[
                "Source"
            ]
            +
            "->"
            +
            exact_only[
                "Target"
            ]
        )

        for node, group in exact_only.groupby(
            "Node"
        ):

            canonical_annotation_gene = (
                resolved_alias_to_canonical.get(
                    node,
                    node,
                )
            )

            annotation = (
                f2_lookup.loc[
                    canonical_annotation_gene
                ]
                if canonical_annotation_gene
                in f2_lookup.index
                else None
            )

            node_rows.append(
                {
                    "Node":
                        node,

                    "Annotation_gene":
                        canonical_annotation_gene,

                    "Node_class":
                        (
                            "F2_internal"
                            if node
                            in seed_graph_genes
                            else
                            "external"
                        ),

                    "Distinct_pairs_exact_required":
                        int(
                            group[
                                "_pair"
                            ]
                            .nunique()
                        ),

                    "Settings_exact_required":
                        int(
                            group[
                                "Setting"
                            ]
                            .nunique()
                        ),

                    "Pair_setting_exact_required":
                        len(
                            group
                        ),

                    "Required_in_E4_X0":
                        bool(
                            (
                                group[
                                    "Setting"
                                ]
                                ==
                                "E4_X0"
                            )
                            .any()
                        ),

                    "Required_in_E4_X1":
                        bool(
                            (
                                group[
                                    "Setting"
                                ]
                                ==
                                "E4_X1"
                            )
                            .any()
                        ),

                    "Required_in_E4_X2":
                        bool(
                            (
                                group[
                                    "Setting"
                                ]
                                ==
                                "E4_X2"
                            )
                            .any()
                        ),

                    "Required_in_E5_X2":
                        bool(
                            (
                                group[
                                    "Setting"
                                ]
                                ==
                                "E5_X2"
                            )
                            .any()
                        ),

                    "Direction":
                        (
                            annotation.get(
                                "Direction"
                            )
                            if annotation
                            is not None
                            else None
                        ),

                    "log2FoldChange":
                        (
                            annotation.get(
                                "log2FoldChange"
                            )
                            if annotation
                            is not None
                            else None
                        ),

                    "F2_tier":
                        (
                            annotation.get(
                                "F2_tier"
                            )
                            if annotation
                            is not None
                            else None
                        ),

                    "CORE_module":
                        core_module_lookup.get(
                            canonical_annotation_gene
                        ),

                    "EXTENDED_module":
                        ext_module_lookup.get(
                            canonical_annotation_gene
                        ),
                }
            )


node_robustness = pd.DataFrame(
    node_rows
)


if not node_robustness.empty:

    node_robustness = (
        node_robustness
        .sort_values(
            by=[
                "Required_in_E4_X0",
                "Settings_exact_required",
                "Distinct_pairs_exact_required",
                "Pair_setting_exact_required",
            ],
            ascending=[
                False,
                False,
                False,
                False,
            ],
        )
        .reset_index(
            drop=True
        )
    )


node_robustness.to_csv(
    OUTPUT_DIR
    / "node_robustness.csv",
    index=False,
)


# ============================================================
# Representative union edges + GraphML
# ============================================================

union_rows = []


for setting_name in SETTING_ORDER:

    usage_dict = (
        representative_edge_usage[
            setting_name
        ]
    )

    for (
        source,
        target
    ), info in usage_dict.items():

        attrs = (
            G[
                source
            ][
                target
            ]
        )

        union_rows.append(
            {
                "Setting":
                    setting_name,

                "Source":
                    source,

                "Target":
                    target,

                "Source_F2_seed":
                    source
                    in seed_graph_genes,

                "Target_F2_seed":
                    target
                    in seed_graph_genes,

                "Effect":
                    attrs.get(
                        "effect"
                    ),

                "Effect_label":
                    attrs.get(
                        "effect_label"
                    ),

                "Layers":
                    attrs.get(
                        "layers"
                    ),

                "Max_curation_effort":
                    attrs.get(
                        "max_curation_effort"
                    ),

                "Supporting_records":
                    attrs.get(
                        "supporting_records"
                    ),

                "Representative_path_usage":
                    info[
                        "path_usage"
                    ],

                "Representative_pair_usage":
                    len(
                        info[
                            "pairs"
                        ]
                    ),
            }
        )

    write_representative_graphml(
        setting=setting_name,

        usage=usage_dict,

        G=G,

        seed_graph_genes=seed_graph_genes,

        output_file=(
            GRAPH_DIR
            /
            f"{setting_name}.graphml"
        ),
    )


union_edges = pd.DataFrame(
    union_rows
)

union_edges.to_csv(
    OUTPUT_DIR
    / "representative_union_edges.csv",
    index=False,
)


# ============================================================
# Setting summary
# ============================================================

setting_rows = []


for setting_name in SETTING_ORDER:

    subset = pair_setting[
        pair_setting[
            "Setting"
        ]
        ==
        setting_name
    ]

    reps = (
        representative_paths_df[
            representative_paths_df[
                "Setting"
            ]
            ==
            setting_name
        ]
        if not representative_paths_df.empty
        else pd.DataFrame()
    )

    exact_required = (
        required_nodes[
            (
                required_nodes[
                    "Setting"
                ]
                ==
                setting_name
            )
            &
            (
                required_nodes[
                    "Exact_required"
                ]
            )
        ]
        if not required_nodes.empty
        else pd.DataFrame()
    )

    setting_rows.append(
        {
            "Setting":
                setting_name,

            "Max_edges":
                SETTINGS[
                    setting_name
                ][
                    "max_edges"
                ],

            "Max_external":
                SETTINGS[
                    setting_name
                ][
                    "max_external"
                ],

            "Configured_pairs":
                len(
                    subset
                ),

            "Testable_pairs":
                int(
                    subset[
                        "Testable"
                    ]
                    .sum()
                ),

            "Exactly_reachable_pairs":
                int(
                    subset[
                        "Reachable_exact"
                    ]
                    .sum()
                ),

            "Representative_paths_retained":
                len(
                    reps
                ),

            "Representative_search_caps_hit":
                int(
                    subset[
                        "Representative_state_cap_hit"
                    ]
                    .sum()
                ),

            "Exact_required_pair_node_hits":
                len(
                    exact_required
                ),
        }
    )


setting_summary = pd.DataFrame(
    setting_rows
)

setting_summary.to_csv(
    OUTPUT_DIR
    / "setting_summary.csv",
    index=False,
)


# ============================================================
# Excel workbook
# ============================================================

with pd.ExcelWriter(
    OUTPUT_XLSX,
    engine="openpyxl",
) as writer:

    setting_summary.to_excel(
        writer,
        sheet_name="Setting summary",
        index=False,
    )

    family_summary.to_excel(
        writer,
        sheet_name="Family summary",
        index=False,
    )

    anchor_audit.to_excel(
        writer,
        sheet_name="Anchor audit",
        index=False,
    )

    pair_robustness.to_excel(
        writer,
        sheet_name="Pair robustness",
        index=False,
    )

    pair_setting.to_excel(
        writer,
        sheet_name="Pair by setting",
        index=False,
    )

    node_robustness.to_excel(
        writer,
        sheet_name="Node robustness",
        index=False,
    )

    required_nodes.to_excel(
        writer,
        sheet_name="Required node tests",
        index=False,
    )

    representative_paths_df.to_excel(
        writer,
        sheet_name="Representative paths",
        index=False,
    )

    union_edges.to_excel(
        writer,
        sheet_name="Representative edges",
        index=False,
    )


format_workbook(
    OUTPUT_XLSX
)


# ============================================================
# Metadata
# ============================================================

metadata = {

    "phase":
        "F5_corrected_omnipath_downstream",

    "organism":
        {
            "name":
                "Mus musculus",

            "ncbi_taxon":
                SPECIES,
        },

    "F2_seed_gene_count":
        len(
            seed_genes
        ),

    "graph_seed_gene_count_after_alias_resolution":
        len(
            seed_graph_genes
        ),

    "source_families":
        SOURCE_FAMILIES,

    "target_families":
        TARGET_FAMILIES,

    "anchor_aliases":
        ANCHOR_ALIASES,

    "settings":
        SETTINGS,

    "support_tiers":
        {
            "A":
                "Reachable under E4_X0",

            "B":
                "First reachable under E4_X1",

            "C":
                "First reachable under E4_X2",

            "D":
                "First reachable under E5_X2",
        },

    "representative_path_limit":
        MAX_REPRESENTATIVE_PATHS,

    "representative_edge_slack":
        REPRESENTATIVE_EDGE_SLACK,

    "important_method_note":
        (
            "Reachability and required-node testing are exact under "
            "the bounded E/X state-space formulation. Representative "
            "paths are capped and are retained only for interpretation "
            "and visualization. Path count is not treated as mechanistic "
            "evidence."
        ),

    "required_node_definition":
        (
            "A candidate intermediate is exact-required if deleting "
            "that node eliminates every source-target route satisfying "
            "the full E/X constraint."
        ),

    "network":
        {
            "nodes":
                G.number_of_nodes(),

            "edges":
                G.number_of_edges(),

            "directed_signaling_rows":
                len(
                    signaling_edges
                ),

            "directed_collectri_rows":
                len(
                    collectri_edges
                ),

            "directed_dorothea_AB_rows":
                len(
                    dorothea_edges
                ),

            "directed_tf_target_rows":
                len(
                    tf_target_edges
                ),
        },

    "raw_metadata":
        raw_metadata,
}


with open(
    OUTPUT_METADATA,
    "w",
    encoding="utf-8",
) as handle:

    json.dump(
        metadata,
        handle,
        indent=2,
        ensure_ascii=False,
        default=str,
    )


# ============================================================
# Console report
# ============================================================

print()
print(
    "=" * 82
)
print(
    "PHASE F5 CORRECTED COMPLETE"
)
print(
    "=" * 82
)

print()
print(
    "Directed prior-knowledge network"
)

print(
    f"    Nodes                         : "
    f"{G.number_of_nodes():,}"
)

print(
    f"    Unique directed edges         : "
    f"{G.number_of_edges():,}"
)

print(
    f"    F2 experimental genes         : "
    f"{len(seed_genes):,}"
)


print()
print(
    "Exact reachability"
)

for _, row in setting_summary.iterrows():

    print(
        f"    "
        f"{row['Setting']:<7} "
        f"reachable="
        f"{int(row['Exactly_reachable_pairs']):>3}/"
        f"{int(row['Testable_pairs']):<3}  "
        f"representative_paths="
        f"{int(row['Representative_paths_retained']):>4}  "
        f"exact_required_hits="
        f"{int(row['Exact_required_pair_node_hits']):>3}"
    )


print()
print(
    "Mechanistic support tiers"
)

for tier in [
    "A",
    "B",
    "C",
    "D",
    "UNREACHABLE",
    "UNTESTABLE",
]:

    count = int(
        (
            pair_robustness[
                "Support_tier"
            ]
            ==
            tier
        )
        .sum()
    )

    print(
        f"    {tier:<11}: "
        f"{count:>3}"
    )


print()
print(
    "Top Tier A pairs "
    "(<=4 edges, no external intermediates)"
)

tier_a = (
    pair_robustness[
        pair_robustness[
            "Support_tier"
        ]
        ==
        "A"
    ]
    .head(
        30
    )
)

if tier_a.empty:

    print(
        "    None."
    )

else:

    for _, row in tier_a.iterrows():

        print(
            f"    "
            f"{row['Source']:<10} "
            f"-> "
            f"{row['Target']:<10} "
            f""
            f"({row['Source_family']} -> "
            f"{row['Target_family']})  "
            f"edges="
            f"{int(row['Shortest_edges_at_min_external'])}"
        )


print()
print(
    "Exact required intermediates"
)

if node_robustness.empty:

    print(
        "    None detected."
    )

else:

    for _, row in (
        node_robustness
        .head(
            25
        )
        .iterrows()
    ):

        print(
            f"    "
            f"{row['Node']:<12} "
            f"class="
            f"{row['Node_class']:<11}  "
            f"pairs="
            f"{int(row['Distinct_pairs_exact_required']):>2}  "
            f"settings="
            f"{int(row['Settings_exact_required'])}  "
            f"E4_X0="
            f"{bool(row['Required_in_E4_X0'])}"
        )


print()
print(
    "Anchor coverage"
)

missing = anchor_audit[
    ~anchor_audit[
        "In_OmniPath_graph"
    ]
]

alias_resolved = anchor_audit[
    anchor_audit[
        "Resolution_status"
    ]
    .isin(
        [
            "alias",
            "multiple_aliases_first_used",
        ]
    )
]


print(
    f"    Alias-resolved anchors        : "
    f"{len(alias_resolved):,}"
)

print(
    f"    Remaining missing anchors     : "
    f"{len(missing):,}"
)


if not alias_resolved.empty:

    print()
    print(
        "    Alias resolutions:"
    )

    for _, row in (
        alias_resolved.iterrows()
    ):

        print(
            f"        "
            f"{row['Configured_gene']:<10} "
            f"-> "
            f"{row['Resolved_graph_gene']}"
        )


if not missing.empty:

    print()
    print(
        "    Still absent from prior-knowledge graph:"
    )

    for _, row in (
        missing.iterrows()
    ):

        print(
            f"        "
            f"{row['Role']:<6} "
            f"{row['Family']:<22} "
            f"{row['Configured_gene']}"
        )


print()
print(
    "Output directory:"
)

print(
    f"    "
    f"{OUTPUT_DIR.relative_to(PROJECT_ROOT)}"
)

print()
print(
    "[NEXT] Interpret support tiers, exact minimum-X routes, "
    "and exact required intermediates before F6 integration."
)