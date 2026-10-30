#!/usr/bin/env python3

"""
F7 — Experimental candidate prioritization
===========================================

This is a reporting / experimental-design step, NOT a new discovery
algorithm and NOT a numerical ranking.

It extracts a deliberately selected mechanistic validation panel from:

    Data/Processed/FongDownstream/F6_Integrated/
        experimental_gene_evidence_matrix.csv

and combines the exact F6 evidence with manually specified biological
questions and discriminating experiments.

Outputs
-------
Data/Processed/FongDownstream/F7_Prioritization/

    F7_candidate_validation_priorities.csv
    F7_candidate_validation_priorities.xlsx
"""

from pathlib import Path

import pandas as pd

from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill, Alignment
from openpyxl.utils import get_column_letter


# ============================================================
# Paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F6_Integrated"
    / "experimental_gene_evidence_matrix.csv"
)

OUTPUT_DIR = (
    PROJECT_ROOT
    / "Data"
    / "Processed"
    / "FongDownstream"
    / "F7_Prioritization"
)

OUTPUT_CSV = (
    OUTPUT_DIR
    / "F7_candidate_validation_priorities.csv"
)

OUTPUT_XLSX = (
    OUTPUT_DIR
    / "F7_candidate_validation_priorities.xlsx"
)


# ============================================================
# Deliberately selected candidates
#
# This ordering is conceptual, not a score.
# ============================================================

CANDIDATES = [

    {
        "Gene": "CLSPN",
        "Validation_group": "Checkpoint / replication-stress",
        "Experimental_role": "Mechanistic lever",
        "Why_prioritize":
            (
                "Represents the newly recovered replication-stress/"
                "checkpoint interface and provides an experimentally "
                "tractable entry point into the proposed stress-to-fate "
                "transition."
            ),
        "Primary_experiment":
            (
                "Reduce CLSPN genetically after Rb-family loss; measure "
                "pCHK1, gamma-H2AX, apoptotic markers, pMLKL and viability."
            ),
        "Discriminating_result":
            (
                "Reduced death would support a death-transmitting checkpoint; "
                "increased death would instead indicate a protective "
                "replication-stress response."
            ),
    },

    {
        "Gene": "TOPBP1",
        "Validation_group": "Checkpoint / replication-stress",
        "Experimental_role": "Mechanistic anchor",
        "Why_prioritize":
            (
                "Replication-stress module member; corrected F5 identified "
                "a strict experimentally anchored TOPBP1-to-PYCARD route "
                "for which E2F1 was an exact required intermediate."
            ),
        "Primary_experiment":
            (
                "Reduce TOPBP1 and measure E2F1, PYCARD, checkpoint activity "
                "and downstream death phenotypes."
            ),
        "Discriminating_result":
            (
                "Loss of downstream E2F1/PYCARD/death activation would "
                "support a replication-stress-to-fate connection."
            ),
    },

    {
        "Gene": "CHEK1",
        "Validation_group": "Checkpoint / replication-stress",
        "Experimental_role": "Mechanistic lever",
        "Why_prioritize":
            (
                "Central checkpoint transducer with multiple short Fong-"
                "internal directed routes toward death-receptor and "
                "regulated-death targets."
            ),
        "Primary_experiment":
            (
                "Genetically reduce CHEK1 and measure DNA-damage burden, "
                "FAS/TNFRSF10B, MLKL and cell death."
            ),
        "Discriminating_result":
            (
                "Direction of the viability effect distinguishes a "
                "protective checkpoint from a checkpoint that contributes "
                "to death signaling."
            ),
    },

    {
        "Gene": "H2AFX",
        "Validation_group": "Checkpoint / replication-stress",
        "Experimental_role": "Temporal marker",
        "Why_prioritize":
            (
                "Provides a direct conceptual bridge to gamma-H2AX/DNA-"
                "damage observations following Rb loss and independently "
                "emerged in the new stress-interface analysis."
            ),
        "Primary_experiment":
            (
                "Perform a time course of gamma-H2AX relative to cell-cycle "
                "entry, pCHK1, downstream fate genes and cell death."
            ),
        "Discriminating_result":
            (
                "DNA damage appearing before fate-effector activation would "
                "support the proposed ordering of the model."
            ),
    },

    {
        "Gene": "E2F1",
        "Validation_group": "Proximal Rb/E2F interface",
        "Experimental_role": "Upstream mechanistic lever",
        "Why_prioritize":
            (
                "Connects the primary Rb/E2F lesion to the newly recovered "
                "downstream architecture and is the exact required node in "
                "the strict TOPBP1-to-PYCARD F5 route."
            ),
        "Primary_experiment":
            (
                "Reduce E2F1 after Rb-family loss and assay MYC, DDIT3, "
                "FAS, TNFRSF10B, PMAIP1, BBC3 and death."
            ),
        "Discriminating_result":
            (
                "Selective loss of downstream branches would establish "
                "which portions of the response are genuinely E2F1-dependent."
            ),
    },

    {
        "Gene": "MYC",
        "Validation_group": "Stress/cycle-to-fate interface",
        "Experimental_role": "Mechanistic lever",
        "Why_prioritize":
            (
                "Strong CORE cross-module connector and one of the very few "
                "genes with direct signed F5 relationships to fate targets."
            ),
        "Primary_experiment":
            (
                "Reduce MYC and measure FAS and TNFRSF10B expression, "
                "followed by receptor-associated death readouts."
            ),
        "Discriminating_result":
            (
                "The prior-knowledge model predicts opposing receptor effects: "
                "loss of MYC should relieve repression of FAS while reducing "
                "positive regulation of TNFRSF10B."
            ),
    },

    {
        "Gene": "DDIT3",
        "Validation_group": "Stress/cycle-to-fate interface",
        "Experimental_role": "Mechanistic lever",
        "Why_prioritize":
            (
                "Stress-response transcription factor repeatedly appearing "
                "in short checkpoint/cycle-to-death-receptor routes."
            ),
        "Primary_experiment":
            (
                "Reduce DDIT3 and measure FAS, TNFRSF10B, caspase-8 "
                "activation and survival."
            ),
        "Discriminating_result":
            (
                "Loss of death-receptor induction and/or rescue of death "
                "would support DDIT3 as part of the missing stress-to-fate "
                "transition."
            ),
    },

    {
        "Gene": "FAS",
        "Validation_group": "Extrinsic / receptor fate",
        "Experimental_role": "Fate discriminator",
        "Why_prioritize":
            (
                "Strong cross-layer candidate with developmental-context "
                "evidence, a stress/death-module position, and extensive "
                "directed connectivity."
            ),
        "Primary_experiment":
            (
                "Reduce or block FAS signaling and measure survival, "
                "caspase-8/caspase-3 activation and alternative death modes."
            ),
        "Discriminating_result":
            (
                "Improved survival supports a pro-death role; worsened "
                "survival would instead support context-dependent protective "
                "FAS signaling."
            ),
    },

    {
        "Gene": "TNFRSF10B",
        "Validation_group": "Extrinsic / receptor fate",
        "Experimental_role": "Fate discriminator",
        "Why_prioritize":
            (
                "Direct positively signed MYC target in corrected F5 and "
                "therefore a comparatively clean extrinsic-death prediction."
            ),
        "Primary_experiment":
            (
                "Reduce TNFRSF10B and assay caspase-8, caspase-3 and "
                "overall survival."
            ),
        "Discriminating_result":
            (
                "Rescue from death would directly support a death-receptor "
                "execution branch."
            ),
    },

    {
        "Gene": "PMAIP1",
        "Validation_group": "Intrinsic apoptosis",
        "Experimental_role": "Fate effector",
        "Why_prioritize":
            (
                "NOXA is strongly altered in Fong and recovered in the "
                "E2F/G1-S apoptotic-effector module despite being absent "
                "from the directed OmniPath target graph."
            ),
        "Primary_experiment":
            (
                "Reduce PMAIP1 alone and together with BBC3; measure "
                "mitochondrial injury, BAX activation, caspase-9 and survival."
            ),
        "Discriminating_result":
            (
                "Rescue—especially after combined PMAIP1/BBC3 reduction—"
                "would support mitochondrial apoptosis as a true execution arm."
            ),
    },

    {
        "Gene": "BBC3",
        "Validation_group": "Intrinsic apoptosis",
        "Experimental_role": "Fate effector",
        "Why_prioritize":
            (
                "PUMA independently emerged in the same E2F/G1-S fate "
                "module and provides a complementary intrinsic-apoptosis test."
            ),
        "Primary_experiment":
            (
                "Reduce BBC3 alone and jointly with PMAIP1; assay "
                "mitochondrial and caspase-9-dependent death."
            ),
        "Discriminating_result":
            (
                "A stronger double-knockdown rescue would suggest partially "
                "redundant or cooperative BH3-only execution."
            ),
    },

    {
        "Gene": "ZBP1",
        "Validation_group": "Regulated inflammatory death",
        "Experimental_role": "Branch initiator",
        "Why_prioritize":
            (
                "Strongly altered Fong gene clustering with PYCARD/MLKL in "
                "the regulated-death module; its absence from the OmniPath "
                "target graph makes experimental validation particularly "
                "informative."
            ),
        "Primary_experiment":
            (
                "Reduce ZBP1 and assay RIPK3/MLKL activation, membrane "
                "integrity and overall survival."
            ),
        "Discriminating_result":
            (
                "Loss of MLKL activation and rescue of death would support "
                "an innate regulated-death branch."
            ),
    },

    {
        "Gene": "MLKL",
        "Validation_group": "Regulated inflammatory death",
        "Experimental_role": "Terminal fate effector",
        "Why_prioritize":
            (
                "Terminal regulated-necrosis effector with multiple Fong-"
                "internal directed routes from the checkpoint/cell-cycle state."
            ),
        "Primary_experiment":
            (
                "Measure phospho-MLKL and genetically reduce MLKL while "
                "monitoring caspase-dependent and caspase-independent death."
            ),
        "Discriminating_result":
            (
                "MLKL activation plus rescue after MLKL reduction would be "
                "direct evidence that regulated necrotic execution occurs."
            ),
    },

    {
        "Gene": "BDNF",
        "Validation_group": "Survival / developmental context",
        "Experimental_role": "Counter-program discriminator",
        "Why_prioritize":
            (
                "Recurrent old/new candidate that represents the trophic "
                "counter-program and preserves developmental-context evidence."
            ),
        "Primary_experiment":
            (
                "Attenuate BDNF/TrkB signaling and, separately, test "
                "recombinant BDNF rescue after Rb-family loss."
            ),
        "Discriminating_result":
            (
                "Worsening after pathway loss and rescue by BDNF would show "
                "that part of the transcriptional response is compensatory."
            ),
    },

    {
        "Gene": "KDR",
        "Validation_group": "Survival / developmental context",
        "Experimental_role": "Counter-program discriminator",
        "Why_prioritize":
            (
                "Represents a second trophic axis and is one of the few "
                "direct positively signed E2F1 target relationships in F5."
            ),
        "Primary_experiment":
            (
                "First test whether E2F1 perturbation changes KDR expression; "
                "then perturb KDR signaling and assay survival."
            ),
        "Discriminating_result":
            (
                "E2F1-dependent KDR expression followed by a survival effect "
                "would connect the proximal Rb/E2F lesion to an adaptive "
                "trophic response."
            ),
    },
]


# ============================================================
# Load exact F6 evidence
# ============================================================

if not INPUT_FILE.exists():

    raise FileNotFoundError(
        f"F6 evidence matrix not found:\n"
        f"{INPUT_FILE}"
    )


OUTPUT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


evidence = pd.read_csv(
    INPUT_FILE,
    low_memory=False,
)


evidence[
    "Gene"
] = (
    evidence[
        "Gene"
    ]
    .astype(str)
    .str.strip()
    .str.upper()
)


candidate_meta = pd.DataFrame(
    CANDIDATES
)


candidate_meta[
    "Gene"
] = (
    candidate_meta[
        "Gene"
    ]
    .str.upper()
)


# ============================================================
# Evidence columns to retain
# ============================================================

requested_evidence_columns = [

    "Gene",

    "Integration_class",
    "Evidence_layers",
    "Candidate_panels",

    "Old_master",
    "Fong_strict_DE",

    "Fong_direction",
    "Fong_log2FC",
    "Fong_padj",

    "F2_tier",
    "F2_categories",

    "Old_adult_embryo_shared",
    "Old_developmental_reversal",
    "Old_Andrusiak_overlap",

    "Julian_E2F4_vs_E2F3_sig_promoter",
    "Julian_binding_expression_overlap",
    "Julian_significant_binding_and_expression",

    "F4_CORE_module",
    "F4_CORE_module_label",
    "F4_CORE_external_degree",
    "F4_CORE_external_module_count",
    "F4_CORE_participation",
    "F4_CORE_betweenness",

    "F4_EXT_module",
    "F4_EXT_module_label",
    "F4_EXT_external_degree",
    "F4_EXT_external_module_count",
    "F4_EXT_participation",
    "F4_EXT_betweenness",

    "F4_interface_CORE_neighbor_count",
    "F4_interface_CORE_modules_touched_count",

    "F5_configured_source",
    "F5_configured_target",
    "F5_target_in_prior_graph",

    "F5_TierA_outgoing_count",
    "F5_TierA_targets",
    "F5_TierA_incoming_count",
    "F5_TierA_sources",

    "F5_direct_outgoing",
    "F5_direct_incoming",

    "F5_exact_required_pair_count",
    "F5_exact_required_pairs",
]


available_columns = [
    column
    for column
    in requested_evidence_columns
    if column
    in evidence.columns
]


selected_evidence = evidence[
    available_columns
].copy()


# ============================================================
# Merge manual interpretation with exact evidence
# ============================================================

final = candidate_meta.merge(
    selected_evidence,
    on="Gene",
    how="left",
)


missing = final[
    final[
        "Integration_class"
    ].isna()
][
    "Gene"
].tolist()


if missing:

    print()
    print(
        "[WARNING] Candidate genes absent from F6 matrix:"
    )

    for gene in missing:
        print(
            f"    {gene}"
        )


# ============================================================
# Reorder columns for readability
# ============================================================

front_columns = [

    "Gene",
    "Validation_group",
    "Experimental_role",

    "Why_prioritize",
    "Primary_experiment",
    "Discriminating_result",

    "Integration_class",
    "Evidence_layers",

    "Fong_direction",
    "Fong_log2FC",
    "Fong_padj",

    "F2_tier",
    "F2_categories",

    "Julian_E2F4_vs_E2F3_sig_promoter",

    "Old_adult_embryo_shared",
    "Old_developmental_reversal",
    "Old_Andrusiak_overlap",

    "F4_CORE_module_label",
    "F4_EXT_module_label",

    "F5_direct_outgoing",
    "F5_direct_incoming",

    "F5_exact_required_pair_count",
    "F5_exact_required_pairs",
]


front_columns = [
    column
    for column
    in front_columns
    if column
    in final.columns
]


remaining_columns = [
    column
    for column
    in final.columns
    if column
    not in front_columns
]


final = final[
    front_columns
    +
    remaining_columns
]


# ============================================================
# Write CSV
# ============================================================

final.to_csv(
    OUTPUT_CSV,
    index=False,
)


# ============================================================
# Experimental-program summary
# ============================================================

program = pd.DataFrame(
    [

        {
            "Stage":
                "1 — Temporal validation",

            "Question":
                (
                    "Does replication/checkpoint stress precede "
                    "fate-effector activation?"
                ),

            "Core_readouts":
                (
                    "Cell-cycle entry; CLSPN/pCHK1; gamma-H2AX; "
                    "MYC/DDIT3; FAS/TNFRSF10B; NOXA/PUMA; "
                    "pMLKL; viability/death"
                ),

            "Decision":
                (
                    "Establish ordering before assigning upstream/"
                    "downstream causal roles."
                ),
        },

        {
            "Stage":
                "2 — Checkpoint perturbation",

            "Question":
                (
                    "Is the replication-stress/checkpoint response "
                    "protective or death-promoting?"
                ),

            "Core_readouts":
                (
                    "CLSPN; CHEK1; TOPBP1 perturbation with "
                    "DNA-damage and fate readouts"
                ),

            "Decision":
                (
                    "Direction of viability change distinguishes "
                    "protective repair from death-transmitting signaling."
                ),
        },

        {
            "Stage":
                "3 — Stress-to-receptor branch",

            "Question":
                (
                    "Does checkpoint/cell-cycle stress feed an "
                    "extrinsic death program?"
                ),

            "Core_readouts":
                (
                    "E2F1; MYC; DDIT3; FAS; TNFRSF10B; "
                    "caspase-8/caspase-3"
                ),

            "Decision":
                (
                    "Perturb upstream factors and receptors to test "
                    "directionality and functional necessity."
                ),
        },

        {
            "Stage":
                "4 — Intrinsic versus regulated-necrotic execution",

            "Question":
                (
                    "Which terminal death mechanism actually accounts "
                    "for loss of viability?"
                ),

            "Core_readouts":
                (
                    "PMAIP1/BBC3 perturbation versus ZBP1/MLKL "
                    "perturbation; mitochondrial markers; caspase-9; "
                    "pMLKL"
                ),

            "Decision":
                (
                    "Rescue patterns discriminate mitochondrial "
                    "apoptosis, regulated necrosis, or mixed execution."
                ),
        },

        {
            "Stage":
                "5 — Survival counter-program",

            "Question":
                (
                    "Are BDNF/KDR responses protective adaptations "
                    "rather than death drivers?"
                ),

            "Core_readouts":
                (
                    "BDNF/TrkB and KDR pathway perturbation/rescue "
                    "with viability and downstream signaling"
                ),

            "Decision":
                (
                    "Loss-of-function worsening plus ligand rescue "
                    "supports compensatory survival signaling."
                ),
        },
    ]
)


# ============================================================
# Write Excel
# ============================================================

with pd.ExcelWriter(
    OUTPUT_XLSX,
    engine="openpyxl",
) as writer:

    final.to_excel(
        writer,
        sheet_name="Candidate priorities",
        index=False,
    )

    program.to_excel(
        writer,
        sheet_name="Experimental program",
        index=False,
    )


# ============================================================
# Format workbook
# ============================================================

wb = load_workbook(
    OUTPUT_XLSX
)


header_fill = PatternFill(
    fill_type="solid",
    fgColor="D9E1F2",
)


for ws in wb.worksheets:

    ws.freeze_panes = "A2"

    if ws.max_row > 1:
        ws.auto_filter.ref = (
            ws.dimensions
        )

    for cell in ws[1]:

        cell.font = Font(
            bold=True
        )

        cell.fill = (
            header_fill
        )

        cell.alignment = Alignment(
            horizontal="center",
            vertical="center",
            wrap_text=True,
        )

    for row in ws.iter_rows():

        for cell in row:

            cell.alignment = Alignment(
                vertical="top",
                wrap_text=True,
            )

    for column_cells in ws.columns:

        letter = get_column_letter(
            column_cells[
                0
            ].column
        )

        max_length = 0

        for cell in column_cells[
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
    OUTPUT_XLSX
)


# ============================================================
# Console summary
# ============================================================

print()
print("=" * 82)
print("F7 — EXPERIMENTAL CANDIDATE PRIORITIZATION")
print("=" * 82)

print(
    f"\nCandidates selected: "
    f"{len(final)}"
)

print()

for group, subset in final.groupby(
    "Validation_group",
    sort=False,
):

    print(
        group
    )

    print(
        "    "
        +
        ", ".join(
            subset[
                "Gene"
            ]
        )
    )

print()
print(
    "Outputs:"
)

print(
    f"    "
    f"{OUTPUT_CSV.relative_to(PROJECT_ROOT)}"
)

print(
    f"    "
    f"{OUTPUT_XLSX.relative_to(PROJECT_ROOT)}"
)

print()
print(
    "[NEXT] Review exact F6 evidence for the 15 candidates "
    "and choose the first experimental perturbation set."
)