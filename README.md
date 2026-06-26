# NG_Analysis

Computational genomics analysis of **Rb-family / E2F-dependent transcriptional regulation** in neural stem and progenitor lineages, with a focus on how promoter binding and expression changes relate to neuronal survival, apoptosis, autophagy, necrosis, neurogenesis, and cell-cycle control.

This repository was developed for the R-Noël project around the question:

> How does loss of Rb-family pocket proteins reshape E2F-regulated transcriptional programs in neural lineages, and which downstream genes or pathways may explain the survival-versus-death phenotype?

---

## Project overview

The project integrates two main types of data:

1. **E2F promoter binding data**
   - ChIP-chip promoter occupancy data from Julian et al. 2016.
   - Main comparison: **E2F4 vs E2F3** promoter binding.
   - Goal: identify genes whose promoter regions show stronger E2F4-associated binding.

2. **Rb/E2F-dependent expression data**
   - Expression datasets from Fong et al. 2022 and related Rb-family / E2F perturbation studies.
   - Main use: compare promoter binding with expression changes after loss of Rb-family regulation.
   - Typical expression fields: gene symbol or Ensembl ID, log2 fold change, adjusted p value.

The analysis asks whether E2F-bound promoters are enriched for biologically meaningful categories and whether those bound genes are transcriptionally altered in Rb-family knockout conditions.

---

## Biological motivation

Rb-family pocket proteins — **Rb/pRb, p107, and p130** — regulate cell-cycle progression largely through E2F transcription factors. In neural stem/progenitor cells and post-mitotic neurons, this pathway is not only a cell-cycle gatekeeper; it is also linked to:

- maintenance of neural stem-cell quiescence,
- activation of adult neural stem cells,
- differentiation and neurogenesis,
- survival of newborn and mature neurons,
- DNA damage responses,
- apoptosis and other cell-death mechanisms.

The working biological idea is that loss of Rb-family repression may release or distort E2F-dependent transcriptional programs. This can produce a mixed phenotype: cell-cycle re-entry, DNA repair activation, altered neurogenic state, and possible engagement of apoptosis, autophagy, or necrosis-related pathways.

---

## Main analysis questions

This repository currently addresses the following questions:

1. **Which promoters are preferentially bound by E2F4 relative to E2F3?**

2. **Which biological processes are enriched among E2F4-biased promoter targets?**

3. **Do E2F4-bound genes overlap with genes differentially expressed after Rb-family loss?**

4. **Among the overlapping genes, which candidates are most relevant to:**
   - apoptosis,
   - autophagy,
   - necrosis,
   - cell cycle,
   - DNA repair,
   - neurogenesis / neural precursor proliferation?

5. **Can promoter binding plus expression direction help prioritize genes for experimental validation?**

---

## Repository structure

A suggested structure is shown below. Folder names can be adapted to the local repository, but large raw data folders should remain ignored by Git.

```text
NG_Analysis/
│
├── README.md
├── .gitignore
│
├── scripts/
│   ├── parse_chipchip.py
│   ├── differential_binding_limma.R
│   ├── go_enrichment.R
│   ├── export_go_results.R
│   ├── integrate_expression.py
│   └── plot_go_categories.py
│
├── data/
│   ├── raw/                 # large raw files; do not commit
│   ├── metadata/            # sample sheets, platform annotations
│   ├── processed/           # cleaned/intermediate tables
│   └── external/            # downloaded reference files, if needed
│
├── results/
│   ├── differential_binding/
│   ├── go_enrichment/
│   ├── expression_integration/
│   └── candidate_tables/
│
├── figures/
│   ├── qc/
│   ├── go_barplots/
│   └── summary_figures/
│
└── docs/
    └── notes/
```

### Important Git note

The raw data can be very large and should **not** be committed. Add these to `.gitignore`:

```gitignore
data/raw/
Data/
*.CEL
*.fastq
*.fq
*.bam
*.sam
*.bedGraph
*.bigWig
*.zip
*.tar.gz
```

If a large data folder was accidentally committed before the first push, undo the commit or remove it from Git history before publishing the branch.

---

## Data inputs

The exact filenames may differ locally, but the pipeline expects the following kinds of inputs.

### 1. ChIP-chip / promoter binding files

Expected information:

- probe ID,
- sample identity,
- IP intensity,
- input/control intensity,
- promoter annotation,
- gene symbol or gene identifier.

The pipeline computes an M-value:

```text
M = log2(IP / Input)
```

This is used as the promoter-binding signal.

### 2. Platform annotation

Used to map probes to genes and to restrict the analysis to bona fide promoter probes.

Expected fields may include:

- probe ID,
- genomic annotation,
- promoter flag or description,
- gene symbol,
- Ensembl or Entrez ID.

### 3. Expression table

Expected fields:

```text
gene_id / symbol
log2FoldChange
padj
```

The expression values are interpreted as:

```text
positive log2FC  = upregulated in the knockout/comparison condition
negative log2FC  = downregulated in the knockout/comparison condition
```

For example:

```text
log2FC = 1      means about 2-fold upregulation
log2FC = -1     means about 2-fold downregulation
log2FC = 0.6    means about 2^0.6 ≈ 1.52-fold upregulation
```

---

## Software requirements

### R packages

Install R and the following packages:

```r
install.packages(c(
  "tidyverse",
  "readr",
  "dplyr",
  "openxlsx",
  "BiocManager"
))

BiocManager::install(c(
  "limma",
  "clusterProfiler",
  "org.Mm.eg.db",
  "GO.db",
  "AnnotationDbi"
))
```

Depending on the script version, additional packages may be needed:

```r
BiocManager::install(c(
  "rtracklayer",
  "GenomicFeatures",
  "TxDb.Mmusculus.UCSC.mm10.knownGene"
))
```

### Python packages

```bash
pip install pandas numpy scipy matplotlib seaborn openpyxl scikit-learn
```

---

## Pipeline summary

### Step 1 — Parse and clean ChIP-chip data

The raw ChIP-chip files are parsed into clean probe-level tables.

Typical operations:

- detect the feature table,
- keep real oligo probes,
- remove control probes,
- remove non-positive IP/input intensities,
- compute M-values,
- save one cleaned table per sample.

Output examples:

```text
processed/chipchip_sample_Mvalues.csv
processed/Mvalue_matrix.csv
```

---

### Step 2 — Restrict to promoter probes

The probe annotation file is used to keep promoter-associated probes only. This avoids mixing promoter binding with intergenic or gene-body signal.

Output examples:

```text
processed/GPL_promoter_probe2gene.csv
processed/promoter_Mvalue_matrix.csv
```

---

### Step 3 — Differential promoter binding

Differential binding is performed with `limma`, usually comparing:

```text
E2F4 vs E2F3
```

The preferred statistical approach is `limma::treat`, which tests for a minimum biologically meaningful fold-change rather than only testing whether the fold-change differs from zero.

Typical threshold:

```text
raw p < 0.05
|log2FC| >= 0.5
```

Output examples:

```text
results/differential_binding/E2F4_vs_E2F3_treat_full.csv
results/differential_binding/E2F4_vs_E2F3_significant_probes.csv
results/differential_binding/E2F4_vs_E2F3_significant_genes.csv
```

---

### Step 4 — GO Biological Process enrichment

Significant promoter-bound genes are tested for Gene Ontology Biological Process enrichment using `clusterProfiler`.

Rather than relying only on keyword search, the preferred approach is to use GO identifiers and/or descendants of chosen GO root terms.

Current major biological categories include:

| Category | GO root term |
|---|---|
| Apoptosis | GO:0006915 |
| Autophagy | GO:0006914 |
| Necrosis | GO:0070265 |
| Cell cycle | GO:0007049 |
| DNA repair | GO:0006281 |
| Neurogenesis | GO:0022008 |
| Neural precursor cell proliferation | GO:0061351 |

Output examples:

```text
results/go_enrichment/GO_results.xlsx
results/go_enrichment/E2F4_cell_cycle_GO.csv
results/go_enrichment/E2F4_DNA_repair_GO.csv
results/go_enrichment/E2F4_apoptosis_GO.csv
results/go_enrichment/E2F4_autophagy_GO.csv
results/go_enrichment/E2F4_necrosis_GO.csv
```

---

### Step 5 — Plot enriched categories

GO terms are visualized as horizontal bar plots, usually ranked by gene count, adjusted p value, or enrichment score.

Output examples:

```text
figures/go_barplots/cell_cycle_GO_barplot.png
figures/go_barplots/DNA_repair_GO_barplot.png
figures/go_barplots/apoptosis_GO_barplot.png
figures/go_barplots/autophagy_GO_barplot.png
figures/go_barplots/neurogenesis_GO_barplot.png
```

---

### Step 6 — Integrate promoter binding with expression

The expression table is collapsed to gene level and merged with the binding table.

The integrated table usually contains:

```text
symbol
binding_logFC
binding_pvalue
expression_log2FC
expression_adjP
GO_category_flags
```

This allows downstream filtering for genes that are both:

1. promoter-bound or differentially bound, and
2. transcriptionally changed in the Rb/E2F perturbation condition.

Output examples:

```text
results/expression_integration/binding_vs_expression.xlsx
results/expression_integration/Significant_both.csv
```

---

### Step 7 — Rank candidate genes

Candidate genes can be ranked using a combined logic such as:

```text
high absolute expression change
+ significant adjusted p value
+ promoter binding evidence
+ membership in relevant GO category
```

For apoptosis/autophagy/necrosis-focused analysis, the ranked table helps select genes for biological interpretation and possible validation.

Output examples:

```text
results/candidate_tables/apoptosis_candidates_ranked.csv
results/candidate_tables/autophagy_candidates_ranked.csv
results/candidate_tables/neurogenesis_candidates_ranked.csv
```

Candidate genes discussed during the project include examples such as:

```text
BDNF
IL6
ZBP1
PMAIP1
XKR5
```

These should be interpreted in the context of the latest pipeline output, not as final validated hits.

---

## Interpreting expression values

Expression values are usually log2 fold changes.

```text
log2FC > 0   upregulated
log2FC < 0   downregulated
log2FC = 0   no change
```

To convert log2FC to fold change:

```text
fold change = 2^(log2FC)
```

Examples:

| log2FC | Fold change | Interpretation |
|---:|---:|---|
| 0.5 | 1.41x | moderately upregulated |
| 0.6 | 1.52x | moderately upregulated |
| 1.0 | 2.00x | two-fold upregulated |
| -0.5 | 0.71x | moderately downregulated |
| -1.0 | 0.50x | two-fold downregulated |

The sign is therefore biologically important. The analysis should not use absolute values alone unless the goal is only to rank by magnitude.

---

## Current biological interpretation

The current analysis supports a broad view of Rb/E2F function in neural lineages:

- Rb-family loss strongly affects genes related to cell-cycle control and DNA repair.
- E2F promoter binding overlaps with genes involved in survival, death, and neural fate programs.
- Apoptosis remains biologically relevant, but it may not be the only or cleanest explanatory category.
- Autophagy-related enrichment appears plausible and should be interpreted alongside stress, survival, and metabolic-state changes.
- Necrosis has not consistently appeared as a strong enriched category in the current pipeline.
- Neurogenesis and neural precursor proliferation terms are prominent and may represent a stronger organizing direction for the project than a purely apoptosis-centered hypothesis.

A careful interpretation is that Rb/E2F dysregulation may not simply “turn on apoptosis.” Instead, it may destabilize neural cell identity and quiescence/activation states, producing cell-cycle re-entry, DNA damage, altered differentiation, and downstream survival defects.

---

## Recommended validation directions

Candidate genes from the computational pipeline can be prioritized for experimental follow-up using:

- qRT-PCR,
- Western blot,
- immunostaining,
- active caspase-3 staining,
- TUNEL,
- cell-cycle markers such as Ki67,
- DNA damage markers such as γH2AX,
- apoptosis/autophagy markers depending on the candidate pathway.

The strongest candidates should ideally satisfy several criteria:

1. significant expression change,
2. promoter-binding evidence,
3. relevant GO/pathway membership,
4. known or plausible neural function,
5. consistency with the Rb/E2F biological model,
6. feasibility of wet-lab validation.

---

## Reproducibility notes

To make the project reproducible:

- keep raw data out of Git,
- keep metadata and sample sheets versioned,
- save every intermediate table,
- document thresholds directly inside scripts,
- avoid manual Excel-only edits for final results,
- use gene symbols consistently,
- preserve Ensembl/Entrez IDs when possible,
- record whether expression comparisons are TKO vs THC, TKO vs DKO, or another contrast.

---

## Suggested citation notes

The project uses or builds on data and biological context from:

- Julian et al. 2016 — E2F promoter occupancy in neural precursor cells.
- Fong et al. 2022 — Rb/E2F axis regulation of adult neural stem-cell quiescence and activation.
- Related Rb-family / pocket-protein studies on neuronal survival, apoptosis, and neurogenesis.

Full bibliographic references should be added to the final report or manuscript.

---

## Project status

This repository contains an active analysis pipeline. The current direction is moving from a narrow apoptosis-only framework toward a broader model of:

```text
Rb-family loss
→ E2F-dependent transcriptional disruption
→ cell-cycle / DNA-repair activation
→ altered neural precursor state and neurogenesis programs
→ survival defects and context-dependent cell-death mechanisms
```

The README should be updated whenever the pipeline thresholds, input datasets, or final candidate-ranking logic change.
