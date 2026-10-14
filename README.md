# SARS-CoV-2 Origins Critical Review & Phylogenetic Reproduction Suite

[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Bioconda](https://img.shields.io/badge/bioconda-mafft%20%7C%20iqtree-green.svg)](https://bioconda.github.io/)

This repository contains the complete computational reproduction pipeline, empirical audit, and visualization suite for our critical review of:

> **Boni, M. F., Lemey, P., Jiang, X., Lam, T. T.-Y., Perry, B. W., Castoe, T. A., Rambaut, A., & Robertson, D. L. (2020).**  
> *Evolutionary origins of the SARS-CoV-2 sarbecovirus lineage responsible for the COVID-19 pandemic.*  
> **Nature Microbiology**, 5(11), 1408–1417.  
> [DOI: 10.1038/s41564-020-0771-4](https://doi.org/10.1038/s41564-020-0771-4)  
> Original code repository: [`plemey/SARSCoV2origins`](https://github.com/plemey/SARSCoV2origins)

The reproduction suite evaluates the validity of **Non-Recombinant Region 1 (NRR1)**, **NRR2**, and **NRA3 / NRR3**, testing whether concatenated "clonal" partitions can legitimately support molecular clock dating in BEAST. It performs *de novo* multiple sequence alignments, maximum-likelihood phylogenetic inference, rigorous topological hypothesis testing (Approximately Unbiased AU, SH, KH tests), Robinson-Foulds distance computation, and Wuhan-Hu-1 sister-clade resolution across both the **original 68-taxon dataset** and the **updated 78-taxon dataset** (incorporating post-publication bat sarbecoviruses `RpYN06`, `RsYN04`, `RmYN02`, `RmYN05`, and `RmYN08`).

---

## Table of Contents

- [Overview & Scientific Findings](#overview--scientific-findings)
- [Repository Structure](#repository-structure)
- [Prerequisites & System Dependencies](#prerequisites--system-dependencies)
- [Installation & Environment Setup](#installation--environment-setup)
  - [Option A: Conda Environment Setup (Recommended)](#option-a-conda-environment-setup-recommended)
  - [Option B: Python Virtual Environment (`venv` + `pip`)](#option-b-python-virtual-environment-venv--pip)
  - [External Bioinformatics Binaries](#external-bioinformatics-binaries)
  - [Upstream Data & Workspace Setup](#upstream-data--workspace-setup)
- [How to Run the Analyses](#how-to-run-the-analyses)
  - [1. Full Master Reproduction Pipeline](#1-full-master-reproduction-pipeline)
  - [2. Step-by-Step Script Execution](#2-step-by-step-script-execution)
  - [3. Unified Breakpoint & Methodological Audit](#3-unified-breakpoint--methodological-audit)
  - [4. Fast Shell-Based Recombination Audit](#4-fast-shell-based-recombination-audit)
- [Detailed Script Documentation](#detailed-script-documentation)
- [Generated Outputs & Results](#generated-outputs--results)
- [Troubleshooting & FAQs](#troubleshooting--faqs)
- [Citation & Contact](#citation--contact)

---

## Overview & Scientific Findings

A central premise of Boni et al. (2020) is that continuous or composite non-recombinant regions can be concatenated to date ancestral sarbecovirus divergence events. Their primary molecular dating backbone—**NRR1 (13.0 kb)**—was constructed by concatenating **Region B** (ORF1a: nsp3–nsp4, 5.5 kb), **Region C** (ORF1a: nsp4–nsp6, 2.5 kb), and **Region A'** (three spliced fragments of ORF1b: nsp12–nsp15, 5.0 kb), while pruning five recombinant taxa (`DQ412042`, `CoVZXC21`, `CoVZC45`, `BM48-31`, and `SC2018`). From this partition, they estimated that SARS-CoV-2 diverged from its bat relative **RaTG13** in **1969** (95% HPD: 1948–1982).

This reproduction suite demonstrates through *de novo* computational re-analysis:

1. **Severe Internal Discordance within NRR1:**  
   Region B, Region C, and Region A' do not share a common evolutionary history. In the 68-taxon dataset, the normalized Robinson-Foulds distance between Region B and Region A' is **0.7077**—identical to the discordance between Region B and the hyper-recombinant **Spike gene** (**0.7077**). Even in the 63-taxon pruned matrix, normalized RF distances remain elevated (**0.5000 to 0.6333**).
2. **Definitive Rejection via AU Topology Tests ($p < 10^{-38}$):**  
   In Approximately Unbiased (AU) topology tests (10,000 multiscale bootstrap replicates in IQ-TREE), Region B statistically rejects Region C ($p = 1.65 \times 10^{-114}$), Region A' ($p = 8.64 \times 10^{-53}$), and NRR1 ($p = 2.04 \times 10^{-41}$). Reciprocally, Region C rejects Region B ($p = 7.17 \times 10^{-39}$) and NRR1 ($p = 2.13 \times 10^{-83}$), while Region A' rejects Region B ($p = 1.20 \times 10^{-90}$), Region C ($p = 3.56 \times 10^{-66}$), and NRR1 ($p = 2.44 \times 10^{-54}$). Concatenating these partitions violates fundamental phylogenetic assumptions.
3. **Displacement of the SARS-CoV-2 / RaTG13 Sister Relationship:**  
   In the 78-taxon dataset, the inclusion of **RpYN06** and **RmYN02** breaks the exclusive sister relationship between SARS-CoV-2 and RaTG13 across NRR1. In Region C, **RpYN06 is the direct sister of SARS-CoV-2 with 100% bootstrap support**, placing RaTG13 in an outer basal position. In Region B and Region A', SARS-CoV-2 clusters with `(RpYN06, RmYN02, PrC31)` with 97–98% support to the exclusion of RaTG13. Only in Spike does RaTG13 remain the immediate sister.
4. **NRA3 Sub-Block Heterogeneity:**  
   Dissecting the 21 kb RDP5 clonal region NRA3 into its four unmasked structural domains (Block 1 ORF1a, Block 2 ORF1b, Block 3 Spike, Block 4 3') reveals significant 4-gamete incompatibilities and reciprocal AU test rejections ($p < 0.05$), demonstrating that NRA3 is also a mosaic of conflicting phylogenetic signals.

---

## Repository Structure

```text
critical-review-boni-et-al-2020/
├── README.md                              # This documentation file
├── environment.yml                        # Conda environment definition (critrev)
├── requirements.txt                       # Python dependencies (Biopython, DendroPy, etc.)
├── run_reproduction_pipeline.py           # Master automated pipeline runner (Scripts 01-08)
├── run_all_analyses.py                    # Unified 7-procedure breakpoint & alignment audit
├── reproduce_breakpoint_analysis.sh       # Fast Bash verification of 3SEQ, GARD, and RAxML logs
├── data/                                  # Extracted unaligned sequences
│   ├── 68_taxa/                           # Original Boni et al. unaligned FASTA files
│   └── 78_taxa/                           # Updated unaligned FASTA files
├── alignments/                            # De novo MAFFT multiple sequence alignments
│   ├── 68_taxa/                           # Region B, C, A', NRR1, NRR2, NRA3, Spike
│   └── 78_taxa/                           # Full & core partition alignments
├── trees/                                 # IQ-TREE phylogenetic trees & AU test results
│   ├── 68_taxa/                           # ML trees (.treefile), UFBoot (.contree), AU tests
│   └── 78_taxa/                           # ML trees, logs, distances, and candidates
├── results/                               # Statistical summaries, CSV tables, markdown reports
│   ├── alignment_metrics_summary.csv      # Lengths, gap counts, alignment comparisons
│   ├── robinson_foulds_distances.csv      # Raw and normalized RF distances
│   ├── topology_tests_summary.csv         # AU, SH, KH, ELW test p-values
│   ├── wuhan_clade_resolution.csv         # Sister lineages, bootstrap support, MRCA stats
│   ├── pairwise_identities_wuhan.csv      # Nucleotide identities to Wuhan-Hu-1
│   ├── nra3_four_gamete_incompatibility.csv # 4-gamete homoplasy scan for NRA3
│   ├── nra3_subblock_topology_tests.csv   # Sub-block AU tests for NRA3
│   ├── nra3_recombination_audit_report.md # Granular NRA3 audit report
│   └── reproduction_summary_report.md     # Comprehensive reproduction report
├── figures/                               # Publication-ready figures (PDF and 300 DPI PNG)
│   ├── hashrf_tree_space_and_heatmaps.png/.pdf       # Figure 1: HashRF distance heatmaps
│   ├── 68_hashrf_comparative_phylogenies.png/.pdf   # Figure 2: 68-taxon comparative trees
│   ├── 78_hashrf_comparative_phylogenies.png/.pdf   # Figure 2: 78-taxon comparative trees
│   ├── pairwise_identities_heatmap.png/.pdf         # Pairwise identity dual-panel heatmap
│   ├── wuhan_clade_resolution_composite_heatmap.png/.pdf # Clade resolution composite
│   ├── wuhan_clade_sister_taxa_heatmap.png/.pdf     # Lineage-resolved sister matrix
│   └── nra3_recombination_audit.png/.pdf           # NRA3 sub-block audit figure
└── scripts/                               # Core pipeline execution scripts
    ├── 01_prepare_unaligned_data.py       # Slices and extracts unaligned sequences
    ├── 02_run_mafft_alignments.py         # De novo MAFFT alignments
    ├── 03_run_phylogenetic_inference.py   # IQ-TREE ML inference (GTR+G4, 1000 UFBoot)
    ├── 04_analyze_nrr1_congruence.py      # RF distances, AU tests, clade resolution
    ├── 05_visualize_hashrf_trees.py       # HashRF tree space & phylogenetic trees
    ├── 06_plot_pairwise_identities_heatmap.py # Pairwise identity heatmaps
    ├── 07_plot_wuhan_clade_support_heatmap.py # Clade support heatmaps
    └── 08_audit_nra3_recombination.py     # NRA3 dissection, AU tests, 4-gamete scan
```

---

## Prerequisites & System Dependencies

- **Operating System:** Linux (x86_64, tested on Ubuntu/Debian/Arch) or macOS.
- **Python:** Version 3.10, 3.11, or 3.12 (Python 3.12 recommended).
- **Package Manager:** `conda` (Miniconda / Anaconda) or `mamba`.
- **System Tools:** `git`, `bash`, `make`, `gcc` / `g++` (if building HashRF from source).

### External Bioinformatics Tools

The reproduction pipeline requires three external binaries:

| Tool | Version Tested | Purpose | Installation Source |
|------|---------------|---------|---------------------|
| **MAFFT** | `7.525` | *De novo* multiple sequence alignments | Bioconda / System package |
| **IQ-TREE 2** | `2.2.6` | Maximum-likelihood trees & AU tests | Bioconda / System package |
| **HashRF** | `6.0.1` | Fast Robinson-Foulds distance matrices | Compiled binary (`hashrf`) |

---

## Installation & Environment Setup

### Option A: Conda Environment Setup (Recommended)

1. **Create the `critrev` conda environment:**
   ```bash
   conda env create -f environment.yml
   ```

2. **Activate the environment:**
   ```bash
   conda activate critrev
   ```

3. **Verify Python dependencies:**
   If `pip` dependencies need manual synchronization from the local directory:
   ```bash
   pip install -r requirements.txt
   ```

### Option B: Python Virtual Environment (`venv` + `pip`)

If you prefer using standard Python without Conda:

```bash
# Create and activate virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install Python packages
pip install --upgrade pip setuptools wheel
pip install -r requirements.txt

# Ensure external binaries are installed via system package manager or direct download:
# Ubuntu/Debian:
# sudo apt-get update && sudo apt-get install -y mafft iqtree
```

### External Bioinformatics Binaries

The scripts automatically detect binaries in your active `PATH`, the active Conda environment (`$CONDA_PREFIX/bin`), and common installation paths. You can also override binary locations using environment variables:

```bash
export MAFFT_BIN="/path/to/mafft"
export IQTREE_BIN="/path/to/iqtree2"
```

#### HashRF Setup (for Script 05)
Script 05 uses HashRF to compute exact, fast Robinson-Foulds distance matrices.
- The pipeline looks for HashRF at `../code-tools/hashrf` or in your system `PATH`.
- To compile HashRF from source if needed:
  ```bash
  git clone https://github.com/sul/hashrf.git /tmp/hashrf_src
  cd /tmp/hashrf_src && make
  # Place the compiled 'hashrf' binary in ../code-tools/ or in your PATH
  ```

### Upstream Data & Workspace Setup

The reproduction suite compares *de novo* findings with Boni et al.'s original data files. The scripts are designed to work in a standard workspace where the upstream repository and the 78-taxon alignment are present:

1. **Clone the original Boni et al. repository (`SARSCoV2origins`):**
   ```bash
   # Clone into the parent directory or create a symlink in the project root:
   git clone https://github.com/plemey/SARSCoV2origins.git ../SARSCoV2origins
   ln -s ../SARSCoV2origins ./SARSCoV2origins
   ```

2. **78-Taxon Whole-Genome Alignment (`mar21update_78vir_wgal_fn.fas`):**
   - Ensure `mar21update_78vir_wgal_fn.fas` is located in the parent directory (`../mar21update_78vir_wgal_fn.fas`) or symlinked.

3. **Workspace Symlink for Reproduction Path Consistency:**
   Scripts `04` and `08` reference the reproduction root path `reproduction/`. When working directly inside `critical-review-boni-et-al-2020`, create a symlink to map `reproduction` to the current folder:
   ```bash
   ln -s . ../reproduction 2>/dev/null || true
   ln -s . ./reproduction 2>/dev/null || true
   ```

---

## How to Run the Analyses

### 1. Full Master Reproduction Pipeline

To execute all 8 reproduction steps end-to-end:

```bash
# Ensure the conda environment is active
conda activate critrev

# Run the master pipeline
python3 run_reproduction_pipeline.py
```

This will run sequentially:
1. `01_prepare_unaligned_data.py` (Sequence slicing & extraction)
2. `02_run_mafft_alignments.py` (De novo MAFFT alignment)
3. `03_run_phylogenetic_inference.py` (IQ-TREE ML tree inference)
4. `04_analyze_nrr1_congruence.py` (RF distances, AU tests, clade resolution)
5. `05_visualize_hashrf_trees.py` (HashRF matrices & tree space plots)
6. `06_plot_pairwise_identities_heatmap.py` (Pairwise identity heatmaps)
7. `07_plot_wuhan_clade_support_heatmap.py` (Wuhan clade resolution heatmaps)
8. `08_audit_nra3_recombination.py` (NRA3 structural audit & 4-gamete test)

> **Note on Runtime:** IQ-TREE phylogenetic inference across 16 partitions with 1,000 bootstrap replicates and AU multiscale tests can take between 20 and 60 minutes depending on CPU core count. By default, scripts utilize 8 CPU threads. To skip recomputing already existing alignments or trees, scripts detect existing files. To force complete re-computation, pass `--redo`.

---

### 2. Step-by-Step Script Execution

Each script can also be run independently:

#### Step 1: Extract Unaligned Sequences
```bash
python3 scripts/01_prepare_unaligned_data.py
```
*Outputs:* Unaligned FASTA files in `data/68_taxa/` and `data/78_taxa/`.

#### Step 2: De Novo MAFFT Alignments
```bash
# Run de novo alignments (uses existing files if already present)
python3 scripts/02_run_mafft_alignments.py

# To force re-alignment from scratch:
python3 scripts/02_run_mafft_alignments.py --redo
```
*Outputs:* Multiple sequence alignments in `alignments/68_taxa/` and `alignments/78_taxa/`.

#### Step 3: Maximum-Likelihood Phylogenies (IQ-TREE)
```bash
# Infer ML trees under GTR+G4 with 1,000 UFBoot
python3 scripts/03_run_phylogenetic_inference.py

# To force re-inference:
python3 scripts/03_run_phylogenetic_inference.py --redo
```
*Outputs:* `.treefile`, `.contree`, `.iqtree`, and `.log` files in `trees/68_taxa/` and `trees/78_taxa/`.

#### Step 4: Topological Congruence, RF Distances & AU Tests
```bash
python3 scripts/04_analyze_nrr1_congruence.py
```
*Outputs:*
- `results/alignment_metrics_summary.csv`
- `results/robinson_foulds_distances.csv`
- `results/topology_tests_summary.csv`
- `results/wuhan_clade_resolution.csv`
- `results/pairwise_identities_wuhan.csv`
- `results/reproduction_summary_report.md`

#### Step 5: HashRF Tree Space & Comparative Topologies
```bash
python3 scripts/05_visualize_hashrf_trees.py
```
*Outputs:*
- `figures/hashrf_tree_space_and_heatmaps.png` and `.pdf` (Figure 1)
- `figures/68_hashrf_comparative_phylogenies.png` and `.pdf` (Figure 2, 68 taxa)
- `figures/78_hashrf_comparative_phylogenies.png` and `.pdf` (Figure 2, 78 taxa)
- `results/hashrf_*.txt` distance matrices

#### Step 6: Pairwise Identity Heatmaps
```bash
python3 scripts/06_plot_pairwise_identities_heatmap.py \
    --input results/pairwise_identities_wuhan.csv \
    --outdir figures
```
*Outputs:*
- `figures/pairwise_identities_heatmap.png` and `.pdf`
- `figures/pairwise_identities_68_taxa_heatmap.png` and `.pdf`
- `figures/pairwise_identities_78_taxa_heatmap.png` and `.pdf`

#### Step 7: Wuhan-Hu-1 Clade Resolution Heatmaps
```bash
python3 scripts/07_plot_wuhan_clade_support_heatmap.py \
    --input results/wuhan_clade_resolution.csv \
    --outdir figures
```
*Outputs:*
- `figures/wuhan_clade_support_heatmap.png` and `.pdf`
- `figures/wuhan_clade_sister_taxa_heatmap.png` and `.pdf`
- `figures/wuhan_clade_resolution_composite_heatmap.png` and `.pdf`

#### Step 8: NRA3 / NRR3 Recombination Audit & 4-Gamete Testing
```bash
python3 scripts/08_audit_nra3_recombination.py
```
*Outputs:*
- `results/nra3_four_gamete_incompatibility.csv`
- `results/nra3_subblock_topology_tests.csv`
- `results/nra3_clade_resolution_summary.csv`
- `results/nra3_recombination_audit_report.md`
- `figures/nra3_recombination_audit.png` and `.pdf`

---

### 3. Unified Breakpoint & Methodological Audit

To audit the raw 3SEQ CSV outputs, GARD breakpoint tables, and RAxML developer logs from the original `SARSCoV2origins` repository in a single execution:

```bash
python3 run_all_analyses.py
```

This runs 7 rigorous verification procedures:
- **Procedure 1:** First 3SEQ Pass — Full genome (30,927 nt) breakpoint extraction & Top 10 BFR identification.
- **Procedure 2:** Second 3SEQ Pass — Auditing Region A internal mosaicism (76 sub-BFRs).
- **Procedure 3:** GARD vs. 3SEQ Discordance Filtering — Quantifying discarded breakpoints within ORF1b.
- **Procedure 4:** Audit of Combined BFRs (NRR1 concatenation logic).
- **Procedure 5:** Alignment Audit — Gap patterns and taxon count verification.
- **Procedure 6:** Breakpoint Density Profile — Identifying genome-wide recombination hotspots.
- **Procedure 7:** RAxML Parameter & Heuristics Review — Evaluating asymmetric bootstrap thresholds and dropped taxa.

*Outputs:* Full console log exported to `analysis_results.md` and CSV profiles in `analysis_results/`.

---

### 4. Fast Shell-Based Recombination Audit

For an immediate command-line verification of the mathematical properties of Boni et al.'s 3SEQ outputs:

```bash
chmod +x reproduce_breakpoint_analysis.sh
./reproduce_breakpoint_analysis.sh
```

---

## Detailed Script Documentation

### `scripts/01_prepare_unaligned_data.py`
Extracts raw nucleotide sequences corresponding to biological partitions from the 68-genome alignment (`nCoVbat.pan6.fst`) and the 78-genome alignment (`mar21update_78vir_wgal_fn.fas`).
- Strips alignment gaps (`-`) and ambiguous pads (`?`).
- Slices Region B (nsp3–4), Region C (nsp4–6), Region A (nsp12–15), Region A' (three spliced fragments: nt 13291–14932, 15405–17162, 18009–19628), NRR2 (11886–21753), and Spike (21563–25384).
- Assembles NRR1 Full (A' + B + C across all taxa) and NRR1 Pruned (excluding 5 recombinant taxa: `DQ412042`, `CoVZXC21`, `CoVZC45`, `NC_014470`, `BM48-31`, `SC2018`).

### `scripts/02_run_mafft_alignments.py`
Performs *de novo* multiple sequence alignment using MAFFT (`--auto --thread 8`) on each extracted unaligned partition. Preserves true biologically informative insertions and deletions without propagating pre-aligned column artifacts.

### `scripts/03_run_phylogenetic_inference.py`
Conducts maximum-likelihood phylogenetic inference using IQ-TREE 2.
- **Substitution Model:** `GTR+G4` (General Time-Reversible model with 4 Gamma rate categories).
- **Branch Support:** 1,000 Ultrafast Bootstrap (UFBoot) replicates (`-B 1000`).
- **Parallelization:** Multi-threaded execution (`-T 8`).

### `scripts/04_analyze_nrr1_congruence.py`
Computes key quantitative metrics of congruence and phylogenetic resolution:
1. **Robinson-Foulds Distances:** Computes raw and normalized RF metrics ($RF_{norm} = RF / (2N - 6)$) across all pairwise partition comparisons for both 68- and 78-taxon datasets.
2. **AU / SH / KH Topology Tests:** Evaluates candidate ML trees against each alignment using IQ-TREE's multiscale bootstrap hypothesis testing (`-zb 10000 -au`).
3. **Clade Resolution:** Identifies immediate sister taxa and Most Recent Common Ancestor (MRCA) bootstrap support for SARS-CoV-2 (Wuhan-Hu-1) across all partitions.

### `scripts/05_visualize_hashrf_trees.py`
Constructs multi-tree files, executes HashRF to generate unweighted and branch-length weighted RF distance matrices, projects tree topologies into 2D space using Classical Multidimensional Scaling (PCoA / MDS), and produces publication-quality comparative phylogenetic trees (Figures 1 and 2).

### `scripts/06_plot_pairwise_identities_heatmap.py`
Generates high-contrast heatmaps showing nucleotide percent identity between Wuhan-Hu-1 and key sarbecovirus reference genomes across partitions.

### `scripts/07_plot_wuhan_clade_support_heatmap.py`
Plots UFBoot support and sister-taxa assignments for Wuhan-Hu-1, visually highlighting the dramatic shift from RaTG13 to RpYN06 / RmYN02 in the updated 78-taxon dataset.

### `scripts/08_audit_nra3_recombination.py`
Dissects the 21 kb RDP5 clonal region NRA3 into four unmasked structural blocks (Block 1 ORF1a, Block 2 ORF1b, Block 3 Spike, Block 4 3'), builds pruned alignments across common core taxa, runs AU topology tests between sub-blocks, and conducts 4-gamete incompatibility scans across NRA3.

---

## Generated Outputs & Results

### Statistical Summary Tables (`results/`)

| File | Description | Key Metric / Result |
|------|-------------|---------------------|
| `robinson_foulds_distances.csv` | Pairwise RF distances between partitions | Region B vs Region A' $RF_{norm} = 0.7077$ (68 taxa) |
| `topology_tests_summary.csv` | AU, SH, KH, and ELW hypothesis test p-values | Region B rejects Region A' tree ($p_{AU} = 8.64 \times 10^{-53}$) |
| `wuhan_clade_resolution.csv` | Sister taxon assignments & UFBoot support | RpYN06 sister to Wuhan-Hu-1 in Region C (100% support) |
| `pairwise_identities_wuhan.csv` | Nucleotide identities to Wuhan-Hu-1 | RpYN06 highest identity in Region C (97.4%) |
| `nra3_four_gamete_incompatibility.csv` | Pairwise homoplasy scan for NRA3 | Widespread 4-gamete violations across structural blocks |
| `nra3_subblock_topology_tests.csv` | AU tests between NRA3 structural blocks | Mutual topology rejection between Blocks 1, 2, 3, and 4 |

### Summary Reports (`results/` & Root)

- `results/reproduction_summary_report.md`: Markdown report detailing alignment metrics, RF discordance tables, AU test p-values, and Wuhan clade support.
- `results/nra3_recombination_audit_report.md`: Audit of the NRA3 RDP5 dataset, documenting sub-block incongruence and 4-gamete tests.
- `analysis_results.md`: Complete console log of `run_all_analyses.py` covering 3SEQ and GARD breakpoint audits.

### Visualizations (`figures/`)

- `hashrf_tree_space_and_heatmaps.{png,pdf}`: Heatmaps of unweighted and branch-length weighted Robinson-Foulds distances.
- `68_hashrf_comparative_phylogenies.{png,pdf}`: Comparative ML trees across Region B, Region C, Region A', NRR1, NRR2, NRA3, and Spike for 68 taxa.
- `78_hashrf_comparative_phylogenies.{png,pdf}`: Comparative ML trees across partitions for 78 taxa with color-coded lineages.
- `pairwise_identities_heatmap.{png,pdf}`: Comparative heatmap of pairwise nucleotide identities to Wuhan-Hu-1.
- `wuhan_clade_resolution_composite_heatmap.{png,pdf}`: Composite heatmap showing clade bootstrap support and sister taxon displacement.
- `nra3_recombination_audit.{png,pdf}`: Multi-panel figure showing NRA3 sub-block topologies, AU test heatmaps, and 4-gamete density.

---

## Troubleshooting & FAQs

### 1. `MAFFT binary not found` or `IQ-TREE binary not found`
If the scripts cannot locate MAFFT or IQ-TREE:
- Verify that your conda environment is activated: `conda activate critrev`
- Ensure the binaries are executable: `which mafft` and `which iqtree2`
- You can explicitly provide their paths via environment variables:
  ```bash
  export MAFFT_BIN=$(which mafft)
  export IQTREE_BIN=$(which iqtree2)
  ```

### 2. `HashRF binary not found`
If `scripts/05_visualize_hashrf_trees.py` reports that HashRF is missing:
- Ensure the compiled binary is placed at `../code-tools/hashrf` or in your system `PATH`.
- Verify execute permissions: `chmod +x ../code-tools/hashrf`

### 3. FileNotFoundError: `../SARSCoV2origins` or `mar21update_78vir_wgal_fn.fas`
If running outside the structured assignment workspace:
- Clone the upstream repository: `git clone https://github.com/plemey/SARSCoV2origins.git ../SARSCoV2origins`
- Symlink the repository to the project root: `ln -s ../SARSCoV2origins ./SARSCoV2origins`
- If you have already executed the pipeline once, all extracted data, alignments, and trees are already stored in `data/`, `alignments/`, and `trees/`.

### 4. Adjusting CPU Core Allocation
By default, MAFFT and IQ-TREE use 8 threads (`--thread 8`, `-T 8`). If running on a system with fewer or more CPU cores:
- Edit the thread argument in `scripts/02_run_mafft_alignments.py` (`--thread <N>`) and `scripts/03_run_phylogenetic_inference.py` (`-T <N>`).

---

## Citation & Contact

If you use this codebase, reproduction pipeline, or audit data in your work, please cite:

```bibtex
@article{Boni2020CriticalReview,
  title={SARS-CoV-2 Origins: Comprehensive Breakpoint & Phylogenetic Incongruence Reproduction},
  author={Critical Review Reproduction Suite},
  year={2026},
  note={Evaluating Boni et al. (2020) Nature Microbiology 5(11): 1408–1417}
}
```
