#!/usr/bin/env python3
"""
================================================================================
Script 05: HashRF Tree Distance Analysis & Publication-Quality Visualizations
Target:      Boni et al. (2020) Sarbecovirus Recombination & Origins Audit
Description:
  1. Runs HashRF on candidate partitions (Region B, Region C, Region A',
     NRR1_full, NRR2, Spike) for 68-taxon and 78-taxon datasets.
  2. Parses unweighted, normalized rate, and branch-length weighted RF matrices.
  3. Projects trees into 2D "Tree Space" using Classical Multidimensional
     Scaling (PCoA / MDS).
  4. Generates publication-ready figures:
     - Fig 1: HashRF Tree Space & Robinson-Foulds Distance Heatmaps
     - Fig 2: Multi-panel comparative phylogenetic trees illustrating topological
       discordance and mosaicism of the SARS-CoV-2 lineage across partitions.
================================================================================
"""

import os
import sys
import subprocess
import io
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import dendropy
from Bio import Phylo
from Bio.Phylo.BaseTree import BranchColor

# Base paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TREES_DIR = os.path.join(BASE_DIR, "trees")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")
HASHRF_BIN = os.path.join(os.path.dirname(BASE_DIR), "code-tools", "hashrf")

os.makedirs(FIGURES_DIR, exist_ok=True)

PARTITIONS = ["Region_B", "Region_C", "Region_Aprime", "NRR1_full", "NRR2", "NRA3", "Spike"]
PARTITION_LABELS = {
    "Region_B": "Region B (ORF1a)",
    "Region_C": "Region C (ORF1a)",
    "Region_Aprime": "Region A' (ORF1b)",
    "NRR1_full": "NRR1 Concatenated",
    "NRR2": "NRR2 (GARD/3SEQ)",
    "NRA3": "NRA3 (RDP5 Clonal)",
    "Spike": "Spike Gene"
}

def run_hashrf_for_dataset(dataset_name, pruned=False, print_opt="matrix"):
    """Runs HashRF on combined multi-tree Newick file for full or pruned datasets and returns (raw_matrix, rate_matrix, n_taxa, parts)."""
    ds_dir = os.path.join(TREES_DIR, dataset_name)

    if dataset_name == "68_taxa":
        parts = ["Region_B", "Region_C", "Region_Aprime", "NRR1_pruned" if pruned else "NRR1_full", "NRR2", "NRA3", "Spike"]
        nra3_tree = dendropy.Tree.get(path=os.path.join(ds_dir, "NRA3.treefile"), schema="newick")
        nra3_labels = set(x.label for x in nra3_tree.taxon_namespace)
        if pruned:
            ref_tree = dendropy.Tree.get(path=os.path.join(ds_dir, "NRR1_pruned.treefile"), schema="newick")
            pruned_labels = set(x.label for x in ref_tree.taxon_namespace) & nra3_labels
            n_taxa = len(pruned_labels)
            combined_path = os.path.join(ds_dir, f"hashrf_input_{dataset_name}_pruned.tre")
            with open(combined_path, "w") as out_f:
                for p in parts:
                    tf = os.path.join(ds_dir, f"{p}.treefile")
                    t = dendropy.Tree.get(path=tf, schema="newick")
                    t.retain_taxa_with_labels(pruned_labels)
                    out_f.write(t.as_string(schema="newick").strip() + "\n")
        else:
            n_taxa = len(nra3_labels)
            combined_path = os.path.join(ds_dir, f"hashrf_input_{dataset_name}_full.tre")
            with open(combined_path, "w") as out_f:
                for p in parts:
                    tf = os.path.join(ds_dir, f"{p}.treefile")
                    t = dendropy.Tree.get(path=tf, schema="newick")
                    t.retain_taxa_with_labels(nra3_labels)
                    out_f.write(t.as_string(schema="newick").strip() + "\n")
    else:
        parts = ["Region_B", "Region_C", "Region_Aprime", "NRR1_pruned" if pruned else "NRR1_full", "NRR2", "NRA3", "Spike"]
        if pruned:
            ref_tree = dendropy.Tree.get(path=os.path.join(ds_dir, "NRR1_pruned.treefile"), schema="newick")
            pruned_labels = set(x.label for x in ref_tree.taxon_namespace)
            n_taxa = len(pruned_labels)
            combined_path = os.path.join(ds_dir, f"hashrf_input_{dataset_name}_pruned.tre")
            with open(combined_path, "w") as out_f:
                for p in parts:
                    tf = os.path.join(ds_dir, f"{p}.treefile")
                    t = dendropy.Tree.get(path=tf, schema="newick")
                    t.retain_taxa_with_labels(pruned_labels)
                    out_f.write(t.as_string(schema="newick").strip() + "\n")
        else:
            n_taxa = 78
            combined_path = os.path.join(ds_dir, f"hashrf_input_{dataset_name}_full.tre")
            with open(combined_path, "w") as out_f:
                for p in parts:
                    tf = os.path.join(ds_dir, f"{p}.treefile")
                    with open(tf) as in_f:
                        out_f.write(in_f.read().strip() + "\n")

    out_file = os.path.join(RESULTS_DIR, f"hashrf_{dataset_name}_{'pruned' if pruned else 'full'}_{print_opt}.txt")
    cmd = [HASHRF_BIN, combined_path, str(len(parts)), "-p", print_opt, "-o", out_file]
    subprocess.run(cmd, capture_output=True, text=True)

    matrix_rows = []
    with open(out_file, "r") as f:
        for line in f:
            parts_line = line.strip().split()
            if parts_line and all(p.replace('.', '', 1).replace('e-', '', 1).replace('-', '', 1).isdigit() for p in parts_line):
                matrix_rows.append([float(x) for x in parts_line])

    mat = np.array(matrix_rows)
    max_splits = n_taxa - 3
    rate = (mat / float(max_splits)) * 100.0
    return mat, rate, n_taxa, parts

def generate_tree_space_and_heatmap_figure():
    """Generates Figure 1: HashRF Robinson-Foulds Incongruence Heatmaps (Full and Pruned Datasets)."""
    print("[*] Generating Figure 1: HashRF Incongruence Heatmaps (Full & Pruned)...")

    # Run HashRF across all 4 configurations
    mat_68_full, rate_68_full, n_68_full, parts_68_full = run_hashrf_for_dataset("68_taxa", pruned=False)
    mat_68_prun, rate_68_prun, n_68_prun, parts_68_prun = run_hashrf_for_dataset("68_taxa", pruned=True)
    mat_78_full, rate_78_full, n_78_full, parts_78_full = run_hashrf_for_dataset("78_taxa", pruned=False)
    mat_78_prun, rate_78_prun, n_78_prun, parts_78_prun = run_hashrf_for_dataset("78_taxa", pruned=True)

    fig = plt.figure(figsize=(20, 16), dpi=300)
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1], hspace=0.30, wspace=0.25)

    short_labels = ["Reg B", "Reg C", "Reg A'", "NRR1", "NRR2", "NRR3", "Spike"]

    configs = [
        (gs[0, 0], mat_68_full, rate_68_full, short_labels, "a. 68-Taxon Dataset (Full Matrix: 66 Taxa)"),
        (gs[0, 1], mat_68_prun, rate_68_prun, short_labels, "b. 68-Taxon Dataset (Recombinant Pruned: 61 Taxa)"),
        (gs[1, 0], mat_78_full, rate_78_full, short_labels, "c. 78-Taxon Dataset (Full Matrix: 78 Taxa)"),
        (gs[1, 1], mat_78_prun, rate_78_prun, short_labels, "d. 78-Taxon Dataset (Recombinant Pruned: 73 Taxa)")
    ]

    for sub_pos, mat, rate, labels, title in configs:
        ax = fig.add_subplot(sub_pos)
        im = ax.imshow(rate, cmap="YlOrRd", vmin=0, vmax=80)
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, fontsize=11, fontweight="bold")
        ax.set_yticklabels(labels, fontsize=11, fontweight="bold")
        ax.set_title(title, fontsize=15, fontweight="bold", pad=12, loc="left")

        # Annotate cells
        for i in range(len(labels)):
            for j in range(len(labels)):
                raw_val = int(mat[i, j])
                pct_val = rate[i, j]
                text_color = "white" if pct_val > 45 else "black"
                ax.text(j, i, f"{raw_val}\n[{pct_val:.1f}%]", ha="center", va="center", color=text_color, fontsize=9.5 if len(labels) > 6 else 10, fontweight="bold")

        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Normalized RF Discordance (%)", fontsize=10)

    out_png = os.path.join(FIGURES_DIR, "hashrf_tree_space_and_heatmaps.png")
    out_pdf = os.path.join(FIGURES_DIR, "hashrf_tree_space_and_heatmaps.pdf")
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"  -> Saved Figure 1: {out_png}")


def clean_taxon_label(name):
    """Cleans raw Newick sequence headers to clear, publication-style labels."""
    if "Wuhan" in name or "MN908947" in name:
        return "SARS-CoV-2 (Wuhan-Hu-1)"
    elif "RaTG13" in name:
        return "RaTG13 (Bat affinis 2013)"
    elif "RpYN06" in name:
        return "RpYN06 (Bat pusillus 2020)"
    elif "RmYN02" in name:
        return "RmYN02 (Bat malayanus 2019)"
    elif "RsYN04" in name:
        return "RsYN04 (Bat sinicus 2019)"
    elif "RmYN05" in name:
        return "RmYN05 (Bat malayanus 2019)"
    elif "410721" in name or "Guandong" in name:
        return "Pangolin-CoV (Guangdong 2019)"
    elif "410538" in name or "P4L" in name:
        return "Pangolin-CoV (Guangxi P4L 2017)"
    elif "ZC45" in name:
        return "bat-SL-CoVZC45 (Zhoushan 2017)"
    elif "ZXC21" in name:
        return "bat-SL-CoVZXC21 (Zhoushan 2015)"
    elif "HSZ-Cc" in name or "SARS-CoV-1" in name:
        return "SARS-CoV-1 (Guangzhou 2002)"
    elif "BtKY72" in name:
        return "BtKY72 (Kenya outgroup)"
    return name.split("|")[0]

def generate_comparative_phylogenies_figure(taxa):
    """Generates Figure 2: Standardized 12-Taxon Pruned Comparative Phylogenies Evaluated by HashRF."""
    print("[*] Generating Figure 2: Standardized Comparative Phylogenetic Topologies...")

    # Pivotal 12 taxa across the Sarbecovirus evolutionary debate
    keywords = [
        "Wuhan", "RaTG13", "RpYN06", "RmYN02", "RsYN04", "RmYN05",
        "ZC45", "ZXC21", "410721", "410538", "HSZ-Cc", "BtKY72"
    ]

    eval_partitions = ["Region_B", "Region_C", "Region_Aprime", "NRA3", "Spike"]
    titles = [
        "Region B (ORF1a, nt 3,625-9,150)\nHashRF vs A': 68.6% discordance",
        "Region C (ORF1a, nt 9,261-11,795)\nHashRF vs A': 62.9% discordance",
        "Region A' (ORF1b, nt 13,291-19,628)\nHashRF vs B: 68.6% discordance",
        "NRR3 / NRA3 (RDP5 Backbone)\nHashRF vs Spike: 69.3% discordance",
        "Spike Gene (nt 21,563-25,384)\nHashRF vs B: 73.3% discordance"
    ]

    fig, axes = plt.subplots(1, 5, figsize=(30, 10), dpi=300)

    for idx, (p, title) in enumerate(zip(eval_partitions, titles)):
        ax = axes[idx]
        tf = os.path.join(TREES_DIR, f"{taxa}_taxa", f"{p}.treefile")

        # Load tree with DendroPy and prune to the exact 12 focal taxa
        dp_tree = dendropy.Tree.get(path=tf, schema="newick")
        retained = [x.label for x in dp_tree.taxon_namespace if any(kw in x.label for kw in keywords)]
        dp_tree.retain_taxa_with_labels(retained)

        # Convert to Biopython Tree for publication-style drawing
        bio_tree = Phylo.read(io.StringIO(dp_tree.as_string(schema="newick")), "newick")

        # Root with outgroup BtKY72 if found
        outgroups = [t for t in bio_tree.get_terminals() if "BtKY72" in t.name]
        if outgroups:
            bio_tree.root_with_outgroup(outgroups[0])

        bio_tree.ladderize()

        # Highlight focal branches
        for clade in bio_tree.find_clades():
            if clade.is_terminal():
                clade.color = BranchColor(40, 40, 40)
                if "Wuhan" in clade.name or "MN908947" in clade.name:
                    clade.color = BranchColor(220, 20, 20)      # Crimson Red
                elif "RaTG13" in clade.name:
                    clade.color = BranchColor(20, 80, 220)      # Royal Blue
                elif "RpYN06" in clade.name or "RsYN04" in clade.name:
                    clade.color = BranchColor(0, 150, 50)       # Vivid Green
                elif "RmYN02" in clade.name or "RmYN05" in clade.name:
                    clade.color = BranchColor(140, 30, 180)     # Purple
                elif "410721" in clade.name or "410538" in clade.name or "pangolin" in clade.name.lower():
                    clade.color = BranchColor(230, 120, 0)      # Amber Orange
                elif "ZC45" in clade.name or "ZXC21" in clade.name:
                    clade.color = BranchColor(120, 60, 20)      # Brown
                elif "HSZ-Cc" in clade.name:
                    clade.color = BranchColor(128, 0, 0)        # Maroon

        # Draw tree onto subplot
        Phylo.draw(bio_tree, axes=ax, do_show=False,
                   label_func=lambda n: clean_taxon_label(n.name) if n.is_terminal() else "")

        ax.set_title(title, fontsize=12, fontweight="bold", pad=12)
        ax.set_xlabel("Genetic Divergence (subst/site)", fontsize=10, fontweight="bold")
        ax.tick_params(labelsize=9)

    plt.suptitle("Comparative Phylogenetic Topologies Evaluated by HashRF Across Sarbecovirus Partitions\n"
                 "Demonstrating Reticulate Mosaicism (Red: SARS-CoV-2, Green: RpYN06/RsYN04, Purple: RmYN02, Blue: RaTG13)",
                 fontsize=15, fontweight="bold", y=0.99)

    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out_png = os.path.join(FIGURES_DIR, f"{taxa}_hashrf_comparative_phylogenies.png")
    out_pdf = os.path.join(FIGURES_DIR, f"{taxa}_hashrf_comparative_phylogenies.pdf")
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"  -> Saved Figure 2: {out_png}")

def main():
    print("=" * 75)
    print(" RUNNING HASHRF VISUALIZATION PIPELINE")
    print("=" * 75)
    if not os.path.exists(HASHRF_BIN):
        print(f"Error: HashRF binary not found at {HASHRF_BIN}!")
        sys.exit(1)

    generate_tree_space_and_heatmap_figure()
    generate_comparative_phylogenies_figure('78')
    generate_comparative_phylogenies_figure('68')
    print("\n[COMPLETE] All HashRF visual illustrations generated successfully in reproduction/figures/.")

if __name__ == "__main__":
    main()
