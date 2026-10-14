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

PARTITIONS = ["Region_B", "Region_C", "Region_Aprime", "NRR1_full", "NRR2", "Spike"]
PARTITION_LABELS = {
    "Region_B": "Region B (ORF1a)",
    "Region_C": "Region C (ORF1a)",
    "Region_Aprime": "Region A' (ORF1b)",
    "NRR1_full": "NRR1 Concatenated",
    "NRR2": "NRR2 (GARD/3SEQ)",
    "Spike": "Spike Gene"
}

def run_hashrf_for_dataset(dataset_name, weighted=False, print_opt="matrix"):
    """Runs HashRF on a combined multi-tree Newick file and returns parsed numpy matrix."""
    ds_dir = os.path.join(TREES_DIR, dataset_name)
    combined_path = os.path.join(ds_dir, f"hashrf_input_{dataset_name}.tre")
    
    # Concatenate trees
    with open(combined_path, "w") as out_f:
        for p in PARTITIONS:
            tf = os.path.join(ds_dir, f"{p}.treefile")
            with open(tf) as in_f:
                out_f.write(in_f.read().strip() + "\n")
                
    out_file = os.path.join(RESULTS_DIR, f"hashrf_{dataset_name}_{'weighted' if weighted else 'unweighted'}_{print_opt}.txt")
    
    cmd = [HASHRF_BIN, combined_path, str(len(PARTITIONS)), "-p", print_opt, "-o", out_file]
    if weighted:
        cmd.append("-w")
        
    subprocess.run(cmd, capture_output=True, text=True)
    
    # Parse matrix file
    matrix_rows = []
    with open(out_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if parts and all(p.replace('.', '', 1).replace('e-', '', 1).replace('-', '', 1).isdigit() for p in parts):
                matrix_rows.append([float(x) for x in parts])
                
    mat = np.array(matrix_rows)
    return mat

def compute_classical_mds(dist_matrix):
    """Classical multidimensional scaling (PCoA) to project tree distances onto 2D."""
    n = dist_matrix.shape[0]
    # Centering matrix
    H = np.eye(n) - np.ones((n, n)) / n
    # Double centering
    B = -0.5 * H.dot(dist_matrix ** 2).dot(H)
    # Eigendecomposition
    eigenvals, eigenvecs = np.linalg.eigh(B)
    # Sort descending
    idx = np.argsort(eigenvals)[::-1]
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]
    # Keep top 2 dimensions
    coords = eigenvecs[:, :2] * np.sqrt(np.maximum(eigenvals[:2], 0))
    var_explained = np.maximum(eigenvals[:2], 0) / np.sum(np.maximum(eigenvals, 0)) * 100
    return coords, var_explained

def generate_tree_space_and_heatmap_figure():
    """Generates Figure 1: HashRF Robinson-Foulds Tree Space & Incongruence Heatmaps."""
    print("[*] Generating Figure 1: HashRF Tree Space and Heatmaps...")
    
    # Get HashRF matrices for 68_taxa and 78_taxa
    mat_68 = run_hashrf_for_dataset("68_taxa", weighted=False, print_opt="matrix")
    mat_78 = run_hashrf_for_dataset("78_taxa", weighted=False, print_opt="matrix")
    mat_68_w = run_hashrf_for_dataset("68_taxa", weighted=True, print_opt="matrix")
    
    # Compute normalized RF rates
    # 68 taxa: max splits = 68 - 3 = 65. 78 taxa: max splits = 78 - 3 = 75.
    rate_68 = (mat_68 / 65.0) * 100
    rate_78 = (mat_78 / 75.0) * 100
    
    # Classical MDS on 68-taxa and 78-taxa
    coords_68, var_68 = compute_classical_mds(rate_68)
    
    fig = plt.figure(figsize=(19, 14), dpi=300)
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[1, 1], hspace=0.32, wspace=0.25)
    
    labels = [PARTITION_LABELS[p] for p in PARTITIONS]
    short_labels = ["Reg B", "Reg C", "Reg A'", "NRR1", "NRR2", "Spike"]
    
    # -------------------------------------------------------------
    # Panel A: Heatmap 68_taxa (Original Dataset)
    # -------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0, 0])
    im_a = ax_a.imshow(rate_68, cmap="YlOrRd", vmin=0, vmax=80)
    ax_a.set_xticks(range(len(short_labels)))
    ax_a.set_yticks(range(len(short_labels)))
    ax_a.set_xticklabels(short_labels, fontsize=11, fontweight="bold")
    ax_a.set_yticklabels(labels, fontsize=11)
    ax_a.set_title("A. HashRF Incongruence Matrix (68-Taxa Original Dataset)\nValues: One-Sided RF Splits [Normalized % Discordance]", fontsize=12, fontweight="bold", pad=12)
    
    # Annotate cells
    for i in range(len(PARTITIONS)):
        for j in range(len(PARTITIONS)):
            raw_val = int(mat_68[i, j])
            pct_val = rate_68[i, j]
            text_color = "white" if pct_val > 45 else "black"
            ax_a.text(j, i, f"{raw_val}\n[{pct_val:.1f}%]", ha="center", va="center", color=text_color, fontsize=10, fontweight="bold")
            
    cbar_a = plt.colorbar(im_a, ax=ax_a, fraction=0.046, pad=0.04)
    cbar_a.set_label("Normalized RF Discordance (%)", fontsize=10)

    # -------------------------------------------------------------
    # Panel B: Heatmap 78_taxa (Updated Dataset with RpYN06 / RmYN02)
    # -------------------------------------------------------------
    ax_b = fig.add_subplot(gs[0, 1])
    im_b = ax_b.imshow(rate_78, cmap="YlOrRd", vmin=0, vmax=80)
    ax_b.set_xticks(range(len(short_labels)))
    ax_b.set_yticks(range(len(short_labels)))
    ax_b.set_xticklabels(short_labels, fontsize=11, fontweight="bold")
    ax_b.set_yticklabels(short_labels, fontsize=11, fontweight="bold")
    ax_b.set_title("B. HashRF Incongruence Matrix (78-Taxa Updated Dataset)\nValues: One-Sided RF Splits [Normalized % Discordance]", fontsize=12, fontweight="bold", pad=12)
    
    for i in range(len(PARTITIONS)):
        for j in range(len(PARTITIONS)):
            raw_val = int(mat_78[i, j])
            pct_val = rate_78[i, j]
            text_color = "white" if pct_val > 45 else "black"
            ax_b.text(j, i, f"{raw_val}\n[{pct_val:.1f}%]", ha="center", va="center", color=text_color, fontsize=10, fontweight="bold")
            
    cbar_b = plt.colorbar(im_b, ax=ax_b, fraction=0.046, pad=0.04)
    cbar_b.set_label("Normalized RF Discordance (%)", fontsize=10)

    # -------------------------------------------------------------
    # Panel C: 2D Tree Space Projection (MDS of 68-Taxa HashRF Matrix)
    # -------------------------------------------------------------
    ax_c = fig.add_subplot(gs[1, 0])
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
    markers = ["o", "s", "^", "D", "v", "X"]
    
    # Custom annotation offsets to eliminate overlaps
    offsets = [
        (1.5, -2.5),   # Reg B
        (0.0, 3.0),    # Reg C
        (-3.0, 2.8),   # Reg A'
        (3.5, 1.5),    # NRR1
        (3.0, -2.5),   # NRR2
        (2.5, 2.0)     # Spike
    ]
    
    for idx, (p, col, mark) in enumerate(zip(PARTITIONS, colors, markers)):
        ax_c.scatter(coords_68[idx, 0], coords_68[idx, 1], color=col, marker=mark, s=190, edgecolors="black", linewidth=1.5, zorder=5, label=PARTITION_LABELS[p])
        dx, dy = offsets[idx]
        ax_c.annotate(short_labels[idx], (coords_68[idx, 0] + dx, coords_68[idx, 1] + dy), fontsize=11, fontweight="bold")
        
    # Draw dashed lines between Region B, Region C, and Region A'
    b_idx, c_idx, a_idx = 0, 1, 2
    triangle = [b_idx, c_idx, a_idx, b_idx]
    ax_c.plot(coords_68[triangle, 0], coords_68[triangle, 1], linestyle="--", color="gray", linewidth=1.5, alpha=0.7, zorder=3, label="NRR1 Sub-partitions (Incongruent)")
    
    ax_c.set_xlabel(f"Tree Space Coordinate 1 ({var_68[0]:.1f}% variance)", fontsize=11, fontweight="bold")
    ax_c.set_ylabel(f"Tree Space Coordinate 2 ({var_68[1]:.1f}% variance)", fontsize=11, fontweight="bold")
    ax_c.set_title("C. 2D Tree Space Projection (MDS of HashRF Distance Matrix)\nDemonstrates Non-Recombinant Discordance & Artificial NRR1 Compromise", fontsize=12, fontweight="bold", pad=12)
    ax_c.grid(True, linestyle=":", alpha=0.6)
    ax_c.legend(loc="lower left", fontsize=9, framealpha=0.95)
    ax_c.set_ylim(coords_68[:, 1].min() - 8, coords_68[:, 1].max() + 8)
    ax_c.set_xlim(coords_68[:, 0].min() - 8, coords_68[:, 0].max() + 8)

    # -------------------------------------------------------------
    # Panel D: Branch-Length Weighted Robinson-Foulds Matrix (HashRF -w)
    # -------------------------------------------------------------
    ax_d = fig.add_subplot(gs[1, 1])
    im_d = ax_d.imshow(mat_68_w, cmap="Purples", vmin=0, vmax=1.5)
    ax_d.set_xticks(range(len(short_labels)))
    ax_d.set_yticks(range(len(short_labels)))
    ax_d.set_xticklabels(short_labels, fontsize=11, fontweight="bold")
    ax_d.set_yticklabels(short_labels, fontsize=11, fontweight="bold")
    ax_d.set_title("D. Branch-Length Weighted Robinson-Foulds Matrix (HashRF -w)\nTopological Discordance Scaled by Evolutionary Branch Lengths", fontsize=12, fontweight="bold", pad=12)
    
    for i in range(len(PARTITIONS)):
        for j in range(len(PARTITIONS)):
            w_val = mat_68_w[i, j]
            text_color = "white" if w_val > 0.8 else "black"
            ax_d.text(j, i, f"{w_val:.2f}", ha="center", va="center", color=text_color, fontsize=10, fontweight="bold")
            
    cbar_d = plt.colorbar(im_d, ax=ax_d, fraction=0.046, pad=0.04)
    cbar_d.set_label("Weighted RF Distance (Substitutions/Site)", fontsize=10)

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

def generate_comparative_phylogenies_figure():
    """Generates Figure 2: Standardized 12-Taxon Pruned Comparative Phylogenies Evaluated by HashRF."""
    print("[*] Generating Figure 2: Standardized Comparative Phylogenetic Topologies...")
    
    # Pivotal 12 taxa across the Sarbecovirus evolutionary debate
    keywords = [
        "Wuhan", "RaTG13", "RpYN06", "RmYN02", "RsYN04", "RmYN05",
        "ZC45", "ZXC21", "410721", "410538", "HSZ-Cc", "BtKY72"
    ]
    
    eval_partitions = ["Region_B", "Region_C", "Region_Aprime", "Spike"]
    titles = [
        "Region B (ORF1a, nt 3,625-9,150)\nHashRF vs A': 68.6% discordance",
        "Region C (ORF1a, nt 9,261-11,795)\nHashRF vs A': 62.9% discordance",
        "Region A' (ORF1b, nt 13,291-19,628)\nHashRF vs B: 68.6% discordance",
        "Spike Gene (nt 21,563-25,384)\nHashRF vs B: 73.3% discordance"
    ]
    
    fig, axes = plt.subplots(1, 4, figsize=(24, 10), dpi=300)
    
    for idx, (p, title) in enumerate(zip(eval_partitions, titles)):
        ax = axes[idx]
        tf = os.path.join(TREES_DIR, "78_taxa", f"{p}.treefile")
        
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
    out_png = os.path.join(FIGURES_DIR, "hashrf_comparative_phylogenies.png")
    out_pdf = os.path.join(FIGURES_DIR, "hashrf_comparative_phylogenies.pdf")
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
    generate_comparative_phylogenies_figure()
    print("\n[COMPLETE] All HashRF visual illustrations generated successfully in reproduction/figures/.")

if __name__ == "__main__":
    main()
