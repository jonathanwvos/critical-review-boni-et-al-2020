#!/usr/bin/env python3
"""
================================================================================
Script 07: Heatmap Visualization of Wuhan-Hu-1 Clade Resolution & Support
Target:      Boni et al. (2020) Sarbecovirus Recombination & Origins Audit
Description: 
  Produces publication-quality heatmaps of the phylogenetic placement, sister
  lineage resolution, and ultrafast bootstrap (UFBoot) support for SARS-CoV-2
  (Wuhan-Hu-1) across all genomic partitions.
  
  Outputs:
    1. reproduction/figures/wuhan_clade_resolution_composite_heatmap.{png,pdf}
       (Dual-panel composite: Panel A = Dataset vs Partition, Panel B = Lineage-resolved matrix)
    2. reproduction/figures/wuhan_clade_support_heatmap.{png,pdf}
       (Primary Dataset vs Partition heatmap with full sister clade annotations)
    3. reproduction/figures/wuhan_clade_sister_taxa_heatmap.{png,pdf}
       (Comprehensive differential analysis: 68-taxa vs. 78-taxa lineage support matrices)
================================================================================
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_DIR = os.path.join(BASE_DIR, "results")
FIGURES_DIR = os.path.join(BASE_DIR, "figures")

DEFAULT_INPUT_CSV = os.path.join(RESULTS_DIR, "wuhan_clade_resolution.csv")

PARTITION_ORDER = [
    "Region_B",
    "Region_C",
    "Region_Aprime",
    "NRR1_pruned",
    "NRR1_full",
    "NRR2",
    "NRA3",
    "Spike"
]

PARTITION_LABELS = [
    "Region B (ORF1a: nsp3–4)",
    "Region C (ORF1a: nsp4–6)",
    "Region A' (ORF1b: nsp12–15)",
    "NRR1 (Pruned: 63/73 taxa)",
    "NRR1 (Full: 68/78 taxa)",
    "NRR2 (GARD/3SEQ candidate)",
    "NRR3 / NRA3 (RDP5 candidate)",
    "Spike Gene (Glycoprotein)"
]

DATASETS = ["68_taxa", "78_taxa"]
DATASET_LABELS = [
    "68-Taxon Dataset\n(Boni et al. Original)",
    "78-Taxon Dataset\n(Updated Lineages)"
]

SISTER_TAXA_78 = [
    "RpYN06",
    "RmYN02",
    "PrC31",
    "RaTG13",
    "Pangolin_GD",
    "Pangolin_GX",
    "CoVZC45",
    "CoVZXC21"
]

SISTER_LABELS_78 = [
    "RpYN06\n(Bat pusillus)",
    "RmYN02\n(Bat malayanus)",
    "PrC31\n(Bat pusillus)",
    "RaTG13\n(Bat affinis)",
    "Pangolin GD\n(Guangdong)",
    "Pangolin GX\n(GX/P1E/2017)",
    "CoVZC45\n(Bat Zhoushan)",
    "CoVZXC21\n(Bat Zhoushan)"
]

SISTER_TAXA_68 = [
    "RaTG13",
    "Pangolin_GD",
    "Pangolin_GX",
    "CoVZC45",
    "CoVZXC21"
]

SISTER_LABELS_68 = [
    "RaTG13\n(Bat affinis)",
    "Pangolin GD\n(Guangdong)",
    "Pangolin GX\n(GX/P1E/2017)",
    "CoVZC45\n(Bat Zhoushan)",
    "CoVZXC21\n(Bat Zhoushan)"
]


def format_sister_text(sister_str):
    """Formats immediate sister clade taxa string cleanly for inside-cell annotation."""
    clean = sister_str.replace("Wuhan-Hu-1, ", "").replace("Wuhan-Hu-1", "Monophyletic")
    if len(clean) > 28:
        parts = clean.split(", ")
        mid = len(parts) // 2
        clean = ", ".join(parts[:mid]) + ",\n" + ", ".join(parts[mid:])
    return clean


def generate_primary_support_heatmap(df, out_dir):
    """Generates standalone primary heatmap: Partitions (Y) vs Datasets (X) colored by UFBoot."""
    mat_support = np.full((len(PARTITION_ORDER), len(DATASETS)), np.nan)
    for i, part in enumerate(PARTITION_ORDER):
        for j, ds in enumerate(DATASETS):
            match = df[(df["Dataset"] == ds) & (df["Partition"] == part)]
            if not match.empty:
                mat_support[i, j] = match.iloc[0]["Clade_UFBoot_Support"]

    fig, ax = plt.subplots(figsize=(11.5, 9.5), dpi=300)
    cmap = plt.cm.YlGnBu
    norm = Normalize(vmin=90.0, vmax=100.0)

    mat_masked = np.ma.masked_invalid(mat_support)
    im = ax.imshow(mat_masked, cmap=cmap, norm=norm, aspect="auto")

    for i, part in enumerate(PARTITION_ORDER):
        for j, ds in enumerate(DATASETS):
            match = df[(df["Dataset"] == ds) & (df["Partition"] == part)]
            if match.empty:
                rect = patches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    linewidth=1,
                    edgecolor="#B0B0B0",
                    facecolor="#F0F0F0",
                    hatch="///"
                )
                ax.add_patch(rect)
                ax.text(
                    j, i,
                    "Not in Dataset\n(N/A*)",
                    ha="center", va="center",
                    color="#666666",
                    fontsize=10,
                    fontweight="bold",
                    fontstyle="italic"
                )
            else:
                row = match.iloc[0]
                supp = int(row["Clade_UFBoot_Support"])
                sister_fmt = format_sister_text(row["Immediate_Sister_Clade_Taxa"])

                text_color = "white" if supp > 96 else "#111111"
                ax.text(
                    j, i,
                    f"{supp}% UFBoot\n(Sister: {sister_fmt})",
                    ha="center", va="center",
                    color=text_color,
                    fontsize=10,
                    fontweight="bold"
                )

    ax.set_xticks(range(len(DATASETS)))
    ax.set_xticklabels(DATASET_LABELS, fontsize=11, fontweight="bold")
    ax.set_yticks(range(len(PARTITION_ORDER)))
    ax.set_yticklabels(PARTITION_LABELS, fontsize=11, fontweight="bold")
    ax.set_ylabel("Genomic Partition (5' to 3')", fontsize=12, fontweight="bold", labelpad=8)
    ax.set_xlabel("Taxa Dataset", fontsize=12, fontweight="bold", labelpad=8)

    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Clade Ultrafast Bootstrap (UFBoot) Support (%)", fontsize=11, fontweight="bold", labelpad=10)
    cbar.ax.tick_params(labelsize=10)

    plt.title(
        "SARS-CoV-2 (Wuhan-Hu-1) Clade Support & Immediate Sister Lineages\n"
        "Across Genomic Partitions (68-Taxon vs. 78-Taxon Datasets)",
        fontsize=13.5, fontweight="bold", pad=15
    )

    plt.tight_layout()
    out_png = os.path.join(out_dir, "wuhan_clade_support_heatmap.png")
    out_pdf = os.path.join(out_dir, "wuhan_clade_support_heatmap.pdf")
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"  [SAVED] Primary clade support heatmap PNG: {out_png}")
    print(f"  [SAVED] Primary clade support heatmap PDF: {out_pdf}")


def generate_sister_taxa_matrix_heatmap(df, out_dir):
    """Generates dual-panel comparative heatmap of Wuhan-Hu-1 MRCA support for 68-taxa and 78-taxa."""
    fig = plt.figure(figsize=(24, 11), dpi=300)
    gs = gridspec.GridSpec(1, 3, width_ratios=[5.8, 8.8, 0.35], wspace=0.18)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    cax = fig.add_subplot(gs[0, 2])

    cmap = plt.cm.YlGnBu
    norm = Normalize(vmin=50.0, vmax=100.0)

    # ---------------- Panel A: 68-Taxon Dataset ----------------
    df_68 = df[df["Dataset"] == "68_taxa"]
    mat_68 = np.full((len(PARTITION_ORDER), len(SISTER_TAXA_68)), np.nan)

    for i, part in enumerate(PARTITION_ORDER):
        match = df_68[df_68["Partition"] == part]
        if match.empty:
            continue
        row = match.iloc[0]
        for j, tax in enumerate(SISTER_TAXA_68):
            supp_val = row.get(f"MRCA_{tax}_Support")
            if pd.notna(supp_val):
                mat_68[i, j] = float(supp_val)

    mat_68_masked = np.ma.masked_invalid(mat_68)
    im1 = ax1.imshow(mat_68_masked, cmap=cmap, norm=norm, aspect="auto")

    for i, part in enumerate(PARTITION_ORDER):
        match = df_68[df_68["Partition"] == part]
        row = match.iloc[0] if not match.empty else None
        for j, tax in enumerate(SISTER_TAXA_68):
            val = mat_68[i, j]
            if np.isnan(val):
                rect = patches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    linewidth=1,
                    edgecolor="#CCCCCC",
                    facecolor="#F4F4F4",
                    hatch="///"
                )
                ax1.add_patch(rect)
                ax1.text(j, i, "Pruned", ha="center", va="center", color="#888888", fontsize=10, fontweight="bold")
            else:
                cnt = int(row.get(f"MRCA_{tax}_Taxa_Count", 0)) if row is not None else 0
                is_sis = bool(row.get(f"Is_{tax}_Sister", False)) if row is not None else False
                text_color = "white" if val > 85 else "#111111"
                if is_sis:
                    label = f"{int(val)}% ★\n(Sister, n={cnt})"
                    ax1.text(j, i, label, ha="center", va="center", color=text_color, fontsize=10.5, fontweight="bold")
                else:
                    label = f"{int(val)}%\n(n={cnt})"
                    ax1.text(j, i, label, ha="center", va="center", color=text_color, fontsize=10, fontweight="bold")

    ax1.set_xticks(range(len(SISTER_TAXA_68)))
    ax1.set_xticklabels(SISTER_LABELS_68, fontsize=10, fontweight="bold")
    ax1.set_yticks(range(len(PARTITION_ORDER)))
    ax1.set_yticklabels(PARTITION_LABELS, fontsize=11, fontweight="bold")
    ax1.set_ylabel("Genomic Partition (5' to 3')", fontsize=12, fontweight="bold", labelpad=8)
    ax1.set_xlabel("Comparative Lineage", fontsize=11, fontweight="bold", labelpad=8)
    ax1.set_title("A. 68-Taxon Dataset (Original Boni et al.)", fontsize=13, fontweight="bold", pad=12)

    # ---------------- Panel B: 78-Taxon Dataset ----------------
    df_78 = df[df["Dataset"] == "78_taxa"]
    mat_78 = np.full((len(PARTITION_ORDER), len(SISTER_TAXA_78)), np.nan)

    for i, part in enumerate(PARTITION_ORDER):
        match = df_78[df_78["Partition"] == part]
        if match.empty:
            continue
        row = match.iloc[0]
        for j, tax in enumerate(SISTER_TAXA_78):
            supp_val = row.get(f"MRCA_{tax}_Support")
            if pd.notna(supp_val):
                mat_78[i, j] = float(supp_val)

    mat_78_masked = np.ma.masked_invalid(mat_78)
    im2 = ax2.imshow(mat_78_masked, cmap=cmap, norm=norm, aspect="auto")

    for i, part in enumerate(PARTITION_ORDER):
        match = df_78[df_78["Partition"] == part]
        row = match.iloc[0] if not match.empty else None
        for j, tax in enumerate(SISTER_TAXA_78):
            val = mat_78[i, j]
            if np.isnan(val):
                rect = patches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    linewidth=1,
                    edgecolor="#CCCCCC",
                    facecolor="#F4F4F4",
                    hatch="///"
                )
                ax2.add_patch(rect)
                ax2.text(j, i, "Pruned", ha="center", va="center", color="#888888", fontsize=10, fontweight="bold")
            else:
                cnt = int(row.get(f"MRCA_{tax}_Taxa_Count", 0)) if row is not None else 0
                is_sis = bool(row.get(f"Is_{tax}_Sister", False)) if row is not None else False
                text_color = "white" if val > 85 else "#111111"
                if is_sis:
                    label = f"{int(val)}% ★\n(Sister, n={cnt})"
                    ax2.text(j, i, label, ha="center", va="center", color=text_color, fontsize=10.5, fontweight="bold")
                else:
                    label = f"{int(val)}%\n(n={cnt})"
                    ax2.text(j, i, label, ha="center", va="center", color=text_color, fontsize=10, fontweight="bold")

    ax2.set_xticks(range(len(SISTER_TAXA_78)))
    ax2.set_xticklabels(SISTER_LABELS_78, fontsize=10, fontweight="bold")
    ax2.tick_params(left=False, labelleft=False)
    ax2.set_xlabel("Comparative Lineage", fontsize=11, fontweight="bold", labelpad=8)
    ax2.set_title("B. 78-Taxon Dataset (Updated with Post-2020 Lineages)", fontsize=13, fontweight="bold", pad=12)

    # Shared colorbar
    cbar = plt.colorbar(im1, cax=cax)
    cbar.set_label("MRCA Ultrafast Bootstrap (UFBoot) Support with SARS-CoV-2 (%)", fontsize=11, fontweight="bold", labelpad=10)
    cbar.ax.tick_params(labelsize=10)

    fig.suptitle(
        "Comparative Differential Analysis of SARS-CoV-2 Clade Support and Node Topology\n"
        "Across Genomic Partitions in the Original (68-Taxon) vs. Updated (78-Taxon) Datasets",
        fontsize=14.5, fontweight="bold", y=0.98
    )

    fig.text(
        0.05, 0.015,
        "★ Denotes immediate sister clade resolution with Wuhan-Hu-1. 'n' indicates the total number of taxa subsumed by the MRCA node (n=2 indicates exclusive sister pairing).\n"
        "Shaded cells denote taxa pruned from that alignment (e.g. CoVZC45/CoVZXC21 pruned in NRR1_pruned). Source: wuhan_clade_resolution.csv.",
        fontsize=9.5, style="italic", color="#333333"
    )

    fig.subplots_adjust(top=0.91, bottom=0.07, left=0.10, right=0.96)
    out_png = os.path.join(out_dir, "wuhan_clade_sister_taxa_heatmap.png")
    out_pdf = os.path.join(out_dir, "wuhan_clade_sister_taxa_heatmap.pdf")
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"  [SAVED] Comprehensive sister taxa support matrix PNG: {out_png}")
    print(f"  [SAVED] Comprehensive sister taxa support matrix PDF: {out_pdf}")


def generate_composite_heatmap(df, out_dir):
    """Generates publication dual-panel composite figure combining Panel A (cohort) and Panel B (78-taxa matrix)."""
    fig = plt.figure(figsize=(24, 11), dpi=300)
    gs = gridspec.GridSpec(1, 3, width_ratios=[7.2, 9.8, 0.35], wspace=0.18)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    cax = fig.add_subplot(gs[0, 2])

    cmap = plt.cm.YlGnBu
    norm = Normalize(vmin=50.0, vmax=100.0)

    # Panel A: Dataset vs Partition
    mat1 = np.full((len(PARTITION_ORDER), len(DATASETS)), np.nan)
    for i, part in enumerate(PARTITION_ORDER):
        for j, ds in enumerate(DATASETS):
            match = df[(df["Dataset"] == ds) & (df["Partition"] == part)]
            if not match.empty:
                mat1[i, j] = match.iloc[0]["Clade_UFBoot_Support"]

    mat1_masked = np.ma.masked_invalid(mat1)
    im1 = ax1.imshow(mat1_masked, cmap=cmap, norm=norm, aspect="auto")

    for i, part in enumerate(PARTITION_ORDER):
        for j, ds in enumerate(DATASETS):
            match = df[(df["Dataset"] == ds) & (df["Partition"] == part)]
            if match.empty:
                rect = patches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    linewidth=1,
                    edgecolor="#B0B0B0",
                    facecolor="#F0F0F0",
                    hatch="///"
                )
                ax1.add_patch(rect)
                ax1.text(
                    j, i,
                    "Not in Dataset\n(N/A*)",
                    ha="center", va="center",
                    color="#666666",
                    fontsize=9.5,
                    fontweight="bold",
                    fontstyle="italic"
                )
            else:
                row = match.iloc[0]
                supp = int(row["Clade_UFBoot_Support"])
                sister_fmt = format_sister_text(row["Immediate_Sister_Clade_Taxa"])

                text_color = "white" if supp > 85 else "#111111"
                ax1.text(
                    j, i,
                    f"{supp}% UFBoot\n({sister_fmt})",
                    ha="center", va="center",
                    color=text_color,
                    fontsize=9.5,
                    fontweight="bold"
                )

    ax1.set_xticks(range(len(DATASETS)))
    ax1.set_xticklabels(DATASET_LABELS, fontsize=10.5, fontweight="bold")
    ax1.set_yticks(range(len(PARTITION_ORDER)))
    ax1.set_yticklabels(PARTITION_LABELS, fontsize=11, fontweight="bold")
    ax1.set_title("A. Immediate Sister Lineage Resolution & Support", fontsize=12.5, fontweight="bold", pad=12)
    ax1.set_ylabel("Genomic Partition (5' to 3')", fontsize=11, fontweight="bold")
    ax1.set_xlabel("Taxa Dataset", fontsize=11, fontweight="bold", labelpad=8)

    # Panel B: Full 78-taxon lineage-resolved matrix
    df_78 = df[df["Dataset"] == "78_taxa"]
    mat2 = np.full((len(PARTITION_ORDER), len(SISTER_TAXA_78)), np.nan)

    for i, part in enumerate(PARTITION_ORDER):
        match = df_78[df_78["Partition"] == part]
        if match.empty:
            continue
        row = match.iloc[0]
        for j, tax in enumerate(SISTER_TAXA_78):
            supp_val = row.get(f"MRCA_{tax}_Support")
            if pd.notna(supp_val):
                mat2[i, j] = float(supp_val)

    mat2_masked = np.ma.masked_invalid(mat2)
    im2 = ax2.imshow(mat2_masked, cmap=cmap, norm=norm, aspect="auto")

    for i, part in enumerate(PARTITION_ORDER):
        match = df_78[df_78["Partition"] == part]
        row = match.iloc[0] if not match.empty else None
        for j, tax in enumerate(SISTER_TAXA_78):
            val = mat2[i, j]
            if np.isnan(val):
                rect = patches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    linewidth=1,
                    edgecolor="#CCCCCC",
                    facecolor="#F4F4F4",
                    hatch="///"
                )
                ax2.add_patch(rect)
                ax2.text(j, i, "Pruned", ha="center", va="center", color="#888888", fontsize=10, fontweight="bold")
            else:
                cnt = int(row.get(f"MRCA_{tax}_Taxa_Count", 0)) if row is not None else 0
                is_sis = bool(row.get(f"Is_{tax}_Sister", False)) if row is not None else False
                text_color = "white" if val > 85 else "#111111"
                if is_sis:
                    label = f"{int(val)}% ★\n(Sister, n={cnt})"
                    ax2.text(j, i, label, ha="center", va="center", color=text_color, fontsize=10.5, fontweight="heavy")
                else:
                    label = f"{int(val)}%\n(n={cnt})"
                    ax2.text(j, i, label, ha="center", va="center", color=text_color, fontsize=10, fontweight="bold")

    ax2.set_xticks(range(len(SISTER_TAXA_78)))
    ax2.set_xticklabels(SISTER_LABELS_78, fontsize=10, fontweight="bold")
    ax2.tick_params(left=False, labelleft=False)
    ax2.set_title("B. 78-Taxon Dataset: Lineage-Resolved MRCA Support & Clade Size", fontsize=12.5, fontweight="bold", pad=12)
    ax2.set_xlabel("Comparative Sarbecovirus Lineage", fontsize=11, fontweight="bold", labelpad=8)

    # Shared colorbar
    cbar = plt.colorbar(im1, cax=cax)
    cbar.set_label("MRCA Ultrafast Bootstrap (UFBoot) Support with SARS-CoV-2 (%)", fontsize=11, fontweight="bold", labelpad=10)
    cbar.ax.tick_params(labelsize=10)

    fig.suptitle(
        "Phylogenetic Resolution & Ultrafast Bootstrap Support for SARS-CoV-2 (Wuhan-Hu-1)\n"
        "Demonstrating the Complete Disruption of the RaTG13 Sister Relationship Across ORF1ab and NRR3 in 78 Taxa",
        fontsize=14, fontweight="bold", y=0.98
    )

    fig.text(
        0.05, 0.015,
        "★ Denotes immediate sister clade resolution with Wuhan-Hu-1. 'n' indicates the total number of taxa subsumed by the MRCA node.\n"
        "Source: reproduction/results/wuhan_clade_resolution.csv (derived from IQ-TREE v2.2.6 ML phylogenies under GTR+G4 with 1,000 UFBoot replicates).",
        fontsize=9.5, style="italic", color="#444444"
    )

    fig.subplots_adjust(top=0.90, bottom=0.08, left=0.08, right=0.95, wspace=0.22)
    out_png = os.path.join(out_dir, "wuhan_clade_resolution_composite_heatmap.png")
    out_pdf = os.path.join(out_dir, "wuhan_clade_resolution_composite_heatmap.pdf")
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"  [SAVED] Composite clade resolution heatmap PNG: {out_png}")
    print(f"  [SAVED] Composite clade resolution heatmap PDF: {out_pdf}")


def main():
    parser = argparse.ArgumentParser(description="Generate heatmaps of Wuhan-Hu-1 clade support and resolution.")
    parser.add_argument("--input", default=DEFAULT_INPUT_CSV, help="Path to wuhan_clade_resolution.csv")
    parser.add_argument("--outdir", default=FIGURES_DIR, help="Directory to save generated figures")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"[ERROR] Input CSV file not found at: {args.input}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.input)
    print(f"\n=======================================================")
    print(f" Generating Wuhan-Hu-1 Clade Resolution Heatmaps")
    print(f" Source: {args.input}")
    print(f" Output Directory: {args.outdir}")
    print(f"=======================================================\n")

    generate_primary_support_heatmap(df, args.outdir)
    generate_sister_taxa_matrix_heatmap(df, args.outdir)
    generate_composite_heatmap(df, args.outdir)

    print("\n[COMPLETE] All Wuhan-Hu-1 clade resolution visualizations generated successfully.")


if __name__ == "__main__":
    main()
