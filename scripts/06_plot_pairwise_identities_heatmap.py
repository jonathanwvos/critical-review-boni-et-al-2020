#!/usr/bin/env python3
"""
================================================================================
Script 06: Heatmap Visualization of Pairwise Sequence Identities to Wuhan-Hu-1
Target:      Boni et al. (2020) Sarbecovirus Recombination & Origins Audit
Description: 
  Produces publication-quality heatmaps of pairwise nucleotide sequence
  identities between SARS-CoV-2 (Wuhan-Hu-1) and key comparative sarbecoviruses.
  
  Layout:
    - Vertical Axis (Y): Genomic partitions (Region B, Region C, Region A',
      NRR1_pruned, NRR1_full, NRR2, Spike Gene)
    - Horizontal Axis (X): Comparative sarbecovirus lineages
  
  Outputs:
    1. reproduction/figures/pairwise_identities_heatmap.{png,pdf}
       (Dual-panel comparative figure: Panel A = 78_taxa, Panel B = 68_taxa)
    2. reproduction/figures/pairwise_identities_78_taxa_heatmap.{png,pdf}
       (Standalone 78-taxon updated dataset heatmap)
    3. reproduction/figures/pairwise_identities_68_taxa_heatmap.{png,pdf}
       (Standalone 68-taxon original dataset heatmap)
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

DEFAULT_INPUT_CSV = os.path.join(RESULTS_DIR, "pairwise_identities_wuhan.csv")

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
    "Region B",
    "Region C",
    "Region A'",
    "NRR1 (Pruned)",
    "NRR1 (Full)",
    "NRR2",
    "NRR3 (NRA3)",
    "Spike"
]

TAXA_ORDER_78 = [
    "RpYN06",
    "RmYN02",
    "RaTG13",
    "Pangolin_GD",
    "Pangolin_GX",
    "CoVZC45",
    "CoVZXC21",
    "RsYN04",
    "RmYN05",
    "RmYN08"
]

TAXA_LABELS_78 = [
    "RpYN06",
    "RmYN02",
    "RaTG13",
    "Pangolin\nGD",
    "Pangolin\nGX (P1E)",
    "CoVZC45",
    "CoVZXC21",
    "RsYN04",
    "RmYN05",
    "RmYN08"
]

TAXA_ORDER_68 = [
    "RaTG13",
    "Pangolin_GD",
    "Pangolin_GX",
    "CoVZC45",
    "CoVZXC21"
]

TAXA_LABELS_68 = [
    "RaTG13",
    "Pangolin\nGD",
    "Pangolin\nGX (P1E)",
    "CoVZC45",
    "CoVZXC21"
]


def load_and_pivot_matrix(df, dataset_name, taxa_order):
    """Filters dataset and returns a 2D numpy array shaped (n_partitions, n_taxa)."""
    df_sub = df[df["Dataset"] == dataset_name]
    pivot = df_sub.pivot(index="Partition", columns="Compared_Taxon", values="Percent_Identity")
    reindexed = pivot.reindex(index=PARTITION_ORDER, columns=taxa_order)
    return reindexed.to_numpy()


def render_heatmap_cells(ax, mat, taxa_labels, cmap, norm, highlight_max=True, font_size=10.5):
    """Draws heatmap raster, annotations, and handles NaN/pruned cells."""
    mat_masked = np.ma.masked_invalid(mat)
    im = ax.imshow(mat_masked, cmap=cmap, norm=norm, aspect="auto")

    n_rows, n_cols = mat.shape
    for i in range(n_rows):
        row_vals = mat[i, :]
        row_max = np.nanmax(row_vals) if not np.all(np.isnan(row_vals)) else None

        for j in range(n_cols):
            val = mat[i, j]
            if np.isnan(val):
                # Pruned / excluded cell
                rect = patches.Rectangle(
                    (j - 0.5, i - 0.5), 1, 1,
                    linewidth=1.2,
                    edgecolor="#AAAAAA",
                    facecolor="#E8E8E8",
                    hatch="///"
                )
                ax.add_patch(rect)
                cell_text = "Not in 78-Taxa\n(N/A*)" if i == PARTITION_ORDER.index("NRA3") else "Pruned\n(N/A*)"
                ax.text(
                    j, i, cell_text,
                    ha="center", va="center",
                    color="#555555",
                    fontsize=font_size - 2.0,
                    fontweight="bold",
                    style="italic"
                )
            else:
                is_max = (val == row_max) and highlight_max
                # Dynamic text color based on colormap luminance
                text_color = "white" if val > 88.0 else "#111111"
                weight = "black" if is_max else "bold"

                ax.text(
                    j, i, f"{val:.2f}%",
                    ha="center", va="center",
                    color=text_color,
                    fontsize=font_size,
                    fontweight=weight
                )

                if is_max:
                    # Red highlight border around highest identity in row
                    rect = patches.Rectangle(
                        (j - 0.48, i - 0.48), 0.96, 0.96,
                        linewidth=2.2,
                        edgecolor="#E63946",
                        facecolor="none"
                    )
                    ax.add_patch(rect)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(taxa_labels, fontsize=10.5, fontweight="bold")
    return im


def generate_dual_panel_heatmap(mat_78, mat_68, out_dir):
    """Generates Figure: Dual-panel comparative heatmap (78-taxa & 68-taxa)."""
    fig = plt.figure(figsize=(21.5, 10.5), dpi=300)
    gs = gridspec.GridSpec(1, 3, width_ratios=[10.0, 5.0, 0.38], wspace=0.18)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
    cax = fig.add_subplot(gs[0, 2])

    cmap = plt.cm.YlGnBu
    norm = Normalize(vmin=65.0, vmax=100.0)

    # Panel A: 78_taxa
    im1 = render_heatmap_cells(ax1, mat_78, TAXA_LABELS_78, cmap, norm, highlight_max=True, font_size=10.5)
    ax1.set_yticks(range(len(PARTITION_ORDER)))
    ax1.set_yticklabels(PARTITION_LABELS, fontsize=11, fontweight="bold")
    ax1.set_title("a.", fontsize=13, fontweight="bold", pad=12, loc="left")
    # ax1.set_xlabel("Compared Comparative Taxa", fontsize=11, fontweight="bold", labelpad=8)
    ax1.set_ylabel("Genomic Partition (5' to 3')", fontsize=11, fontweight="bold", labelpad=8)

    # Panel B: 68_taxa
    im2 = render_heatmap_cells(ax2, mat_68, TAXA_LABELS_68, cmap, norm, highlight_max=True, font_size=10.5)
    ax2.tick_params(left=False, labelleft=False)
    ax2.set_title("b.", fontsize=13, fontweight="bold", pad=12, loc="left")

    # ax2.set_xlabel("Compared Comparative Taxa", fontsize=11, fontweight="bold", labelpad=8)

    # Colorbar
    cbar = plt.colorbar(im1, cax=cax)
    cbar.set_label("Pairwise Nucleotide Identity to SARS-CoV-2 / Wuhan-Hu-1 (%)", fontsize=11, fontweight="bold", labelpad=12)
    cbar.ax.tick_params(labelsize=10)

    # Title & Notes
    # fig.suptitle(
    #     "Pairwise Sequence Identity to SARS-CoV-2 (Wuhan-Hu-1) Across Genomic Partitions\n"
    #     "Red bounding boxes highlight the closest relative (highest pairwise identity) per genomic region",
    #     fontsize=14.5, fontweight="bold", y=0.98
    # )

    # fig.text(
    #     0.05, 0.02,
    #     "* Note: CoVZC45 and CoVZXC21 were pruned out of the NRR1_pruned partition to enforce tree clonality (Boni et al. 2020). "
    #     "Their identities across the full concatenated NRR1 partition are shown under NRR1 (Full).\n"
    #     "Data source: reproduction/results/pairwise_identities_wuhan.csv (calculated over non-gap pairwise sites).",
    #     fontsize=9, style="italic", color="#444444"
    # )

    out_png = os.path.join(out_dir, "pairwise_identities_heatmap.png")
    out_pdf = os.path.join(out_dir, "pairwise_identities_heatmap.pdf")
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"  [SAVED] Combined dual-panel heatmap PNG: {out_png}")
    print(f"  [SAVED] Combined dual-panel heatmap PDF: {out_pdf}")


def generate_single_dataset_heatmap(mat, taxa_labels, dataset_title, filename_prefix, out_dir):
    """Generates standalone publication heatmap for an individual dataset."""
    n_cols = mat.shape[1]
    fig_width = 13.5 if n_cols > 5 else 9.0
    fig, ax = plt.subplots(figsize=(fig_width, 9.5), dpi=300)

    cmap = plt.cm.YlGnBu
    norm = Normalize(vmin=65.0, vmax=100.0)

    im = render_heatmap_cells(ax, mat, taxa_labels, cmap, norm, highlight_max=True, font_size=11)
    ax.set_yticks(range(len(PARTITION_ORDER)))
    ax.set_yticklabels(PARTITION_LABELS, fontsize=11, fontweight="bold")
    ax.set_ylabel("Genomic Partition (5' to 3')", fontsize=12, fontweight="bold", labelpad=8)
    # ax.set_xlabel("Compared Comparative Taxa", fontsize=12, fontweight="bold", labelpad=8)

    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    cbar.set_label("Pairwise Sequence Identity to Wuhan-Hu-1 (%)", fontsize=11, fontweight="bold", labelpad=10)
    cbar.ax.tick_params(labelsize=10)

    # plt.title(
    #     f"{dataset_title}\nPairwise Sequence Identity (%) to SARS-CoV-2 (Wuhan-Hu-1)",
    #     fontsize=13.5, fontweight="bold", pad=15
    # )

    plt.tight_layout()
    out_png = os.path.join(out_dir, f"{filename_prefix}.png")
    out_pdf = os.path.join(out_dir, f"{filename_prefix}.pdf")
    plt.savefig(out_png, bbox_inches="tight", dpi=300)
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"  [SAVED] Standalone heatmap PNG: {out_png}")
    print(f"  [SAVED] Standalone heatmap PDF: {out_pdf}")


def main():
    parser = argparse.ArgumentParser(description="Generate pairwise sequence identity heatmaps to SARS-CoV-2 (Wuhan-Hu-1).")
    parser.add_argument("--input", default=DEFAULT_INPUT_CSV, help="Path to pairwise_identities_wuhan.csv")
    parser.add_argument("--outdir", default=FIGURES_DIR, help="Directory to save generated figures")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"[ERROR] Input CSV file not found at: {args.input}")
        sys.exit(1)

    os.makedirs(args.outdir, exist_ok=True)
    df = pd.read_csv(args.input)
    print(f"\n=======================================================")
    print(f" Generating Pairwise Identity Heatmaps")
    print(f" Source: {args.input}")
    print(f" Output Directory: {args.outdir}")
    print(f"=======================================================")

    mat_78 = load_and_pivot_matrix(df, "78_taxa", TAXA_ORDER_78)
    mat_68 = load_and_pivot_matrix(df, "68_taxa", TAXA_ORDER_68)

    # 1. Dual-panel combined figure
    generate_dual_panel_heatmap(mat_78, mat_68, args.outdir)

    # 2. Standalone 78-taxon figure
    generate_single_dataset_heatmap(
        mat_78, TAXA_LABELS_78,
        "78-Taxon Dataset (Updated Sarbecovirus Lineages)",
        "pairwise_identities_78_taxa_heatmap",
        args.outdir
    )

    # 3. Standalone 68-taxon figure
    generate_single_dataset_heatmap(
        mat_68, TAXA_LABELS_68,
        "68-Taxon Dataset (Boni et al. 2020 Original)",
        "pairwise_identities_68_taxa_heatmap",
        args.outdir
    )

    print("\n[COMPLETE] All heatmap visualizations generated successfully.")


if __name__ == "__main__":
    main()
