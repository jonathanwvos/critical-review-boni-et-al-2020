#!/usr/bin/env python3
"""
Master Reproduction Pipeline for Boni et al. (2020) Critical Review.

Executes all reproduction steps in sequential order:
1. 01_prepare_unaligned_data.py - Extracts unaligned sequences for 68-taxon and 78-taxon sets.
2. 02_run_mafft_alignments.py    - Performs de novo MAFFT multiple sequence alignments.
3. 03_run_phylogenetic_inference.py - Infers ML trees using IQ-TREE (GTR+G4, 1000 UFBoot).
4. 04_analyze_nrr1_congruence.py - Calculates metrics, RF distances, AU tests, and clade resolution.
5. 05_visualize_hashrf_trees.py  - Generates HashRF distance matrices, heatmaps, and tree space plots.
6. 06_plot_pairwise_identities_heatmap.py - Generates pairwise sequence identity heatmaps.
7. 07_plot_wuhan_clade_support_heatmap.py - Generates Wuhan clade resolution composite heatmap.
8. 08_audit_nra3_recombination.py - Conducts 4-gamete tests, sub-block AU tests, and NRA3 audit.
"""

import os
import sys
import subprocess
import time

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = [
    os.path.join(SCRIPT_DIR, "scripts", "01_prepare_unaligned_data.py"),
    os.path.join(SCRIPT_DIR, "scripts", "02_run_mafft_alignments.py"),
    os.path.join(SCRIPT_DIR, "scripts", "03_run_phylogenetic_inference.py"),
    os.path.join(SCRIPT_DIR, "scripts", "04_analyze_nrr1_congruence.py"),
    os.path.join(SCRIPT_DIR, "scripts", "05_visualize_hashrf_trees.py"),
    os.path.join(SCRIPT_DIR, "scripts", "06_plot_pairwise_identities_heatmap.py"),
    os.path.join(SCRIPT_DIR, "scripts", "07_plot_wuhan_clade_support_heatmap.py"),
    os.path.join(SCRIPT_DIR, "scripts", "08_audit_nra3_recombination.py"),
]

def main():
    print("="*75)
    print(" SARS-CoV-2 ORIGINS (Boni et al. 2020) COMPLETE REPRODUCTION PIPELINE")
    print("="*75)
    total_start = time.time()

    for idx, script_path in enumerate(SCRIPTS, 1):
        script_name = os.path.basename(script_path)
        print(f"\n[STEP {idx}/{len(SCRIPTS)}] Executing: {script_name} ...")
        t0 = time.time()
        res = subprocess.run([sys.executable, script_path], text=True)
        elapsed = time.time() - t0
        if res.returncode != 0:
            print(f"\n[ERROR] Pipeline aborted: {script_name} failed with exit code {res.returncode}")
            sys.exit(res.returncode)
        print(f"[STEP {idx}/{len(SCRIPTS)} SUCCESS] {script_name} finished in {elapsed:.1f}s")

    total_time = time.time() - total_start
    print("\n" + "="*75)
    print(f" ALL REPRODUCTION PROCEDURES COMPLETED SUCCESSFULLY in {total_time/60:.2f} minutes")
    print(" Output files generated:")
    print("  - Alignments: reproduction/alignments/{68_taxa, 78_taxa}/")
    print("  - Trees:       reproduction/trees/{68_taxa, 78_taxa}/")
    print("  - Results:     reproduction/results/")
    print("="*75)

if __name__ == "__main__":
    main()
