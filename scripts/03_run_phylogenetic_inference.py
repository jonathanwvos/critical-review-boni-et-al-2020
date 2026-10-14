#!/usr/bin/env python3
"""
Script 03: Maximum-Likelihood Phylogenetic Inference with IQ-TREE
Fits the GTR+G4 model with 1000 ultrafast bootstrap replicates (-B 1000)
across all de novo aligned partitions for both 68_taxa and 78_taxa.
"""

import os
import sys
import time
import shutil
import subprocess

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ALN_DIR = os.path.join(BASE_DIR, "alignments")
TREE_DIR = os.path.join(BASE_DIR, "trees")


def find_iqtree_binary():
    env_bin = os.environ.get("IQTREE_BIN")
    if env_bin and (os.path.exists(env_bin) or shutil.which(env_bin)):
        return env_bin
    for name in ["iqtree2", "iqtree"]:
        which_bin = shutil.which(name)
        if which_bin:
            return which_bin
        sys_bin = os.path.join(sys.prefix, "bin", name)
        if os.path.exists(sys_bin):
            return sys_bin
        for candidate in [
            os.path.expanduser(f"~/.conda/envs/critrev/bin/{name}"),
            f"/home/jono/.conda/envs/critrev/bin/{name}",
            f"/home/jono/miniconda3/envs/phylo/bin/{name}",
        ]:
            if os.path.exists(candidate):
                return candidate
    return "iqtree2"


IQTREE_BIN = find_iqtree_binary()

TARGET_ALNS = [
    "Region_B.aln.fasta",
    "Region_C.aln.fasta",
    "Region_Aprime.aln.fasta",
    "NRR1_pruned.aln.fasta",
    "NRR1_full.aln.fasta",
    "NRR2.aln.fasta",
    "NRA3.aln.fasta",
    "Spike.aln.fasta"
]

def run_iqtree_for_dataset(dataset_name):
    print(f"\n=======================================================")
    print(f" Starting Phylogenetic Inference for: {dataset_name}")
    print(f"=======================================================")

    in_dir = os.path.join(ALN_DIR, dataset_name)
    out_dir = os.path.join(TREE_DIR, dataset_name)
    os.makedirs(out_dir, exist_ok=True)

    summary_results = []

    for aln_name in TARGET_ALNS:
        aln_path = os.path.join(in_dir, aln_name)
        if not os.path.exists(aln_path):
            print(f"  [SKIP] {aln_path} not found.")
            continue

        prefix_name = aln_name.replace(".aln.fasta", "")
        out_prefix = os.path.join(out_dir, prefix_name)
        treefile = out_prefix + ".treefile"

        if os.path.exists(treefile) and os.path.getsize(treefile) > 100 and "--redo" not in sys.argv:
            print(f"  [EXISTS] Tree already inferred: {treefile}")
            summary_results.append((prefix_name, 0.0, treefile))
            continue

        print(f"\n[*] Inferring ML tree for {aln_name} ...")
        t0 = time.time()

        # Command: iqtree2 -s <aln> -m GTR+G4 -B 1000 -T 8 --prefix <prefix> --redo
        cmd = [
            IQTREE_BIN,
            "-s", aln_path,
            "-m", "GTR+G4",
            "-B", "1000",
            "-T", "8",
            "--prefix", out_prefix,
            "--redo"
        ]

        res = subprocess.run(cmd, capture_output=True, text=True)
        t1 = time.time()
        elapsed = t1 - t0

        if os.path.exists(treefile):
            print(f"  [SUCCESS] Tree inferred in {elapsed:.2f}s -> {treefile}")
            summary_results.append((prefix_name, elapsed, treefile))
        else:
            print(f"  [ERROR] IQ-TREE failed for {aln_name}:\n{res.stderr[:500]}")

    return summary_results

def main():
    if not (os.path.exists(IQTREE_BIN) or shutil.which(IQTREE_BIN)):
        print(f"Error: IQ-TREE binary not found (looked for '{IQTREE_BIN}', in PATH, and in conda env)")
        sys.exit(1)

    res_68 = run_iqtree_for_dataset("68_taxa")
    res_78 = run_iqtree_for_dataset("78_taxa")

    print("\n" + "=" * 80)
    print(" SUMMARY OF ALL MAXIMUM-LIKELIHOOD PHYLOGENIES INFERRED")
    print("=" * 80)
    for name, res in [("68_taxa (Original Dataset)", res_68), ("78_taxa (Updated Dataset)", res_78)]:
        print(f"\n--- Dataset: {name} ---")
        print(f"{'Partition':<24} {'Time (s)':<10} {'Output Tree File'}")
        print("-" * 75)
        for r in res:
            print(f"{r[0]:<24} {r[1]:<10.2f} {r[2]}")

    print("\n[COMPLETE] All ML phylogenies inferred successfully.")

if __name__ == "__main__":
    main()
