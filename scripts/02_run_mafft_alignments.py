#!/usr/bin/env python3
"""
Script 02: De Novo Multiple Sequence Alignment with MAFFT
Executes MAFFT across all unaligned regional partitions for both:
  - 68_taxa (Original Boni et al. dataset)
  - 78_taxa (Updated dataset with RpYN06, RsYN04, RmYN02, RmYN05, RmYN08)
"""

import os
import sys
import time
import shutil
import subprocess

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, "data")
ALN_DIR = os.path.join(BASE_DIR, "alignments")


def find_mafft_binary():
    env_bin = os.environ.get("MAFFT_BIN")
    if env_bin and (os.path.exists(env_bin) or shutil.which(env_bin)):
        return env_bin
    which_bin = shutil.which("mafft")
    if which_bin:
        return which_bin
    sys_bin = os.path.join(sys.prefix, "bin", "mafft")
    if os.path.exists(sys_bin):
        return sys_bin
    for candidate in [
        os.path.expanduser("~/.conda/envs/critrev/bin/mafft"),
        "/home/jono/.conda/envs/critrev/bin/mafft",
        "/home/jono/miniconda3/envs/phylo/bin/mafft",
    ]:
        if os.path.exists(candidate):
            return candidate
    return "mafft"


MAFFT_BIN = find_mafft_binary()

TARGET_REGIONS = [
    ("Region_B_unaligned.fasta", "Region_B.aln.fasta"),
    ("Region_C_unaligned.fasta", "Region_C.aln.fasta"),
    ("Region_Aprime_unaligned.fasta", "Region_Aprime.aln.fasta"),
    ("NRR1_pruned_unaligned.fasta", "NRR1_pruned.aln.fasta"),
    ("NRR1_full_unaligned.fasta", "NRR1_full.aln.fasta"),
    ("NRR2_unaligned.fasta", "NRR2.aln.fasta"),
    ("Spike_unaligned.fasta", "Spike.aln.fasta")
]

def check_alignment_stats(aln_file):
    seqs = {}
    with open(aln_file, "r") as f:
        cur = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                cur = line[1:]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    lengths = {len("".join(v)) for v in seqs.values()}
    return len(seqs), list(lengths)[0] if lengths else 0

def run_mafft_for_dataset(dataset_name):
    print(f"\n=======================================================")
    print(f" Starting De Novo MAFFT Alignments for: {dataset_name}")
    print(f"=======================================================")
    
    in_dir = os.path.join(DATA_DIR, dataset_name)
    out_dir = os.path.join(ALN_DIR, dataset_name)
    os.makedirs(out_dir, exist_ok=True)

    results_summary = []

    for in_fname, out_fname in TARGET_REGIONS:
        in_path = os.path.join(in_dir, in_fname)
        out_path = os.path.join(out_dir, out_fname)

        if not os.path.exists(in_path):
            print(f"  [SKIP] {in_path} does not exist.")
            continue

        if os.path.exists(out_path) and os.path.getsize(out_path) > 100 and "--redo" not in sys.argv:
            taxa_count, aln_len = check_alignment_stats(out_path)
            print(f"  [EXISTS] {out_fname} already aligned | Taxa: {taxa_count} | Length: {aln_len} nt")
            results_summary.append((out_fname, taxa_count, aln_len, 0.0))
            continue

        print(f"\n[*] Aligning {in_fname} -> {out_fname} ...")
        t0 = time.time()
        
        # Run MAFFT with 8 threads, auto algorithm
        cmd = [MAFFT_BIN, "--auto", "--thread", "8", in_path]
        with open(out_path, "w", encoding="utf-8") as out_f:
            res = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True)

        t1 = time.time()
        elapsed = t1 - t0

        if res.returncode != 0:
            print(f"  [ERROR] MAFFT failed for {in_fname}:\n{res.stderr}")
        else:
            taxa_count, aln_len = check_alignment_stats(out_path)
            print(f"  [SUCCESS] Completed in {elapsed:.2f}s | Taxa: {taxa_count} | Aligned Length: {aln_len} nt")
            results_summary.append((out_fname, taxa_count, aln_len, elapsed))

    return results_summary

def main():
    if not (os.path.exists(MAFFT_BIN) or shutil.which(MAFFT_BIN)):
        print(f"Error: MAFFT binary not found (looked for '{MAFFT_BIN}', in PATH, and in conda env)")
        sys.exit(1)

    res_68 = run_mafft_for_dataset("68_taxa")
    res_78 = run_mafft_for_dataset("78_taxa")

    print("\n" + "=" * 80)
    print(" SUMMARY OF ALL DE NOVO MAFFT ALIGNMENTS")
    print("=" * 80)
    for name, res in [("68_taxa (Original)", res_68), ("78_taxa (Updated)", res_78)]:
        print(f"\n--- Dataset: {name} ---")
        print(f"{'Alignment Target':<28} {'Taxa':<8} {'Length (nt)':<14} {'Time (s)':<10}")
        print("-" * 65)
        for r in res:
            print(f"{r[0]:<28} {r[1]:<8} {r[2]:<14} {r[3]:<10.2f}")

    print("\n[COMPLETE] All de novo multiple sequence alignments finished.")

if __name__ == "__main__":
    main()
