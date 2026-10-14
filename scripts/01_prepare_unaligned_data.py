#!/usr/bin/env python3
"""
Script 01: Extract Raw Unaligned Sequences for Regional Alignments
Datasets:
  1. Original 68 sarbecovirus genomes (Boni et al.) from nCoVbat.pan6.fst
  2. Updated 78 sarbecovirus genomes from mar21update_78vir_wgal_fn.fas
Regions Extracted:
  - Region B (ORF1a: nsp3-nsp4)
  - Region C (ORF1a: nsp4-nsp6)
  - Region A (ORF1b: nsp12-nsp15)
  - Region A' (three spliced non-recombinant fragments: A'_1, A'_2, A'_3)
  - NRR1 Full (A' + B + C across all taxa)
  - NRR1 Pruned (A' + B + C with 5 recombinant taxa removed)
  - NRR2 (ORF1b/GARD longest candidate segment)
  - Spike Gene (Full Spike glycoprotein)
  - Spike RBD (Receptor-binding domain / variable loop)
"""

import os
import sys

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT_ROOT = os.path.dirname(BASE_DIR)
DATA_DIR = os.path.join(BASE_DIR, "data")
os.makedirs(DATA_DIR, exist_ok=True)

# 5 taxa removed by Boni et al. to enforce NRR1 clonality
PRUNED_TAXA = {
    "DQ412042",
    "CoVZXC21",
    "CoVZC45",
    "NC_014470",
    "BM48-31",
    "SC2018",
    "MK211374"
}

def is_pruned(header):
    for p in PRUNED_TAXA:
        if p in header:
            return True
    return False

def read_fasta_alignment(path):
    seqs = {}
    cur = None
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                cur = line[1:]
                seqs[cur] = []
            elif cur:
                seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}

def slice_and_strip(aln_seq, start_1based, end_1based):
    sub = aln_seq[start_1based - 1 : end_1based]
    return sub.replace("-", "").replace("?", "").upper()

def write_fasta(filepath, seq_dict):
    with open(filepath, "w", encoding="utf-8") as f:
        for h, s in seq_dict.items():
            f.write(f">{h}\n{s}\n")
    print(f"  [SAVED] {filepath} ({len(seq_dict)} sequences)")

def process_dataset(name, fasta_path, coord_dict):
    print(f"\n=======================================================")
    print(f" Processing Dataset: {name}")
    print(f" Source: {fasta_path}")
    print(f"=======================================================")
    
    aln = read_fasta_alignment(fasta_path)
    print(f" Total genomes loaded: {len(aln)}")

    out_prefix = os.path.join(DATA_DIR, name)
    os.makedirs(out_prefix, exist_ok=True)

    # 1. Whole genome unaligned
    wg_unaligned = {k: v.replace("-", "").replace("?", "").upper() for k, v in aln.items()}
    write_fasta(os.path.join(out_prefix, "whole_genome_unaligned.fasta"), wg_unaligned)

    # 2. Region B
    s_b, e_b = coord_dict["Region_B"]
    b_unaligned = {k: slice_and_strip(v, s_b, e_b) for k, v in aln.items()}
    write_fasta(os.path.join(out_prefix, "Region_B_unaligned.fasta"), b_unaligned)

    # 3. Region C
    s_c, e_c = coord_dict["Region_C"]
    c_unaligned = {k: slice_and_strip(v, s_c, e_c) for k, v in aln.items()}
    write_fasta(os.path.join(out_prefix, "Region_C_unaligned.fasta"), c_unaligned)

    # 4. Region A (full)
    s_a, e_a = coord_dict["Region_A"]
    a_unaligned = {k: slice_and_strip(v, s_a, e_a) for k, v in aln.items()}
    write_fasta(os.path.join(out_prefix, "Region_A_unaligned.fasta"), a_unaligned)

    # 5. Region A' (three spliced fragments)
    s_a1, e_a1 = coord_dict["A_prime_1"]
    s_a2, e_a2 = coord_dict["A_prime_2"]
    s_a3, e_a3 = coord_dict["A_prime_3"]
    aprime_unaligned = {}
    for k, v in aln.items():
        f1 = slice_and_strip(v, s_a1, e_a1)
        f2 = slice_and_strip(v, s_a2, e_a2)
        f3 = slice_and_strip(v, s_a3, e_a3)
        aprime_unaligned[k] = f1 + f2 + f3
    write_fasta(os.path.join(out_prefix, "Region_Aprime_unaligned.fasta"), aprime_unaligned)

    # 6. NRR1 Full (A' + B + C)
    nrr1_full = {}
    for k in aln.keys():
        nrr1_full[k] = aprime_unaligned[k] + b_unaligned[k] + c_unaligned[k]
    write_fasta(os.path.join(out_prefix, "NRR1_full_unaligned.fasta"), nrr1_full)

    # 7. NRR1 Pruned (5 taxa removed)
    nrr1_pruned = {k: v for k, v in nrr1_full.items() if not is_pruned(k)}
    write_fasta(os.path.join(out_prefix, "NRR1_pruned_unaligned.fasta"), nrr1_pruned)

    # 8. NRR2 (ORF1b/GARD longest segment)
    s_n2, e_n2 = coord_dict["NRR2"]
    nrr2_unaligned = {k: slice_and_strip(v, s_n2, e_n2) for k, v in aln.items()}
    write_fasta(os.path.join(out_prefix, "NRR2_unaligned.fasta"), nrr2_unaligned)

    # 9. Spike Gene
    s_sp, e_sp = coord_dict["Spike"]
    spike_unaligned = {k: slice_and_strip(v, s_sp, e_sp) for k, v in aln.items()}
    write_fasta(os.path.join(out_prefix, "Spike_unaligned.fasta"), spike_unaligned)

def main():
    # Coordinates in 68-genome alignment (nCoVbat.pan6.fst)
    coords_68 = {
        "Region_B": (3625, 9150),
        "Region_C": (9261, 11795),
        "Region_A": (13291, 19628),
        "A_prime_1": (13291, 14932),
        "A_prime_2": (15405, 17162),
        "A_prime_3": (18009, 19628),
        "NRR2": (11886, 21753),
        "Spike": (21563, 25384)
    }

    # Corresponding coordinates in 78-genome alignment (mar21update_78vir_wgal_fn.fas)
    coords_78 = {
        "Region_B": (3648, 9205),
        "Region_C": (9316, 11848),
        "Region_A": (13343, 19671),
        "A_prime_1": (13343, 14984),
        "A_prime_2": (15457, 17212),
        "A_prime_3": (18053, 19671),
        "NRR2": (11939, 21794),
        "Spike": (21604, 25307)
    }

    path_68 = os.path.join(PROJECT_ROOT, "SARSCoV2origins", "recombinationAnalyses", "nCoVbat.pan6.fst")
    path_78 = os.path.join(PROJECT_ROOT, "mar21update_78vir_wgal_fn.fas")

    process_dataset("68_taxa", path_68, coords_68)
    process_dataset("78_taxa", path_78, coords_78)

    # Ensure NRA3 alignments are in place for both cohorts
    aln_dir_68 = os.path.join(BASE_DIR, "alignments", "68_taxa")
    aln_dir_78 = os.path.join(BASE_DIR, "alignments", "78_taxa")
    os.makedirs(aln_dir_68, exist_ok=True)
    os.makedirs(aln_dir_78, exist_ok=True)

    nra3_68 = os.path.join(aln_dir_68, "NRA3.aln.fasta")
    orig_nra3 = os.path.join(PROJECT_ROOT, "SARSCoV2origins", "recombinationAnalyses", "3 Clonal Regions to be used in main analysis", "NRA3", "NAR3.fas")
    if os.path.exists(orig_nra3) and (not os.path.exists(nra3_68) or os.path.getsize(nra3_68) < 100):
        import shutil
        shutil.copyfile(orig_nra3, nra3_68)
        print(f"  [SAVED] Copied 68-taxon NRA3 alignment: {nra3_68}")
    else:
        print(f"  [EXISTS] 68-taxon NRA3 alignment: {nra3_68}")

    nra3_78 = os.path.join(aln_dir_78, "NRA3.aln.fasta")
    if not os.path.exists(nra3_78) or os.path.getsize(nra3_78) < 100:
        from importlib import import_module
        mod = import_module("08_audit_nra3_recombination")
        mod.construct_78_taxon_nra3()
    else:
        print(f"  [EXISTS] 78-taxon NRA3 alignment: {nra3_78}")

    print("\n[COMPLETE] All unaligned regional datasets and NRA3 alignments prepared.")

if __name__ == "__main__":
    main()
