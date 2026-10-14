#!/usr/bin/env python3
"""
Step 4: Analyze NRR1 Congruence, Topological Discordance, and Clade Relationships.

This script performs comprehensive comparative analyses on the de novo alignments
and inferred ML trees for both the 68-taxon (original) and 78-taxon (updated) datasets:
1. Alignment statistics (de novo MAFFT vs Boni et al. sliced alignments).
2. Pairwise sequence identities to SARS-CoV-2 (Wuhan-Hu-1).
3. Topological incongruence via raw and normalized Robinson-Foulds distances (full and pruned sets).
4. Rigorous statistical hypothesis testing of tree topologies via IQ-TREE (AU, SH, KH, ELW tests).
5. Detailed clade resolution and bootstrap support around SARS-CoV-2 (Wuhan-Hu-1).
"""

import os
import sys
import shutil
import subprocess
import itertools
import pandas as pd
import numpy as np
from Bio import SeqIO, Phylo
import dendropy

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
REPRO_DIR = os.path.join(BASE_DIR, "reproduction")
ALIGN_DIR = os.path.join(REPRO_DIR, "alignments")
TREE_DIR = os.path.join(REPRO_DIR, "trees")
RESULTS_DIR = os.path.join(REPRO_DIR, "results")
ORIG_DIR = os.path.join(BASE_DIR, "SARSCoV2origins")


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

os.makedirs(RESULTS_DIR, exist_ok=True)

PARTITIONS = [
    "Region_B",
    "Region_C",
    "Region_Aprime",
    "NRR1_pruned",
    "NRR1_full",
    "NRR2",
    "NRA3",
    "Spike"
]

DATASETS = ["68_taxa", "78_taxa"]

def find_wuhan_id(seq_records):
    for rec_id in seq_records.keys():
        if "wuhan" in rec_id.lower():
            return rec_id
    return None

def compute_alignment_stats():
    print("\n" + "="*70)
    print(" 1. COMPUTING ALIGNMENT METRICS & PAIRWISE IDENTITIES")
    print("="*70)

    align_records = []
    identity_records = []

    # Focal comparative taxa
    focal_targets = [
        "RaTG13",
        "Pangolin_GD",
        "Pangolin_GX",
        "RpYN06",
        "RmYN02",
        "RsYN04",
        "RmYN05",
        "RmYN08",
        "CoVZC45",
        "CoVZXC21"
    ]

    def match_target_header(target_name, headers):
        for h in headers:
            hl = h.lower()
            if target_name == "Pangolin_GD":
                # Specifically match Guangdong pangolin coronavirus (EPI_ISL_410721)
                # Avoid false-positive match to bat virus HKU3-7 from Guangdong
                if "410721" in h or ("pangolin" in hl and "guandong" in hl):
                    return h
            elif target_name == "Pangolin_GX":
                # Specifically match Guangxi pangolin isolate P1E (EPI_ISL_410539) used in Boni et al. Fig 1b/1c
                if "410539" in h or "p1e" in hl:
                    return h
            elif target_name.lower() in hl:
                return h
        return None

    for ds in DATASETS:
        ds_dir = os.path.join(ALIGN_DIR, ds)
        for part in PARTITIONS:
            aln_file = os.path.join(ds_dir, f"{part}.aln.fasta")
            if not os.path.exists(aln_file):
                continue
            
            recs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(aln_file, "fasta")}
            num_seqs = len(recs)
            aln_len = len(next(iter(recs.values())))

            # Check gap-only columns
            matrix = [list(s) for s in recs.values()]
            gap_cols = 0
            for col_idx in range(aln_len):
                col_chars = set(matrix[r][col_idx] for r in range(num_seqs))
                if col_chars == {'-'} or col_chars == {'-', 'N'}:
                    gap_cols += 1

            total_chars = num_seqs * aln_len
            total_gaps = sum(s.count('-') for s in recs.values())
            gap_pct = (total_gaps / total_chars) * 100.0 if total_chars > 0 else 0

            # Compare to original if 68_taxa
            orig_len = None
            if ds == "68_taxa":
                if part == "Region_B":
                    orig_f = os.path.join(ORIG_DIR, "recombinationAnalyses/1 First 3SEQ Iteration/1.1 Phylogenies/Region_B_3625_9150.fasta")
                elif part == "Region_C":
                    orig_f = os.path.join(ORIG_DIR, "recombinationAnalyses/1 First 3SEQ Iteration/1.1 Phylogenies/Region_C_9261_11795.fasta")
                elif part == "Region_Aprime":
                    orig_f = None
                elif part == "NRR1_pruned":
                    orig_f = os.path.join(ORIG_DIR, "alignments/sarbecovirus/NRR1/NRR1.fas")
                elif part == "NRR2":
                    orig_f = os.path.join(ORIG_DIR, "alignments/sarbecovirus/NRR2/NRR2.fas")
                elif part == "NRA3":
                    orig_f = os.path.join(ORIG_DIR, "alignments/sarbecovirus/NAR3/NAR3.fas")
                else:
                    orig_f = None
                
                if orig_f and os.path.exists(orig_f):
                    first_seq = next(SeqIO.parse(orig_f, "fasta"))
                    orig_len = len(first_seq.seq)

            align_records.append({
                "Dataset": ds,
                "Partition": part,
                "Num_Taxa": num_seqs,
                "De_Novo_Length_nt": aln_len,
                "Boni_Original_Length_nt": orig_len if orig_len is not None else "N/A",
                "Gap_Only_Columns": gap_cols,
                "Total_Gap_Percent": round(gap_pct, 2)
            })

            # Pairwise identities to Wuhan-Hu-1
            wuhan_id = find_wuhan_id(recs)
            if not wuhan_id:
                continue
            wuhan_seq = recs[wuhan_id]

            for target_name in focal_targets:
                matched_id = match_target_header(target_name, recs.keys())
                if matched_id and matched_id != wuhan_id:
                    target_seq = recs[matched_id]
                    # Calculate identity over non-gap positions
                    valid_pairs = 0
                    matches = 0
                    for c1, c2 in zip(wuhan_seq, target_seq):
                        if c1 not in ['-', 'N'] and c2 not in ['-', 'N']:
                            valid_pairs += 1
                            if c1 == c2:
                                matches += 1
                    pct_ident = (matches / valid_pairs * 100.0) if valid_pairs > 0 else 0.0
                    identity_records.append({
                        "Dataset": ds,
                        "Partition": part,
                        "Focal_Sequence": "Wuhan-Hu-1",
                        "Compared_Taxon": target_name,
                        "Matched_Header": matched_id,
                        "Valid_Aligned_Sites": valid_pairs,
                        "Identical_Sites": matches,
                        "Percent_Identity": round(pct_ident, 2)
                    })

    df_align = pd.DataFrame(align_records)
    df_ident = pd.DataFrame(identity_records)

    align_csv = os.path.join(RESULTS_DIR, "alignment_metrics_summary.csv")
    ident_csv = os.path.join(RESULTS_DIR, "pairwise_identities_wuhan.csv")
    df_align.to_csv(align_csv, index=False)
    df_ident.to_csv(ident_csv, index=False)

    print(f"[*] Alignment metrics saved to: {align_csv}")
    print(f"[*] Wuhan pairwise identities saved to: {ident_csv}")
    print("\n--- Summary of De Novo vs Original Alignment Lengths ---")
    print(df_align[["Dataset", "Partition", "Num_Taxa", "De_Novo_Length_nt", "Boni_Original_Length_nt", "Total_Gap_Percent"]].to_string(index=False))

def compute_robinson_foulds():
    print("\n" + "="*70)
    print(" 2. COMPUTING ROBINSON-FOULDS TOPOLOGICAL DISTANCES")
    print("="*70)

    rf_records = []

    for ds in DATASETS:
        tree_dir = os.path.join(TREE_DIR, ds)
        
        # Test A: Full Taxon Set
        full_parts = ["Region_B", "Region_C", "Region_Aprime", "NRR1_full", "NRR2", "NRA3", "Spike"]
            
        for p1, p2 in itertools.combinations(full_parts, 2):
            t1_file = os.path.join(tree_dir, f"{p1}.treefile")
            t2_file = os.path.join(tree_dir, f"{p2}.treefile")
            if not os.path.exists(t1_file) or not os.path.exists(t2_file):
                continue
            
            t1 = dendropy.Tree.get(path=t1_file, schema="newick")
            t2 = dendropy.Tree.get(path=t2_file, schema="newick")
            
            labels1 = set(x.label for x in t1.taxon_namespace)
            labels2 = set(x.label for x in t2.taxon_namespace)
            common_labels = labels1 & labels2
            
            if len(labels1) != len(common_labels):
                t1.retain_taxa_with_labels(common_labels)
            if len(labels2) != len(common_labels):
                t2.retain_taxa_with_labels(common_labels)
                
            tns = dendropy.TaxonNamespace(common_labels)
            t1_reloaded = dendropy.Tree.get(data=t1.as_string(schema="newick"), schema="newick", taxon_namespace=tns)
            t2_reloaded = dendropy.Tree.get(data=t2.as_string(schema="newick"), schema="newick", taxon_namespace=tns)
            
            t1_reloaded.encode_bipartitions()
            t2_reloaded.encode_bipartitions()
            
            rf = dendropy.calculate.treecompare.symmetric_difference(t1_reloaded, t2_reloaded)
            n_taxa = len(common_labels)
            max_rf = 2 * (n_taxa - 3)
            norm_rf = rf / max_rf if max_rf > 0 else 0.0

            rf_records.append({
                "Dataset": ds,
                "Comparison_Type": "Full_Taxa",
                "Partition_1": p1,
                "Partition_2": p2,
                "Num_Taxa": n_taxa,
                "Raw_RF_Distance": rf,
                "Max_Possible_RF": max_rf,
                "Normalized_RF": round(norm_rf, 4)
            })

        # Test B: Pruned Taxon Set
        t_pruned_file = os.path.join(tree_dir, "NRR1_pruned.treefile")
        if os.path.exists(t_pruned_file):
            t_ref = dendropy.Tree.get(path=t_pruned_file, schema="newick")
            pruned_labels = set(x.label for x in t_ref.taxon_namespace)
            pruned_parts = ["Region_B", "Region_C", "Region_Aprime", "NRR1_pruned", "NRR2", "NRA3"]
            for p1, p2 in itertools.combinations(pruned_parts, 2):
                t1_file = os.path.join(tree_dir, f"{p1}.treefile")
                t2_file = os.path.join(tree_dir, f"{p2}.treefile")
                if not os.path.exists(t1_file) or not os.path.exists(t2_file):
                    continue
                t1 = dendropy.Tree.get(path=t1_file, schema="newick")
                t2 = dendropy.Tree.get(path=t2_file, schema="newick")
                
                common_labels = pruned_labels & set(x.label for x in t1.taxon_namespace) & set(x.label for x in t2.taxon_namespace)
                t1.retain_taxa_with_labels(common_labels)
                t2.retain_taxa_with_labels(common_labels)
                
                tns = dendropy.TaxonNamespace(common_labels)
                t1_reloaded = dendropy.Tree.get(data=t1.as_string(schema="newick"), schema="newick", taxon_namespace=tns)
                t2_reloaded = dendropy.Tree.get(data=t2.as_string(schema="newick"), schema="newick", taxon_namespace=tns)
                
                t1_reloaded.encode_bipartitions()
                t2_reloaded.encode_bipartitions()

                rf = dendropy.calculate.treecompare.symmetric_difference(t1_reloaded, t2_reloaded)
                n_taxa = len(common_labels)
                max_rf = 2 * (n_taxa - 3)
                norm_rf = rf / max_rf if max_rf > 0 else 0.0

                rf_records.append({
                    "Dataset": ds,
                    "Comparison_Type": "Pruned_Taxa",
                    "Partition_1": p1,
                    "Partition_2": p2,
                    "Num_Taxa": n_taxa,
                    "Raw_RF_Distance": rf,
                    "Max_Possible_RF": max_rf,
                    "Normalized_RF": round(norm_rf, 4)
                })

    df_rf = pd.DataFrame(rf_records)
    rf_csv = os.path.join(RESULTS_DIR, "robinson_foulds_distances.csv")
    df_rf.to_csv(rf_csv, index=False)
    print(f"[*] Robinson-Foulds distances saved to: {rf_csv}")
    print("\n--- Key NRR1 Intra-Partition Topological Discordance ---")
    sub = df_rf[(df_rf["Partition_1"].isin(["Region_B", "Region_C"])) & 
                (df_rf["Partition_2"].isin(["Region_C", "Region_Aprime", "NRR1_pruned", "NRR1_full"]))]
    print(sub[["Dataset", "Comparison_Type", "Partition_1", "Partition_2", "Num_Taxa", "Raw_RF_Distance", "Normalized_RF"]].to_string(index=False))

def run_topology_tests():
    print("\n" + "="*70)
    print(" 3. EXECUTING STATISTICAL TOPOLOGY TESTS (AU, SH, KH, ELW)")
    print("="*70)

    test_records = []

    for ds in DATASETS:
        tree_dir = os.path.join(TREE_DIR, ds)
        align_dir = os.path.join(ALIGN_DIR, ds)
        
        # Prepare candidate trees file for Region B, C, A', NRR1_full, and NRR2
        cand_trees_file = os.path.join(tree_dir, "candidates_5trees.tre")
        tree_names = ["Region_B", "Region_C", "Region_Aprime", "NRR1_full", "NRR2"]
        
        with open(cand_trees_file, "w") as out:
            for tname in tree_names:
                tfile = os.path.join(tree_dir, f"{tname}.treefile")
                with open(tfile) as inf:
                    out.write(inf.read().strip() + "\n")
                    
        # Test each alignment against all candidate trees
        test_parts = ["Region_B", "Region_C", "Region_Aprime", "NRR1_full", "NRR2"]
        for test_part in test_parts:
            aln_file = os.path.join(align_dir, f"{test_part}.aln.fasta")
            if not os.path.exists(aln_file):
                continue
            prefix = os.path.join(tree_dir, f"au_test_5trees_{test_part}")
            iqtree_out = f"{prefix}.iqtree"
            
            if not os.path.exists(iqtree_out) or "--redo" in sys.argv:
                cmd = [
                    IQTREE_BIN,
                    "-s", aln_file,
                    "-m", "GTR+G4",
                    "--trees", cand_trees_file,
                    "--test", "10000",
                    "--test-au",
                    "--prefix", prefix,
                    "-T", "4",
                    "-redo"
                ]
                print(f"[*] Running AU test: {ds} {test_part} ...")
                res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                if res.returncode != 0:
                    print(f"  [ERROR] AU test failed for {ds} {test_part}: {res.stderr[:200]}")
                    continue
            else:
                print(f"[*] Using existing AU test report: {iqtree_out}")
                
            # Parse .iqtree report table
            if os.path.exists(iqtree_out):
                with open(iqtree_out) as f:
                    lines = f.readlines()
                in_table = False
                for line in lines:
                    if "USER TREES" in line:
                        in_table = True
                        continue
                    if in_table and line.strip().startswith("Tree"):
                        continue
                    if in_table and line.strip().startswith("------"):
                        continue
                    if in_table and line.strip().startswith("deltaL"):
                        break
                    if in_table and line.strip():
                        parts = line.split()
                        if len(parts) >= 8 and parts[0].isdigit():
                            t_idx = int(parts[0]) - 1
                            if t_idx < len(tree_names):
                                evaluated_tree_name = tree_names[t_idx]
                                logL = float(parts[1])
                                deltaL = float(parts[2])
                                bp_rell = float(parts[3])
                                p_kh = float(parts[5]) if parts[4] in ['+', '-'] else float(parts[4])
                                p_sh = float(parts[7]) if parts[6] in ['+', '-'] else float(parts[6])
                                c_elw = parts[-3] if parts[-2] in ['+', '-'] else parts[-2]
                                p_au = float(parts[-1]) if parts[-1] not in ['+', '-'] else float(parts[-2])
                                sig_excl = "-" in line

                                test_records.append({
                                    "Dataset": ds,
                                    "Test_Alignment": test_part,
                                    "Evaluated_Tree": evaluated_tree_name,
                                    "LogL": logL,
                                    "deltaL": deltaL,
                                    "p_KH": p_kh,
                                    "p_SH": p_sh,
                                    "c_ELW": c_elw,
                                    "p_AU": p_au,
                                    "Statistically_Excluded": sig_excl
                                })

        # For 68_taxa, also parse NRA3 AU test report if present
        if ds == "68_taxa":
            nra3_iqtree = os.path.join(tree_dir, "au_test_NRA3.iqtree")
            if os.path.exists(nra3_iqtree):
                nra3_tree_names = ["Region_B", "Region_C", "Region_Aprime", "NRR1_full", "NRR2", "NRA3"]
                with open(nra3_iqtree) as f:
                    lines = f.readlines()
                in_table = False
                for line in lines:
                    if "USER TREES" in line:
                        in_table = True
                        continue
                    if in_table and line.strip().startswith("Tree"):
                        continue
                    if in_table and line.strip().startswith("------"):
                        continue
                    if in_table and line.strip().startswith("deltaL"):
                        break
                    if in_table and line.strip():
                        parts = line.split()
                        if len(parts) >= 8 and parts[0].isdigit():
                            t_idx = int(parts[0]) - 1
                            if t_idx < len(nra3_tree_names):
                                evaluated_tree_name = nra3_tree_names[t_idx]
                                logL = float(parts[1])
                                deltaL = float(parts[2])
                                bp_rell = float(parts[3])
                                p_kh = float(parts[5]) if parts[4] in ['+', '-'] else float(parts[4])
                                p_sh = float(parts[7]) if parts[6] in ['+', '-'] else float(parts[6])
                                c_elw = parts[-3] if parts[-2] in ['+', '-'] else parts[-2]
                                p_au = float(parts[-1]) if parts[-1] not in ['+', '-'] else float(parts[-2])
                                sig_excl = "-" in line

                                test_records.append({
                                    "Dataset": ds,
                                    "Test_Alignment": "NRA3",
                                    "Evaluated_Tree": evaluated_tree_name,
                                    "LogL": logL,
                                    "deltaL": deltaL,
                                    "p_KH": p_kh,
                                    "p_SH": p_sh,
                                    "c_ELW": c_elw,
                                    "p_AU": p_au,
                                    "Statistically_Excluded": sig_excl
                                })

    df_tests = pd.DataFrame(test_records)
    test_csv = os.path.join(RESULTS_DIR, "topology_tests_summary.csv")
    df_tests.to_csv(test_csv, index=False)
    print(f"[*] Statistical topology test results saved to: {test_csv}")
    print("\n--- Summary of Hypothesis Testing (Approximately Unbiased p-AU) ---")
    print(df_tests[["Dataset", "Test_Alignment", "Evaluated_Tree", "deltaL", "p_SH", "p_AU", "Statistically_Excluded"]].to_string(index=False))

def analyze_wuhan_clade_resolution():
    print("\n" + "="*70)
    print(" 4. RESOLVING SARS-CoV-2 SISTER TAXA & CLADE SUPPORT")
    print("="*70)

    clade_records = []

    for ds in DATASETS:
        tree_dir = os.path.join(TREE_DIR, ds)
        for part in PARTITIONS:
            tree_file = os.path.join(tree_dir, f"{part}.treefile")
            if not os.path.exists(tree_file):
                continue
            
            tree = Phylo.read(tree_file, "newick")
            
            wuhan = None
            for leaf in tree.get_terminals():
                if "wuhan" in leaf.name.lower():
                    wuhan = leaf
                    break
            if not wuhan:
                continue

            # Trace lineage to root
            # Trace lineage to root
            path = tree.get_path(wuhan)
            if len(path) < 2:
                continue

            # In Bio.Phylo, path[-1] is wuhan itself; immediate parent node is path[-2]
            parent = path[-2]
            wuhan_tax = wuhan.name.split('|')[0]
            other_terms = sorted([t.name.split('|')[0] for t in parent.get_terminals() if t != wuhan])
            parent_terms = [wuhan_tax] + other_terms
            sister_terms = other_terms
            parent_support = int(parent.confidence) if parent.confidence is not None else None

            # Key comparative lineages to resolve against Wuhan-Hu-1
            taxa_targets = {
                "RaTG13": "ratg13",
                "RpYN06": "rpyn06",
                "RmYN02": "rmyn02",
                "PrC31": "prc31",
                "Pangolin_GD": "410721",
                "Pangolin_GX": "410539",
                "CoVZC45": "zc45",
                "CoVZXC21": "zxc21"
            }

            record = {
                "Dataset": ds,
                "Partition": part,
                "Immediate_Sister_Clade_Taxa": ", ".join(parent_terms),
                "Immediate_Sister_Taxa": ", ".join(sister_terms),
                "Clade_UFBoot_Support": parent_support,
            }

            for tname, pattern in taxa_targets.items():
                matches = [l for l in tree.get_terminals() if pattern in l.name.lower() or (tname == "Pangolin_GD" and "guandong" in l.name.lower() and "pangolin" in l.name.lower())]
                if matches:
                    target = matches[0]
                    mrca = tree.common_ancestor(wuhan, target)
                    mrca_terms = [t.name.split('|')[0] for t in mrca.get_terminals()]
                    mrca_supp = int(mrca.confidence) if mrca.confidence is not None else 100
                    is_sister = (set(mrca_terms) == {wuhan_tax, target.name.split('|')[0]})
                    record[f"MRCA_{tname}_Support"] = mrca_supp
                    record[f"MRCA_{tname}_Taxa_Count"] = len(mrca_terms)
                    record[f"Is_{tname}_Sister"] = is_sister
                else:
                    record[f"MRCA_{tname}_Support"] = None
                    record[f"MRCA_{tname}_Taxa_Count"] = 0
                    record[f"Is_{tname}_Sister"] = False

            # Backwards compatibility fields
            record["Is_RaTG13_Sister"] = record["Is_RaTG13_Sister"]
            record["MRCA_Wuhan_RaTG13_Support"] = record["MRCA_RaTG13_Support"]
            record["MRCA_Wuhan_RaTG13_Taxa_Count"] = record["MRCA_RaTG13_Taxa_Count"]
            record["Is_RpYN06_Sister"] = record["Is_RpYN06_Sister"]
            record["MRCA_Wuhan_RpYN06_Support"] = record["MRCA_RpYN06_Support"]
            record["MRCA_Wuhan_RpYN06_Taxa_Count"] = record["MRCA_RpYN06_Taxa_Count"]

            clade_records.append(record)

    df_clade = pd.DataFrame(clade_records)
    clade_csv = os.path.join(RESULTS_DIR, "wuhan_clade_resolution.csv")
    df_clade.to_csv(clade_csv, index=False)
    print(f"[*] Clade resolution saved to: {clade_csv}")
    print("\n--- Summary of SARS-CoV-2 Sister Relationships ---")
    print(df_clade[["Dataset", "Partition", "Is_RaTG13_Sister", "Is_RpYN06_Sister", "MRCA_Wuhan_RaTG13_Taxa_Count", "Immediate_Sister_Clade_Taxa"]].to_string(index=False))

def generate_text_report():
    report_file = os.path.join(RESULTS_DIR, "reproduction_summary_report.md")
    df_align = pd.read_csv(os.path.join(RESULTS_DIR, "alignment_metrics_summary.csv"))
    df_rf = pd.read_csv(os.path.join(RESULTS_DIR, "robinson_foulds_distances.csv"))
    df_tests = pd.read_csv(os.path.join(RESULTS_DIR, "topology_tests_summary.csv"))
    df_clade = pd.read_csv(os.path.join(RESULTS_DIR, "wuhan_clade_resolution.csv"))
    df_ident = pd.read_csv(os.path.join(RESULTS_DIR, "pairwise_identities_wuhan.csv"))

    with open(report_file, "w") as out:
        out.write("# Empirical Reproduction & Topological Congruence Report\n\n")
        out.write("## 1. Alignment Metrics & Gap Elimination\n\n")
        out.write(df_align.to_markdown(index=False) + "\n\n")
        
        out.write("## 2. Wuhan-Hu-1 Pairwise Identity (%) Across Partitions\n\n")
        out.write(df_ident.to_markdown(index=False) + "\n\n")
        
        out.write("## 3. Robinson-Foulds Incongruence Distances\n\n")
        out.write(df_rf.to_markdown(index=False) + "\n\n")
        
        out.write("## 4. Statistical Topology Hypothesis Testing (AU / SH / KH Tests)\n\n")
        out.write(df_tests.to_markdown(index=False) + "\n\n")
        
        out.write("## 5. Clade Resolution Around SARS-CoV-2 (Wuhan-Hu-1)\n\n")
        out.write(df_clade.to_markdown(index=False) + "\n\n")

    print(f"[*] Comprehensive reproduction summary report saved to: {report_file}")

if __name__ == "__main__":
    compute_alignment_stats()
    compute_robinson_foulds()
    run_topology_tests()
    analyze_wuhan_clade_resolution()
    generate_text_report()
    print("\n[COMPLETE] Step 4 analyses completed successfully.")
