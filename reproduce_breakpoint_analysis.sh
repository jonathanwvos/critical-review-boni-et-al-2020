#!/usr/bin/env bash
# ==============================================================================
# Script: reproduce_breakpoint_analysis.sh
# Purpose: Reproduce and audit the breakpoint and phylogenetic analyses from
#          Boni et al. (2020) using the SARSCoV2origins repository data.
# Note:    All Python snippets use Python 3 syntax to avoid Python 2 errors.
# ==============================================================================

set -euo pipefail

# Ensure we are operating in the repository/assignment root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "======================================================================"
echo "SARS-CoV-2 Origins: Breakpoint & Phylogenetic Incongruence Reproduction"
echo "Working directory: $(pwd)"
echo "======================================================================"

# ------------------------------------------------------------------------------
# STEP 1: First 3SEQ Pass - Full Genome (30,927 nt) Breakpoint & BFR Extraction
# ------------------------------------------------------------------------------
echo -e "\n[STEP 1] Extracting Breakpoints & Top 10 BFRs (Regions A-J) from First 3SEQ Pass..."

python3 -c "
import operator

rec_file = 'SARSCoV2origins/recombinationAnalyses/1 First 3SEQ Iteration/att01.3s.rec.csv'
with open(rec_file) as f:
    next(f)  # skip header
    bps = []
    for line in f:
        parts = line.strip().split(',')
        for p in parts[12:]:
            pairs = p.split('&')
            if len(pairs) == 2:
                r1, r2 = pairs[0].split('-'), pairs[1].split('-')
                bps.extend([int((int(r1[0]) + int(r1[1])) / 2.0), int((int(r2[0]) + int(r2[1])) / 2.0)])

sorted_bps = sorted(list(set(bps)))
alignment_len = 30927
LL = [[1, sorted_bps[0], sorted_bps[0] - 1]]
for i in range(len(sorted_bps) - 1):
    LL.append([sorted_bps[i], sorted_bps[i+1], sorted_bps[i+1] - sorted_bps[i]])
LL.append([sorted_bps[-1], alignment_len, alignment_len - sorted_bps[-1]])

print(f'Total Unique Breakpoints Detected: {len(sorted_bps)}')
print(f'Total Candidate BFR Intervals: {len(LL)}')
print('\nTop 10 Longest Breakpoint-Free Regions (Regions A-J):')
print('----------------------------------------------------------------------')
print(f'{\"ID\":<5} {\"Coordinates (nt)\":<25} {\"Script Len\":<12} {\"Inclusive Len\":<15}')
print('----------------------------------------------------------------------')
for i, l in enumerate(sorted(LL, key=operator.itemgetter(2), reverse=True)[:10]):
    inclusive_len = l[1] - l[0] + 1
    print(f'Reg {chr(65+i):<2} nt {l[0]:>5} - {l[1]:<6}       {l[2]:<12} {inclusive_len:<15}')
print('----------------------------------------------------------------------')
"

# ------------------------------------------------------------------------------
# STEP 2: Second 3SEQ Pass - Proving Region A is Riddled with Recombination
# ------------------------------------------------------------------------------
echo -e "\n[STEP 2] Auditing Second 3SEQ Pass on Region A (ORF1b, 6338 nt)..."

python3 -c "
import operator

rec_file = 'SARSCoV2origins/recombinationAnalyses/2 Second 3SEQ Iteration/RegA.3s.rec.csv'
with open(rec_file) as f:
    next(f)
    bps = []
    for line in f:
        parts = line.strip().split(',')
        for p in parts[12:]:
            pairs = p.split('&')
            if len(pairs) == 2:
                r1, r2 = pairs[0].split('-'), pairs[1].split('-')
                bps.extend([int((int(r1[0]) + int(r1[1])) / 2.0), int((int(r2[0]) + int(r2[1])) / 2.0)])

sorted_bps = sorted(list(set(bps)))
regA_len = 6338
LL = [[1, sorted_bps[0], sorted_bps[0] - 1]]
for i in range(len(sorted_bps) - 1):
    LL.append([sorted_bps[i], sorted_bps[i+1], sorted_bps[i+1] - sorted_bps[i]])
LL.append([sorted_bps[-1], regA_len, regA_len - sorted_bps[-1]])

print(f'Internal Breakpoints Detected INSIDE Region A: {len(sorted_bps)}')
print(f'Internal Sub-BFRs Created: {len(LL)}')
print('\nTop 5 Sub-regions in Region A:')
for i, l in enumerate(sorted(LL, key=operator.itemgetter(2), reverse=True)[:5]):
    print(f'  {i+1}. Relative nt {l[0]} - {l[1]} (length: {l[2]} nt)')
print('\n* Note: Subregions A1 (3872-4716) and A4 (1642-2113) caused phylogenetic conflict')
print('  and were excised, creating Region A-prime (A\').')
"

# ------------------------------------------------------------------------------
# STEP 3: GARD vs. 3SEQ Concordance Filter (NRR2 Discarded Breakpoints)
# ------------------------------------------------------------------------------
echo -e "\n[STEP 3] Evaluating GARD vs. 3SEQ Discordance Filter (Tolerance = 50 nt)..."

python3 -c "
rec_file = 'SARSCoV2origins/recombinationAnalyses/3 Clonal Regions to be used in main analysis/NRR2 - GARD and 3SEQ/att01.3s.rec.csv'
with open(rec_file) as infile:
    lst_bp_to_check = [1684, 3046, 4448, 9237, 11894, 11885, 14415, 15590, 16466, 20015, 21753, 22773, 22778, 24598, 24628, 29238]
    lst_bp_confirmed = []
    tol = 50
    next(infile)
    for line in infile:
        l = line.split(',')
        for item in l[12:]:
            pairs = item.split('&')
            if len(pairs) == 2:
                r1 = [int(x) for x in pairs[0].split('-')]
                r2 = [int(x) for x in pairs[1].split('-')]
                for bp in lst_bp_to_check:
                    if (r1[0] - tol <= bp <= r1[1] + tol) or (r2[0] - tol <= bp <= r2[1] + tol):
                        if bp not in lst_bp_confirmed:
                            lst_bp_confirmed.append(bp)

confirmed = sorted(lst_bp_confirmed)
discarded = sorted(list(set(lst_bp_to_check) - set(lst_bp_confirmed)))

print(f'GARD Candidate Breakpoints: {sorted(list(set(lst_bp_to_check)))}')
print(f'\nConfirmed Breakpoints (within 50 nt of 3SEQ): {confirmed}')
print(f'Discarded Breakpoints (>50 nt from 3SEQ): {discarded}')
print('\n* Methodological Flaw: GARD breakpoints 14415, 15590, 16466, and 20015 inside ORF1b')
print('  were eliminated because 3SEQ called them at 13291 and 19628 (>50 nt away).')
print('  This falsely created a 9.9 kb \"non-recombinant\" region (NRR2: 11886-21753).')
"

# ------------------------------------------------------------------------------
# STEP 4: Breakpoint Support Density Across Genome (Replicating Histogram)
# ------------------------------------------------------------------------------
echo -e "\n[STEP 4] Scanning Genomic Breakpoint Support Density Across 31 kb Alignment..."

python3 -c "
import numpy as np

SEQLEN = 30927
X = np.zeros(SEQLEN + 1)
rec_file = 'SARSCoV2origins/recombinationAnalyses/1 First 3SEQ Iteration/att01.3s.rec.csv'

with open(rec_file) as f:
    next(f)
    for line in f:
        l = line.strip().split(',')
        for p in l[12:]:
            pairs = p.split('&')
            if len(pairs) == 2:
                r1 = [int(x) for x in pairs[0].split('-')]
                r2 = [int(x) for x in pairs[1].split('-')]
                X[r1[0]:r1[1]+1] += 1
                X[r2[0]:r2[1]+1] += 1

print('Top 5 Peak Recombination Breakpoint Hotspots:')
peak_indices = np.argsort(X)[::-1][:5]
for idx in peak_indices:
    print(f'  Genomic position {idx:>5}: supported by {int(X[idx]):>4} recombinant triplets')
"

# ------------------------------------------------------------------------------
# STEP 5: Verifying Molecular Dating Matrices (NRR1, NRR2, NRA3)
# ------------------------------------------------------------------------------
echo -e "\n[STEP 5] Checking Dimensions of Final Non-Recombinant Alignments..."

python3 -c "
import glob

for p in sorted(glob.glob('SARSCoV2origins/alignments/sarbecovirus/**/*.fas', recursive=True)):
    with open(p) as f:
        seqs = {}
        cur = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                cur = line[1:]
                seqs[cur] = []
            else:
                seqs[cur].append(line)
        lens = {len(''.join(v)) for v in seqs.values()}
        name = p.split('/')[-2]
        print(f'  Dataset {name:<5}: {len(seqs):>2} taxa | Alignment Length: {list(lens)[0]:>5} nt')
"

# ------------------------------------------------------------------------------
# STEP 6: Reviewing RAxML Inference Parameters & Warning Outputs
# ------------------------------------------------------------------------------
echo -e "\n[STEP 6] Inspecting RAxML Parameters & Taxon Redundancies for Region A..."

head -n 25 "SARSCoV2origins/recombinationAnalyses/1 First 3SEQ Iteration/1.1 Phylogenies/RAxML_info.RegA" | grep -E "WARNING|distinct" || true
grep -A 2 "RAxML was called as follows:" "SARSCoV2origins/recombinationAnalyses/1 First 3SEQ Iteration/1.1 Phylogenies/RAxML_info.RegA"

echo -e "\n======================================================================"
echo "Reproduction & verification run complete."
echo "======================================================================"
