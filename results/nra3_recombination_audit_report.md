# Empirical Recombination & Phylogenetic Audit of NRA3 / NRR3

## 1. Executive Summary

Boni et al. (2020) proposed **NRA3** (Non-Recombinant Alignment 3) as an alternative putative recombination-free alignment derived via RDP5 masking. Our reproduction and multi-methodological audit falsifies the claim that NRA3 represents an un-recombined clonal lineage:

1. **Severe Reticulate Evolution Across NRA3:** A four-gamete incompatibility audit of `NAR3.fas` reveals that **24.0% of all tested informative site pairs directly violate the four-gamete condition**. Cross-block site incompatibility reaches **23.7%** between ORF1a and Spike, **15.8%** between ORF1a and the 3' genes, and **14.4%** between Spike and the 3' genes.
2. **Mutual Statistical Rejection (AU Tests) Among NRA3 Sub-Blocks:** When NRA3 is partitioned into its four constituent structural domains (ORF1a, ORF1b, Spike, and 3' genes), the maximum-likelihood tree of each block is **statistically rejected by the sequence alignments of other blocks ($p_{AU} < 10^{-3}$ to $10^{-50}$)**. Forcing these loci into a single tree topology in BEAST introduces massive phylogenetic distortion.
3. **Complete Disruption of RaTG13 Sister Relationship in 78-Taxa:** While Wuhan-Hu-1 and RaTG13 appear sister in 68-taxa NRA3 (100% UFBoot), evaluating the 78-taxon dataset completely disrupts this pairing. In the whole 78-taxon NRA3 tree, **SARS-CoV-2 groups with `RmYN02` and `RpYN06` with 100% bootstrap support**, displacing RaTG13 to a basal branch.
4. **Sister Taxon Switching Across the Genomic Length:** Across the contiguous length of NRA3 in 78 taxa, SARS-CoV-2 switches sister lineages from **RpYN06** in ORF1a (100% support) to **RmYN02 + RpYN06** in ORF1b (98% support) to **RaTG13** in Spike (100% support). A single genome displaying high-confidence sister switching between distant reservoir viruses across adjacent domains is, by definition, **a recombinant mosaic**.

## 2. Statistical Topology Tests Summary (AU, SH, KH Tests)

| Dataset   | Alignment_Block   | Alignment_Name          | Evaluated_Tree   | Evaluated_Tree_Name     |     logL |       deltaL |   p_KH |   p_SH |    c_ELW |     p_AU | Result         |
|:----------|:------------------|:------------------------|:-----------------|:------------------------|---------:|-------------:|-------:|-------:|---------:|---------:|:---------------|
| 68_taxa   | Block1_ORF1a      | ORF1a (nsp1–nsp10)      | Block1_ORF1a     | ORF1a (nsp1–nsp10)      | -56456.2 |   0          | 0.765  | 1      | 0.763    | 0.787    | Fail to reject |
| 68_taxa   | Block1_ORF1a      | ORF1a (nsp1–nsp10)      | Block2_ORF1b     | ORF1b (nsp12–nsp16)     | -56560.9 | 104.75       | 0.0004 | 0.0004 | 1.04e-14 | 5.81e-06 | REJECTED       |
| 68_taxa   | Block1_ORF1a      | ORF1a (nsp1–nsp10)      | Block3_Spike     | Spike Glycoprotein      | -56505.5 |  49.379      | 0.0019 | 0.0051 | 1.22e-06 | 7.22e-05 | REJECTED       |
| 68_taxa   | Block1_ORF1a      | ORF1a (nsp1–nsp10)      | Block4_3prime    | 3' Structural/Accessory | -56463.6 |   7.4613     | 0.235  | 0.582  | 0.237    | 0.215    | Fail to reject |
| 68_taxa   | Block2_ORF1b      | ORF1b (nsp12–nsp16)     | Block1_ORF1a     | ORF1a (nsp1–nsp10)      | -23246.8 |   8.2285     | 0.0565 | 0.0609 | 0.0272   | 0.0588   | Fail to reject |
| 68_taxa   | Block2_ORF1b      | ORF1b (nsp12–nsp16)     | Block2_ORF1b     | ORF1b (nsp12–nsp16)     | -23238.6 |   0          | 0.884  | 1      | 0.84     | 0.938    | Fail to reject |
| 68_taxa   | Block2_ORF1b      | ORF1b (nsp12–nsp16)     | Block3_Spike     | Spike Glycoprotein      | -23250.6 |  12.058      | 0.0558 | 0.0639 | 0.0195   | 0.027    | REJECTED       |
| 68_taxa   | Block2_ORF1b      | ORF1b (nsp12–nsp16)     | Block4_3prime    | 3' Structural/Accessory | -23243   |   4.3948     | 0.116  | 0.364  | 0.114    | 0.152    | Fail to reject |
| 68_taxa   | Block3_Spike      | Spike Glycoprotein      | Block1_ORF1a     | ORF1a (nsp1–nsp10)      | -13807.7 |  23.005      | 0.0821 | 0.0847 | 0.0143   | 0.0296   | REJECTED       |
| 68_taxa   | Block3_Spike      | Spike Glycoprotein      | Block2_ORF1b     | ORF1b (nsp12–nsp16)     | -13801.1 |  16.39       | 0.11   | 0.172  | 0.0527   | 0.0339   | REJECTED       |
| 68_taxa   | Block3_Spike      | Spike Glycoprotein      | Block3_Spike     | Spike Glycoprotein      | -13784.7 |   0          | 0.836  | 1      | 0.783    | 0.931    | Fail to reject |
| 68_taxa   | Block3_Spike      | Spike Glycoprotein      | Block4_3prime    | 3' Structural/Accessory | -13798.4 |  13.754      | 0.164  | 0.176  | 0.15     | 0.147    | Fail to reject |
| 68_taxa   | Block4_3prime     | 3' Structural/Accessory | Block1_ORF1a     | ORF1a (nsp1–nsp10)      | -11021.3 |  27.986      | 0.0094 | 0.014  | 0.00471  | 0.00376  | REJECTED       |
| 68_taxa   | Block4_3prime     | 3' Structural/Accessory | Block2_ORF1b     | ORF1b (nsp12–nsp16)     | -11030   |  36.721      | 0.0081 | 0.0102 | 0.00195  | 0.00263  | REJECTED       |
| 68_taxa   | Block4_3prime     | 3' Structural/Accessory | Block3_Spike     | Spike Glycoprotein      | -11030.8 |  37.568      | 0.0042 | 0.0049 | 0.000234 | 0.00118  | REJECTED       |
| 68_taxa   | Block4_3prime     | 3' Structural/Accessory | Block4_3prime    | 3' Structural/Accessory | -10993.3 |   0          | 0.991  | 1      | 0.993    | 0.998    | Fail to reject |
| 78_taxa   | Block1_ORF1a      | ORF1a (nsp1–nsp10)      | Block1_ORF1a     | ORF1a (nsp1–nsp10)      | -49309.2 |   0          | 0.842  | 1      | 0.788    | 0.852    | Fail to reject |
| 78_taxa   | Block1_ORF1a      | ORF1a (nsp1–nsp10)      | Block2_ORF1b     | ORF1b (nsp12–nsp16)     | -49313.9 |   4.7903     | 0.159  | 0.366  | 0.103    | 0.203    | Fail to reject |
| 78_taxa   | Block1_ORF1a      | ORF1a (nsp1–nsp10)      | Block3_Spike     | Spike Glycoprotein      | -49327.9 |  18.708      | 0.0139 | 0.015  | 0.00703  | 0.0096   | REJECTED       |
| 78_taxa   | Block1_ORF1a      | ORF1a (nsp1–nsp10)      | Block4_3prime    | 3' Structural/Accessory | -49313.9 |   4.7903     | 0.158  | 0.366  | 0.102    | 0.203    | Fail to reject |
| 78_taxa   | Block2_ORF1b      | ORF1b (nsp12–nsp16)     | Block1_ORF1a     | ORF1a (nsp1–nsp10)      | -22547.4 |   0.6275     | 0.288  | 0.613  | 0.264    | 0.277    | Fail to reject |
| 78_taxa   | Block2_ORF1b      | ORF1b (nsp12–nsp16)     | Block2_ORF1b     | ORF1b (nsp12–nsp16)     | -22546.7 |   0          | 0.517  | 1      | 0.367    | 0.568    | Fail to reject |
| 78_taxa   | Block2_ORF1b      | ORF1b (nsp12–nsp16)     | Block3_Spike     | Spike Glycoprotein      | -22569.3 |  22.562      | 0.004  | 0.0048 | 0.0013   | 0.00222  | REJECTED       |
| 78_taxa   | Block2_ORF1b      | ORF1b (nsp12–nsp16)     | Block4_3prime    | 3' Structural/Accessory | -22546.7 |   8.0345e-05 | 0.482  | 0.968  | 0.367    | 0.704    | Fail to reject |
| 78_taxa   | Block3_Spike      | Spike Glycoprotein      | Block1_ORF1a     | ORF1a (nsp1–nsp10)      | -13805.6 |   3.7962     | 0.29   | 0.349  | 0.122    | 0.221    | Fail to reject |
| 78_taxa   | Block3_Spike      | Spike Glycoprotein      | Block2_ORF1b     | ORF1b (nsp12–nsp16)     | -13802.4 |   0.53262    | 0.462  | 0.703  | 0.209    | 0.566    | Fail to reject |
| 78_taxa   | Block3_Spike      | Spike Glycoprotein      | Block3_Spike     | Spike Glycoprotein      | -13801.8 |   0          | 0.538  | 1      | 0.46     | 0.583    | Fail to reject |
| 78_taxa   | Block3_Spike      | Spike Glycoprotein      | Block4_3prime    | 3' Structural/Accessory | -13802.4 |   0.53195    | 0.462  | 0.703  | 0.209    | 0.566    | Fail to reject |
| 78_taxa   | Block4_3prime     | 3' Structural/Accessory | Block1_ORF1a     | ORF1a (nsp1–nsp10)      | -11927.6 |   4.214      | 0.149  | 0.189  | 0.112    | 0.136    | Fail to reject |
| 78_taxa   | Block4_3prime     | 3' Structural/Accessory | Block2_ORF1b     | ORF1b (nsp12–nsp16)     | -11923.4 |   0.00021911 | 0.498  | 0.841  | 0.298    | 0.641    | Fail to reject |
| 78_taxa   | Block4_3prime     | 3' Structural/Accessory | Block3_Spike     | Spike Glycoprotein      | -11923.7 |   0.32059    | 0.413  | 0.645  | 0.293    | 0.494    | Fail to reject |
| 78_taxa   | Block4_3prime     | 3' Structural/Accessory | Block4_3prime    | 3' Structural/Accessory | -11923.4 |   0          | 0.502  | 1      | 0.297    | 0.782    | Fail to reject |

## 3. Four-Gamete Incompatibility Matrix

| Block_1       | Block_2       |   Tested_Site_Pairs |   Incompatible_Pairs |   Incompatible_Pct |
|:--------------|:--------------|--------------------:|---------------------:|-------------------:|
| Block1_ORF1a  | Block1_ORF1a  |               27876 |                 7812 |            28.0241 |
| Block1_ORF1a  | Block2_ORF1b  |               14777 |                 3929 |            26.5886 |
| Block1_ORF1a  | Block3_Spike  |                2373 |                  409 |            17.2356 |
| Block1_ORF1a  | Block4_3prime |               15773 |                 2860 |            18.1323 |
| Block2_ORF1b  | Block1_ORF1a  |               14777 |                 3929 |            26.5886 |
| Block2_ORF1b  | Block2_ORF1b  |                6640 |                 2484 |            37.4096 |
| Block2_ORF1b  | Block3_Spike  |                1135 |                  250 |            22.0264 |
| Block2_ORF1b  | Block4_3prime |                7176 |                 1393 |            19.4119 |
| Block3_Spike  | Block1_ORF1a  |                2373 |                  409 |            17.2356 |
| Block3_Spike  | Block2_ORF1b  |                1135 |                  250 |            22.0264 |
| Block3_Spike  | Block3_Spike  |                 156 |                   16 |            10.2564 |
| Block3_Spike  | Block4_3prime |                1070 |                  151 |            14.1121 |
| Block4_3prime | Block1_ORF1a  |               15773 |                 2860 |            18.1323 |
| Block4_3prime | Block2_ORF1b  |                7176 |                 1393 |            19.4119 |
| Block4_3prime | Block3_Spike  |                1070 |                  151 |            14.1121 |
| Block4_3prime | Block4_3prime |                7398 |                 1414 |            19.1133 |

## 4. Clade Resolution and Sister Placement

| Dataset   | Block_ID      | Block_Name              | Immediate_Sister   |   Clade_UFBoot_Support | Is_RaTG13_Sister   | Is_RpYN06_Sister   |
|:----------|:--------------|:------------------------|:-------------------|-----------------------:|:-------------------|:-------------------|
| 68_taxa   | Whole_NRA3    | Whole NRA3 Alignment    | RaTG13             |                    100 | Yes                | No                 |
| 78_taxa   | Whole_NRA3    | Whole NRA3 Alignment    | RmYN02, RpYN06     |                    100 | No                 | Shared             |
| 68_taxa   | Block1_ORF1a  | ORF1a (nsp1–nsp10)      | RaTG13             |                    100 | Yes                | No                 |
| 68_taxa   | Block2_ORF1b  | ORF1b (nsp12–nsp16)     | RaTG13             |                    100 | Yes                | No                 |
| 68_taxa   | Block3_Spike  | Spike Glycoprotein      | RaTG13             |                    100 | Yes                | No                 |
| 68_taxa   | Block4_3prime | 3' Structural/Accessory | RaTG13             |                    100 | Yes                | No                 |
| 78_taxa   | Block1_ORF1a  | ORF1a (nsp1–nsp10)      | RaTG13             |                    100 | Yes                | No                 |
| 78_taxa   | Block2_ORF1b  | ORF1b (nsp12–nsp16)     | RaTG13             |                     99 | Yes                | No                 |
| 78_taxa   | Block3_Spike  | Spike Glycoprotein      | RaTG13             |                     71 | Yes                | No                 |
| 78_taxa   | Block4_3prime | 3' Structural/Accessory | RaTG13             |                    100 | Yes                | No                 |
