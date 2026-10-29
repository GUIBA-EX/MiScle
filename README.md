# MiScle

[中文](README.zh-CN.md)

**MiScle, a set of PCR primer candidates for metabarcoding environmental DNA from scleractinian corals**

## Abstract

MiScle compares two homologous mitochondrial 16S windows for scleractinian coral eDNA. S retains ordinary and D pairs. M now uses one shared forward primer and two reverse primers prioritizing Refertina and Vacatina: MR+MV for the ordinary formulation and MRD+MVD for the degenerate formulation. R/V denote design emphasis and D denotes degeneracy. Both suborders are targets; cross-coverage is allowed. Each M formulation contains three distinct oligos, with the shared F added once.

Among 304 evaluable species at M, the ordinary combination covers 282 (92.76%) with 78 exact paired matches (25.66%); the D combination covers 284 (93.42%) with 98 exact matches (32.24%). The original M and MD each covered 124 species, with 1 and 81 exact matches. Strict WC/GT background unions increase from 0/0 to 103/111 records. These are exploratory coverage-oriented candidates with a substantive background and exact-match tradeoff.

Under the shared-F constraint, the new M-D union exactly matches only 12/154 Vacatina (7.79%), versus 75/154 (48.70%) for the original MD. Combining S with the new M-D formulation gives 160/305 exact matches, versus 206/305 for the original SD+MD. The wider M coverage therefore does not resolve exact matching in both suborders simultaneously. Tissue and mixed-template validation should retain the original M/MD as controls.

## Primers and names

All target-specific sequences are written 5′→3′, without adapters. Y=C/T; R=A/G. S/M label homologous windows; product lengths vary among lineages. Degeneracy is capped at eight per oligo; each of the three M-D oligos selected here has degeneracy two.

| Pair | Forward, 5′→3′ | Reverse, 5′→3′ | F/R degeneracy |
|---|---|---|---:|
| MiScle-S | `TTTGATAAGTGGGACAGTTTGG` | `TTAACGGGTCACTTAAACCC` | 1 / 1 |
| MiScle-SD | `YTKGATAAGTGRGACAGTTTGG` | `YYAACGGGTCACTTWAACCC` | 8 / 8 |
| MiScle-MR | `AAGTCCCCATGGAACTTTAC` | `CTTAAGCCCGCCTTCGTTAC` | 1 / 1 |
| MiScle-MV | `AAGTCCCCATGGAACTTTAC` | `TCATAGTTCGTTTTCGTTAC` | 1 / 1 |
| MiScle-MRD | `AAGTCCCYATGGAACTTTAC` | `CTTAAGYCCGCCTTCGTTAC` | 2 / 2 |
| MiScle-MVD | `AAGTCCCYATGGAACTTTAC` | `CCATAGTTCRTTTTCGTTAC` | 2 / 2 |

MR and MV share exactly the same F; MRD and MVD share exactly the same F. All four M pairs use F1561–1580 and R1741–1760. Corresponding ordinary/D endpoints are identical. The two degeneracy values describe individual oligos, not the total complexity of the reaction.

[Primer table](primer_pairs.csv) · [FASTA](primers.fasta) · [Parameters](selected.json) · [Original comparator sequences](results/baseline/primer_pairs.csv)

## Materials and methods

### References and taxonomy

Hisata et al. (2025) developed a **12S** system. Here, annotated mitochondrial 16S was extracted from eligible Table 1 accessions and supplemented with public Scleractinia records, including non-reef and deep-water taxa. Of 110 audited paper accessions, 65 contained eligible 16S and 45 did not; the latter are not counted as amplification failures. A target query returned 1,974 records. Combined curation retained 879 eligible 16S records. One representative per source taxid and binomial species label was selected using gene completeness, length and low ambiguity, without selecting for primer compatibility.

The frozen panel contains 307 species labels, 131 genera and 29 source family labels. Source lineages assign 148 to Refertina, 157 to Vacatina and two as unassigned, using the framework recognized by WoRMS and Vaga et al. (2025). Guynia and Stenocyathus retain their unassigned source status. Family labels include older placements and do not assert 29 currently accepted families; genus and species results are also supplied. Only correctly oriented annotated 16S/mitochondrial large-subunit rRNA features were extracted. Records shorter than 150 nt or with more than 5% ambiguity were excluded. MAFFT was used for alignment (Katoh & Standley, 2013).

Coordinates are 1-based positions in the extracted, ungapped 16S gene of Acropora tenuis NC_003522.1, not whole-mitogenome positions. S is F1696–1717/R1826–1845 (150-bp reference span); M is F1561–1580/R1741–1760 (200 bp). Actual predicted primer-inclusive products among covered references are 147–209 bp for S, 186–298 bp for MR/MRD and 99–186 bp for MV/MVD.

![Reference binding positions](figures/reference_windows.svg)

### Fixed-coordinate suborder search

The initial design adopted alignment, conserved-anchor, lineage-compensation and paired-evaluation principles from MiFish (Miya et al., 2015). This round fixed all four M endpoints and generated candidates from independently verified binding templates and homologous columns of missed references. Separate pools were built for the shared F, Refertina R and Vacatina R. Ordinary candidates used observed A/C/G/T templates and consensus. D candidates iteratively merged variation supported by homologous target columns, following the variation-compensation logic of Anth-28S-eDNA (McCartin et al., 2024), with degeneracy ≤8 per oligo. Ambiguity was not restricted to the 3′ half.

The retained pools contained 40 F, 13 Refertina R and 32 Vacatina R candidates, comprising 79 unique oligos. All underwent raw-sequence site scans in both orientations, followed by paired orientation, spacing and terminal checks. Conservative and relaxed background screens were compared; the original low-background gate did not yield substantial Refertina compensation. The selected exploratory gate was ≤150 strict WC, ≤150 strict GT and ≤1,200 WC3 background records per pair; passing this gate does not imply negligible cross-amplification. Within a loss of at most two species from the best union coverage, exact matches were prioritized, followed by coverage and coverage of the weaker suborder. This finite heuristic search does not establish a global optimum.

### Thermodynamics and mixtures

Primer3 used 50 mM monovalent salt, 1.5 mM Mg²⁺ and 0.6 mM dNTP (Untergasser et al., 2012). Two scenarios were checked: 50 nM per concrete variant, and 200 nM total shared F with 100 nM total for each R. With twofold degeneracy for all three D oligos, the latter gives 100 nM per F variant and 50 nM per R variant. All expanded sequences were checked for hairpins, self-dimers and cross-dimers, including R–R interactions. Required limits were maximum structure Tm <55°C, F-to-either-R duplex-Tm midpoint difference ≤5°C and minimum duplex Tm ≥45°C. Maximum predicted structure Tm is 44.46°C for both M formulations. The final three positions must include a definite G/C base, and potential five-base homopolymers were excluded. Adapter-bearing full-length oligos were not evaluated.

### Coverage, background and independent verification

Coverage allows ≤2 mismatches per primer with IUPAC-compatible final three bases. Exact matching requires zero mismatches at both ends. Both sequence orientations are scanned; target products must be 70–500 bp with ≥20 internal bases, and background internal distances may be 0–1,000 bp. WC uses Watson–Crick compatibility; GT computationally permits G–T pairing. WC3, GT2 and WC4 are relaxed sensitivity screens, not laboratory detection probabilities.

The frozen background contains 9,145 records: 1,938 stratified cnidarian and 7,207 echinoderm references. Cnidarian queries were capped; bacteria, symbiotic algae and all marine biota are not systematically represented. M has 304 evaluable species and three unknowns. At least one of S/M is evaluable for 305 species. Missing binding-window information is not counted as failure.

The six released pairs and eight additional mixture pairs (F/F, R/R, V/V and R/V for each M formulation) were independently scanned by exhaustive C++ and bit-parallel Python implementations. All 84 target and strict-background coverage/exact set comparisons agree. The eight additional pairs have no target or strict WC/GT background hits. Mixture metrics are computational compatibility unions plus interaction checks; they do not establish absence of competition or amplification bias under a common PCR protocol. Cross-window S/M combinations describe separate assays followed by taxonomic aggregation.

## Results

| Pair | Evaluable | Covered | Exact | Strict WC/GT | Relaxed WC3/GT2/WC4 |
|---|---:|---:|---:|---:|---:|
| MiScle-S | 304 | 153 (50.33%) | 64 (21.05%) | 0 / 0 | 0 / 0 / 0 |
| MiScle-SD | 304 | 155 (50.99%) | 125 (41.12%) | 0 / 0 | 6 / 0 / 211 |
| MiScle-MR | 304 | 160 (52.63%) | 78 (25.66%) | 103 / 111 | 542 / 111 / 1152 |
| MiScle-MV | 304 | 122 (40.13%) | 0 (0.00%) | 0 / 0 | 0 / 0 / 0 |
| MiScle-MRD | 304 | 160 (52.63%) | 93 (30.59%) | 103 / 111 | 645 / 111 / 1168 |
| MiScle-MVD | 304 | 124 (40.79%) | 5 (1.64%) | 0 / 0 | 0 / 0 / 4 |

### M combinations versus original assays

| Scheme | Evaluable | Covered | Exact | Strict WC/GT background |
|---|---:|---:|---:|---:|
| Original M | 304 | 124 (40.79%) | 1 (0.33%) | 0 / 0 |
| MR+MV | 304 | 282 (92.76%) | 78 (25.66%) | 103 / 111 |
| Original MD | 304 | 124 (40.79%) | 81 (26.64%) | 0 / 0 |
| MRD+MVD | 304 | 284 (93.42%) | 98 (32.24%) | 103 / 111 |

The M-D union covers 145/148 Refertina (97.97%) and 137/154 Vacatina (88.96%), with 85/148 (57.43%) and 12/154 (7.79%) exact matches, respectively. Both unassigned species are covered, with one exact match. Twenty evaluable species remain missed at M-D. Across windows, SD+MRD+MVD covers 285/305 (93.44%) with 160/305 (52.46%) exact matches; the original SD+MD gave 279/305 (91.48%) and 206/305 (67.54%).

All 103 WC background hits are Actiniaria (70) and Corallimorpharia (33); GT gives 72 and 39. These represent 1.13% and 1.21% of all background records, or 5.31% and 5.73% of the cnidarian subset. The Refertina reverse primer introduces these hits. MRD has 645 WC3 and 1,168 WC4 hits, requiring close non-target mixed-template challenges.

[Per species](by_species.csv) · [Per suborder](by_suborder.csv) · [Combined suborders](combination_by_suborder.csv) · [Combinations](two_pair_comparison.csv) · [Background hit list](background_hits.csv) · [M-D misses](missed_species.csv) · [Original results](results/baseline/comparison.csv)

## Base matching and barcode resolution

Figures group references by suborder and source family. Both ends are shown 5′→3′, with R templates reverse-complemented. E denotes exact paired matching, C allowed-mismatch coverage, – no predicted hit and ? insufficient reference information. Blue indicates compatibility, orange mismatch and grey gap/unknown; the last three bases are boxed. Covered rows use independently verified raw binding templates; other rows show homologous reference columns. R/V names indicate design emphasis rather than exclusive amplification of that suborder.

[Base data](base_matching.csv) · [Verified templates and products](verified_templates.csv) · [Barcode resolution](barcode_resolution.csv) · [Shared barcodes](shared_barcodes.csv)

For MR/MRD, 48/160 covered species have an internal barcode not shared identically with another covered species; the counts are 24/122 for MV and 26/124 for MVD. Shared barcodes, including cross-genus sharing, remain. These are observations within the frozen library, not identification-success estimates. Improving binding compatibility does not automatically improve internal-barcode species resolution.

![MiScle-S](figures/base_matching_MiScle-S.svg)

![MiScle-SD](figures/base_matching_MiScle-SD.svg)

![MiScle-MR](figures/base_matching_MiScle-MR.svg)

![MiScle-MV](figures/base_matching_MiScle-MV.svg)

![MiScle-MRD](figures/base_matching_MiScle-MRD.svg)

![MiScle-MVD](figures/base_matching_MiScle-MVD.svg)

## Reproduction and files

- Run `python3 scripts/check_release.py` to check names, sequences, coordinates, shared F, degeneracy, figures and results.
- Run `python3 scripts/verify.py` for independent rescanning against frozen inputs using Python's standard library, without network access.
- [Verification summary](verification.json) · [84 comparisons](independent_checks.csv) · [Search and tradeoffs](results/suborder_search/README.md) · [Input manifest](manifest.json) · [Accession audit](data/paper_accession_audit.csv).
- New generation and selection scripts are documented in the search notes. Other early design scripts retain workspace-specific paths as provenance. Final release hashes are recorded in `SHA256SUMS`.

Priority experiments compare original M/MD and the R/V formulations on tissue DNA and balanced/unbalanced mixtures containing sea anemones and corallimorpharians. Further improvement of Vacatina exact matching should compare complete suborder-specific F/R pairs with the shared-F constraint relaxed, rather than simply adding more R ambiguity.

## References

1. Hisata, K., Nagata, T., Kanai, M., et al. (2025). An eDNA metabarcoding system for detecting scleractinian corals to the generic level along the Japanese coast. *Galaxea, Journal of Coral Reef Studies*, 27, 13–29. https://doi.org/10.3755/galaxea.G27D-5
2. Miya, M., Sato, Y., Fukunaga, T., et al. (2015). MiFish, a set of universal PCR primers for metabarcoding environmental DNA from fishes: detection of more than 230 subtropical marine species. *Royal Society Open Science*, 2, 150088. https://doi.org/10.1098/rsos.150088
3. McCartin, L., Saso, E., Vohsen, S. A., et al. (2024). Nuclear eDNA metabarcoding primers for anthozoan coral biodiversity assessment. *PeerJ*, 12, e18607. https://doi.org/10.7717/peerj.18607
4. Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution*, 30, 772–780. https://doi.org/10.1093/molbev/mst010
5. Untergasser, A., Cutcutache, I., Koressaar, T., et al. (2012). Primer3—new capabilities and interfaces. *Nucleic Acids Research*, 40, e115. https://doi.org/10.1093/nar/gks596
6. Vaga, C. F., Quattrini, A. M., Galvão de Lossio e Seiblitz, I., et al. (2025). A global coral phylogeny reveals resilience and vulnerability through deep time. *Nature*, 648, 377–382. https://doi.org/10.1038/s41586-025-09615-6
7. World List of Scleractinia. [Taxonomic classification](https://www.marinespecies.org/scleractinia/aphia.php?p=stats). Accessed 2026-10-30.

## Copyright

Copyright © 2026 Fei Xia (GUIBA-EX). All rights reserved. Republication and redistribution of protected original materials require prior written permission. Public database records retain their existing rights and conditions. See [Copyright and permissions](COPYRIGHT.md).
