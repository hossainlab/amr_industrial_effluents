# Genomic Profiling of Extended-Spectrum ARGs and Mobilizable Resistance Determinants in Industrial Effluents

## Shotgun Metagenomic Analysis Report

**Project:** AMR in Industrial Effluents (Bangladesh)
**Samples:** IE1, IE2, IE3, IE4 (paired-end Illumina, 150 bp)
**Date:** 2026-09-03
**Analysis pipeline:** fastp QC/trimming → MEGAHIT assembly → ABRicate ARG/MGE/virulence annotation → ISEScan IS detection → ARG-MGE co-localization + Mobility Index

---

## 1. Quality Control and Trimming

All four samples were processed with fastp (v1.3.6) for adapter detection, quality trimming (sliding window Q20, cut_front/cut_tail), and length filtering (≥50 bp).

| Sample | Raw reads | Raw bases | Trimmed reads | Trimmed bases | Q30 (%) | GC (%) | Retained (%) |
|--------|-----------|-----------|---------------|---------------|---------|--------|--------------|
| IE1 | 82,637,880 | 12.48 Gb | 80,836,740 | 12.18 Gb | 93.1 | 58.4 | 97.8 |
| IE2 | 79,744,950 | 12.04 Gb | 78,675,110 | 11.86 Gb | 94.7 | 50.6 | 98.7 |
| IE3 | 68,004,508 | 10.27 Gb | 67,168,486 | 10.13 Gb | 95.3 | 46.3 | 98.8 |
| IE4 | 76,832,668 | 11.60 Gb | 75,698,726 | 11.42 Gb | 95.0 | 48.4 | 98.5 |

**Total:** ~302M read pairs retained across 4 samples. Q30 rates of 93–95% indicate high-quality sequencing. IE1 has an elevated GC content (58.4%) compared to IE2–IE4 (46–51%), which may reflect a distinct microbial community composition and warrants investigation with metadata.

QC reports (fastp HTML per sample + MultiQC summary) are saved alongside this report.

---

## 2. Metagenomic Assembly

MEGAHIT (v1.2.9, meta mode, --min-contig-len 1000) was used for per-sample assembly and a co-assembly of all four samples.

| Assembly | Contigs | Total bp | N50 (bp) | Max contig (bp) | Runtime |
|----------|---------|----------|----------|-----------------|---------|
| IE1 | 93,406 | 322.6 Mb | 5,414 | 1,300,145 | ~84 min |
| IE2 | 102,355 | 335.8 Mb | 4,830 | 619,544 | ~126 min |
| IE3 | 127,765 | 351.3 Mb | 3,410 | 486,665 | ~129 min |
| IE4 | 112,884 | 306.1 Mb | 3,300 | 222,504 | ~117 min |
| Co-assembly | 419,531 | 1,262.9 Mb | 4,060 | 1,300,139 | ~583 min |

The modest N50 (3.3–5.4 kb) reflects the inherent limitation of short-read metagenomic assembly — the "assembly bottleneck" described by Kerkvliet et al. (2024) [44]. ARG-MGE co-localization results must be interpreted as **putative** in this context: short contigs may fragment what are in reality single mobile genetic units. Long-read sequencing would be needed to resolve complete mobilization contexts [15, 17].

---

## 3. Resistome Catalog

ARG annotation was performed with ABRicate (v1.4.0) against three databases (CARD, ResFinder, NCBI AMR). Hits were merged and deduplicated by sample/contig/gene/position. A total of **1,469 unique ARG hits** were identified across 15 drug classes.

### 3.1 ARG burden per sample

| Sample | Unique ARG genes | ARG-bearing contigs |
|--------|-----------------|---------------------|
| IE1 | 130 | 68 |
| IE2 | 18 | 11 |
| IE3 | 200 | 127 |
| IE4 | 228 | 130 |
| Co-assembly | 345 | 234 |

IE2 is a striking outlier with only 18 unique ARG genes (11 contigs), while IE4 carries the highest ARG diversity (228 genes, 130 contigs). IE3 and IE4 together account for the majority of resistance determinants.

### 3.2 ARG by drug class (unique genes)

| Drug class | IE1 | IE2 | IE3 | IE4 | Co-asm |
|-----------|-----|-----|-----|-----|--------|
| Beta-lactam | 50 | 6 | 70 | 56 | 109 |
| Aminoglycoside | 17 | 0 | 27 | 48 | 61 |
| Macrolide/Lincosamide/Streptogramin | 27 | 9 | 50 | 53 | 75 |
| Tetracycline/Glycylcycline | 15 | 1 | 42 | 29 | 55 |
| Fluoroquinolone/Quinolone | 8 | 0 | 23 | 14 | 40 |
| Phenicol | 10 | 0 | 19 | 19 | 37 |
| Sulfonamide/Diaminopyrimidine | 7 | 2 | 12 | 19 | 23 |
| Colistin/Polymyxin | 0 | 0 | 5 | 2 | 12 |
| Disinfectant/QAC | 2 | 0 | 13 | 4 | 19 |
| Fosfomycin | 2 | 0 | 2 | 0 | 4 |
| Oxazolidinone | 4 | 0 | 0 | 5 | 7 |
| Rifamycin | 2 | 1 | 9 | 1 | 16 |
| Nitroimidazole | 0 | 0 | 1 | 0 | 2 |
| Nucleoside | 0 | 0 | 3 | 2 | 5 |

### 3.3 Critical resistance determinants

The following high-priority ARGs were detected:

**Carbapenemases:**
- **blaNDM-1** — detected in IE4 only (contig k141_275398, 1,178 bp). NDM-1 is a metallo-β-lactamase conferring resistance to virtually all β-lactams including carbapenems. Its presence in industrial effluent is a serious public health concern, as NDM-1 is highly transmissible via mobile genetic elements. However, the 1,178 bp contig is too short to resolve its mobilization context — this is a direct consequence of the short-read assembly bottleneck [44].
- **GES-5** (carbapenem-hydrolyzing) — detected in IE1. GES-type enzymes with carbapenemase activity are emerging in Pseudomonas and Acinetobacter.
- **IDC-1** (integron-encoded cephalosporinase with carbapenem activity) — detected in IE1.

**ESBLs and other β-lactamases:**
- **blaVEB-1, VEB-9, VEB-15** — present in all four samples. VEB-type ESBLs are common in Pseudomonas and Acinetobacter in South and Southeast Asia.
- **blaGES-2, GES-14** — IE1, IE3.
- **blaPER-1** — IE4. PER-1 is an ESBL frequently found in Acinetobacter baumannii.
- **blaTEM-1** — IE1, IE3, IE4 (broad-spectrum β-lactamase, not ESBL).
- **blaCMY-4, CMY-32, CMY-59** — IE1, IE4 (AmpC cephalosporinases).
- **blaDHA-1** — IE4 (plasmid-mediated AmpC).
- **53 unique OXA-type variants** — including OXA-1, OXA-10, OXA-1108, OXA-1136, etc. OXA enzymes are particularly diverse in Acinetobacter.

**Notable absence:** No CTX-M or SHV β-lactamases were detected in any sample. CTX-M (especially CTX-M-15) is the dominant ESBL in clinical Enterobacteriaceae in South Asia [25]. The absence of CTX-M and the predominance of VEB/GES/PER/OXA suggests these effluent microbial communities are dominated by non-Enterobacteriaceae (likely Pseudomonas/Acinetobacter), which carry different β-lactamase repertoires. This is consistent with environmental/industrial effluent communities differing from clinical gut flora.

**Colistin resistance:**
- **mcr-5.1** — detected in IE4 (contig k141_6907, 7,558 bp) and the co-assembly. MCR-5.1 is a phosphoethanolamine transferase conferring colistin resistance. Unlike mcr-1 (the most globally widespread variant), mcr-5 was originally described in Salmonella and Enterococcus. Its detection in industrial effluent extends the known environmental reservoir of mobilizable colistin resistance [32].

**16S rRNA methyltransferases (high-level aminoglycoside resistance):**
- **rmtF, rmtF2** — IE4 (contig k141_183049, 1,341 bp). These confer high-level resistance to all clinically important aminoglycosides including amikacin and gentamicin. rmtF is increasingly reported in South Asia.
- **rmtD2** — IE1.

**Fluoroquinolone resistance:**
- **18 unique qnr variants** (QnrA1, QnrS1, QnrS2, QnrVC1, QnrVC4, QnrVC6) — across IE1, IE3, IE4. QnrVC variants are integron-mediated, consistent with their co-localization patterns (see Section 5).

**Other notable determinants:**
- **sul1, sul2, sul3, sul4** — sulfonamide resistance (IE2–IE4). sul4 is a recently described sulfonamide resistance gene.
- **dfrA1, dfrA14, dfrA15, dfrA16, dfrA19** — trimethoprim resistance (all samples).
- **52 unique tet variants** — including tet(A), tet(M), tet(32), tet(39), tet(40) (all samples).
- **qacL, qacEdelta1** — biocide/disinfectant resistance genes (IE1, IE3, IE4). IE3 has 13 disinfectant/QAC resistance genes, the highest of any sample — potentially relevant for industrial settings where biocides are used.
- **floR, cmlA5, cmlB1** — phenicol resistance.
- **ere(A), ere(B), mph(A–G), lnu(B–I)** — macrolide resistance.

---

## 4. Mobile Genetic Elements

### 4.1 Insertion sequences (ISEScan)

| Sample | IS elements | Contigs with IS | IS families |
|--------|-------------|-----------------|-------------|
| IE1 | 1,783 | 1,559 | 27 |
| IE2 | 1,011 | 927 | 26 |
| IE3 | 3,827 | 3,542 | 27 |
| IE4 | 2,853 | 2,647 | 28 |

IS elements are abundant across all samples, with IE3 carrying the highest burden (3,827). IS families detected include IS3, IS5, IS21, IS91, IS110, IS1380, IS1634, IS200/IS605, and others. IS91 is of particular interest as it uses rolling-circle transposition and is known to mobilize adjacent resistance genes.

### 4.2 Plasmid replicons (PlasmidFinder)

| Sample | Replicon types |
|--------|---------------|
| IE1 | IncA/C2, IncN5, IncP1, IncQ1, IncQ2, IncW, Col440II |
| IE2 | None detected |
| IE3 | IncA/C, IncP1, IncP6, IncQ1, IncQ2, Col440II, Col(pHAD28), ColKP3, repUS2, repUS43 |
| IE4 | IncA/C, IncA/C2, IncQ1, IncQ2, IncU, Col(pHAD28), ColE10, rep33, repUS43, repUS50 |
| Co-asm | All above + IncI1, IncW |

**Broad-host-range plasmids** (IncA/C, IncA/C2, IncP1, IncW) are present in IE1, IE3, and IE4. These incompatibility groups are associated with multidrug resistance and can transfer across diverse bacterial taxa — their presence in industrial effluent indicates potential for inter-species ARG dissemination. IE2 has no detectable plasmid replicons, consistent with its low overall ARG/MGE profile.

### 4.3 Integron-associated ARGs

Sixty-five ARG hits have product descriptions indicating integron association (e.g., "integron-encoded aminoglycoside nucleotidyltransferase", "integron-mediated quinolone resistance"). These include aadA variants, dfrA variants, ereA, qnrVC, and OXA enzymes. This serves as proxy evidence for class 1 integron-mediated resistance gene capture, though precise integron boundaries were not determined (IntegronFinder was unavailable; see Limitations).

---

## 5. ARG-MGE Co-localization and Mobility Index

### 5.1 Mobility Index

The Mobility Index (MI) was computed as the fraction of ARG-bearing contigs (≥1 kb) that carry at least one MGE (IS, plasmid replicon, or integron-associated gene) within 5 kb. Both raw MI and overlap-adjusted MI are reported — the latter excludes MGE hits overlapping the ARG by >50% of ARG length, correcting for cargo-sequence homology artifacts as recommended by Slizovskiy et al. (2020) [46].

| Sample | ARG contigs | Contigs with MGE | MI (raw) | MI (overlap-adj.) |
|--------|-------------|------------------|----------|-------------------|
| IE1 | 68 | 21 | 0.309 | 0.250 |
| IE2 | 11 | 3 | 0.273 | 0.182 |
| IE3 | 127 | 21 | 0.165 | 0.118 |
| IE4 | 130 | 38 | 0.292 | 0.185 |
| Co-assembly | 234 | 25 | 0.107 | 0.034 |

**Interpretation:** IE1 and IE4 show the highest mobilization potential (~29–31% of ARG contigs co-localized with MGEs). The overlap-adjusted MI is consistently lower, confirming that a portion of raw co-localizations represent annotation overlap rather than true physical linkage [46]. The co-assembly MI is lower partly because IS element data was not available for the co-assembly (see Limitations).

### 5.2 Key co-localization events

| ARG | Sample | Contig | Contig length | MGE | Distance (bp) | Significance |
|-----|--------|--------|--------------|-----|---------------|--------------|
| rmtF/rmtF2 | IE4 | k141_183049 | 1,341 bp | IS91 | 36 | IS91 uses rolling-circle transposition; 36 bp proximity strongly suggests IS91-mediated mobilization of this 16S rRNA methyltransferase |
| blaTEM-1 | IE3 | k141_144258 | 1,430 bp | IS3 | 646 | IS3 family element near TEM-1 β-lactamase |
| OXA-1108 | IE4 | k141_188303 | 2,064 bp | IS3 | 850 | IS3 near OXA-1108 carbapenemase |
| OXA-1136 | IE3 | k141_202590 | 7,259 bp | Integron cassette (aadA + qacL + OXA-21) | 810 | Classic integron cassette arrangement on a well-assembled contig |
| qnrVC1/4/6 | IE1, IE3, IE4 | multiple | 1.3–1.9 kb | Integron-associated | 0 | qnrVC variants are integron-encoded by definition |
| blaVEB-1/9 | Co-asm | k141_412928 | — | Integron-associated | 0 | VEB ESBLs are typically integron-encoded |

### 5.3 Critical determinants without detectable MGE co-localization

- **blaNDM-1** (IE4, 1,178 bp contig) — no MGE within 5 kb. The contig is too short to resolve context; NDM-1 is almost always associated with Tn125 transposon and IncX3 plasmids in clinical isolates, but this cannot be confirmed from the assembly.
- **mcr-5.1** (IE4, 7,558 bp contig) — no IS, plasmid replicon, or integron element within 5 kb. mcr-5.1 has been reported on small ColE-type plasmids and in chromosomal contexts; the absence of plasmid replicon on this contig may indicate chromosomal integration or an unresolved plasmid context.

---

## 6. Virulence Factor Co-occurrence

VFDB annotation identified **509 virulence factor hits** across all samples, with **16 hits on ARG-bearing contigs**. The virulence genes co-occurring with ARGs include:

- **adeF, adeG** — AdeFGH efflux pump (Acinetobacter baumannii); biofilm formation
- **ompA** — outer membrane protein A (E. coli); invasion
- **acrA, acrB** — AcrAB efflux pump (Klebsiella pneumoniae); antimicrobial activity/competitive advantage
- **pmrA** — PmrAB two-component regulator (Salmonella); regulates polymyxin resistance modifications
- **lpsB, lpxL** — LPS biosynthesis (Acinetobacter baumannii); immune modulation

The co-occurrence of efflux pump genes (adeFGH, acrAB) with ARGs on the same contigs is notable, as these systems can confer both virulence and multidrug resistance phenotypes simultaneously. The pmrA regulator on an ARG-bearing contig is particularly relevant given the detection of mcr-5.1 in the same sample (IE4), as PmrAB and MCR systems both modify lipid A to confer colistin resistance.

---

## 7. Sample Comparison Summary

| Feature | IE1 | IE2 | IE3 | IE4 |
|---------|-----|-----|-----|-----|
| Read pairs (M) | 40.4 | 39.3 | 33.6 | 37.8 |
| GC (%) | 58.4 | 50.6 | 46.3 | 48.4 |
| Assembly N50 (kb) | 5.4 | 4.8 | 3.4 | 3.3 |
| Unique ARG genes | 130 | 18 | 200 | 228 |
| IS elements | 1,783 | 1,011 | 3,827 | 2,853 |
| Plasmid replicon types | 7 | 0 | 10 | 10 |
| MI (raw) | 0.309 | 0.273 | 0.165 | 0.292 |
| MI (adj.) | 0.250 | 0.182 | 0.118 | 0.185 |
| Virulence hits | 22 | 9 | 173 | 47 |
| NDM-1 | — | — | — | **Present** |
| mcr-5.1 | — | — | — | **Present** |
| VEB ESBL | **Present** | **Present** | **Present** | **Present** |
| GES (incl. carbapenemase) | **Present** | — | **Present** | — |
| PER-1 | — | — | — | **Present** |
| rmtF/rmtF2 | — | — | — | **Present** |
| rmtD2 | **Present** | — | — | — |
| IDC-1 | **Present** | — | — | — |

**IE4** is the highest-risk sample, carrying both blaNDM-1 and mcr-5.1, the highest ARG diversity (228 genes), rmtF/rmtF2, PER-1, DHA-1, and the highest number of ARG contigs with MGE co-localization (38 contigs, MI=0.292).

**IE2** is the lowest-burden sample with only 18 ARG genes, no plasmid replicons, and no critical determinants (NDM, mcr, rmt, GES, PER). Its resistome is limited to VEB-1, sul1/sul2, mph, and blaPAU-1.

**IE3** has the highest IS element burden (3,827) and virulence factor hits (173), with a well-assembled integron cassette (7,259 bp contig carrying aadA + qacL + OXA-1136 + OXA-21).

**IE1** has the highest GC content (58.4%), the highest overlap-adjusted MI (0.250), and carries GES-5 (carbapenemase), IDC-1, and rmtD2.

---

## 8. Limitations

1. **Sample size (n=4):** This is a descriptive/cataloguing study. No between-group statistical comparisons are possible. Findings are hypothesis-generating.

2. **Short-read assembly bottleneck:** N50 of 3.3–5.4 kb limits resolution of ARG-MGE co-localization. Critical genes (NDM-1 on 1,178 bp contig) cannot have their mobilization context resolved. Long-read sequencing (Nanopore/PacBio) would be needed to reconstruct complete mobile genetic units [15, 17, 44]. All co-localization results are **putative**.

3. **Co-assembly IS data incomplete:** ISEScan on the co-assembly did not complete due to compute time constraints. Co-assembly MI is therefore based on plasmid replicon and integron evidence only, and is likely underestimated.

4. **Plasmid classification:** PlasClass and PlasFlow (per-contig plasmid probability prediction) could not be installed due to unresolved legacy Python dependency conflicts. Plasmid evidence is limited to replicon detection (PlasmidFinder), which misses plasmids lacking known replicon markers and cannot classify contigs without replicon hits.

5. **Integron boundary mapping:** IntegronFinder was unavailable. Integron evidence is proxy-based — ARGs with "integron-encoded" product descriptions were classified as integron-associated. This means integron co-localization is partially self-referential (the ARG itself is the evidence). Precise integron boundary definition (attC sites, intI integrase) was not performed. IS and plasmid replicon co-localizations are independent evidence and more robust.

6. **Host/human read removal:** Not performed. Human reads in environmental metagenomes are unlikely to assemble into bacterial ARG-bearing contigs, but this step is recommended before data deposition to public repositories.

7. **Taxonomic profiling:** Kraken2/Bracken taxonomic profiling was not completed (database download was time-prohibitive within compute constraints). The bacterial community composition underlying the resistome is therefore inferred indirectly (e.g., from ARG types and virulence factor species associations) rather than directly profiled.

8. **Normalization:** ARG abundance is reported as hit counts, not normalized by read depth, contig length, or 16S rRNA gene copies. Cross-sample comparisons of absolute abundance should be interpreted cautiously. Read-recruitment (mapping trimmed reads back to assemblies) would enable RPKM/TPM normalization.

9. **Sample metadata:** Correct sample metadata (site description, effluent type, location, collection date, physicochemical parameters) was not available at analysis time. Metadata-dependent interpretation (e.g., correlating ARG profiles with industrial process type or antibiotic usage) is deferred.

10. **AMRFinderPlus and RGI:** These tools could not be installed due to conda dependency conflicts. ABRicate with CARD, ResFinder, and NCBI databases provides robust ARG detection, but AMRFinderPlus would add protein-level variant identification and RGI would provide CARD-specific curated resistance models.

---

## 9. Key Conclusions

1. **Industrial effluents in Bangladesh harbor a diverse and clinically significant resistome**, with 1,469 unique ARG hits across 15 drug classes. The detection of **blaNDM-1** (carbapenemase) and **mcr-5.1** (colistin resistance) in IE4 is particularly alarming, as these target last-resort antibiotics.

2. **The β-lactamase profile is dominated by VEB/GES/PER/OXA rather than CTX-M/SHV**, suggesting non-Enterobacteriaceae (likely Pseudomonas/Acinetobacter) dominate these effluent communities. This differs from the typical clinical ESBL landscape in South Asia and highlights the importance of environmental surveillance beyond clinical settings.

3. **ARG mobilization potential is substantial**: 16–31% of ARG-bearing contigs co-localize with MGEs (raw MI), with IS elements (particularly IS91 near rmtF at 36 bp) providing the strongest evidence for active mobilization. Broad-host-range plasmids (IncA/C, IncP1, IncW) are present in 3 of 4 samples.

4. **IE4 represents the highest-risk sample**, combining critical determinants (NDM-1, mcr-5.1, rmtF, PER-1), high ARG diversity, and significant MGE co-localization.

5. **The assembly bottleneck fundamentally limits mobilization assessment** — the NDM-1 gene sits on a 1,178 bp contig with no resolvable MGE context. This underscores the need for long-read metagenomics in environmental AMR surveillance, as emphasized by Kerkvliet et al. (2024) [44] and demonstrated by Dai et al. (2022) [15] and Che et al. (2019) [17].

6. **Integron-mediated resistance is widespread**, with qnrVC, aadA, dfrA, and OXA variants showing integron association across samples. A well-assembled integron cassette (7,259 bp) in IE3 carries aadA + qacL + OXA-1136 + OXA-21 — a classic multi-resistance integron arrangement.

---

## 10. Deliverables

| File | Description |
|------|-------------|
| `arg_all_hits.tsv` | All 1,469 deduplicated ARG hits (sample, contig, gene, position, database, product, resistance) |
| `arg_by_drug_class.tsv` | Unique ARG gene counts by drug class per sample |
| `mge_all_hits.tsv` | All MGE hits (IS, plasmid replicons, integron-associated) |
| `colocalization_table.tsv` | 392 ARG-MGE co-localization events (gene, MGE type, distance) |
| `mobility_index.tsv` | Raw and overlap-adjusted Mobility Index per sample |
| `virulence_all_hits.tsv` | All 509 VFDB virulence factor hits |
| `virulence_on_arg_contigs.tsv` | 16 virulence hits co-occurring on ARG-bearing contigs |
| `qc_summary.tsv` | QC statistics per sample |
| `figures/arg_heatmap.png` | Heatmap of top 30 ARG gene families per sample |
| `figures/key_genes_bar.png` | Bar plot of key ARG determinants per sample |
| `figures/mobility_index_bar.png` | Raw vs overlap-adjusted Mobility Index per sample |
| `qc/` | fastp HTML reports + MultiQC report |

---

## References

Key citations from the accompanying literature review:
- [9] Sharif et al. (2025) — Bangladesh environmental water PCR surveillance
- [11] McInnes et al. (2021) — Bangladesh surface-water metagenome
- [15] Dai et al. (2022) — long-read metagenomics for ARG resolution
- [17] Che et al. (2019) — nanopore plasmid resolution
- [25] Mohsin et al. (2021) — tet(X4) in South Asia
- [32] Ben Khedher et al. (2020) — mcr environmental reservoirs
- [34] Bengtsson-Palme et al. (2019) — industrial WWTP resistome
- [41] Ibitoye et al. (2025) — systematic review of environmental AMR
- [44] Kerkvliet et al. (2024) — assembly bottleneck in metagenomic ARG-MGE analysis
- [46] Slizovskiy et al. (2020) — co-localization methodology critique
