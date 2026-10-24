# Literature Review, Gap Analysis & Contribution Framing

**Project:** Genomic Profiling of Extended-Spectrum Antimicrobial Resistance Genes (ARGs) and Mobilizable Resistance Determinants in Industrial Effluents

**Scope:** Bangladesh-centered, framed against global industrial-effluent resistome evidence.

**Review depth:** Rigorous narrative synthesis with structured gap/evidence tables (manuscript-introduction grade).

---

## 1. Background and Rationale

Antimicrobial resistance (AMR) is a leading global health threat, associated with an estimated 4.95 million deaths in 2019, of which ~1.27 million were directly attributable to drug-resistant infections [9]. The environment is now recognized as a central compartment in the One Health AMR framework: wastewater, industrial effluents, and receiving water bodies act as reservoirs where antibiotic selection pressure, dense microbial communities, and mobile genetic elements (MGEs) converge to promote horizontal gene transfer (HGT) of resistance determinants [34, 41].

Industrial effluents are of particular concern because they can carry high concentrations of antibiotics, heavy metals, biocides, and other co-selective agents that impose strong selective pressure on microbial communities far exceeding that of municipal sewage [34, 36]. Textile, tannery, and pharmaceutical manufacturing discharges have repeatedly been shown to enrich antibiotic resistance genes (ARGs) and alter microbial community structure in receiving waters [36, 39, 41]. Within this landscape, four resistance targets are especially consequential for clinical and public health:

- **ESBL genes (blaCTX-M, blaTEM, blaSHV) and metallo-β-lactamases (MBLs):** confer resistance to third-generation cephalosporins and carbapenems—critical antibiotics on the WHO AWaRe "Watch/Reserve" list.
- **mcr colistin-resistance genes:** plasmid-mediated determinants compromising colistin, a last-resort polymyxin; mcr-1 is now globally disseminated and harbors an environmental/water reservoir [32].
- **tet(X) family (tet(X4)):** flavin-dependent monooxygenases that inactivate tetracyclines including tigecycline, a last-resort glycylcycline; emerging on transferable plasmids in South Asia [25].
- **Class-1 integrons:** the primary capture/dissemination platform for multiple ARG cassettes and a widely used marker of anthropogenic ARG pollution and HGT potential [41, 46].

Characterizing the prevalence, genomic context, and mobilizable potential of these determinants in industrial effluents—particularly in low- and middle-income countries (LMICs) with limited wastewater treatment infrastructure—is therefore a research priority.

---

## 2. Literature Review

### 2.1 Global industrial (textile / tannery / pharmaceutical) effluent resistomes and MGEs

A growing body of metagenomic evidence confirms that industrial effluents enrich ARGs and reshape microbial communities. A recent systematic review (PRISMA 2020; 15 included studies, 2016–2024) concluded that industrial wastewater from textile, tannery, and pharmaceutical sources markedly increases the abundance and diversity of ARGs and MGEs in receiving water bodies, with β-lactamase, macrolide, and sulfonamide resistance genes and associated integrons/plasmids frequently detected, and wastewater treatment plants (WWTPs) showing limited ARG removal efficiency [41]. The review also flagged that the evidence base is geographically skewed (Europe and East Asia dominant) and sparse in less-industrialized and rural regions [41].

Pharmaceutical manufacturing effluents are the most intensively studied industrial source. Bengtsson-Palme et al. showed that an industrial WWTP receiving azithromycin-manufacturing waste in Croatia enriched ARGs and altered microbial community structure, with treated effluent still discharging resistance determinants into receiving sediments and prompting calls for regulatory limits on antibiotic emissions [34]. Subsequent work along antibiotic-concentration gradients from pharmaceutical WWTPs to receiving rivers demonstrated that ARG and MGE abundance tracks antibiotic concentrations, reinforcing the selection-pressure mechanism [35]. A 2026 preprint on a pharmaceutical manufacturing hub in Baddi, India, found ARGs and MGEs more diverse and abundant than in non-industrial areas, with clinically relevant ESBL and carbapenemase genes enriched in industrial effluent [40]. Earlier pharmaceutical-wastewater surveys detected 100% ARG prevalence across sampled sites with integrons (intl1, intl2) indicating active HGT potential [33].

Tannery effluents have been characterized metagenomically since Wang et al. (2013), who found abundant ARGs in tannery WWTP sludge (sul1 >20% of ARG reads; tet genes enriched anaerobically) alongside >70 insertion-sequence types and prevalent class-1 integrase genes [39]. A 2026 metagenomic snapshot of treated tannery effluent in Paraguay confirmed retention of tetracycline, β-lactam, and aminoglycoside ARGs despite treatment [42]. Industrial effluents have also been shown to boost resistome risk in coastal sediments, with higher ARG diversity and abundance in sediments than the water column [36]. Pesticide WWTPs similarly harbor ARGs whose mobility and host range are shaped by MGEs [16].

Across municipal and hospital WWTPs (the most-studied effluent types), metagenomic studies converge on a consistent message: treatment reduces overall ARG abundance but residual ARGs linked to MGEs persist in effluent, sustaining dissemination risk [13, 20, 22, 23, 37]. A global analysis of human wastewater treatment systems mapped ARG diversity and distribution across continents, identifying major carrier phyla and distinct WWTP-associated resistome composition [38]. Genome-resolved metagenomics in Wales found ~13.6% of wastewater bacteria carry ARGs [28], and a 165-metagenome genome-centric analysis underscored that MGEs are crucial for ARG transmission to pathogens in activated sludge and wastewater [18]. Long-read metagenomics further revealed that most WWTP ARGs are plasmid-borne [17], and that ARG–MGE associations shift from sewage to activated sludge [15], while HGT events and microdiversity are statistically linked to fluctuating antibiotic levels [19].

### 2.2 ARG–MGE co-localization and mobility inference methods

Inferring mobilizable ARGs requires linking resistance genes to MGEs (plasmids, insertion sequences, transposons, integrons, phages) on the same genomic context. The dominant short-read workflow is: metagenomic assembly → contig-level ARG annotation (AMRFinderPlus, RGI/CARD, ResFinder via ABRicate, DeepARG) → MGE/plasmid prediction (PlasClass, PlasFlow, ISEScan, MobileElementFinder, IntegronFinder) → co-localization by ARG–MGE co-occurrence on a contig [44, 46, 50].

Two methodological caveats are critical for any study using this approach:

1. **Assembly is the primary bottleneck.** Kerkvliet et al. (2024) benchmarked the MetaMobilePicker pipeline (metaSPAdes + PlasClass + ISEScan + DeepVirFinder + ABRicate/ResFinder) on a simulated seven-species metagenome and showed that the metagenomic assembly step—not the annotation tools—is the main cause of missed mobile ARGs: repetitive IS elements and multi-copy ARGs are collapsed during de Bruijn graph assembly, so mARGs and IS elements are frequently not recovered in correct copy numbers (IS sensitivity 0.585 overall but 0.932 among assembled elements; ARG sensitivity 0.691 overall but 0.982 among assembled genes) [44]. Plasmid prediction precision was moderate (0.57) and lowering the contig length cutoff from 1 kb to 500 bp increased false positives without improving recall [44]. The authors conclude long-read and hybrid sequencing are the most viable route to reliable MGE reconstruction, though plasmids remain challenging even with long reads [44].

2. **Co-localization results are approach-dependent and error-prone.** Slizovskiy et al. (2020) compared alignment-based and assembly-based resistome–mobilome colocalization on the same metagenomes and found the two yielded "contradictory and incomplete" results [46]. Major sources of error include (i) inconsistent identity/coverage thresholds across pipelines, (ii) substantial sequence homology between AMR and MGE database accessions (MGE databases contain "cargo" ARG sequences, inflating apparent co-localization), and (iii) short contig N50 limiting the sequence space needed to detect genuine co-localization [46]. They introduced an overlap-adjusted Mobility Index (MI) showing raw vs adjusted MI differed by 1–2 orders of magnitude, and urged investigators to triangulate multiple approaches and integrate database-validation tools (DeepARG, ARGDIT) [46].

Complementary advances include MobileElementFinder, a web tool for detecting MGEs and their insertion sites in bacterial genomes (developed for isolates, not metagenomes) [43]; Sarand, which extracts ARG neighborhoods directly from metagenomic assembly graphs with higher precision/sensitivity than contig-only methods [49]; and hybrid assembly as the most comprehensive strategy for ARG–MGE association in isolates [45]. A recent framework paper proposes an evidence ladder for moving from resistome detection to mechanistic mobility inference [48]. The practical implication: short-read metagenomic co-localization results should be reported with explicit assembly-quality metrics (N50, contig lengths) and interpreted as *putative* mobility, not confirmed transfer.

### 2.3 Last-resort resistance determinants (mcr, tet(X)) in wastewater and environmental spread

**Colistin resistance (mcr).** A massive analysis of 64,628 bacterial genomes identified 6,651 significant matches to nine MCR variants across 39 genera and >1,050 species, with MCR-1 most abundant in *Escherichia coli* (862 hits) and many variants found in environmental/water-associated bacteria—evidence of a broad water reservoir and possible non-antibiotic ecological roles for these enzymes [32]. mcr-1 is now routinely detected in wastewater globally, including a conjugative hybrid plasmid co-harboring mcr-1.1, mcr-3.39, and tet(X4) in *E. coli* from Seoul wastewater—a direct co-occurrence of colistin and tigecycline resistance on a transferable plasmid [24]. Hospital wastewater is a recognized hotspot for last-resort ARGs [26].

**Tigecycline resistance (tet(X)/tet(X4)).** The plasmid-mediated tet(X4) gene has emerged as a serious threat, with documented environmental spread. Mohsin et al. reported the emergence of plasmid-mediated tet(X4) in *E. coli* from poultry, food, and the environment in South Asia (Pakistan), highlighting regional surveillance needs [25]. tet(X4) has since been found in wastewater in Turkey (ST609 *E. coli*) [27], Norway (first report, ST167, co-occurring with blaCTX-M-14 and qnrS1) [30], hospital sewage in China (on transferable IncC/IncF/IncHI1 hybrid plasmids, alongside tmexCD-toprJ1) [29], and swine wastewater in China (stable, conjugative IncFII-like plasmid) [31]. These studies collectively show tet(X4) is environmentally mobile, often co-located with other ARGs on hybrid plasmids, and present in wastewater across both high- and low-income settings.

### 2.4 Bangladesh / South Asia water AMR evidence

Bangladesh-specific evidence is growing but remains methodologically skewed toward culture/PCR and toward hospital, surface, and slaughterhouse waters rather than industrial effluents.

- **Metagenomic studies (scarce, non-industrial):** McInnes et al. (2021) performed shotgun metagenomics on 24 rural and urban surface-water/sediment sites across Mymensingh, Shariatpur, and Dhaka, identifying 114 ARGs across 16 classes with ARG abundance varying 1,525-fold between sites and highest in urban surface water; ARG abundance correlated strongly (R² = 0.73) with human gut-bacteria contribution, implicating untreated sewage as the urban driver [11]. Notably, urban ARGs were more often plasmid-associated (predicted via PlasFlow), indicating higher HGT potential [11]. This is the most comprehensive Bangladeshi metagenomic resistome study to date, but it targeted surface waters/sediments, not industrial effluents, and did not perform ARG–MGE co-localization on assembled contigs.
- **Culture/PCR studies (dominant):** Sharif et al. (2025) screened 155 enteric isolates from 50 water samples across 24 sites in Savar (including industrial effluent, hospital, municipal, river, canal, lake, pond, irrigation) by culture and PCR, detecting blaTEM, blaCTX-M-15, blaOXA, qnrA/qnrB, tetA, sul-1, and mcr-1, with β-lactam resistance genotype prevalence of 30–37%, colistin (mcr-1) 3–10%, and significant genotype–phenotype associations; isolates from industrial effluent clustered closely with clinical and other environmental isolates, indicating cross-source circulation [9]. Mou et al. (2023) found ~69% of *E. coli* from surface water adjacent to Bangladeshi pharmaceutical industries were ESBL-positive (blaCTX-M group 60.8%, blaCTX-M-15 54.9%, blaTEM 70.6%), with 100% multidrug resistance [7]. Manik et al. (2023) reported blaCTX-M-1 group in 71.7% and blaTEM in 76.1% of hospital wastewater isolates [10]. Mannan et al. (2024) found 89/184 β-lactamase-positive hospital wastewater isolates were ESBL producers, with blaNDM-1 the most prevalent carbapenemase [12]. Chowdhury et al. (2025) characterized antimicrobial-resistant, biofilm-producing *E. coli* in Dhaka slaughterhouse wastewater [3]. Asaduzzaman et al. (2022) showed ESBL-*E. coli* ubiquitous in urban/rural Bangladesh, with colistin-resistant *E. coli* more prevalent in urban wastewater and rural rivers [6]. Siddique et al. (2024) detected multiple ARB in Dhaka municipal wastewater [8].
- **Regional context:** South Asia more broadly is an established hotspot for ESBL (especially CTX-M-15) and for emerging tet(X4) [25]. Tunisian industrial/biological WWTP studies detected ESBL-producing Enterobacteriaceae with mcr-1 in industrial effluent [4] and CTX-M-15 dominance in wastewater/river isolates [5]. South African wastewater MDR *E. coli* carried ARGs on mobile support [2].

The consistent pattern: Bangladesh has strong culture/PCR evidence of ESBL, mcr-1, tet, and class-1 integron circulation across water sources including industrial effluent [9], and one high-quality surface-water metagenomic study [11]—but **no published shotgun metagenomic study focused specifically on Bangladeshi industrial effluents with ARG–MGE co-localization analysis.**

---

## 3. Gap Analysis

| # | Gap theme | Current evidence | What is missing | Relevance to this study |
|---|-----------|------------------|-----------------|--------------------------|
| G1 | **Industrial-effluent shotgun metagenomics in Bangladesh** | Global industrial effluent resistomes characterized (pharma [34,40], tannery [39,42], textile [41]); BD metagenomics limited to surface waters/sediments [11] | No BD study performs shotgun metagenomics directly on industrial effluents; most BD work is culture/PCR [7,9,10,12] | IE1–IE4 directly address this: first BD industrial-effluent shotgun resistome catalog |
| G2 | **ARG–MGE co-localization on contigs in BD industrial effluents** | Co-localization done globally [15,17,18,44]; BD metagenomic study used PlasFlow plasmid prediction only, no contig-level ARG–MGE linkage [11] | No BD industrial-effluent contig-level co-localization of ESBL/mcr/tet(X) with IS/transposons/integrons/plasmids | Core objective 4; PlasClass + MobileElementFinder + IntegronFinder pipeline fills this |
| G3 | **Integrated last-resort gene surveillance (mcr + tet(X)) in BD industrial settings** | mcr-1 and tetA detected by PCR in BD industrial effluent [9]; tet(X4) emerging in South Asia [25] but not surveyed metagenomically in BD industry | No BD metagenomic quantification of mcr variants and tet(X) family in industrial effluents; no co-occurrence with ESBLs | Shotgun reads enable variant-level resolution (mcr-1..10, tet(X/X4)) and co-occurrence mapping |
| G4 | **Methodological rigor in short-read co-localization** | Assembly is the bottleneck [44]; co-localization is approach-dependent and error-prone [46] | Many studies report co-localization without assembly-quality metrics or overlap-adjustment | This study can model best practice: report N50/contig lengths, overlap-adjusted MI, triangulate tools |
| G5 | **Geographic/evidence skew** | Systematic review notes sparse data in less-industrialized/rural regions; heterogeneity limits synthesis [41] | LMIC industrial effluents under-represented globally | Bangladesh industrial effluents add a high-burden LMIC data point |
| G6 | **Virulence–resistance co-occurrence in effluents** | BD slaughterhouse/clinical work links AMR + virulence [3]; rarely integrated with resistome context | Limited metagenomic co-occurrence of ARGs and virulence markers on industrial effluent contigs | Virulence screening (VFDB/ABRicate) on ARG-bearing contigs addresses this |

---

## 4. Contribution Framing

### 4.1 How IE1–IE4 address the gaps

The four shotgun metagenomic samples (IE1–IE4) from Bangladeshi industrial effluents position this study to make several distinct contributions:

1. **First BD industrial-effluent shotgun resistome catalog (G1, G5).** Moving beyond the culture/PCR evidence [7,9] and the surface-water metagenomics of McInnes et al. [11], this study will provide a direct, variant-level metagenomic profile of ARGs in industrial effluent—quantifying ESBL variants (blaCTX-M families, blaTEM, blaSHV), MBLs, mcr variants, and tet(X) family genes with normalized abundance (RPKM/TPM, copies per 16S). This adds a high-burden LMIC industrial data point to a globally skewed evidence base [41].

2. **Contig-level ARG–MGE co-localization in BD industrial effluents (G2).** Using metagenomic assembly followed by PlasClass (plasmid probability), MobileElementFinder/ISEScan (IS and transposons), and IntegronFinder (class-1 integrons), the study will quantify how often ESBL/mcr/tet(X) genes co-localize with mobilizable elements on the same contig—a direct test of horizontal-transfer potential that no BD industrial-effluent study has performed.

3. **Integrated last-resort surveillance (G3).** Shotgun data enable simultaneous, variant-resolved detection of mcr (colistin) and tet(X)/tet(X4) (tigecycline) alongside ESBLs—testing whether last-resort determinants co-occur with each other and with MGEs in industrial effluents, paralleling the hybrid-plasmid co-occurrence reported elsewhere [24,29].

4. **Methodological best-practice demonstration (G4).** By reporting assembly-quality metrics (N50, contig length distribution), applying an overlap-adjusted Mobility Index [46], and triangulating multiple MGE predictors, the study can model rigorous short-read co-localization interpretation—explicitly framing results as *putative* mobility given the known assembly bottleneck [44].

5. **Virulence–resistance co-occurrence (G6).** Screening ARG-bearing contigs against VFDB adds a virulence dimension rarely integrated in effluent resistome studies.

### 4.2 Novelty statement

This will be among the first studies to combine shotgun metagenomic resistome profiling with contig-level ARG–MGE co-localization analysis focused specifically on industrial effluents in Bangladesh, integrating ESBL, mcr, and tet(X) surveillance under a unified mobility-risk framework. It bridges a gap between Bangladesh's extensive culture/PCR AMR evidence [7,9,10,12] and the global metagenomic industrial-effluent literature [34,36,39,41], and applies current best-practice caveats for short-read co-localization [44,46].

### 4.3 Limitations to acknowledge explicitly

- **Sample size (n = 4):** suitable for a prevalence/cataloguing and co-localization-feasibility study, not for robust between-site statistical comparison or temporal dynamics. Framing should be descriptive and hypothesis-generating.
- **No sample-level metadata:** industry type, location, season, physicochemical parameters, and antibiotic concentrations are not in the uploaded data; contribution framing stays at cohort/study level. Obtaining metadata would materially strengthen interpretation (see follow-up).
- **Short-read only:** co-localization is bounded by assembly contiguity; ARG–MGE linkage will be under-resolved relative to long-read/hybrid approaches [15,17,44]. Results must be reported as putative mobility with assembly-quality context.
- **No culture/phenotypic data uploaded:** the PCR and phenotypic-resistance objectives cannot be addressed from the uploaded data alone; the contribution centers on the metagenomic arm.

---

## 5. Concrete Metagenomic Analysis Plan (IE1–IE4)

### 5.1 Pipeline overview

```
Raw FASTQ (IE1–IE4, PE)
  → QC + trimming + host removal
  → Taxonomic profiling
  → Metagenomic assembly (per-sample + co-assembly)
  → ARG annotation (variant-level)
  → MGE / integron / plasmid annotation
  → ARG–MGE co-localization on contigs
  → Virulence co-occurrence
  → (Optional) MAGs + host-ARG linkage
  → Abundance tables, co-localization table, figures, report
```

### 5.2 Steps, tools, and parameters

**Step 1 — QC and preprocessing**
- `FastQC` + `MultiQC`: per-sample read quality, adapter content, duplication.
- `fastp` (or Trimmomatic SLIDINGWINDOW:4:15 MINLEN:50): adapter/quality trimming.
- Host/human read removal: `Bowtie2` against GRCh38 (or `KneadData`); retain non-host reads for resistome analysis and separately account for human reads for ethics/privacy. This is important because effluent metagenomes can contain human-associated DNA.

**Step 2 — Taxonomic profiling**
- `Kraken2` + `Bracken` (custom DB including bacteria, archaea, viral, human) OR `MetaPhlAn` for composition. Report phylum/genus-level community structure to contextualize resistome (cf. Proteobacteria/Cyanobacteria dominance patterns [11]).

**Step 3 — Metagenomic assembly**
- Per-sample assembly: `MEGAHIT` (default, meta mode) and/or `metaSPAdes` (meta mode). Given n = 4 and likely moderate complexity, also run a **co-assembly** of all four samples to improve contiguity for low-abundance ARGs/MGEs.
- Report N50, number of contigs >1 kb, max contig length, total assembled bp. Use ≥1 kb contig cutoff for co-localization (lowering to 500 bp increases false positives [44]).
- Pilot on IE1 first to estimate runtime/memory before scaling.

**Step 4 — ARG annotation (variant-level)**
- Primary: `AMRFinderPlus` (NCBI) and `RGI` (CARD) on assembled contigs; cross-validate with `ABRicate` (ResFinder). Optionally `DeepARG` for divergent/novel genes [46].
- Target gene groups for explicit reporting: blaCTX-M (family + variants e.g. CTX-M-15, -55, -14), blaTEM, blaSHV, MBLs (blaNDM, blaVIM, blaIMP), mcr-1..10, tet(X)/tet(X4), plus class-1 integron-associated cassettes.
- Abundance normalization: map reads back to contigs (`Bowtie2` + `samtools`), compute RPKM/TPM per ARG and copies per 16S rRNA gene (normalize to bacterial load).

**Step 5 — MGE / integron / plasmid annotation**
- Plasmid prediction: `PlasClass` (per-contig plasmid probability) + `PlasFlow`; `PlasmidFinder` via ABRicate for replicon types.
- IS/transposons: `ISEScan` (IS) and `MobileElementFinder` (composite transposons, insertion sites—note: isolate-oriented, use on contigs with caution) [43].
- Integrons: `IntegronFinder` (class-1, -2, -3; complete vs CALIN vs In0).
- Phages (optional): `DeepVirFinder` or `VirSorter2`.

**Step 6 — ARG–MGE co-localization**
- Define co-localization: ARG and MGE/integron/plasmid feature on the **same contig**, optionally within a distance window (e.g., ≤5 kb) to capture transposon-scale linkage.
- Compute a **Mobility Index (MI)**: fraction of ARG-bearing contigs that also carry ≥1 MGE feature, per sample [46]. Apply **overlap-adjusted MI** to correct for cargo-sequence homology between AMR and MGE databases [46].
- Output a co-localization table: contig ID | length | PlasClass score | ARG(s) | MGE type(s) | integron | distance | MI category. Annotate ESBL/mcr/tet(X) co-localization events specifically.
- Report assembly-quality context alongside every co-localization claim (per G4 best practice).

**Step 7 — Virulence co-occurrence**
- `ABRicate` against `VFDB` on ARG-bearing contigs; flag contigs carrying both ARG and virulence factors.

**Step 8 — (Optional) MAGs and host-ARG linkage**
- If coverage permits: bin contigs (`MetaBAT2`/`CONCOCT`), assess with `CheckM2`, assign taxonomy (GTDB-Tk), and link ARGs to MAG hosts—strengthening host-context risk interpretation [18,28].

**Step 9 — Outputs and figures**
- Per-sample ARG abundance table (normalized); MGE/integron summary; co-localization table; MI table.
- Figures: heatmap of ARG classes by sample; bar chart of ESBL/mcr/tet(X) prevalence; co-localization network (ARG–MGE bipartite) or bar of MI per sample; community composition stacked bar.
- Summary report integrating findings with the gap analysis.

### 5.3 Compute/resource estimate

- Input: ~21 GB compressed FASTQ → ~60–90 GB decompressed; per-sample assembly of metagenomes of this size typically needs 30–60 GB RAM for metaSPAdes (MEGAHIT more memory-frugal, ~16–32 GB); co-assembly of all four may exceed 64 GB and should be run on a right-sized machine or HPC.
- Plan: pilot IE1 on a 16 CPU / 64 GB sandbox machine to measure assembly runtime/memory and ARG annotation throughput; scale to remaining samples and co-assembly accordingly, using HPC tools (MEGAHIT/SPAdes available) if co-assembly exceeds sandbox memory.
- Disk: allocate ≥200 GB working space in `/workspace/` for intermediates; checkpoint final tables/figures to `/mnt/results/`.

### 5.4 Quality criteria / acceptance

- All four samples pass QC with acceptable read metrics (post-trim Q30 ≥ 80%, host fraction reported).
- ARG annotation reported with ≥2 concordant tools for key genes (ESBL, mcr, tet(X)).
- Every co-localization claim accompanied by contig length and PlasClass/MGE evidence; MI reported raw and overlap-adjusted.
- Limitations section explicitly states short-read co-localization is putative.

---

## 6. References

Inline numbered citations are used throughout (e.g., [9], [11], [34], [41], [44], [46]); the full reference list is rendered by the platform from the session reference registry. Key sources grounding this review include: McInnes et al. 2021 (mSystems) [11]; Sharif et al. 2025 (npj Clean Water) [9]; Ibitoye et al. 2025 (F1000Research) [41]; Kerkvliet et al. 2024 (PeerJ) [44]; Slizovskiy et al. 2020 (Front Microbiol) [46]; Bengtsson-Palme et al. 2019 (Water Res) [34]; Johansson et al. 2020 (JAC) [43]; Dai et al. 2022 (Microbiome) [15]; Che et al. 2019 (Microbiome) [17]; Mohsin et al. 2021 (STOTEN) [25]; Ben Khedher et al. 2020 (Sci Rep) [32]; Wang et al. 2013 (PLoS ONE) [39]; Su et al. 2023 (Environ Int) [36]; Mou et al. 2023 [7]; Manik et al. 2023 [10]; Mannan et al. 2024 [12]; Asaduzzaman et al. 2022 [6]; and the additional sources cited above.
