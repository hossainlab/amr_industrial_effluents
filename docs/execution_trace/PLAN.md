# Execution Plan: Shotgun Metagenomic ARG/MGE Analysis of IE1–IE4

## Summary
Run the approved metagenomic pipeline on the 4 industrial-effluent paired-end samples (IE1–IE4, Illumina 150 bp, ~21 GB compressed) to produce a variant-level resistome catalog, MGE/integron/plasmid annotation, contig-level ARG–MGE co-localization with an overlap-adjusted Mobility Index, virulence co-occurrence, and figures + final report. The pipeline is metadata-independent; the correct sample metadata (pending from user) will be integrated for labeling/interpretation once provided. The wrong metadata file (`Metagenomics sample information.docx`, farm/hospital samples) is excluded.

## Data confirmed
- FASTQs: IE1–IE4, R1/R2, Illumina 150 bp PE (instrument LH00129). ~2.3–3 GB/file, ~21 GB total compressed.
- Metadata: pending. Current uploaded docx is from a different project (CW/CH/PW/PH/HW farm/hospital samples) — not used.

## Execution target & resources
- Provision a right-sized sandbox machine via ManageMachine: 16 CPU / 64 GB RAM, ample local disk (/workspace). Default worker-0 is minimal-footprint and insufficient for assembly.
- Workload: metagenomic assembly (MEGAHIT, memory-frugal ~16–32 GB per sample; metaSPAdes alternative needs more). Co-assembly of all 4 may approach/exceed 64 GB → fallback to HPC MEGAHIT tool if OOM.
- Strategy: **pilot IE1 end-to-end through assembly** to measure runtime/memory/throughput, then scale to IE2–IE4 and co-assembly. Chunk with checkpoints to /mnt/shared-workspace; final deliverables to /mnt/results.
- Disk: stage decompressed reads and intermediates in /workspace (≥200 GB needed — verify available); checkpoint assemblies + annotation tables to /mnt/shared-workspace.

## Pipeline steps (behavior-level)
1. **Stage & QC**: extract FASTQs from zip to /workspace; FastQC + MultiQC per sample; fastp trim (adapter auto-detect, Q15 sliding window, min len 50). Report read counts pre/post, %Q30, adapter content.
2. **Host/human removal**: Bowtie2 vs GRCh38 (or KneadData); separate human reads (ethics/privacy), retain non-host for resistome.
3. **Taxonomic profiling**: Kraken2 + Bracken (standard+minusbacteria DB or available prebuilt) on non-host reads; report phylum/genus composition.
4. **Assembly**: per-sample MEGAHIT (meta mode, --min-contig-len 1000); plus a co-assembly of all 4. Report N50, # contigs ≥1 kb, max/total length. Use ≥1 kb contigs for co-localization (per Kerkvliet 2024: 500 bp cutoff raises false positives).
5. **ARG annotation (variant-level)**: AMRFinderPlus + RGI/CARD on contigs; cross-validate key genes with ABRicate/ResFinder. Explicitly report blaCTX-M families/variants, blaTEM, blaSHV, MBLs (blaNDM/VIM/IMP), mcr-1..10, tet(X)/tet(X4), class-1 integron cassettes. Map reads back (Bowtie2) → RPKM/TPM and copies-per-16S normalization.
6. **MGE/integron/plasmid**: PlasClass + PlasFlow (plasmid probability), ISEScan (IS), MobileElementFinder (transposons/insertion sites, contig-level with caution), IntegronFinder (class-1/2/3), PlasmidFinder (replicons via ABRicate).
7. **ARG–MGE co-localization**: ARG + MGE/integron/plasmid on same contig (and ≤5 kb window variant). Compute Mobility Index = ARG-bearing contigs with ≥1 MGE / all ARG-bearing contigs, per sample; report **raw and overlap-adjusted MI** (corrects cargo-sequence homology, per Slizovskiy 2020). Output table: contig | length | PlasClass score | ARG(s) | MGE type(s) | integron | distance | MI. Flag ESBL/mcr/tet(X) co-localizations.
8. **Virulence co-occurrence**: ABRicate vs VFDB on ARG-bearing contigs.
9. **(Optional) MAGs**: if coverage permits — MetaBAT2/CONCOCT binning, CheckM2, GTDB-Tk taxonomy, host-ARG linkage.
10. **Figures**: heatmap of ARG classes by sample; bar of ESBL/mcr/tet(X) normalized abundance; co-localization network (ARG–MGE bipartite) + MI bar; community composition stacked bar. Save .png + .svg; run media output check on each.
11. **Final report**: markdown `/mnt/results/report_amr_industrial_effluents_analysis.md` with methods, per-sample results, co-localization findings (with assembly-quality context and putative-mobility caveat), figures, and limitations (n=4, short-read only, metadata pending).

## Deliverables (in /mnt/results)
- `arg_abundance_table.tsv` (per-sample, normalized), `mge_integron_summary.tsv`, `colocalization_table.tsv`, `mobility_index.tsv`, `taxonomy_summary.tsv`.
- Figures in `figures/` (png+svg).
- `report_amr_industrial_effluents_analysis.md`.
- Checkpointed assemblies/annotations in /mnt/shared-workspace.

## Assumptions & decisions
- Assembler: MEGAHIT primary (memory-frugal); metaSPAdes only if higher contiguity needed and memory allows.
- Co-localization reported as **putative** mobility with N50/contig-length context (short-read limitation, per Kerkvliet 2024 & Slizovskiy 2020).
- n=4 → descriptive/cataloguing study, no between-group statistical testing unless metadata defines groups.
- Begin metadata-independent steps (1–8) now; integrate metadata for labels/grouping when user provides correct file.
- Wrong metadata docx excluded from analysis.
- Follow-up questions block included in final answer per house style.
