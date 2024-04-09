# Pseudogene identification in paleo-polyploid species

These scripts obtain the required data from PLAZA v5.0 (`S1_fetch_data.sh`), run PseudoPipe (`S2_run_PseudoPipe_array.sh`) and do some post-processing and filtering of the PseudoPipe output (`S3_postprocessing_PseudoPipe.sh` and `S4_filter_pseudogenes.R`).

`S1_fetch_data.sh`:
First, data was obtained from PLAZA v5.0 and the genomes of the paleo-polyploid species were repeat-masked using RepeatMasker v4.1.1 (Chen, 2004). 

`S2_run_PseudoPipe_array.sh`:
Putative pseudogenes were identified using the publicly available PseudoPipe pipeline (Zhang et al., 2006). Briefly, PseudoPipe attempts to identify pseudogenes by searching for homologous sequences of functional genes in the non-coding part of the genome. The proteome of a species is blasted (with tblastn) against the repeat- and gene-masked genome of that species. The hits are then further processed: redundant pseudogenes are eliminated, neighboring hits are merged, and a unique parent is determined for each pseudogene. Hits are then re-aligned using tfasty of the FASTA v36.3.8d suite (Pear-son, 2016). Furthermore, pseudogenes are classified as “processed” (PSSD; presence of poly-A tail, no introns, > 40% protein percent identity with functional gene, E-value < 1e-10 and ≥ 70% aligned to the functional gene), “duplicated” (DUP; no poly-A tail, multiple exons, > 40% percent identity with functional gene and E-value < 1e-10) and “fragmented” (FRAG; sequence similarity with functional gene, but below the cut-offs of DUP and PSSD pseudogenes).

`S3_postprocessing_PseudoPipe.sh`, `S4_filter_pseudogenes.R` and `S5_filter_overlapping_pseudogenes.R`: After running the pipeline, some further post-processing was performed: putative pseudogenes with more than 30 bps overlap with an exon were removed, pseudogenes were filtered (≥ 5% aligned to the functional gene, ≥ 20% protein percent identity and E-value ≤ 1e-5) and only the best pseudogene is retained if there is overlap.
