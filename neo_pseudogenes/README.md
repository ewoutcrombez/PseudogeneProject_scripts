# Detection of post-polyploid pseudogenization and translocation in neo-polyploid species

## i-ADHoRe run

i-ADHoRe was run using the same parameters as for the paleopolyploids:

```
cluster_type=colinear
prob_cutoff=0.01
multiple_hypothesis_correction=FDR
gap_size=35
cluster_gap=40
q_value=0.75
anchor_points=3
alignment_method=gg2
max_gaps_in_alignment=40
output_path=output
number_of_threads=4
verbose_output=true
level_2_only=true
```

`genelist_iadhore_creation.sh` and `settingsfile_iadhore_creation.sh` to get input files for i-AHoRe
`run_blastp.sh` to run blastp

## Pipeline requirements

To run the different steps of the pipeline, a folder is required which has a data and results folder.

The data folder contains:
- GFF folder containing GFF files for each sub-genome (named as ${subgenome}.protein-coding.gene.gff3)
- CDS folder containing CDS fasta files (named as ${subgenome}.CDS.fa)
- prot folder containing protein fasta files (named as ) 
- genome folder containing fasta files for each sub-genome (named as ${subgenome}_haplotype_genome.fasta)
- iadhore_input folder containing the input that was used to run i-ADHoRe, specifically it should contain a genelists folder which contains folders of the gene lists of the sub-genomes

The results folder contains:
- i-ADHoRe-run result folder containing the output of i-ADHoRe
- blastp folder containing the output of an all-vs-all BLASTP run (this BLASTP file may be already filtered based on cut-offs, see methods, named allvsallblastp_filtered_wlen.tsv.gz)

It is assumed that gene names follow the structure "${subgenome}-${name}". The subgenome names should just be one letter, e.g. A, B, C, D. Gene and protein names should also be the same.

## Collinearity retention analysis

### Confirm i-ADHoRe output and look for other homologous in collinear segment

See `process_i-ADHoRe_wrap.sh` (to run `process_i-ADHoRe_multiplicon.sh` for each collinear segment), `process_i-ADHoRe_multiplicon.sh`.
This runs the scripts:
- `obtain_multiplicon_gene_positions.py`: Get gene positions on chromosomes from i-ADHoRe gene lists
- `find_tandem_or_matched.R`: search for other homologous gene in the collinear segment (i.e. tandem or genes with high similarity)

### Collinear retention level assessment and missing homologue search
After running the above for all collinear segments, the results of all collinear segments are combined into one data frame and then collinear segments are grouped together across sub-genomes to assess whether there is collinearity retention across all sub-genomes or only for a subset of sub-genomes, i.e. infer collinear retention level over the sub-genomes. Then check for levels 3, 2 and 1 whether we can find translocated genes and/or collinear segments between the sub-genome where a homologue is missing.
See `collinear_retention_analysis.sh` and `collinear_retention_analysis.R` and `search_missing_homologues.py`. To assess whether genes are translocated we use Syri which infers translocated regions (see `syri_translocation_detection.sh`).

## Pseudogene search

For each collinear segment, search for pseudogenes based on the missing homologues (`search_pseudogenes.sh`; `obtain_bed_for_regions.R` to obtain BED file for the two regions of a collinear block).
Find hits/remnants in the non-coding region that match with lonely genes (`find_remnants_segment1.sh` and `find_remnants_segment2.sh`)

In this step, for all lonely genes a hit is searched in the non-coding region of the other region of the collinear block using exonerate. The best hits for each gene are then selected based on alignment coverage and percent identity. Additionally, hits that overlap with each other should get the same hit name (otherwise during the plotting it will visualize multiple pseudogenes while it is actually the same one that matches multiple lonely genes). Hits that are within 5000 bp and match to the same lonely gene and have the same orientation are merged. After these initial merging and filtering steps, exonerate is run again but using a different mode whereby the "CDS" sequences of the pseudogenes are aligned against the CDS sequences of the genes. Based on the results of this second alignment we filter a final time based on query alignment length and percent identity. Finally, we then align the putative gene-pseudogene pairs using MACSE to take a look at the specific alignments and summarize the alignments.

- `process_exonerate_result.py`: parse exonerate result file to obtain GFF and result table separately.
- `combine_and_get_best_hits.R`: get a table whereby hits that match to the same region of a gene are indicated. For these we will only select one (the one with the highest percent identity and coverage).
- `process_gff.py`: hits that match to the same gene and are within 5000 bp are merged. Here, new hit names are also given, thus hit names of files before do not correspond with the new hit names! A hit will have the same name if it overlaps with another hit and has the same orientation. It does not matter what gene the hit is from.
    - output: GFF file with exons
- `get_hit_table.R`: Get table with summary stats for merged hits
    - input: GFF file with exons from process_gff.py
    - output: summary stats hit table
- `combine_exons.py`: join fasta sequences if they are from the same hit
- `filter_hits.R`: filter hits a final time by percent identity (> 20%) and alignment length (> 20% gene length)
- `get_result_table_segment.R`: get exonerate result table
- `get_g_psg_alignments.py`: align pseudogene and gene with MACSE and summarize alignment statistics in MACSE result table

The results of the pseudogene search are then combined using the `get_final_pseudogene_df.sh` and `get_final_pseudogene_df.R` scripts.

## Get final result


