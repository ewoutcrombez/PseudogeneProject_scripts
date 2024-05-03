# Overview pseudogene/remnant identification

All scripts can be found in the `autopol_pseudogene_scripts` folder. The pipeline consists of three major steps:

0. **Run i-ADHoRe** (`genelist_iadhore_creation.sh` and `settingsfile_iadhore_creation.sh` to get input files for i-AHoRe)

1. **Processing of the i-ADHoRe output** (`process_i-ADHoRe_multiplicon.sh`)

   In this first step, i-ADHoRe output files are parsed to get specific files that are required for visualization. In addition, genes that did not have a counter part in the collinear block (thus a "lonely" gene, no homologous gene on the other region of the collinear block) were confirmed by searching for tandem repeat genes, or appear to still have a match in the other region of the collinear block or on the same or different chromosomes (between the same haplotypes). The lonely genes that still have a match are removed from the "lonely genes" list. The remaining lonely genes may still have a counterpart in the other collinear segment, and thus have a "pseudogene", or at least a counterpart in the non-coding region. Finally, the genic regions of both regions of the collinear block are masked to (in the next step) search for hits of the lonely genes in the non-coding regions.

   - `obtain_i-ADHoRe_chrom_for_visualization_between_species.py`: Create a file with the required columns for creating a plot (only gene-gene pairs).
     - input: gene list (also used as input for i-ADHoRe) of the chromosome and haplotype of the selected multiplicon.
     - output: file ending with tab.txt containing all genes of the specific chromosome and species. In `process_i-ADHoRe_multiplicon.sh` the specific region of the multiplicon is then selected.
   - `find_tandem_or_matched.R`: Identify tandem repeats and lonely genes that still match with a gene in the region, chromosome, or haplotype. Matches with other genes is based on BLAST results (60% percent identity and alignment length >= 0.7xquery length) and tandem repeats are parsed from i-ADHoRe output file. Finally color genes according to identity (i.e. lonely gene, tandem repeat, matches between lonely genes thus not really lonely, matches between lonely gene and other chromosome, matches with gene of different chromosome). These are always between the two haplotypes of that multiplicon.
     - input: list of "lonely" genes, blast result
     - output: many files of different gene identities
   - `obtain_bed_for_regions.R`: extract the regions that are part of a collinear region and mask the genic regions in this region using the right coordinate system

2. **Searching for "pseudogenes"**, or at least hits/remnants in the non-coding region that match with lonely genes (`find_remnants_segment1.sh` and `find_remnants_segment2.sh`)

   In this step, for all lonely genes a hit is searched in the non-coding region of the other region of the collinear block using exonerate. The best hits for each gene are then selected based on alignment coverage and percent identity. Additionally, hits that overlap with each other should get the same hit name (otherwise during the plotting it will visualize multiple pseudogenes while it is actually the same one that matches multiple lonely genes). Hits that are within 5000 bp and match to the same lonely gene and have the same orientation are merged. After these initial merging and filtering steps, exonerate is run again but using a different mode whereby the "CDS" sequences of the pseudogenes are aligned against the CDS sequences of the genes. Based on the results of this second alignment we filter a final time based on query alignment length and percent identity. Finally, we can then also align the putative gene-pseudogene pairs using another alignment tool FSA to then also take a look at the specific alignments.

   - `process_exonerate_result.py`: parse exonerate result file to obtain GFF and result table separately.
   - `combine_and_get_best_hits.R`: get a table whereby hits that match to the same region of a gene are indicated. For these we will only select one (the one with the highest percent identity and coverage).
   - `process_gff.py`: hits that match to the same gene and are within 5000 bp are merged. Here, new hit names are also given, thus hit names of files before do not correspond with the new hit names! A hit will have the same name if it overlaps with another hit and has the same orientation. It does not matter what gene the hit is from.
     - output: GFF file with exons
   - `get_hit_table.R`: Get table with summary stats for merged hits 
     - input: GFF file with exons from `process_gff.py`
     - output: summary stats hit table
   - `combine_exons.py`: join fasta sequences if they are from the same hit
   - `filter_hits.R`: filter hits a final time by percent identity (> 20%) and alignment length (> 20% gene length)

3) **Visualizing the multiplicon** with all identified pairs shown on it (`visualize_with_remnants.sh`)

   In this step, I visualize the multiplicons with all gene-gene, gene-tandem gene, gene-pseudogene pairs visualized...
   
   - `add_orientation_gene.R`: get the orientation of the gene and add it to BED file
   - `bed_to_table.py`: convert BED file with both genes and pseudogenes to a table that can be used as input for GenoPlotR.
   - `visualize_i-ADHoRe_multiplicon.R`: Visualize multiplicons. If "YES" is given this makes a plot with pseudogenes