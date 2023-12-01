# Collinearity search

These scripts prepare the input for running i-ADHoRe v3.0.
`iadhore_preparation.sh` is a wrapper script of the other preparation scripts.

i-ADHoRe v3.0 (Proost et al., 2012) was executed to identify collinear blocks, i.e. regions with-in the genome with similar gene content and order. As input for i-ADHoRe, gene-pseudogene pairs were extracted from the DUP and FRAG classes of PseudoPipe, and gene-gene pairs were extracted from an all-vs-all BLASTP run (E-value ≤ 1e-5, ≥ 30% alignment coverage, and length difference ≤ length smallest protein) (see `run_blastp.sh`). 

Parameters that were used for i-ADHoRe were: (1) A minimum number of 3 anchorpairs, i.e. pairs of putative homologous (pseudo)genes in a collinear block, (2) a maximum gap size that should exist between two points of 35, (3) a max-imum cluster gap size of 40, (4) maximum number of gaps in alignment of 40, (5) a probability cut-off of 0.01, and (6) a q-value cut-off of 0.75 (see `settingsiADHoRe.txt`). 

Collinear blocks were visually inspected using GenoPlotR (Guy et al., 2010) (see `visualization` folder).
