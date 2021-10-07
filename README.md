## Aim

This repository contains the scripts and code used in the DISCERNS paper. It is based on a snakemake pipeline that includes all code to simulate an RNA-seq dataset and to generate GTF files with reduced annotation. Besides, the workflow compares STAR and hisat2 using the simulated data. It also predicts novel exons with DISCERNS and StringTie and generates Precision-Recall curves.

Forked from https://github.com/khembach/novel_exon_pipeline (1st Sept 2021)

## Design changes

Projected/implemented

- i (1 Sept): get rid of conda envs, just a virtenv and base R (4.1.0). Reason why: envs contains multiple R3.6 dependencies, which are tedious to update. Rather, use a virtenv with the most uptodate soft, then generate the env files (if any)
- i 6 Oct: go back to conda, but with mamba
- p (.): no downsampling to chr19 and 22, whole genome instead.
- i 6 Oct: Most recent GTF and reference FASTA, gencode
- i 6 Oct: add data download steps
- i (.) : make paths relative to the working directory WD, update 6 Oct: partial, to deal with hardoded paths
- p (.): estimate theta for RSEM using real data?
- p (.), remove unused soft stacks, e.g. is tophat used anyway? this goes for featurecounts, tophat, salmon etc

## Run

0. install mini- or anaconda
1. install deps using `soft_installs.sh`
2. `snakemake -s Snakefile -j NUMCORES -n -p` (but without the `-n -p`), potentially with --use-conda
