# Python-Bioinformatics-scripts

This repository contains a sort of Python scripts for daily bioinformatics tasks.

## N50 fasta statistics

The script `N50_stats.py` will give the stats of a draft genome:

### Usage

`N50_stats.py −f arrow/cgun_LG12.arrow.fa`

### Output

`###### N50 stats #####`

`Assembly file name: pilon/cgun_LG12.pilon.fasta`

`Assembly size (bp): 22,504,182`

`Number of sequences: 8`

`Largest sequence: 9,135,553`

`N50 (bp): 6,574,633`

`L50: 2`

## Genomic circle plot with PyCairo

The script `Circle_plot_pyCairo` will create a genomic Circle plot using two tsv files:

* A expression fold change data with genomic location of each gene
* A Fst tsv file for each SNP with the genomic coordinates

### Usage

`./plot.py −e Gacu_FoldChange_GenomCoords.tsv −f batch_1.phistats_fw2−oc.tsv`

where `-e` refers to the expression data and `-f` refers to the Fst data

### Output
![image](https://user-images.githubusercontent.com/32884929/114288508-22f00980-9a36-11eb-9d71-90b346053c71.png)

