#!/bin/bash

#### THIS SCRIPT WILL RUN STAR ALIGNMENT TOOL ###

### GENERATE BAM FILES FROM FASTQ FILES AFTER ALIGNED TO A GENOME ####

### GENOME NEED TO BE INDEX FIRST, FOR WATERHEMP AND PALMER INDEXES ARE AVAILABLE AT ~/genomes #####


# ARGUMENTS

## FIRST ARGUMENTS - FILE PATH
file_path=$1
## SECOND ARGUMENT - OUTPUT PATH
out_path=$2
## THIRD ARGUMENT - NUMBER OF THREADS
threads=$3
## FOURTH ARGUMENT - INDEX PATH
genome_index=$4

### RUN SCRIPT ####

for file in $(ls ${file_path}/*.fastq.gz | rev | cut -c 10- | rev | uniq); do \
    echo "############ Mapping ${file} ################"; \
    STAR --runThreadN $threads \
         --genomeDir $genome_index \
         --readFilesCommand zcat \
         --readFilesIn "${file_path}"/"${file}".fastq.gz \
         --outSAMtype BAM SortedByCoordinate \
         --alignIntronMax 50000 \
         --outFileNamePrefix "${out_path}"/"${file}". \
         --outTmpDir ./
done

echo "######## Finished STAR mapping #############";

### CREATE INDEX FOR ALL BAM FILES ####

echo "############## Creating bam files indexes ###########";
parallel samtools index ::: "${out_path}"/*.bam
