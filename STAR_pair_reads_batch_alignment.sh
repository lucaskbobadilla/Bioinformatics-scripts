#!/bin/bash

#### THIS SCRIPT WILL RUN STAR ALIGNMENT TOOL ###

### GENERATE BAM FILES FROM FASTQ FILES AFTER ALIGNED TO A GENOME ####

### GENOME NEED TO BE INDEX FIRST, FOR WATERHEMP AND PALMER INDEXES ARE AVAILABLE AT ~/genomes #####

### Test the ls command inside the for loop before running change the parameter number after cut -c according to
# your file names


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

cd ${file_path}



for file in $(ls *.fastq.gz | rev | cut -c 21- | rev | uniq); do \
    echo "############ Mapping ${file} ################"; \
    STAR --runThreadN $threads \
         --genomeDir $genome_index \
         --readFilesCommand zcat \
         --readFilesIn "${file}"_R1_trimmed.fastq.gz "${file}"_R2_trimmed.fastq.gz\
         --outSAMtype BAM SortedByCoordinate \
         --alignIntronMax 50000 \
         --outFileNamePrefix "${out_path}"/"${file}". \
         --outTmpDir "${out_path}"/temp_dir \
	 --limitBAMsortRAM 6643038087
done

echo "######## Finished STAR mapping #############";

### CREATE INDEX FOR ALL BAM FILES ####

echo "############## Creating bam files indexes ###########";
parallel samtools index ::: "${out_path}"/*.bam

echo "################ Move bam files to folder ############";

mkdir "${out_path}"/bam_files
mv "${out_path}"/*.bam "${out_path}"/bam_files
mv "${out_path}"/*.bai "${out_path}"/bam_files
