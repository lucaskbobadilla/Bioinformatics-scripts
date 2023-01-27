#!/bin/bash

#
## THIS SCRIPT RUNS KALLISTO ON MULTIPLE FILES WITHIN A DIRECTORY
#

# ARGUMENTS

kallisto_index=$1 # GIVE THE FULL PATHWAY FOR IT

fastq_folder=$2

out_dir=$3


# MOVE TO FOLDER
cd $fastq_folder

#RUN KALLISTO
for sample in $(ls *.fastq.gz | rev | cut -c 13- | rev | uniq); do \
    echo "Running Kallisto for sample" ${sample}
    kallisto quant \
             -i $kallisto_index \
             -o ${out_dir}/${sample} \
             -b 100 \
             -t 30 \
             ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz;
done
