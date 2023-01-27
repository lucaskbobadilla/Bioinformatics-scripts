#!/bin/bash


# this script will convert bam files into fastq files

echo Program to convert bam files into fastq files
echo by Lucas K Bobadilla

# get variables

echo "Do you have pair or single read fastq files mapped to the bam? (answer pair or single)"
read read_type

bam_input_dir=$1

fastq_out=$2


# run program

cd $bam_input_dir

for bam_file in *.bam; do
    if samtools view -H $bam_file | head -n 1 | grep -q "queryname"
    then
        echo "${bam_file} already sorted by name. Converting to fastq: "
        if [ $read_type == "pair"  ]
        then
            bedtools bamtofastq -i $bam_file -fq ${fastq_out}/${bam_file::-4}_R1.fastq -fq2 ${fastq_out}/${bam_file::-4}_R2.fastq
            pigz -p 25 ${fastq_out}/${bam_file::-4}_R1.fastq
            pigz -p 25 ${fastq_out}/${bam_file::-4}_R2.fastq
        else
            bedtools bamtofastq -i $bam_file -fq ${fastq_out}/${bam_file::-4}.fastq
            pigz -p 25 ${fastq_out}/${bam_file::-4}.fastq
        fi
        echo "${bam_file::-4} converted"
    else
        echo "${bam_file} not sorted by name. Sorting:"
