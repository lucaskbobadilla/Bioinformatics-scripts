#!/bin/python

### SAM PARSE PROGRAM - LUCAS KOPECKY BOBADILLA #####

import argparse

p = argparse.ArgumentParser()

#add options

p.add_argument('--input', '-i', required=True) # file path

args = p.parse_args()


input_file = args.input  # open the fastq file

sam = open(input_file,'r') # open user specified file  

# count the number of unmapped reads
unmapped =  4 # flag for unmapped
secondary = 256 # flag for secondary flag
supp = 2048 # flag for supplementary alignment
counter_unmapped = 0 # conter for unmapped
counter_primary = 0 # Counter for primary reads
counter_record = 0 # counter all reads
counter_supp = 0
for line in sam: # for each line in the  file
    if line[0] == "@":  # skip  all  headers
        continue # if it is not a header continue
    else:
        record = line.split("\t") # transform each line into an array
        flag = int(record[1])  # get the flag value
        counter_record += 1 # add one to the record counter
        if ((flag & unmapped) == unmapped): # check bitwise comparison for unmapped
            counter_unmapped += 1 # if it is unmapped add one to the counter
        elif ((flag & secondary) == secondary): # check if it is a secondary alignment
            counter_primary += 1 # if it is secondary count it
        elif ((flag & supp == supp)):
            counter_supp += 1
not_sec  = counter_record - (counter_primary + counter_unmapped + counter_supp)# subtract the secondary alignements from the total
print "%s reads unmapped" % (counter_unmapped) # print unmapped
print "%s reads were primary alignments" %  (not_sec) # print primary alignments
