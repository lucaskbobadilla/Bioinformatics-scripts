#!/bin/python

##### N50 and assembly stats calculation #####
### Lucas Kopecky Bobadilla

# import modules
import argparse # for arguments
import gzip # to open gzip files if necessary
import sys # access to variables for SDOUT

#
# expecify arguments
#
p = argparse.ArgumentParser() #parse arguments
# specify file argument path
p.add_argument('-f', type=str, required=True, help = "This argument is the file path")
args = p.parse_args() # get arguments

# open file - first check if it is gzipped or not
if args.f[-3:] == ".gz":
    fasta = gzip.open(args.f,'rb')
else:
    fasta = open(args.f,'rb')


# get the legnths of each contig/scaffold
lengths = [] # empty array to store lengths of each sequences
length = 0 # counter for each length
n_seq = 0 # counter of number of sequences
for line in fasta: # loop over fasta lines
    line = line.strip("\n") # remove \n
    if line[0] == ">": # check if the line is a header
        n_seq += 1 # count sequence
        if length == 0: # if this is the first sequence continue
            continue
        else: # otherwise, append the length from the previous sequence into the array
            lengths.append(length) # append length
            length = 0 # restart counter for next sequence
    else:
        length += len(line) # increment the length of the sequence with each line

# sort lengths
len_sort = sorted(lengths, reverse = True) # sort array from largest to smallest

# get total and max values
total = sum(len_sort)
large = max(len_sort)

# get n50 and L50
mid = 0 # counter for mid point
n50 = 0 # counter for n50
L50 = 0 # counter for L50
for i in range(len(len_sort)): # loop over len_sort index
    mid += len_sort[i] # increment mid point
    if mid >= total/2: # check if n50 was achieved
        n50 = len_sort[i] # if it was, pass current value to n50
        L50 = i + 1 # if it was, save current index + 1 as L50
        break # stop loop once n50 is achieved

# print results into STDOUT
sys.stdout.write("###### N50 stats #####\nAssembly file name: {}\nAssembly size (bp): {:,}\nNumber of sequences: {}\nLargest sequence: {:,}\nN50 (bp): {:,}\nL50: {}\n".format(args.f,total,n_seq,large,n50,L50) )

fasta.close()
