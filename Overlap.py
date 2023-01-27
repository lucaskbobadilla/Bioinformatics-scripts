#!/bin/python3

######  Overlap Graphs ########

## LIBRARY
import argparse
import itertools 

# GET ARGUMENT FROM USER
p = argparse.ArgumentParser()

## add options
p.add_argument('--input','-i',required=True, help= "Pathway to the query fasta file") # add pathway value
p.add_argument('--kval', '-k',required=True,type=int, help= "k value for graph edges") # add pathway value

## parse arguments
args = p.parse_args()

# Define fasta generator function - from Biopython
def input_fasta(fasta_file):
    name, seq = None, []
    for line in fasta_file:
        line = line.rstrip()
        if line[0] == ">":
            if name: 
                yield (name, ''.join(seq))
            name, seq = line, []    
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))
    

## get fasta file path
fp = args.input # open the fasta file

# parse fasta file into a dictonary
seqs = {}
with open(fp, 'r') as fasta:
    for name, seq in input_fasta(fasta):
        seqs[name[1:]] = seq


k = args.kval # get k value specified by user
# using combination for itertools to make all possible combinations in the dict
edges = []
for name1, name2 in itertools.combinations(seqs.keys(), 2):
    if seqs[name1][-k:] == seqs[name2][:k]:
        edges.append(name1 + " " + name2)
    elif seqs[name2][-k:] == seqs[name1][:k]:
        edges.append(name1 + " " + name2)

for i in edges:
    print(i)



    