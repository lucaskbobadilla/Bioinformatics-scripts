#define function for motif finding
import sys

DNA = str(sys.argv[1])
DNA = open(DNA, "r").read().rstrip("\n")

motif = str(sys.argv[2])
motif = open(motif,"r").read().rstrip("\n")


def finding_motif(DNA, motif):
    result = [] # create an empy vector
    mlen = len(motif) # get length of motif
    slen = len(DNA) # get length of DNA
    for i in range(0, slen - mlen + 1): # for each letter in the sequence checking every length of the motif
        if DNA[i:i+mlen] == motif: # if the sequence is equal to the motif
            result.append(i+1)     # save the i value into the results 
    return result

print(finding_motif(DNA,motif))





