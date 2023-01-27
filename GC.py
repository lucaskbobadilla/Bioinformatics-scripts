# count the GC content in a fasta file

from Bio import SeqIO
import argparse


# get all fasta files created and read it as a dictionary list
import os
   
files = {}
 
for filename in os.listdir():
    if os.path.isfile(filename) \
       and filename.startswith("Rosalind") \
       and not filename in files:
        with open(filename, "r") as file:
            files[filename] = file.read().rstrip("\n")
            
#remove the first line of each DNA file and print it
sequences = {}
for filename, text in files.items():
	text = text.split("\n",1)[1]
	sequences[filename] = text


GC = {} # create an empty dictionary to store GC values
#create a new dictionary with the proportions for each GC 
for filename,text in sequences.items():
	total = int(len(text))
	c = text.count("C")
	g = text.count("G")
	gc_total = c + g
	gc_content = gc_total/total
	GC[filename] = round(gc_content*100,6)

#get larger value of GC content

max_value = max(GC.values())  # maximum value
max_keys = [k for k, v in GC.items() if v == max_value] # getting all keys containing the `maximum`

print(GC)
print(*max_keys,'\n', max_value) # print answer
