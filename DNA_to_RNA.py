# change a DNA to RNA

file_name = "rosalind_rna.txt" # file name

file_content = open(file_name).read().rstrip("\n") # read file

dna = file_content # change file name

print(dna.replace("T", "U")) # replace T for U