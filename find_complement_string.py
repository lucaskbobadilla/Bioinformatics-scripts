# find the complement string of a DNA string needs to be reverse since the reading
# direction (3-5) is oposite in the complement dna

file_name = "rosalind_revc.txt" # file name

file_content = open(file_name).read().rstrip("\n") # read file

dna = file_content # change file name


import string as str

trans = {ord('A'): 'T', ord('T'): 'A', ord('G'): 'C', ord('C'): 'G'} # specify the complements
dna_trans = dna.translate(trans) # Translate the DNA


def reverse(dna_trans): 
    if len(dna_trans) == 0: 
        return dna_trans 
    else: 
        return reverse(dna_trans[1:]) + dna_trans[0]  # create a reverse function

dna_translate_reverse = reverse(dna_trans) # reverse the string
print(dna_translate_reverse) #print the string
