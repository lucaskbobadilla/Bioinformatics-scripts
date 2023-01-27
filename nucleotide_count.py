# Nucleotide count - rosalind exercise

file = "rosalind_dna.txt" # file name

my_file = open(file) # open file

my_file_contents = my_file.read() # read file

dna = my_file_contents.rstrip("\n")


count_A = dna.count("A")

count_G = dna.count("G")

count_T = dna.count("T")

count_C = dna.count("C")

print(str(count_A) + " " + str(count_C) + " " + str(count_G) + " " + str(count_T))


