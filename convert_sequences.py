#
# PYTHON PROGRAM TO CONVERT cDNA SEQUENCES TO PROTEIN SEQUENCES
#


from Bio import SeqIO
import sys

dna_file = sys.argv[1]

from Bio.SeqRecord import SeqRecord

with open("AAseq.fasta", 'w') as aa_fa:
    for dna_record in SeqIO.parse(dna_file, 'fasta'):

        # use both fwd and rev sequences
        dna_seqs = [dna_record.seq, dna_record.seq.reverse_complement()]

        # generate all translation frames
        aa_seqs = (s[i:].translate(to_stop=True) for i in range(3) for s in dna_seqs)

        # select the longest one
        max_aa = max(aa_seqs, key=len)

        # write new record
        aa_record = SeqRecord(max_aa, id=dna_record.id, description="translated sequence")
        SeqIO.write(aa_record, aa_fa, 'fasta')
