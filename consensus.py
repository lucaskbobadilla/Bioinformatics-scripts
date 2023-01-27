#!/usr/bin/env python

from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio import AlignIO
import os
import numpy as np
import sys


#find the consensus alignment
alignment = AlignIO.read(sys.argv[1], 'fasta')
summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus(.5)

# score matrix
my_pssm = summary_align.pos_specific_score_matrix(consensus) 

print(consensus)
print(my_pssm)