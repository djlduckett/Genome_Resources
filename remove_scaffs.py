from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy
import sys
import re


###Definitions###

infile = sys.argv[1]
r_file = sys.argv[2]

###Functions###

def get_scaff_names(file):
    names = []
    for seq in SeqIO.parse(file, 'fasta'):
        names.append(seq.id)
    return names
    

###Main###

scaff_names = get_scaff_names(infile) # get scaffold lengths
scaff_index = SeqIO.index(infile, 'fasta')

f = open(r_file, 'r')
r_list = f.readlines()
r_list = [i.strip() for i in r_list]
clean_names = [x for x in scaff_names if x not in r_list]

outfile = infile.split('.')[0] + '_cleaned.fasta'
with open(outfile, 'w') as out:
    for scaff in clean_names:
    	print(scaff_index[scaff].format('fasta'), file = out)

