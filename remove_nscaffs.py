from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy
import sys
import re


###Definitions###

infile = sys.argv[1]

###Functions###

def get_scaff_names(file):
    names = []
    for seq in SeqIO.parse(file, 'fasta'):
        names.append(seq.id)
    return names

def get_scaffold_stats(scaff_index, scaff_name):
    seq = re.sub(">.+\n", "", scaff_index.get_raw(scaff_name).decode()) # Get sequence record andremove metadata
    seq = re.sub("\s+", "", seq).upper() # remove whitespace and capitalize
    num_A = len(re.findall('A', seq))
    num_T = len(re.findall('T', seq))
    num_G = len(re.findall('G', seq))
    num_C = len(re.findall('C', seq))
    num_N = len(re.findall('N', seq))
    length = num_A + num_T + num_G + num_C + num_N
    if num_N != length:
    	return scaff_name
    

###Main###

scaff_names = get_scaff_names(infile) # get scaffold lengths
scaff_index = SeqIO.index(infile, 'fasta')

n_names = [get_scaffold_stats(scaff_index, name) for name in scaff_names]
clean_names = [x for x in scaff_names if x]

outfile = infile.split('.')[0] + '2.fasta'
with open(outfile, 'w') as out:
    for scaff in clean_names:
    	print(scaff_index[scaff].format('fasta'), file = out)

