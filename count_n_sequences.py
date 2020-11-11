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
    length = len(seq)
    num_A = len(re.findall('A', seq))
    num_T = len(re.findall('T', seq))
    num_G = len(re.findall('G', seq))
    num_C = len(re.findall('C', seq))
    num_N = len(re.findall('N', seq))
    len_non_N = num_A + num_T + num_G + num_C
    if num_N == length:
    	return 1
    else:
    	return 0
    

###Main###

scaff_names = get_scaff_names(infile) # get scaffold lengths
scaff_index = SeqIO.index(infile, 'fasta')

n_sequences = [get_scaffold_stats(scaff_index, name) for name in scaff_names]

print(numpy.sum(n_sequences))
