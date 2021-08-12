#!/usr/bin/env python


###Imports###

import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy


###Definitions###

infile = sys.argv[1]

###Functions###

def get_seq_names(file):
    names = []
    for seq in SeqIO.parse(file, 'fasta'):
        names.append(seq.id)
    return names

def get_seq_stats(seq_index, seq_name):
    seq = re.sub(">.+\n", "", seq_index.get_raw(seq_name).decode()) # Get sequence record andremove metadata
    seq = re.sub("\s+", "", seq).upper() # remove whitespace and capitalize
    length = len(seq)
    num_N = len(re.findall('N', seq))
    return length, num_N, seq


###Main###

seq_names = get_seq_names(infile)  # get sequence names
seq_index = SeqIO.index(infile, 'fasta')
#seq_dict = SeqIO.to_dict(SeqIO.parse(infile, "fasta"))

seq_lengths, Ns, seqs = map(list, zip(*[get_seq_stats(seq_index, name) for name in seq_names])) # get lists of lengths, number of ns, and sequences
seq_dict = dict(zip(seq_names, seqs)) # create dictionary with sequence names as keys and sequences as values

percent_N = [N / s for N, s in zip(Ns, seq_lengths)] # calculate percent N for each sequence
N_above_ten = [ind for ind,val in enumerate(percent_N) if val > 0.1] # get indices of sequences with a large %N

too_long = [ind for ind,val in enumerate(seq_lengths) if val > 1000] # get indices of sequences longer than 1000bp

inds_to_remove = N_above_ten.append(too_long) # combine indices of sequences too long or with too many Ns

if len(inds_to_remove) > 0: # if list not empty, remove duplicate indices
    inds_to_remove = list(set(inds_to_remove))
    keys_to_remove = [list(seq_dict.keys())[q] for q in inds_to_remove] # match indices to remove to samples (dict keys)

    for k in keys_to_remove: # remove bad keys/values from dictionary
        seq_dict.pop(k, None)

    seq_lengths = [len for ind,len in enumerate(seq_lengths) if ind not in inds_to_remove] # remove sequence lengths of bad indices
    

