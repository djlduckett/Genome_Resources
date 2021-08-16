#!/usr/bin/env python


###Imports###

import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import random


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

#dup_inds = [idx for idx,val in enumerate(seq_names) if val in seq_names[:idx]] # get indices of any duplicates
#seq_names = [idx for val,idx in enumerate(seq_names) if val not in dup_inds] # remove duplicate sequence names
#seqs = [idx for val,idx in enumerate(seqs) if val not in dup_inds] # remove duplicate sequences
seq_dict = dict(zip(seq_names, seqs)) # create dictionary with sequence names as keys and sequences as values

percent_N = [N / s for N, s in zip(Ns, seq_lengths)] # calculate percent N for each sequence
N_above_ten = [ind for ind,val in enumerate(percent_N) if val > 0.1] # get indices of sequences with a large %N
#print(N_above_ten)

too_long = [ind for ind,val in enumerate(seq_lengths) if val > 1000] # get indices of sequences longer than 1000bp
#print(too_long)

inds_to_remove = N_above_ten + too_long # combine indices of sequences too long or with too many Ns
#print(inds_to_remove)
if bool(inds_to_remove) == True: # if list not empty, remove duplicate indices
    inds_to_remove = list(set(inds_to_remove))
    keys_to_remove = [list(seq_dict.keys())[q] for q in inds_to_remove] # match indices to remove to samples (dict keys)
    for k in keys_to_remove: # remove bad keys/values from dictionary
        black_hole = seq_dict.pop(k, None)
#print(seq_dict.keys())
seq_lengths_final = [len(b) for b in list(seq_dict.values())] # get lengths of remaining sequences
#seq_names_final = list(seq_dict.keys()) # get names of remaining sequences

longest_value = sorted(seq_lengths_final, reverse = True)[0] # get longest length
longest_name = [x for x in list(seq_dict.keys()) if len(seq_dict[x]) == longest_value] # get sequence name of longest length
#print(longest_name)
if len(longest_name) > 1: # if there is more than 1 name (multiple sequences have same length), randomly choose one
    longest_name = random.choice(longest_name)

longest_name = re.search("uce.+_contigs", str(longest_name)).group(0) # reformat sequence name

print(">" + longest_name + "\n" + seq_dict[longest_name]) # print sequence in fasta format

    

