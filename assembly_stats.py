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
    num_gaps = len(re.findall('N{10,}', seq))
    return length, num_A, num_T, num_G, num_C, num_N, len_non_N, num_gaps

def get_N50_L50(scaff_lengths, assembly_length):
    half_assembly = int(numpy.ceil(assembly_length // 2))
    sorted_lengths = sorted(scaff_lengths, reverse = True)
    bp_count = 0
    scaff_count = 0
    while bp_count < half_assembly:
        bp_count += sorted_lengths[scaff_count]
        scaff_count += 1
    return sorted_lengths[scaff_count - 1], scaff_count


###Main###

scaff_names = get_scaff_names(infile) # get scaffold lengths
scaff_index = SeqIO.index(infile, 'fasta')

scaff_lengths, As, Ts, Gs, Cs, Ns, non_Ns, gaps = map(list, zip(*[get_scaffold_stats(scaff_index, name) for name in scaff_names]))

num_scaffs = len(scaff_lengths)
assembly_length = sum(scaff_lengths)
assembly_percent_A = round(sum(As) / assembly_length * 100, 3)
assembly_percent_T = round(sum(Ts) / assembly_length * 100, 3)
assembly_percent_G = round(sum(Gs) / assembly_length * 100, 3)
assembly_percent_C = round(sum(Cs) / assembly_length * 100, 3)
assembly_Ns = sum(Ns)
assembly_non_Ns = sum(non_Ns)
assembly_gaps = sum(gaps)

max_scaff_length = max(scaff_lengths)
min_scaff_length = min(scaff_lengths)

N50, L50 = get_N50_L50(scaff_lengths, assembly_length)

assembly_percent_N = round(assembly_Ns / assembly_length * 100, 3)

percent_GC = round((sum(Gs) + sum(Cs)) / assembly_non_Ns * 100, 3)

outfile = infile.split('.')[0] + '.stats'
with open(outfile, 'w') as out:
    print('Assembly Length = ' + str(assembly_length), file = out)
    print('# Scaffolds = ' + str(num_scaffs), file = out)
    print('Longest Scaffold = ' + str(max_scaff_length), file = out)
    print('Shortest Scaffold = ' + str(min_scaff_length), file = out)
    print('N50 = ' + str(N50), file = out)
    print('L50 = ' + str(L50), file = out)
    print('# Gaps = ' + str(assembly_gaps), file = out)
    print('%N = ' + str(assembly_percent_N), file = out)
    print('%A = ' + str(assembly_percent_A), file = out)
    print('%T = ' + str(assembly_percent_T), file = out)
    print('%G = ' + str(assembly_percent_G), file = out)
    print('%C = ' + str(assembly_percent_C), file = out)
    print('%GC = ' + str(percent_GC), file = out)
