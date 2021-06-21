#!/usr/bin/env python


###Imports###

import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy as np
import textwrap


###Definitions###

infile = sys.argv[1]
r_file = sys.argv[2]


###Functions###

def get_scaff_names(file):
    names = []
    for seq in SeqIO.parse(file, 'fasta'):
        names.append(seq.id)
    return names


def split_scaff(adapter_line, outfile):
	scaff_name = adapter_line.split('\t')[0]
	coords = adapter_line.split('\t')[2:-1]
	scaff_length = adapter_line.split('\t')[1]
	coords = coords[0].replace('..', ',').split(',')
	coords.insert(0, '0')
	coords.append(scaff_length)
	num_pairs = int(len(coords) / 2)
	sub_list = [0,1] * num_pairs
	coords = [int(x) - int(y) for x,y in zip(coords, sub_list)]
	coords[len(coords) - 1] += 1
	coords_pairs = np.array_split(coords, len(coords) / 2)
	whole_scaff = scaff_index.get_raw(scaff_name).decode()
	clean_scaff = re.sub("^.+\n", "", whole_scaff) # get rid of scaffold name
	clean_scaff = re.sub("\n", "", clean_scaff) # get rid of line breaks
	out2 = open(outfile, 'a')
	for i in range(0,num_pairs):
		new_name = scaff_name + '_' + str(i)
		new_scaff = clean_scaff[coords_pairs[i].tolist()[0]: coords_pairs[i].tolist()[1]]
		new_scaff = textwrap.fill(new_scaff, width = 60)
		if len(new_scaff) < 200 or len(re.findall('N', new_scaff.upper())) == len(new_scaff):
			continue
		out2.writelines('>' + scaff_name + '_' + str(i) + '\n' + new_scaff + '\n\n')
	out2.close()
	

    

###Main###

scaff_names = get_scaff_names(infile) # get scaffold lengths
scaff_index = SeqIO.index(infile, 'fasta')

f = open(r_file, 'r')
r_list = f.readlines()
f.close()

bad_names = [x.split('\t')[0] for x in r_list]
clean_names = [x for x in scaff_names if x not in bad_names]


outfile = infile.split('.')[0] + '_cleaned.fasta'
with open(outfile, 'w') as out:
    for scaff in clean_names:
    	print(scaff_index[scaff].format('fasta'), file = out)
    	
    	
for line in r_list:
	split_scaff(line, outfile)





