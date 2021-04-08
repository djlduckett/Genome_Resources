#!/usr/bin/env python


###Imports###
import sys

###Definitions###

fasta_file = sys.argv[1]

###Functions###
def get_scaff_names(file): # get scaffold names in fasta file
    names = []
    for seq in SeqIO.parse(file, 'fasta'):
        names.append(seq.id)
    return names

###Main###
scaff_names = get_scaff_names(fasta_file)
scaff_index = SeqIO.index(fasta_file, 'fasta')


out2 = open(outfile, 'a')
for name in scaff_names:
    whole_scaff = scaff_index.get_raw(scaff_name).decode()
	clean_scaff = re.sub("[^A-Z]+", "", whole_scaff)
	
	new_scaff = clean_scaff
	new_scaff = textwrap.fill(new_scaff, width = 60)
	if len(new_scaff) < 200 or len(re.findall('N', new_scaff.upper())) == len(new_scaff):
		continue
	out2.writelines('>' + scaff_name + '_' + str(i) + '\n' + new_scaff + '\n\n')
	out2.close()