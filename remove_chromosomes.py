#!/usr/bin/env python

###Imports###
import numpy as np
import pandas as pd
import re
import sys
import collections

###Definitions###
vcf_file = sys.argv[1]
chr_file = sys.argv[2]
out_file = sys.argv[3]

###Functions###

def get_header_lines(vcf_file):
    lines = []
    comment = True    
    with open(vcf_file, 'r') as f:
        while comment == True: # while line begins with '##'
            line = f.readline()
            comment = bool(re.search('##', line))
            if comment == True:
                lines.append(line)
    return(lines)


###Main###
header_lines = get_header_lines(vcf_file) # get vcf header lines with '##'

f = open(chr_file, 'r')
chr_list = f.readlines()
chr_list = [x.strip('\n') for x in chr_list]

vcf_df = pd.read_table(vcf_file, delimiter = '\t', header = 1, skiprows = len(header_lines) - 1) # read in vcf
total_snps = len(vcf_df.index) # get starting number of snps
print("Total SNPs: %s" % total_snps)

chrom = vcf_df['#CHROM'] # get chromosomes

snp_dict = collections.defaultdict(list)
for key, value in zip(chrom.values, chrom.index): # create dictionary with locus ids as keys and indices as values
    snp_dict[key].append(value)

for k in chr_list: # remove chromosomes from dictionary
    snp_dict.pop(k, None)
 
keep_indices = [item for sublist in snp_dict.values() for item in sublist] # get snp indices to keep
snp_mod_df = vcf_df.iloc[keep_indices] # use single snp indices to subsample vcf dataframe

single_snps = len(snp_mod_df) # get number of subsampled snps
print("Subsampled SNPs: %s" % single_snps)

with open(out_file, 'a') as out:
    out.write(''.join(header_lines)) # write vcf header lines
    snp_mod_df.to_csv(out, sep = '\t', header = True, index = False)