#!/usr/bin/env python

###Imports###
import numpy as np
import pandas as pd
import re
import sys
import collections

###Definitions###
vcf_file = sys.argv[1]

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
vcf_df = pd.read_table(vcf_file, delimiter = '\t', header = 1, skiprows = len(header_lines) - 1) # read in vcf
total_snps = len(vcf_df.index) # get starting number of snps
print("Total SNPs: %s" % total_snps)

id = vcf_df['ID'] # get IDs
loci = id.apply(lambda x: x.split('_')[0]) # get locus names

d = collections.defaultdict(list)
for key, value in zip(loci.values, loci.index): # create dictionary with locus ids as keys and indices as values
    d[key].append(value)

single_snp_inds = [int(np.random.choice(ind_list, 1)) for ind_list in d.values()] # sample a single snp index for each locus 

single_snp_df = vcf_df.iloc[single_snp_inds] # use single snp indices to subsample vcf dataframe

single_snps = len(single_snp_df) # get number of subsampled snps
print("Subsampled SNPs: %s" % single_snps)

out_file = ''.join([vcf_file.split('vcf')[0], '1snp.vcf'])

with open(out_file, 'a') as out:
    out.write(''.join(header_lines)) # write vcf header lines
    single_snp_df.to_csv(out, sep = '\t', header = True, index = False)