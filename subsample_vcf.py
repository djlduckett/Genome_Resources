#!/usr/bin/env python

###Imports###
import numpy as np
import pandas as pd
import re
import sys
import random

###Definitions###
vcf_file = sys.argv[1]
n_snps = int(sys.argv[2])
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

vcf_df = pd.read_table(vcf_file, delimiter = '\t', header = 1, skiprows = len(header_lines) - 1) # read in vcf
total_snps = len(vcf_df.index) # get starting number of snps
print("Total SNPs: %s" % total_snps)

if n_snps > total_snps:
    print("Number to subsample is more than the total number of SNPs in the VCF!")

keep_indices = random.sample(range(0,total_snps - 1), n_snps) # get snp indices to keep
snp_mod_df = vcf_df.iloc[keep_indices] # use subsampled snp indices to subsample vcf dataframe

single_snps = len(snp_mod_df) # get number of subsampled snps
print("Subsampled SNPs: %s" % single_snps)

with open(out_file, 'a') as out:
    out.write(''.join(header_lines)) # write vcf header lines
    snp_mod_df.to_csv(out, sep = '\t', header = True, index = False)