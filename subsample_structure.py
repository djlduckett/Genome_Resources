#!/usr/bin/env python

###Imports###
import numpy as np
import pandas as pd
import re
import sys
import random

###Definitions###
str_file = sys.argv[1]
n_snps = int(sys.argv[2])
out_file = sys.argv[3]

###Main###
str_df = pd.read_table(str_file, delimiter = '\s+', header = 0) # read in structure file
total_snps = len(str_df.columns) - 2 # get starting number of snps
print("Total SNPs: %s" % total_snps)

if n_snps > total_snps:
    print("Number to subsample is more than the total number of SNPs in the STRUCTURE file!")

keep_indices = random.sample(range(0,total_snps - 1), n_snps) # get snp indices to keep
snp_mod_df = str_df.iloc[:, keep_indices] # use subsampled snp indices to subsample structure dataframe

single_snps = len(snp_mod_df.columns) # get number of subsampled snps
print("Subsampled SNPs: %s" % single_snps)

with open(out_file, 'a') as out:
    snp_mod_df.to_csv(out, sep = ' ', header = True, index = True)