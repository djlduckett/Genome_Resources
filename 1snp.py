#!/usr/bin/env python

###Imports###
import numpy as np
import pandas as pd

###Definitions###
vcf_file = sys.argv[1]

###Functions###


id = vcf_df['ID'] # get IDs
loci = id.apply(lambda x: x.split('_')[0]) # get locus names

d = collections.defaultdict(list)
for key, value in zip(loci.values, loci.index): # create dictionary with locus ids as keys and indices as values
    d[key].append(value)

single_snp_inds = [int(np.random.choice(ind_list, 1)) for ind_list in d.values()] # sample a single snp index for each locus 

single_snp_df = vcf_df.iloc[single_snp_inds] # use single snp indices to subsample vcf dataframe


###Main###
vcf_df = pd.read_table(vcf_file, delimiter = '\t', header = 1, skiprows = 10) # read in vcf
total_snps = len(vcf_df.index) # get starting number of snps