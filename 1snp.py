#!/usr/bin/env python

###Imports###
import numpy as np
import pandas as pd

###Definitions###
vcf_file = sys.argv[1]

###Functions###


id = vcf_df['ID'] # get IDs
loci = id.apply(lambda x: x.split('_')[0]) # get locus names


###Main###
vcf_df = pd.read_table(vcf_file, delimiter = '\t', header = 1, skiprows = 10) # read in vcf
total_snps = len(vcf_df.index) # get starting number of snps