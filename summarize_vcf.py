#!/usr/bin/env python

###Imports###
import numpy as np
import pandas as pd
import re
import sys

###Definitions###
vcf_file = sys.argv[1]

###Functions###
def make_sample_dict(sample_name):
    allele1 = vcf_df[sample_name].apply(lambda x: x.split('/')[0])
    allele2 = vcf_df[sample_name].apply(lambda x: x.split('/')[1].split(':')[0])
    depth = vcf_df[sample_name].apply(lambda x: x.split(':')[1])



###Main###
vcf_df = pd.read_table(vcf_file, delimiter = '\t', header = 1, skiprows = 9) # read in vcf
num_snps = len(vcf_df) # get number of SNPs
num_samples = len(vcf_df.columns[9:]) # get number of samples
