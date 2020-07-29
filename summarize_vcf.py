#!/usr/bin/env python

###Imports###
import numpy as np
import pandas as pd
import re
import sys

###Definitions###
vcf_file = sys.argv[1]

###Functions###
def get_sample_stats(vcf_df, sample_name):
    allele1 = vcf_df[sample_name].apply(lambda x: x.split('/')[0]) # get allele 1
    allele2 = vcf_df[sample_name].apply(lambda x: x.split('/')[1].split(':')[0]) # get allele 2
    snp_bin = pd.Series(map(make_binary, allele1, allele2)) # convert genotype to binary
    n_missing = allele1.value_counts()['.'] + allele2.value_counts()['.'] # count number of missing sites
    n_hets = snp_bin.value_counts()['1'] # get number of heterozygous sites
    depths = vcf_df[sample_name].apply(lambda x: x.split(':')[1]) # get depth for each site
    mean_depth = round(pd.to_numeric(depths).mean(), 2) # calculate mean depth
    #sample_df = pd.DataFrame([allele1, allele2, depth], index = np.arange(len(allele1)), columns = [sample_name])
    return(list([n_missing, n_hets, mean_depth]))

def make_binary(a1, a2):
    if a1 == "0" and a2 == "0": # homozygote 1
        return("0")
    elif a1 == "1" and a2 == "1": # homozygote 2
        return("2")
    elif a1 == "0" and a2 == "1": # heterozygote 1
        return("1")
    elif a1 == "1" and a2 == "0": # heterozygote 2
        return("1")
    else:
        return(np.NaN) 



###Main###
vcf_df = pd.read_table(vcf_file, delimiter = '\t', header = 1, skiprows = 9) # read in vcf
num_snps = len(vcf_df) # get number of SNPs
num_samples = len(vcf_df.columns[9:]) # get number of samples
sample_df = pd.DataFrame([get_sample_stats(vcf_df, s) for s in vcf_df.columns[9:]]) # get df of stats over samples
sample_df.index = vcf_df.columns[9:] # set rownames to sample names
sample_df.columns = ['N_count', 'het_count', 'mean_depth'] # set column names to statistic names
sample_df['N_percent'] = round(sample_df['N_count'] / num_snps * 2 * 100, 2) # calculate percent Ns per sample
sample_df['het_percent'] = round(sample_df['het_count'] / num_snps * 100, 2) # calculate percent hets per sample
sample_df_sorted = sample_df[['N_count', 'N_percent', 'het_count', 'het_percent', 'mean_depth']]

out_file = 'vcf_summary.txt'

with open(out_file, 'a') as out:
    out.write('VCF file contains %s SNPs and %s samples\n\n' % (num_snps, num_samples)) # write vcf header lines
    sample_df_sorted.to_csv(out, sep = '\t', header = True, index = True)
