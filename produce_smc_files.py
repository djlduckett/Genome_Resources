#!/usr/bin/env python

###Imports###
import numpy as np
import pandas as pd
import re
import sys
import collections
import subprocess

###Definitions###
vcf_file = sys.argv[1]
population = sys.argv[2]
samples = sys.argv[3] # separated by commas
path = sys.argv[4] # path to smc++

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

chrom = vcf_df['#CHROM'] # get chromosomes
chrom = chrom.astype('str')
chrom_list = list(set(chrom))
print("# chromosomes/contigs: %s" % len(chrom_list))

for x in chrom_list:
    smc_file = x + '.smc.gz'
    command = path + 'smc++ vcf2smc ' + vcf_file + ' ' + smc_file + ' ' + x + ' ' + population + ':' + samples
    #print(command)
    subprocess.call(command, shell = True)
