#!/usr/bin/env python

###Imports###
import numpy as np
import pandas as pd
import re
import sys
import collections
import csv

###Definitions###
vcf_file = sys.argv[1]
out_file = 'lengths.txt'

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
contig_lines = [x for x in header_lines if '##contig=' in x] # get contig header lines
contig_names = [re.search('ID=[0-9a-zA-Z]+', i).group().split('=')[1] for i in contig_lines]
contig_lengths = [int(re.search('length=\d+', j).group().split('=')[1]) for j in contig_lines]

with open(out_file, 'w') as out:
    writer = csv.writer(out)
    writer.writerows(zip(contig_names, contig_lengths)) # write contig name, length as csv

print("# contigs: %s" % len(contig_names))
print("Mean length: %s" % np.mean(contig_lengths))
print("Median length: %s" % np.median(contig_lengths))
print("Contig min,max: %s,%s" % (min(contig_lengths), max(contig_lengths)))
