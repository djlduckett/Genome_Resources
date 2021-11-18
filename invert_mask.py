#!/usr/bin/env python

### This script will invert a bed file, listing regions to be masked instead of regions to be unmasked.
### Note: This will not work if needing to go from masked to unmasked regions


###Imports###
import numpy as np
import pandas as pd
import sys


###Definitions###
bed_file = sys.argv[1]
chr_name = sys.argv[2]
chr_length = sys.argv[3]


###Main###
pre_mask = pd.read_csv(bed_file, sep = '\t', header = None) # no header in bed file
pre_mask = pre_mask.rename({0:'chr', 1:'start', 2:'end'}, axis = 1) # change column names

start = pre_mask['start'] # get starting positions
start2 = np.append(start, chr_length) # add ending coordinate
end = pre_mask['end'] # get ending positions
end2 = np.append(np.array([[1]]), end) # add 1 to beginning
chr = pre_mask['chr'] # get chr names
chr2 = np.append(chr, chr_name) # add an additional entry to match start2 and end2

final_mask = pd.concat([pd.Series(chr2), pd.Series(end2), pd.Series(start2)], axis = 1) # starting and ending coordinates switch places to output masked sites instead of unmasked sites
final_mask_len = len(final_mask)
print(chr_name + ', ' + str(final_mask_len))
final_mask.to_csv(chr_name + '_fixed_mask.bed', header = None, index = None, sep = '\t')