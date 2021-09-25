#!/usr/bin/env python

###Imports###
import re
import sys
import os
import glob

###Definitions###
out_file = sys.argv[1]

###Main###
out = open(out_file, 'w+')
dirs = glob.glob('*/') # get directories
for dir in dirs:
    os.chdir(dir)
    taxon = os.getcwd().split('/')[-1] # get taxon name
    o_files = glob.glob('*_single.txt') # get output files in folder
    if 'COI_single.txt' in o_files:
        with open('COI_single.txt', 'r') as z:
            line = z.readlines()[0]
            threshold = line.split('threshold time:\t')[-1].strip()
        out_line = taxon + ',COI,' + threshold + '\n'
        out.write(out_line)
    if 'cytb_single.txt' in o_files:
        with open('cytb_single.txt', 'r') as z:
            line = z.readlines()[0]
            threshold = line.split('threshold time:\t')[-1].strip()
        out_line = taxon + ',cytb,'+ threshold + '\n'
        out.write(out_line)
    os.chdir('../') # go back a directory
out.close()