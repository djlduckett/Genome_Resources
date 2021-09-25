#!/usr/bin/env python

###Imports###
import glob
import os
import re
import sys

###Definitions###
out_file = sys.argv[1]

###Main###
out = open(out_file, 'w+')
dirs = glob.glob('*/') # get directories
for dir in dirs:
    os.chdir(dir)
    taxon = os.getcwd().split('/')[-1] # get taxon name
    o_files = glob.glob('*.xml.o*') # get output files in folder
    for file in o_files:
        marker = file.split('.xml')[0] # get marker name
        with open(file, 'r') as z:
            for line in z:
                if bool(re.search('  \d+ taxa', line)) == True:
                    num_taxa = line.split(' ')[2] # get number of taxa
                elif bool(re.search('  \d+ sites', line)) == True:
                    num_sites = line.split(' ')[2] # get number of sites
        out_line = taxon + ',' + marker + ',' + num_taxa + ',' + num_sites + '\n'
        out.write(out_line)
    os.chdir('../') # go back a directory
out.close()
