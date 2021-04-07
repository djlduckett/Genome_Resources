#!/usr/bin/env python


###Imports###
import pandas as pd


###Definitions###

gff_file = sys.argv[1]
gofeat_file = sys.argv[2]


###Functions###
def merge_annotations(locus_tag, length, gff_squeeze):
    locs = gff_squeeze.loc[gff_squeeze['ID'] == locus_tag]

###Main###
gff = pd.read_table(gff_file, header = None, sep = '\t')
gff_squeeze = gff[[0,3,4,8]]
gff_squeeze.columns = ['Scaffold', 'Start', 'Stop', 'ID']
gff_squeeze['ID'] = gff_squeeze['ID'].apply(lambda x: x.str.replace('^ID=', '').str.replace(':.*', ''))
gff_squeeze['Start', 'Stop'] = gff_squeeze['Start', 'Stop'].apply(pd.to_numeric)
gff_squeeze['Length'] = gff_squeeze['Stop'] - gff_squeeze['Start']


gofeat = pd.read_table(gofeat_file, sep = ';')
gofeat_squeeze = gofeat[['Locus tag', 'Length', 'Product']]
