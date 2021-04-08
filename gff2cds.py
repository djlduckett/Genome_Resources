#!/usr/bin/env python


###Imports###
import sys
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import numpy as np
import textwrap
import re
pd.options.mode.chained_assignment = None # turn off warning for chained assignment

###Definitions###

gff_file = sys.argv[1]
gofeat_file = sys.argv[2]
scaff_file = sys.argv[3]


###Functions###
def merge_annotations(locus_tag, product, gff_squeeze): # combine info from gff and gofeat: scaffold name, start, stop, length, ID, product
    locs = gff_squeeze.loc[gff_squeeze['ID'] == locus_tag] # get gff ID that matches gofeat locus tag
    locs['Product'] = product
    return locs

def get_scaff_names(file): # get scaffold names in fasta file
    names = []
    for seq in SeqIO.parse(file, 'fasta'):
        names.append(seq.id)
    return names

def extract_cds(scaff_name, start, stop, product, scaff_index, outfile): # extract sequences corresponding to gff/gofeat annotations and write them to a new fasta file
    whole_scaff = scaff_index.get_raw(str(scaff_name)).decode()
    clean_scaff = re.sub("[^A-Z]+", "", whole_scaff)   
    out = open(outfile, 'a')
    new_name = str(scaff_name) + '_' + str(product)
    new_scaff = clean_scaff[int(start)-1:int(stop)-1]
    new_scaff = textwrap.fill(new_scaff, width = 60)
    out.writelines('>' + new_name + '\n' + new_scaff + '\n')
    out.close()

###Main###
gff = pd.read_table(gff_file, header = None, sep = '\t')
gff_squeeze = gff[[0,3,4,8]]
gff_squeeze.columns = ['Scaffold', 'Start', 'Stop', 'ID']
gff_squeeze['ID'] = pd.Series(gff_squeeze['ID'], dtype="string")
gff_squeeze['ID'] = gff_squeeze['ID'].str.replace('^ID=', '').str.replace(';.*', '') # to match gofeat locus tags
gff_squeeze['Start'] = pd.Series(gff_squeeze['Start'], dtype="int64")
gff_squeeze['Stop'] = pd.Series(gff_squeeze['Stop'], dtype="int64")
gff_squeeze['Length'] = gff_squeeze['Stop'] - gff_squeeze['Start']


gofeat = pd.read_table(gofeat_file, sep = ';')
gofeat_squeeze = gofeat[['Locus tag', 'Length', 'Product']]

db_list = [merge_annotations(x, y, gff_squeeze) for x,y in zip(gofeat_squeeze['Locus tag'], gofeat_squeeze['Product'])] # create list of dfs - 1 df per id
db_df = pd.concat(db_list) # combine to single df
db_df.to_csv('annotations.db', sep = '\t')


scaff_names = get_scaff_names(scaff_file)
scaff_index = SeqIO.index(scaff_file, 'fasta')

outfile = 'cds.fasta'
for row in range(0, len(db_df)):
    try:
        extract_cds(db_df.loc[row, 'Scaffold'], db_df.loc[row, 'Start'], db_df.loc[row, 'Stop'], db_df.loc[row, 'Product'], scaff_index, outfile)
    except:
        print('Error, scaffold ' + str(db_df.loc[row, 'Scaffold']))
        continue