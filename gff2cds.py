#!/usr/bin/env python


###Imports###
import sys
import pandas as pd


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

def extract_cds(scaff_name, start, stop, product, scaff_index, outfile) # extract sequences corresponding to gff/gofeat annotations and write them to a new fasta file
    whole_scaff = scaff_index.get_raw(row['ID'].astype('string')).decode()
	clean_scaff = re.sub("[^A-Z]+", "", whole_scaff)
    
    out = open(outfile, 'a')
	new_name = scaff_name + '_' + row['Product'].astype('string')
	new_scaff = clean_scaff[row['Start'].astype('int64')-1:row['Stop'].astype('int64')-1]
	new_scaff = textwrap.fill(new_scaff, width = 60)
	out.writelines('>' + scaff_name + '\n' + new_scaff + '\n\n')
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

db_list = [merge_annotations(x, y, gff_squeeze) for x in gofeat_squeeze['Locus tag'] for y in gofeat_squeeze['Product']] # create list of dfs - 1 df per id
db_df = pd.concat(db_list) # combine to single df
db_df.to_csv('annotations.db', sep = '\t')


scaff_names = get_scaff_names(scaff_file)
scaff_index = SeqIO.index(scaff_file, 'fasta')

for row in len(db_df):
    extract_cds(row, outfile)