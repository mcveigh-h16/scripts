# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 08:14:58 2022

@author: mcveigh
"""
#script to count gene names and protein names from a gbk file and tally the results
#accession prefixes in rejectacc are excluded
#output to a tab delimited file


import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
import sys
from datetime import datetime
pd.set_option('display.max_rows', 100)
genecount = 0
#feature_count = 0
cdscount = 0
#proteincount = 0
orgname = []
genename = []
genes = []
protname = []
proteins = []
#RNAs = []
#table = []
phash = []
ghash = []

inputfile = sys.argv[1]
lastgfile = sys.argv[2]
#fh = open(lastgfile)
rejectacc = ["AP", "BS", "AL", "BX", "CR", "CT", "CU", "FP", "FQ", "FR", "AE", "CP", "CY"]

startTime = datetime.now()
print("Start time is ", startTime)

for seq_record in SeqIO.parse(inputfile, "genbank"): 
    seq_record.description = seq_record.annotations["organism"]
    orgname = seq_record.annotations["organism"]
    str_id = str(seq_record.id)
    #print("found acc", str_id)
    two = str_id[:2]
    if two not in rejectacc:   
        #print("test2", str_id)
        if orgname != "Severe acute respiratory syndrome coronavirus 2":
            for feature in seq_record.features:
                if feature.type == 'CDS':
                    protname = str(feature.qualifiers.get("product")) 
                    #protname = feature.qualifiers.get("product")
                    cdscount += 1
                    phash = seq_record.id, protname
                    phash = [sub.replace('[\'', '') for sub in phash]
                    phash = [sub.replace('\']', '') for sub in phash]
                    proteins.append(phash)
                if feature.type == 'gene':
                    genename = str(feature.qualifiers.get("gene")) 
                    #genename = feature.qualifiers.get("gene")
                    genecount += 1
                    ghash = seq_record.id, genename
                    ghash = [sub.replace('[\'', '') for sub in ghash]
                    ghash = [sub.replace('\']', '') for sub in ghash]
                    genes.append(ghash)
    elif two in rejectacc: 
        print("reject accession found", str_id)
        #else:
            #print(orgname)
            
#print(genes)
#print(proteins)

print("cds count ", cdscount)
print("gene count ", genecount)

protein_df = pd.DataFrame(proteins, columns=['accession', 'protein_name'])
#print(protein_df.head(100))
proteincount = protein_df.groupby('protein_name').size().sort_values(ascending=False).reset_index()
proteincount.rename( columns={0 :'count'}, inplace=True )
gene_df = pd.DataFrame(genes, columns=['accession', 'gene_name'])

genecount_df = gene_df.groupby('gene_name').size().sort_values(ascending=False).reset_index()
genecount_df.rename( columns={0 :'count1'}, inplace=True )
print("current genecount df")
print(genecount_df.head(20))

#lastgcount_df = pd.read_csv(lastgfile, sep='\t', index_col=None, usecols=[1,2,3], na_values=['-'], names=["gene", "count1x","count1"])
lastgcount_df = pd.read_csv(lastgfile, sep='\t', index_col=None, usecols=[1,2,3], na_values=['-'])
lastgcount_df.rename(columns={"count1_x": "lastcount", "count1_y": "count1"}, inplace=True)
print("last count df")
print(lastgcount_df.head(20))

totalgene_df = pd.concat([genecount_df, lastgcount_df]).groupby('gene_name')['count1'].sum().sort_values(ascending=False).reset_index()

print("total gene df")
print(totalgene_df.head(20))

run_total_df = pd.merge(genecount_df, totalgene_df, left_on='gene_name', right_on='gene_name')
print("run total")
print(run_total_df.head(20))

#print(proteincount)
proteincount.to_csv('protein_count', index=True, header=True, sep ='\t')
#print(genecount)
genecount_df.to_csv('gene_count', index=True, header=True, sep ='\t')
#print("total count df")
#print(totalgene_df)
run_total_df.to_csv('gene_count_total', index=True, header=True, sep ='\t')




stopTime = datetime.now()
print("Stop time is ", stopTime)