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
definition = []
total_record_count = 0
mito_chloro_genome = []
organelle_genome_count = 0
records_counted = 0

inputfile = sys.argv[1]
lastgenefile = sys.argv[2]
lastproteinfile = sys.argv[3]
#fh = open(lastgfile)
rejectacc = ["AP", "BS", "AL", "BX", "CR", "CT", "CU", "FP", "FQ", "FR", "AE", "CP", "CY"]
organelle = ["mitochondrion, complete genome", "mitochondrion, partial genome", "chloroplast, complete genome", "chloroplast, partial genome"]

starttime = datetime.now()
print("Start time is ", starttime)

for seq_record in SeqIO.parse(inputfile, "genbank"): 
    total_record_count += 1
    definition = str(seq_record.description)
    seq_record.description = seq_record.annotations["organism"]
    orgname = seq_record.annotations["organism"]
    str_id = str(seq_record.id)
    #print("found acc", str_id)
    accprefix = str_id[:2]
    if accprefix not in rejectacc:   
        #print("test2", str_id)
        #Insert organelle exclusion
        if any(ext in definition for ext in organelle):
            #print(seq_record.id, definition)
            organelle_genome_count += 1
            #mito_chloro_genome.append(seq_record)
        else:
            if orgname != "Severe acute respiratory syndrome coronavirus 2":
                #print("counting loop", seq_record.id, definition)
                records_counted += 1
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
    elif accprefix in rejectacc: 
        print("reject accession found", str_id)
        #else:
            #print(orgname)
            
#print(genes)
#print(proteins)

print("cds count ", cdscount)
print("gene count ", genecount)
print("total records ", total_record_count)
print("organelle genome count ", organelle_genome_count)
print("records we counted ", records_counted)

protein_df = pd.DataFrame(proteins, columns=['accession', 'protein_name'])
#print(protein_df.head(100))
proteincount_df = protein_df.groupby('protein_name').size().sort_values(ascending=False).reset_index()
proteincount_df.rename( columns={0 :'count1'}, inplace=True )

gene_df = pd.DataFrame(genes, columns=['accession', 'gene_name'])
genecount_df = gene_df.groupby('gene_name').size().sort_values(ascending=False).reset_index()
genecount_df.rename( columns={0 :'count1'}, inplace=True )
print("current genecount df")
print(genecount_df.head(20))

#determine gene counts current and cummulative
#lastgenecount_df = pd.read_csv(lastgenefile, sep='\t', index_col=None, usecols=[1,2,3], na_values=['-'], names=["gene", "count1x","count1"])
lastgenecount_df = pd.read_csv(lastgenefile, sep='\t', index_col=None, usecols=[1,2,3], na_values=['-'])
lastgenecount_df.rename(columns={"count1_x": "lastcount", "count1_y": "count1"}, inplace=True)
print("last gene count df")
print(lastgenecount_df.head(20))
totalgene_df = pd.concat([genecount_df, lastgenecount_df]).groupby('gene_name')['count1'].sum().sort_values(ascending=False).reset_index()

print("total gene df")
print(totalgene_df.head(20))
run_genetotal_df = pd.merge(genecount_df, totalgene_df, left_on='gene_name', right_on='gene_name')
print("run gene total")
print(run_genetotal_df.head(20))

#determine protein counts current and cummulative
lastprotcount_df = pd.read_csv(lastproteinfile, sep='\t', index_col=None, usecols=[1,2,3], na_values=['-'])
lastprotcount_df.rename(columns={"count1_x": "lastcount", "count1_y": "count1"}, inplace=True)
print("last protein count df")
print(lastprotcount_df.head(20))
totalprotein_df = pd.concat([proteincount_df, lastprotcount_df]).groupby('protein_name')['count1'].sum().sort_values(ascending=False).reset_index()

run_proteintotal_df = pd.merge(proteincount_df, totalprotein_df, left_on='protein_name', right_on='protein_name')
print("run protein total")
print(run_proteintotal_df.head(20))

#print(proteincount_df)
proteincount_df.to_csv('protein_count', index=True, header=True, sep ='\t')
#print(genecount)
genecount_df.to_csv('gene_count', index=True, header=True, sep ='\t')
#print("total count df")
#print(totalgene_df)
run_genetotal_df.to_csv('gene_count_total', index=True, header=True, sep ='\t')
run_proteintotal_df.to_csv('protein_count_total', index=True, header=True, sep ='\t')



stoptime = datetime.now()
print("Stop time is ", stoptime)