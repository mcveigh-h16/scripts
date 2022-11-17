# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 08:14:58 2022

@author: mcveigh
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqFeature
import sys
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

for seq_record in SeqIO.parse(inputfile, "genbank"): 
    seq_record.description = seq_record.annotations["organism"]
    orgname = seq_record.annotations["organism"]
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
            
#print(genes)
#print(proteins)

print("cds count ", cdscount)
print("gene count ", genecount)

protein_df = pd.DataFrame(proteins, columns=['accession', 'protein_name'])
#print(protein_df.head(100))
proteincount = protein_df.groupby('protein_name').size().sort_values(ascending=False)
print(proteincount)
proteincount.to_csv('protein_count', index=True, header=True, sep ='\t')

gene_df = pd.DataFrame(genes, columns=['accession', 'gene_name'])
#print(gene_df.head)
genecount = gene_df.groupby('gene_name').size().sort_values(ascending=False)
print(genecount)
genecount.to_csv('gene_count', index=True, header=True, sep ='\t')