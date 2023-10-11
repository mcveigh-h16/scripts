# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 14:22:38 2020

@author: mcveigh

Takes a fasta file as input, looks up the taxname using srcchk an rewrites a new
fasta file adding a hard coded line of text to the definition lines
23S ribosomal RNA sequence in this example
replacement for rna_fixer.pl

"""



import Bio
import os
import sys
import pandas as pd

inputfile = sys.argv[1]
outputfile = sys.argv[2]

from Bio import SeqIO
sequences = []  
accessions = []

for seq_record in SeqIO.parse(inputfile, "fasta"): 
    s = seq_record
    str_id = seq_record.id
    accessions.append(str_id)
with open('acclist', 'w') as filehandle:
    for listitem in accessions:
        filehandle.write('%s\n' % listitem)
os.system("/netopt/ncbi_tools64/bin/srcchk -i acclist -f taxname,taxid -o acclist.taxdata")
 
taxdata_file_name = (r'acclist.taxdata')    
taxdata_df = pd.read_csv(taxdata_file_name, sep='\t', index_col=None, low_memory=False)
header = taxdata_df.iloc[0]
taxdata_df.rename(columns = header)

taxdata_df.to_csv('acclist.taxmap', sep='\t', index=False, header=False, columns=['accession','taxid'])
#print(taxdata_df)
    
for seq_record in SeqIO.parse(inputfile, "fasta"): 
    s = seq_record
    tmp_df = taxdata_df[taxdata_df['accession'] == seq_record.id]
    if not tmp_df.empty:
        taxname = tmp_df['organism'].iloc[0]
        print("I found ", taxname)
        seq_record.description = (taxname + " 23S ribosomal RNA sequence")
    sequences.append(s)   
        
        
        #s = seq_record
   # str_id = seq_record.id
   # print(str_id)
    #if str_id in taxdata_df['accession'].tolist():
    #if str_id == taxdata_df['accession'].tolist():
        #taxname = taxdata_df['organism'].iloc[0]
        #print("I found this name ", taxname)
        #seq_record.description = (taxname + " 23S ribosomal RNA sequence")
    #sequences.append(s)
    
SeqIO.write(sequences, outputfile, "fasta")  

fh = open(outputfile)
n = 0
for line in fh:
    if line.startswith(">"):
        n += 1
fh.close()
print("The number of sequences in the final output file is: ", n)
print ("Sourceupdater script is done and output is saved in:", outputfile)