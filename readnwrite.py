# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 08:39:05 2020

@author: mcveigh
"""

#
# Selectively sort fasta file removing selected seqeunces and saving them to a new file
#
import pandas as pd
import Bio
import subprocess

  
file_name_string = (r'ITS_reject_seqs')
rejectlist_df = pd.read_csv(file_name_string, sep='\t', index_col=None, low_memory=False, header=None, names=["accession", "type", "reason"])
rejectlist = rejectlist_df['accession']
print (rejectlist) 
   
from Bio import SeqIO
sequences = [] 
found = []
#print(rejectlist_df['accession'])
for seq_record in SeqIO.parse("shortfasta.fsa", "fasta"): 
    s = seq_record    
    if seq_record.id.find('.') != -1:        
        seq_record.id = seq_record.id[:seq_record.id.find('.')]        
    if seq_record.id not in rejectlist_df['accession'].tolist():                       
        print("I didn't find these seq_records: ", seq_record.id)
        sequences.append(s)
    else:
        found.append(s)
    
    SeqIO.write(sequences, "stripped.fsa", "fasta")  
    SeqIO.write(found, "found.fsa", "fasta") 
    
subprocess.run(["ls", "-l"])