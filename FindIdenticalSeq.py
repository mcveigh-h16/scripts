# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 09:58:43 2022

@author: mcveigh
"""

#Working trial to FeatureSlicer16S-2.py. 
#Find and identify identical and unique sequences from feature splice of 16S rRNA from genome 

unique = []

from Bio import SeqIO
records = list(SeqIO.parse("genome.out", "fasta"))
d = dict()
for record in records:
    if record.seq in d:
        d[record.seq].append(record)
    else:
        d[record.seq] = [record]
        unique.append(record)
for record.seq, record_set in d.items():
    print (record.seq + ': (' + str(len(record_set)) + ' sequence found)')
    for record in record_set:
        print ('identical seqs are: ' + record.description)
#print (unique.id)
for record in unique:
    print ('uniques are: ' + record.description)