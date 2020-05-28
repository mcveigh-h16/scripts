# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 10:07:23 2020

#Remove sequences with duplicate seqids and write a fasta file with the unique sequences only

@author: mcveigh
"""
import sys
inputfile = sys.argv[1]
outputfile = sys.argv[2]

from Bio import SeqIO

with open(outputfile, 'a') as outFile:
    record_ids = list()
    for record in SeqIO.parse(inputfile, 'fasta'):
        if record.id not in record_ids:
            record_ids.append( record.id)
            SeqIO.write(record, outFile, 'fasta')