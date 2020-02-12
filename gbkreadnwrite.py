# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:22:23 2020

@author: mcveigh
"""

from Bio import SeqIO

for seq_record in SeqIO.parse("test.gbk", "genbank"):
    print(seq_record.id)
    seq_record.description = seq_record.description + " hiho"
    print(seq_record.description)
    print(seq_record)
    #print(repr(seq_record.seq))
    #print(len(seq_record))