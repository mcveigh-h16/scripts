# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:32:45 2020

@author: mcveigh
"""
#Feature Slicer to extract 23S rRNA sequences from a larger file

import Bio
import sys
#from Bio import SeqIO, SeqFeature
import os

inputfile = sys.argv[1]
outputfile1 = sys.argv[2]

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
feature_count = 0
sequence = []  
for seq_record in SeqIO.parse(inputfile, "genbank"):    
    str_id = seq_record.id
    prot_id = []
    proteinacc = []
    counter = 0
    for feature in seq_record.features:
        if feature.type == 'CDS':
            seq_record.description = seq_record.annotations["organism"]
            protname = str(feature.qualifiers.get("product"))  
            if protname == "['cytochrome b']":                              
                counter += 1
                if counter == 1:
                    sequence_of_interest = feature.location.extract(seq_record).seq
                    #sub_record.id = str_id
                    defline = (seq_record.description + " cytochrome b")
                    sequence.append(SeqRecord(sequence_of_interest,description=defline,id=str_id))
                    feature_count += 1
                elif counter > 1:
                    print(str(str_id)+" has "+str(counter)+" cytochrome b features")
                else:
                    print(str(str_id)+" has no cytochrome b features")
#print(sequence)

SeqIO.write(sequence, outputfile1, "fasta") 

#for seq_record in SeqIO.parse(outputfile, "fasta"): 
#    seq = seq_record.seq.translate(table=5)
    #seq_record.id = "trans_"+ seq_record.id 
#    protein.append(seq)
    #seq_record.id = "trans_"+ seq_record.id
    #protein.append(seq)
#SeqIO.write(protein, "cytb_ex.aa", "fasta")
    
fh = open(outputfile1)
n = 0
for line in fh:
    if line.startswith(">"):
        n += 1
fh.close()
print("The number of sequences in the final output file is: ", n)
print ("FeatureSlicer script is done and output is saved in:", outputfile1)   
    
    

    