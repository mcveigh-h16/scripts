# -*- coding: utf-8 -*-
"""
Created on Wed May 13 14:15:50 2020

@author: mcveigh
"""
#splits a fasta file into multiple files with 1000 sequences each, assigning each file a new name


def batch_iterator(iterator, batch_size):
 
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                 entry = iterator.__next__()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

from Bio import SeqIO

record_iter = SeqIO.parse(open("LatestFungiITSb.out"),"fasta")
for i, batch in enumerate(batch_iterator(record_iter, 1000)):
    filename = "group_%i.fasta" % (i + 1)
    with open(filename, "w") as handle:
        count = SeqIO.write(batch, handle, "fasta")
    print("Wrote %i records to %s" % (count, filename))
    
