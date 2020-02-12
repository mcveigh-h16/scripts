import pandas as pd
import Bio
import os
from os import path


def main():
  # Open the files for writing and create it if it doesn't exist
  f1 = open("found.fsa","w+")
  origfasta = open("shortfasta.txt", "r")
  rejectlist = open("rejectlist2.txt", "r")
  f4 = open("stripped.fsa", "w+")

  # close the files when done
  f1.close()
  origfasta.close()
  rejectlist.close()
  f4.close()

  rejectlist = (r'/net/snowman/vol/export2/mcveigh/python/rejectlist2.txt')
  rejectlist_df = pd.read_csv(rejectlist, sep='\t', index_col=None, low_memory=False, header=None, names=["accession", "type", "reason"])
  rejectlist = rejectlist_df['accession']
  print (rejectlist) 

  from Bio import SeqIO
  sequences = [] 
  found = []
  #print(rejectlist_df['accession'])
  for seq_record in SeqIO.parse("shortfasta.txt", "fasta"): 
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
    
if __name__ == "__main__":
    main()