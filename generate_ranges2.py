# -*- coding: utf-8 -*-
"""
Created on Mon May 23 14:21:36 2022

@author: mcveigh
"""

import sys
import pandas as pd
import itertools

inputfile = sys.argv[1]
#outputfile = sys.argv[2]

acclist = pd.read_csv(inputfile, sep='\t', index_col=None, low_memory=False, header=None, usecols=[0], names=["accession"])
acclist.drop_duplicates(keep = 'first', inplace = True)
#acclist['accession'] = acclist['accession'].astype(str)
acclist.sort_values('accession', ascending=True, inplace=True)

acclist['letter'] = acclist['accession'].str.extract('([A-Za-z]+)')
acclist['number'] = acclist['accession'].str.replace('([A-Za-z]+)', '')
acclist['number'] = pd.to_numeric(acclist['number'])
#print(acclist)

prefixlist = acclist['letter']
prefixlist.drop_duplicates(keep = 'first', inplace = True)
list = prefixlist.tolist()
#print (list)


##else statement for the single is never true so singles are never printed
##logic error in dealing with the AF vs. multiple ranges in MT
            
first = []
last = []
single = []
number = []
prefix = []
count = int(1)
            
for x, row in acclist.iterrows():
    prefix = acclist.at[x, 'letter'] 
    single = acclist.at[x, 'accession']
    number = acclist.at[x, 'number']
    print("next row ", prefix, number, count)
        
    if count == 1:
        prefix = acclist.at[x, 'letter'] 
        single = acclist.at[x, 'accession']
        number = acclist.at[x, 'number']
        first = acclist.at[x, 'accession']
        count += 1
    if count > 1:
        if prefix == acclist.at[x, 'letter']:
            
        
    
    
    
    elif prefix != 0:
        if prefix == acclist.at[x, 'letter']:
            nextnumber = acclist.at[x, 'number']
            if nextnumber == number + 1:
                count += 1
                last == acclist.at[x, 'accession']
            elif nextnumber != number + 1:
                if count == 1:
                    print(single)
                elif count > 1:
                    print (first, last)
        elif prefix != acclist.at[x, 'letter']:    
            prefix = acclist.at[x, 'letter'] 
            single = acclist.at[x, 'accession']
            first = acclist.at[x, 'accession']
            number = acclist.at[x, 'number']
    else:
          print(single)




