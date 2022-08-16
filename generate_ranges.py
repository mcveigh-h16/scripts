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

count = []

for i in list:
    first = []
    last = []
    single = []
    number = []
    count = 0
    accset = acclist[acclist['letter'] == i]
    for x, row in accset.iterrows():
        single = accset.at[x, 'accession']
        print(i, single, count)
        if count == 0:
            first = accset.at[x, 'accession']    
            number = accset.at[x, 'number']
            print("zero count is ", number, count)
            count = 1
            
        elif count == 1:
            nextnumber = accset.at[x, 'number'] 
            print("one count is ", count, number, nextnumber)
            if nextnumber == number + 1:
                last = accset.at[x, 'accession']
                number = nextnumber
                print("lap", count, first, last, number)
            elif nextnumber != number + 1:
                count = 0
                print("exit", first, last)
                
            else:
                print("single", single)
    print("range", first, "-", last)


##else statement for the single is never true so singles are never printed
##logic error in dealing with the AF vs. multiple ranges in MT
            
            





