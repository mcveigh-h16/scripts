# -*- coding: utf-8 -*-
"""
Created on Fri May 27 09:46:50 2022

@author: mcveigh
"""

##Script takes an itemized list of accessions and compacts the list into hyphenated
##ranges wherever possible. Singles are listed and duplicates are dropped. 
##Versions are stripped from any accession.version
##RefSeq accesssion that are all zeros will also cause and error as these are not integers
##RefSeq or WGS require editting the :2 or 2: fro the elif statement


#all_accs = {'AF123456', 'AF123457', 'MT082704', 'MT082705', 'MT082706',
#    'MT082707', 'MT082708', 'MT082709', 'MT082711', 'MT082712', 'MT082713',
#    'MT082714', 'MT082715', 'MT082716', 'MT082717', 'MT082718', 'MT082719',
#    'MT082720', 'MT157270', 'U12345', 'U12346', 'MT082716',}

import sys
inputfile = sys.argv[1]

with open(inputfile) as f:
    #all_accs = f.readlines.strip()
    readline=f.read().splitlines()
all_accs = readline
final_list = []
prev_acc = None 
for acc in sorted(all_accs):
    acc = acc.split('.')[0]
    if not prev_acc:
        new_accs = [acc]
    elif prev_acc[:2] == acc[:2] and int(prev_acc[2:]) + 1 == int(acc[2:]):
        new_accs.append(acc)
    else:
        final_list.append(f'{new_accs[0]}-{new_accs[-1]}') if len(new_accs)>1 else final_list.append(new_accs[0])
        new_accs = [acc]
    prev_acc = acc 
final_list.append(f'{new_accs[0]}-{new_accs[-1]}') if len(new_accs)>1 else final_list.append(new_accs[0])    
print(final_list)