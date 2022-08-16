# -*- coding: utf-8 -*-
"""
Created on Fri May 27 09:46:50 2022

@author: mcveigh
"""

##Script takes an itemized list of accessions and compacts the list into hyphenated
##ranges wherever possible. Singles are listed and duplicates are dropped. 
##Accessions with .version will generate an error. Add function to strip .version if needed
##RefSeq accesssion that are all zeros will also cause and error as these are not integers


#all_accs = {'AF123456', 'AF123457', 'MT082704', 'MT082705', 'MT082706',
#    'MT082707', 'MT082708', 'MT082709', 'MT082711', 'MT082712', 'MT082713',
#    'MT082714', 'MT082715', 'MT082716', 'MT082717', 'MT082718', 'MT082719',
#    'MT082720', 'MT157270', 'U12345', 'U12346', 'MT082716',}

import sys
import re
inputfile = sys.argv[1]

with open(inputfile) as f:
    #all_accs = f.readlines.strip()
    readline=f.read().splitlines()
all_accs = readline
final_list = []
acc_prefix = []
acc_digit = []
prev_acc_prefix = []
prev_acc_digit = []

prev_acc = None 
for acc in sorted(all_accs):
    if not prev_acc:
        new_accs = [acc]
    #elif prev_acc[:2] == acc[:2] and int(prev_acc[2:]) + 1 == int(acc[2:]):
    #elif re.findall(r'[A-Za-z]+', prev_acc) == re.findall(r'[A-Za-z]+', acc) and int(re.findall(r'(\d+)', prev_acc)) + 1 == int(re.findall(r'(\d+)', acc)):
    else:
        acc_prefix = re.findall(r'[A-Za-z]+', acc)
        acc_digit = re.findall(r'(\d+)', acc)
        prev_acc_prefix = re.findall(r'[A-Za-z]+', prev_acc)
        prev_acc_digit = re.findall(r'(\d+)', prev_acc)
        acc_digit = list(map(int, acc_digit))
        prev_acc_digit = list(map(int, prev_acc_digit))
        if prev_acc_prefix == acc_prefix and prev_acc_digit + 1 == acc_digit:
            new_accs.append(acc)
 
        else:
            final_list.append(f'{new_accs[0]}-{new_accs[-1]}') if len(new_accs)>1 else final_list.append(new_accs[0])
            new_accs = [acc]
    prev_acc = acc 
final_list.append(f'{new_accs[0]}-{new_accs[-1]}') if len(new_accs)>1 else final_list.append(new_accs[0])    
print(final_list)