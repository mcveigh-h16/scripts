# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:43:02 2022

@author: mcveigh
"""
import sys
import re

p = re.compile(r'([a-zA-Z]+)[_]*(\d+)')
accs = ['EAT123456.7', 'EAT123456', 'EA_123456.7']
for acc in accs:
    acc = acc.split('.')[0]
    x = p.search(acc)
    #for i in [1, 2]:
        #print(x.group(i))
    print(x[1])
    print(x[2])