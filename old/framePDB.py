#!\usr\bin\env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 07:59:14 2024

@author: allen
"""
import pandas as pd

source = 'pdb/5h9g.pdb'

columns = {'record': (1,4),
           'number': (7,11),
           'name': (13,16),
           'alt': (17,17),
           'res': (18,20),
           'chain': (22,22),
           'seqnum': (23,26),
           'insert': (27,27),
           'x': (31,38),
           'y': (39,46),
           'z': (47,54),
           'occ': (55,60),
           'B': (61,66),
           'segid': (73,76),
           'element': (77,78),
           'charge': (79,80)
            }
data = dict( zip( columns.keys(), [ [] for i in range(len(columns)) ] ) ) 
lines=[]
with open(source,'r') as file:
    for line in file:
        lineSplit = line.split()
        if lineSplit[0].strip()=='ATOM':
            lines.append(line)

for line in lines:
    for item in columns:
        data[item].append( line[ columns[item][0]-1: columns[item][1]].strip() )

df = pd.DataFrame( data )
        
        
