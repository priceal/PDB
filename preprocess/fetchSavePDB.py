#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:44:14 2025

downloads PDBs and creates (if needed) directory to store

@author: allen
"""

pdbCodeFile = 'initialCodes10.csv'
dataDirectory = 'pdb10'
bioassembly = '1'
######################################################################
import os
import numpy as np
import pandas as pd
import pdbTools as pt
######################################################################
with open(pdbCodeFile) as f:
    pdbCodes = f.read().split('\n')

pdbCodes = [ cd.lower() for cd in pdbCodes ]

if not os.path.exists(dataDirectory):
    os.mkdir(dataDirectory)
    
for code in pdbCodes:
    data = pt.load(code,ba=bioassembly)
    if data:
        datalines = [ line+'\n' for line in data ]
        path=os.path.join(dataDirectory,code+'.pdb')
        with open(path,'w') as f:
            f.writelines(datalines)
   