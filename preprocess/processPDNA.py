#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 14:44:14 2025

downloads PDBs and creates (if needed) directory to store

@author: allen
"""

pdbCodeFile = 'initialCodes.csv'
dataDirectory = 'pdna'

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
    data = pt.load(code)
    if data:
        df = pt.makeDataFrame( data )
        path=os.path.join(dataDirectory,code+'.csv')
        df.to_csv(path)
   