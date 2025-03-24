#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 12:50:36 2024

@author: allen

loads and plots protein CA ribbon and DNA backbone trace for proteib-DNA
complex. each chain is a different color, and DNA phosphates are represented
by circles.

"""
import os
import pandas as pd

structureDirectory = '../DATA/db/assemblies'
pdbCodeFile = '../db/filtered.csv'       # file containing pdb ids
maxNumber = 1000   # limit to first maxNumber ids   
sortedFile = 'filteredSort_20250324.csv'   # leave ''  if you do not want to save                  

###########################################################################
# load in pdb ids, use below if csv file with one column labeled 'pdbid'
if os.path.splitext(pdbCodeFile)[-1] == '.csv':  
    df = pd.read_csv(pdbCodeFile)
    pdbCodes=list(df['pdbid'])
else:    # or this if a simple whitespace separated list of ids
    with open(pdbCodeFile) as f:
        fileRead=f.read()
    pdbCodes = fileRead.strip().split()
    
pdbCodes = pdbCodes[:maxNumber]   # limit the number
pdbCodes = [p.lower() for p in pdbCodes]

# create list of files paths
dirList = os.listdir(structureDirectory)
fileList=[ p for p in dirList if p[:4] in pdbCodes]
pdbids = [ f[:4] for f in fileList]
pathList = [ os.path.join(structureDirectory,f) for f in fileList ]

dataDf = pd.DataFrame({'pdbid':pdbids, 
                       'path':pathList, 
                       'sort':['u']*len(pdbids)}
                      )
    
if sortedFile:
    dataDf.to_csv(sortedFile)
    