#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 12:50:36 2024

@author: allen

creates a dataframe including pdbid and path to mmCIF file. also includes
a sort column to be used by a sorting script.

"""
import os
import pandas as pd

# inputs
structureDirectory = '../DATA/db/assemblies'
pdbCodeFile = ''  # can limit which processed with this, leave '' for not
maxNumber = 100000   # limit to first maxNumber ids   

# outputs
sortedFile = './csv/sort715.csv'   # leave ''  if you do not want to save                  

###########################################################################

# create list of file names
dirList = os.listdir(structureDirectory)

# load in pdb ids, use below if csv file with one column labeled 'pdbid'
if not pdbCodeFile:  # in not filtering - use all files
    pathList = [ os.path.join(structureDirectory,f) for f in dirList ]
    
else:    # read in pdbids to filter and apply
    if os.path.splitext(pdbCodeFile)[-1] == '.csv':  
        df = pd.read_csv(pdbCodeFile)
        pdbCodes=list(df['pdbid'])
    else:    # or this if a simple whitespace separated list of ids
        with open(pdbCodeFile) as f:
            fileRead=f.read()
        pdbCodes = fileRead.strip().split()    
        
    pdbCodes = [p.lower() for p in pdbCodes]
    fileList=[ p for p in dirList if p[:4] in pdbCodes]
    pathList = [ os.path.join(structureDirectory,f) for f in fileList ]

if maxNumber:
    pathList = pathList[:maxNumber]

pdbids = [ os.path.basename(p)[:4] for p in pathList ]

dataDf = pd.DataFrame({'pdbid':pdbids, 
                       'path':pathList, 
                       'sort':['u']*len(pdbids)}
                      )
    
if sortedFile:
    dataDf.to_csv(sortedFile, index=False)
    