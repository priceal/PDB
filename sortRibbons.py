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
import numpy as  np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Bio.PDB import MMCIFParser, is_aa

###############################################################################
def extractCAP( model ):
    caList=[]; phList=[]
    for chain in model:
        ca=[]; ph=[]
        for res in chain:
            if is_aa(res):
                try:
                    ca.append(res['CA'].get_coord())
                except:
                    pass
            else:
                try:
                    ph.append(res['P'].get_coord())
                except:
                    pass
        if ca:
            caList.append( np.array(ca) )
        if ph:
            phList.append( np.array(ph) )
    return caList, phList
###############################################################################

def plotCAP( caList, phList ):
    
    fig=plt.figure(figsize=(10,10))
    ax=fig.add_subplot(projection='3d')
    for ca in caList:
        ca=ca.T
        ax.plot3D(ca[0],ca[1],ca[2],'-')
        
    for ph in phList:
        ph=ph.T
        ax.plot3D(ph[0],ph[1],ph[2],'o-')
    return 

###############################################################################

###############################################################################
############################# main ############################################
###############################################################################

if __name__ == "__main__":
    
    structureDirectory = '../DATA/db/assemblies'
    pdbCodeFile = '../db/pdbListAll.csv'       # file containing pdb ids
    maxNumber = 20   # limit to first maxNumber ids   
    sortedFile = 'sorted.csv'   # leave ''  if you do not want to save                  

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
    
    # create list of files paths
    dirList = os.listdir(structureDirectory)
    fileList=[ p for p in dirList if p[:4] in pdbCodes]
    pdbids = [ f[:4] for f in fileList]
    pathList = [ os.path.join(structureDirectory,f) for f in fileList ]
    
    dataDf = pd.DataFrame({'pdbid':pdbids, 
                           'path':pathList, 
                           'sort':['u']*len(pdbids)}
                          )
    parser = MMCIFParser(QUIET=True)
    for i in dataDf.index:
        
        structure = parser.get_structure('structure', dataDf.at[i,'path'] )
        
        # list number of models, and number of chains in model 0
        print(dataDf.at[i,'pdbid'],'-',len(structure),'model(s)')
        chains = []
        model = structure[0]
        for chain in model:
            chains.append( chain.get_id() )
        print(len(model),'chain(s) in model 0:',chains)
        cas, phs = extractCAP( model )
        plotCAP(cas,phs)
        plt.show(block=True)
        selection=input( 'select:  [s]ave  [d]iscard  [f]or latter review  :')
        if selection:
            dataDf.at[i,'sort']=selection
        
    if sortedFile:
        dataDf.to_csv(sortedFile)
        