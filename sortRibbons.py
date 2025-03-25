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
    ax.set_aspect('equal')
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
    
    sortedFile = 'filteredSort_20250324.csv'   # leave ''  if you do not want to save                  

    ###########################################################################

    dataDf = pd.read_csv( sortedFile )
    parser = MMCIFParser(QUIET=True)
    for i in dataDf.index:
        
        if dataDf.at[i,'sort'] == 'u': # only sort if unsorted
            structure = parser.get_structure('structure', dataDf.at[i,'path'] )
            
            # list number of models, and number of chains in model 0
            print(f'{i}/{len(dataDf)}')
            print(dataDf.at[i,'pdbid'],'-',len(structure),'model(s)', end='\n\t\t')
            chains = []
            model = structure[0]
            for chain in model:
                chains.append( chain.get_id() )
            print(len(model),'chain(s) in model 0:',chains, end='\n\t\t' )
            cas, phs = extractCAP( model )
            plotCAP(cas,phs)
            plt.show(block=True)
            selection=input( '[s]ave, [d]iscard, [f]or latter review, or [E]xit: ')
            if selection:   #skip if no letter is input
                if selection=='E':   # break on 'E'xit
                    break
                dataDf.at[i,'sort']=selection
        
    if sortedFile:
        dataDf.to_csv(sortedFile)