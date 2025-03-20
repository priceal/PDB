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
    
    fig=plt.figure()
    ax=fig.add_subplot(projection='3d')
    for ca in caList:
        ca=ca.T
        ax.plot3D(ca[0],ca[1],ca[2],'-')
        
    for ph in phList:
        ph=ph.T
        ax.plot3D(ph[0],ph[1],ph[2],'o-')

###############################################################################

###############################################################################
############################# main ############################################
###############################################################################

if __name__ == "__main__":
    
    structureDirectory = '../DATA/db/assemblies'
    structureFile = '1y1w-assembly1.cif'
    
    ###########################################################################
      
    # assume pdb id code given by first 4 characters of file name
    code = structureFile[:4]
    
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(code, \
                            os.path.join(structureDirectory,structureFile))
    
    
    # list number of models, and number of chains in model 0
    print(len(structure),'model(s)')
    chains = []
    model = structure[0]
    for chain in model:
        chains.append( chain.get_id() )
    print(len(model),'chain(s) in model 0:',chains)
    
    cas, phs = extractCAP( model )
    plotCAP( cas, phs )
    
    