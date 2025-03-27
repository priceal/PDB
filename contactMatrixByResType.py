# -*- coding: utf-8 -*-
"""
Spyder Editor

Calculates and plots as a heat map the contact matrix for two structures,
with contacts aggregated by residue type.
input files must be .npz files, with keys that match those in inputs below.

coordinate arrays must have shape (N,maxResSize,3) where N is number of 
residues, maxResSize is number of atoms/res. Therefore, indices indicate 
(residue #, atom #, coordinate).

"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from distanceMatrixBySeq_func import loadStructureAndSequence, \
                                     distanceMatrix, \
                                     displayHeatMap

'''
           proteins:   bb = backbone, sc = sidechain
           dna:        ph = phosphate, rb = ribose, ba = base         
'''



# inputs
structureListFile = './csv/sort713.csv'
fileDirectory = '../DATA/PDB/npz'
groups = { 'protein': {'bb'} ,'dna': {'ba'} }
maxNumber = 800

# parameters
cutoff = 4.0  # required for contact map
mapTitle = 'contact map by residue type'
colorMap = 'OrRd'

###############################################################################
###############################################################################
# vocabulary for residue types --- assumes usual dna/protein residues
dnaResidues = 'ACGT'
proteinResidues = 'ARNDCEQGHILKMFPSTWYV'

structureDf = pd.read_csv(structureListFile)

resCM = np.zeros((20,4))
for entry in structureDf[:maxNumber].itertuples():
    
    chains= entry.protein.split(',')
    files=[entry.pdbid+'_'+c+'.npz' for c in chains]
    structureProtein,seqProtein = loadStructureAndSequence(
                                files, groups['protein'], 
                                Directory=fileDirectory
                                )
    chains= entry.dna.split(',')
    files=[entry.pdbid+'_'+c+'.npz' for c in chains]
    structureDNA,seqDNA = loadStructureAndSequence(
                                files, groups['dna'], 
                                Directory=fileDirectory
                                )
    #    print(structure,sequence)
    # the vector displacement matrix: shape = (N1,N2,R1,R2,3)
    DM = distanceMatrix(structureProtein, structureDNA)    
    CM = np.less(DM,cutoff)             # applies cutoff
    
    # perform "similarity" transform on sequence based contact map to 
    # generate the residue type based map
    resCM += np.linalg.multi_dot([seqProtein.T,CM,seqDNA])  # the magic step
    
displayHeatMap( {'variable': resCM,
                 'title': mapTitle,
                 'xlabel': 'dna '+','.join(groups['dna']),
                 'ylabel': 'protein '+','.join(groups['protein']),
                 'colormap': colorMap,
                 'logarithmic': False
                 }, xticks=dnaResidues, yticks=proteinResidues   
              )
   
