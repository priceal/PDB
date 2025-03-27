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
import matplotlib.pyplot as plt
from distanceMatrixBySeq_func import loadStructureAndSequence, \
                                     distanceMatrix, \
                                     displayHeatMap

'''
structure file information --- format of dictionary entries:
    label : [ [files], {groups} ]

label = used to label axes
[files] = list of .npz files for this set of structure
{groups} = set of subgroups of atoms to consider
allowed groups ...
           proteins:   bb = backbone, sc = sidechain
           dna:        ph = phosphate, rb = ribose, ba = base

N.B. all .npz files under same label MUST have structure arrays of same shape!           
'''

inputStructures = { 'protein':  [ ['1le5_A.npz','1le5_B.npz'],    \
                                 {'sc'} ],
                    'dna':      [ ['1le5_C.npz','1le5_D.npz'],    \
                                 {'ba'} ]
                  }
 
# optional structure file directory, can leave undefined '' or '.'
fileDirectory = '../DATA/PDB/npz'

cutoff = 7  # required for contact map
mapTitle = 'contact map by residue type: sidechain-base'
colorMap = 'OrRd'

###############################################################################
###############################################################################
# vocabulary for residue types --- assumes usual dna/protein residues
dnaResidues = 'ACGT'
proteinResidues = 'ARNDCEQGHILKMFPSTWYV'

   
structure = {}; sequence = {}; polymerType = {}    
for k,(files,group) in inputStructures.items():
    if group&{'sc','bb'}:
        polymerType[k] = 'protein'
    elif group&{'ph','rb','ba'}:
        polymerType[k] = 'dna'    
#    print(k,files,group)
    structure[k], sequence[k] = loadStructureAndSequence(
                                files, group, 
                                Directory=fileDirectory
                                )
#    print(structure,sequence)
# the vector displacement matrix: shape = (N1,N2,R1,R2,3)
label2, structure2 = structure.popitem() # pull out coordinates LIFO
label1, structure1 = structure.popitem()
DM = distanceMatrix(structure1, structure2)    
CM = np.less(DM,cutoff)             # applies cutoff

# perform "similarity" transform on sequence based contact map to 
# generate the residue type based map
label2, seq2 = sequence.popitem() # pull out sequences LIFO
label1, seq1 = sequence.popitem()
resCM = np.linalg.multi_dot([seq1.T,CM,seq2])  # the magic step

label2, type2 = polymerType.popitem() # pull out types LIFO
label1, type1 = polymerType.popitem()

if type1=='dna':
    yticks=dnaResidues
else:
    yticks=proteinResidues
if type2=='dna':
    xticks=dnaResidues
else:
    xticks=proteinResidues

displayHeatMap( {'variable': resCM,
                 'title': mapTitle,
                 'xlabel': label2,
                 'ylabel': label1,
                 'colormap': colorMap,
                 'logarithmic': False
                 }, xticks=xticks, yticks=yticks
              )
   
