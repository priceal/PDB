# -*- coding: utf-8 -*-
"""
Spyder Editor

Calculates and plots as a heat map the contact matrix for input coordinates
from a protein structure.

input file must be a .npy file, containing the coordinates of the CA atoms.
input array must have have shape (N,3) where N=number of residues

"""
import numpy as np
import matplotlib.pyplot as plt

inputFile = 'data/2dpu_A.npz'
mapTitle = '2dpu_A backbone distance map'
axisLabel = 'residue'
colorMap = 'afmhot_r'
logorithmic = False

group = 0 # 0 for backbone and 1 for sidechain

###############################################################################
###############################################################################

structure=np.load(inputFile)
xyzRaw = structure['sc'] if group == 1 else structure['bb']
#xyz=xyzRaw.reshape((-1,3))
xyz= xyzRaw[ ~np.isnan(xyzRaw) ].reshape((-1,3))
'''
first calculate the distance map
use broadcasting to calculate map efficiently. in following, indices
are numbered from (1) fastest to slowest.
the shape of the N residue xyz array is A=(N(2),3(1))
the shape of the higher rank created array is B=(N(3),:(2),3(1))
then DMV = B-A. explanation: B's index(2) is stretched to N, 
so B=(N(3),N(2),3(1)). and then A is expanded to rank 3 and stretched to
A=(N(3), N(2), 3(1)). Note that if we consider this a matrix of 3-vectors
that is NxN, we have A has one column per coordinate (all same along
the new stretched slow 3-axis), and B has one row per coordinate (all same
along the new stretched slower 2-axis).
'''
DMV = xyz[:,np.newaxis,:] - xyz  # broadcasting creates matrix displacement 
DM = np.sqrt( (DMV*DMV).sum(axis=2) )  # sum along fast index

# create amninoacid sequence map
heatMapDict = {'variable': DM,
               'name': axisLabel,
               'title': mapTitle,
               'colormap': colorMap
              }
if logorithmic: heatMapDict[ 'variable' ] = np.log(heatMapDict[ 'variable' ])
        
###############################################################################
fig,ax=plt.subplots()
ax.imshow(heatMapDict['variable'])
ax.set_xlabel( heatMapDict['name'] +' 1' )
ax.set_ylabel( heatMapDict['name'] +' 2' )
ax.xaxis.set_label_position("top")
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
if heatMapDict['title']:
    ax.set_title( heatMapDict['title'] )
else:
    ax.set_title('distance matrix for '+ heatMapDict['name'])
ax.imshow( heatMapDict['variable'], cmap=heatMapDict['colormap'])
