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

inputFile1 = 'data/1bdi_A.npz'
inputFile2 = 'data/1bdi_B.npz'
mapTitle = '1bdi A-B all atom distance map'
axisLabel = 'residue'
colorMap = 'afmhot_r'
logorithmic = False

group1 = 'sc' # 0 for backbone and 1 for sidechain
group2 = 'sc' # 0 for backbone and 1 for sidechain
cutoff = 10    # non-zero for a contact map with cutoff value

###############################################################################
###############################################################################

structure1=np.load(inputFile1)[group1]
structure2=np.load(inputFile2)[group2]
xyz1 = structure1[ ~np.isnan(structure1) ].reshape((-1,3))
xyz2 = structure2[ ~np.isnan(structure2) ].reshape((-1,3))

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

DMV = xyz1[:,np.newaxis,:] - xyz2[np.newaxis,:,:]   
DM = np.sqrt( (DMV*DMV).sum(axis=2) )  # sum along fast index

# now create the contact ap
if cutoff:
    DM = np.less(DM,cutoff)

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
