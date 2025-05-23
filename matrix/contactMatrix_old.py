# -*- coding: utf-8 -*-
"""
Spyder Editor

Calculates and plots as a heat map the distance matrix for the CA coordinates
of a protein structure.

input file must be a .npy file, containing the coordinates of the CA atoms.
input array must have have shape (N,3) where N=number of residues

"""

inputFile = '5MBA_ca.npy'
mapTitle = '5MBA contact map'
axisLabel = 'residue'
colorMap = 'afmhot_r'
logorithmic = False
contactDistance = 7.0

###############################################################################
###############################################################################
import numpy as np
import matplotlib.pyplot as plt

xyz=np.load(inputFile)

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

# now create the contact ap
CM = np.less(DM,contactDistance)

#code snippet to plot a heat map
heatMapDict = {'variable': CM,
               'name': axisLabel,
               'title': mapTitle,
               'colormap': colorMap,
               'logorithmic': False
              }
###############################################################################
if heatMapDict['logorithmic']: 
    heatMapDict['variable'] = np.log(heatMapDict['variable'])
fig,ax=plt.subplots()
ax.imshow(heatMapDict['variable'])
ax.set_xlabel( heatMapDict['name'] +' 1' )
ax.set_ylabel( heatMapDict['name'] +' 2' )
ax.xaxis.set_label_position("top")
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
if heatMapDict['title']:
    ax.set_title( heatMapDict['title'] )
else:
    ax.set_title('heat map for '+ heatMapDict['name'])
ax.imshow( heatMapDict['variable'], cmap=heatMapDict['colormap'])
