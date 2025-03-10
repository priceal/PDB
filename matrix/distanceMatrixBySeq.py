# -*- coding: utf-8 -*-
"""
Spyder Editor

Calculates and plots as a heat map the distance matrix for two structures.
input files must be .npz files, with keys that match those in inputs below.

coordinate arrays must have shape (N,maxResSize,3) where N is number of 
residues, maxResSize is number of atoms/res. Therefore, indices indicate 
(residue #, atom #, coordinate).

"""
import os
import numpy as np
import matplotlib.pyplot as plt

# input structure files, format of dictionary entries:
# label: [ file, group (set) ]
# keys:   proteins:   bb = backbone, sc = sidechain
#         dna:        ph = phosphate, rb = ribose, ba = base
# label is used in labeling axes 
inputStructures = { 'A':  ['1b72_A.npz', {'sc'}],
                    'D':  ['1b72_D.npz', {'ba'}]
                  }
# optional structure file directory, can leave undefined '' or '.'
fileDirectory = 'data'

cutoff = 7  # non-zero for a contact map with cutoff value
mapTitle = 'contact map'
colorMap = 'OrRd'
logorithmic = False

###############################################################################
###############################################################################
# load and parse structure files. store in dictionary
structure = {}
for k,(file,group) in inputStructures.items():
    if fileDirectory: file=os.path.join(fileDirectory,file)
    tempStructure = np.load(file)
    coords = []
    for g in group:
        coords.append(tempStructure[g])
    structure[k] = np.concatenate(coords,axis=1)

'''
first calculate the distance map. use broadcasting to calculate map efficiently. 
the shape of structure1 array is S1=(N1,R1,3) and that of structure2 is
S2=(N2,R2,3). New axes are created to increase the ranks to that of
S1*=(N1,1,R1,1,3) and S2*=(1,N2,1,R2,3). Broadcasting will stretch the axes 
of dim=1 thus: S1* -> (N1,(N2),R1,(R2),3) and S2* -> ((N1),N2,(R1),R2,3).
The shape of the resulting difference will be 
S1*-S2* = (N1(0),N2(1),R1(2),R2(3),3(4)). where the indices meanings are:
    
    (0) residue from structure1
    (1) residue from structure2
    (2) atom from structure1 residue
    (3) atom from structure2 residue
    (4) coordinate x, y, or z

To get distance between atoms, take the sqrt of sums along axis 4.
To get minimum distance between residues,take minimum along axes 2 and 3.
'''

# the vector displacement matrix: shape = (N1,N2,R1,R2,3)
label2, structure2 = structure.popitem() # pull out coordinates LIFO
label1, structure1 = structure.popitem()
vDM = structure1[          :, np.newaxis,           :, np.newaxis, : ] - \
      structure2[ np.newaxis,           :, np.newaxis,          :, : ]
      
# the distance matrix: shape = (N1,N2,R1,R2)
aDM = np.sqrt( (vDM*vDM).sum( axis = 4 ) )  # 4 is axis of coordinates (xyz)
DM = np.nanmin( aDM, axis=(2,3))    # (2,3) are axes of atom # w/in residues

# now create the contact map if asked
if cutoff: DM = np.less(DM,cutoff)

# create amninoacid sequence map
heatMapDict = {'variable': DM,
               'title': mapTitle,
               'xlabel': label2,
               'ylabel': label1,
               'colormap': colorMap
              }
if logorithmic: heatMapDict[ 'variable' ] = np.log(heatMapDict[ 'variable' ])
        
###############################################################################
fig,ax=plt.subplots()
ax.imshow(heatMapDict['variable'])
ax.set_xlabel( heatMapDict['xlabel'] )
ax.set_ylabel( heatMapDict['ylabel'] )
ax.xaxis.set_label_position("top")
ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
if heatMapDict['title']:
    ax.set_title( heatMapDict['title'] )
else:
    ax.set_title('distance matrix for '+ heatMapDict['name'])
ax.imshow( heatMapDict['variable'], cmap=heatMapDict['colormap'])
