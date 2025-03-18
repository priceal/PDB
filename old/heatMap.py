# -*- coding: utf-8 -*-
"""
Spyder Editor

Plots a heat map that has been stored in a text file. the input file must
be a 2D array with white space separated values

"""
inputFile = 'f/1gzm.dM'
mapTitle = 'Fortran'
axisLabel = 'residue'
colorMap = 'afmhot_r'
logorithmic = False

###############################################################################
###############################################################################
import numpy as np
import matplotlib.pyplot as plt

heatMap = np.loadtxt(inputFile)

# create amninoacid sequence map
heatMapDict = {'variable': heatMap,
               'name': axisLabel,
               'title': mapTitle,
               'colormap': colorMap
              }
if logorithmic: heatMapDict['variable'] = np.log(heatMapDict['variable'])

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
