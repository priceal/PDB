#!/home/allen/anaconda3/bin/python
# -*- coding: utf-8 -*-
 
"""
Spyder Editor

Plots a heat map that has been stored in a text file. the input file must
be a 2D array with white space separated values

"""

mapTitle = 'Heat Map'
axisLabel = 'residue'
colorMap = 'afmhot_r'
logorithmic = False

###############################################################################
###############################################################################
import numpy as np
import matplotlib.pyplot as plt
import sys

# parameter and data input from standard input
heatMap = []
for line in sys.stdin:
    data = line.split()
    heatMap.append( [ float(s) for s in data ] )
heatMap = np.array(heatMap)         

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
plt.show()
print('end of program')
