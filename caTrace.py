#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 12:50:36 2024

@author: allen
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

ca=loadCoords('1aki.pdb',ca=True)
cas = np.array(ca).T

fig=plt.figure()
ax=fig.add_subplot(projection='3d')
ax.plot3D(cas[0],cas[1],cas[2],'-')