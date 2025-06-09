#Plot showing the luminosity distance vs confidence area relationship for simulated BBH events, 
#color-coded based off of the nymber of detectors used to observe each event.

#Uses locally accessed gw simulation data: combined_bbh_O5.csv, 
#and NERSC-retrieved detector data: detector_data6.npy

############################################################################### Imports, Data, and Variables to Plot

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as matpatch
import math

events = np.loadtxt(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/combined_bbh_O5.csv", skiprows=1)
detectors = np.loadtxt(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/detector_data6.npy", dtype=str)

D = events[:,1]
A_70 = events[:,11]

simID = detectors[:, 2].astype(int)

############################################################################### Groups and Color-coding

sharedID = np.intersect1d(simID, events[:,0], return_indices=True)
colorMatch = detectors[sharedID[1], 1].astype(int)
colorMap = {1:'gold', 2:'green', 3:'blue', 4:'purple'}
opaqueMap = {1:1, 2:1, 3:1, 4:1}

############################################################################### Plotting

plt.figure(figsize=(5,5))
plt.scatter(D, A_70, c=[colorMap[n] for n in colorMatch], alpha=[opaqueMap[n] for n in colorMatch])

plt.xscale('log')
plt.yscale('log')
plt.xlabel('dL [Mpc]')
plt.ylabel('area (70% CL) [deg^2]')
plt.title('area (70% CL) vs dL')

patches = [ matpatch.Patch(color=colorMap[n], label=f'{n} Detectors') for n in colorMap ]
plt.legend(handles=patches, loc='lower right')

plt.show()