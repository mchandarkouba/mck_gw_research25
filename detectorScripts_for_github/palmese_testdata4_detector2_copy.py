#Plot showing the luminosity distance vs confidence area relationship for simulated BBH events, 
#color-coded based off of the nymber of detectors used to observe each event.
#includes a searching method for "detectors of interest", 
#given that ALL detectors of interest (at least) were used to observe the event.

#Uses locally accessed gw simulation data: combined_bbh_O5.csv, 
#and NERSC-retrieved detector data: detector_data6.npy

############################################################################### Imports, data, and Plotting variables

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as matpatch
import math

events = np.loadtxt(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/combined_bbh_O5.csv", skiprows=1)
detectors = np.loadtxt(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/detector_data6.npy", dtype=str)

D = events[:,1]
A_70 = events[:,11]

############################################################################### Search criteria and color-coding

detectorsInterest = {'L1', 'H1'}
colorMap = {True:{1:'gold', 2:'green', 3:'blue', 4:'purple'}, 
            False:{1:'gray', 2:'gray', 3:'gray', 4:'gray'}}
opaqueMap = {True:1, False:0.25}

simID = detectors[:, 2].astype(int)
sharedID = np.intersect1d(simID, events[:,0], return_indices=True)

colorMatch = detectors[sharedID[1], 1].astype(int)
detectorMatch = detectors[sharedID[1], 0].astype(str)
detectorMatchBools = [all((d in elem.split(',')) for d in detectorsInterest) for elem in detectorMatch]

colors = [ colorMap[detectorMatchBools[n]][colorMatch[n]] for n in range(np.size(colorMatch))]
opaques = [ opaqueMap[detectorMatchBools[n]] for n in range(np.size(colorMatch))]

############################################################################### Plotting

plt.figure(figsize=(5,5))
plt.scatter(D, A_70, c=colors, alpha=opaques)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('dL [Mpc]')
plt.ylabel('area (70% CL) [deg^2]')
plt.title(f'area (70% CL) vs dL for {detectorsInterest} detectors')

patches = [ matpatch.Patch(color=colorMap[True][n], label=f'{n} Detectors') for n in colorMap[True] ]
plt.legend(handles=patches, loc='lower right')

plt.show()