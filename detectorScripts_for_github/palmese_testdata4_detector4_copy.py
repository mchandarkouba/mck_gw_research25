#Plot showing the luminosity distance vs confidence area relationship for simulated BBH events, 
#color-coded based off of the nymber of detectors used to observe each event.
#includes a searching method for "detectors of interest" given various conditions:
    #any: at least detector within "detectors of interest" must be used
    #all: all detectors within "detectors of interest" must be used
    #none: no detectors within "detectors of interest" must be used
    #unique: exactly one detector from "detectors of interest" must be used
#given the number of detectors used to observe an event is within the "number of interest" set


#Uses locally accessed gw simulation data: combined_bbh_O5.csv, 
#and NERSC-retrieved detector data: detector_data6.npy

############################################################################### Imports, data, scatterplot variables

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as matpatch
import math

events = np.loadtxt(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/combined_bbh_O5.csv", skiprows=1)
detectors = np.loadtxt(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/detector_data6.npy", dtype=str)

D = events[:,1]
A_70 = events[:,11]

############################################################################### Search conditions and color-coding

detectorsInterest = set(['H1', 'L1']) # H1, L1 and K1, V1
conditionInterest = 'any' #'any', 'all', 'none', 'uniques'
numberInterest = set([3])

colorMap = {True:{1:'gold', 2:'green', 3:'blue', 4:'purple'}, 
            False:{1:'gray', 2:'gray', 3:'gray', 4:'gray'}}
opaqueMap = {True:1, False:0.25}

simID = detectors[:, 2].astype(int)
sharedID = np.intersect1d(simID, events[:,0], return_indices=True)

colorMatch = detectors[sharedID[1], 1].astype(int)
detectorMatch = detectors[sharedID[1], 0].astype(str)

detectorMatchBools = []
for m in detectorMatch:
    ofInterest = True
    m = set(m.split(','))
    n = detectorsInterest
 
    if conditionInterest=='any': ofInterest = n&m!=set()
    if conditionInterest=='all': ofInterest = n|m==m
    if conditionInterest=='none': ofInterest = n&m==set()
    if conditionInterest=='uniques': ofInterest = len(n&m)==1    
    
    ofInterest &= len(m) in numberInterest
    detectorMatchBools.append(ofInterest)

colors = [ colorMap[detectorMatchBools[n]][colorMatch[n]] for n in range(np.size(colorMatch))]
opaques = [ opaqueMap[detectorMatchBools[n]] for n in range(np.size(colorMatch))]

############################################################################### Plotting

plt.figure(figsize=(5,5))
plt.scatter(D, A_70, c=colors, alpha=opaques)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('dL [Mpc]')
plt.ylabel('area (70% CL) [deg^2]')
plt.title(f"""area (70% CL) vs dL: 
          events observed using {numberInterest} detectors including {conditionInterest} of {detectorsInterest} highlighted""")

patches = [ matpatch.Patch(color=colorMap[True][n], label=f'{n} Detectors') for n in colorMap[True] ]
plt.legend(handles=patches, loc='lower right')

plt.show()