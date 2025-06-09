#Plot showing the luminosity distance vs confidence area relationship for simulated BBH events, 

#COLOR CODING UPDATED: now color-codes based off of categories defined using search parameters (below) 

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
detectorData = np.loadtxt(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/detector_data6.npy", dtype=str)

D = events[:,1]
A_70 = events[:,11]

colors = ['gray' for _ in events]
opaques = [1 for _ in events]

simID = detectorData[:, 2].astype(int)
sharedID = np.intersect1d(simID, events[:,0], return_indices=True)

colorMatch = detectorData[sharedID[1], 1].astype(int)
detectorMatch = detectorData[sharedID[1], 0].astype(str)

############################################################################### search criteria

def colorPop(detectors, condition, number, color, opaque):
    detectorsInterest = detectors
    conditionInterest = condition
    numberInterest = number
    color = color
    opaque = opaque

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
    
    for n in range(len(detectorMatchBools)):
        if detectorMatchBools[n]:
            colors[n] = color
            opaques[n] = opaque
        
############################################################################### Population definition and color-coding
param = [   [set(['H1','L1']), 'all', set([2,3,4]), 'purple', 0.5],
            [set(['H1','L1','K1','V1']), 'any', set([1]), 'gold', 0.5],
            [set(['H1','L1']), 'uniques', set([2,3]), 'green', 0.5],
            [set(['K1','V1']), 'all', set([2]), 'blue', 0.5]  ]
        
for pop in param: colorPop(*pop)

############################################################################### Plotting

plt.figure(figsize=(5,5))
plt.scatter(D, A_70, c=colors, alpha=opaques)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('dL [Mpc]')
plt.ylabel('area (70% CL) [deg^2]')
plt.title(f"""area (70% CL) vs dL""")

patches = [ matpatch.Patch(color=param[n][3], label=f'{param[n][2]} detectors including {param[n][1]} of {param[n][0]}') for n in range(len(param)) ]
plt.legend(handles=patches, loc='lower right')

plt.show()