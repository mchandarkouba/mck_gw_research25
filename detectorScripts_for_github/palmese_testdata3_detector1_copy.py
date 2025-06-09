#Test script analyzing the number of detectors used to measure each event,
#though plotting the wrong relationship...see later versions

#Uses locally accessed gw simulation data: combined_bbh_O5.csv, 
#and NERSC-retrieved detector data: detector_data6.npy

############################################################################### Imports, constants, and scatterplot variables

import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.patches as matpatch

events = np.loadtxt(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/combined_bbh_O5.csv", skiprows=1)
detectors = np.loadtxt(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/detector_data6.npy", dtype=str)

rho_agn = 10**(-4.75) # 1/MpC^3 #volume density of AGN in space
rho_event = 1e-4 #expected proportion of AGN flares associated with BBH mergers

V_90 = events[:,18]
M = events[:,5] + events[:,6]
#log_AGN = np.log(V_90 * rho_agn * rho_event)

simID = detectors[:, 2].astype(int)

############################################################################### Color-coding

sharedID = np.intersect1d(simID, events[:,0], return_indices=True)
colorMatch = detectors[sharedID[1], 1].astype(int)
colorMap = {1:'gold', 2:'green', 3:'blue', 4:'purple'}
opaqueMap = {1:1, 2:1, 3:1, 4:1}

############################################################################### Plotting

plt.figure(figsize=(5,5))
plt.scatter(V_90, M, c=[colorMap[n] for n in colorMatch], alpha=[opaqueMap[n] for n in colorMatch])

#plt.colorbar(label='log10(N_AGN)')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('vol 90% CL [Mpc^3]')
plt.ylabel('m1+m2 [M0]')
plt.title('M_total vs vol (90% CL)')

patches = [ matpatch.Patch(color=colorMap[n], label=f'{n} Detectors') for n in colorMap ]
plt.legend(handles=patches, loc='lower right')

plt.show()