#Scatterplot illustrating the correlation between luminosity distance 
#and 70% confidence level area of the sky-location for BBH events.
#The color gradient indicates the number of AGNs associated with each event.

import numpy as np
import matplotlib.pyplot as plt
import math

events = np.loadtxt(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/combined_bbh_O5.csv", skiprows=1)

rho_agn = 10**(-4.75) # 1/MpC^3 #volume density of AGN in space
rho_event = 1e-4 #expected proportion of AGN flares associated with BBH mergers

############################################################################### Scatter plot variables

D = events[:,1]
A_70 = events[:,11]
M = events[:,5] + events[:,6]
V_90 = events[:,18]
log_AGN = np.log(V_90 * rho_agn * rho_event)

############################################################################### Plotting

plt.figure(figsize=(10,10))
plt.scatter(D, A_70, c=log_AGN)

plt.colorbar(label='log10(N_AGN)')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('dL [Mpc]')
plt.ylabel('area (70% CL) [deg^2]')
plt.title('area (70% CL) vs dL')

plt.show()