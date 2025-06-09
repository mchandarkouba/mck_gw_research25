#Plots the relationship between total mass and volume for BBH events, color coded by the number of AGNs expected in each event.

#Uses locally accessed gw simulation data: combined_bbh_O5.csv, 

############################################################################### Imports and Constants

import numpy as np
import matplotlib.pyplot as plt
import math

events = np.loadtxt(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/combined_bbh_O5.csv", skiprows=1)

rho_agn = 10**(-4.75) # 1/MpC^3 # volume density of AGN in space
rho_event = 1e-4 # proportion of AGN flares associated with BBH mergers

############################################################################### Variables and Plotting

V_90 = events[:,18]
M = events[:,5] + events[:,6]
log_AGN = np.log(V_90 * rho_agn * rho_event)


plt.figure(figsize=(5,5))
plt.scatter(V_90, M, c=log_AGN)

plt.colorbar(label='log10(N_AGN)')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('vol 90% CL [Mpc^3]')
plt.ylabel('m1+m2 [M0]')
plt.title('M_total vs vol (90% CL)')
plt.show()