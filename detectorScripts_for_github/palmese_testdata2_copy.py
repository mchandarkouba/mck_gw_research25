#Scatterplot illustrating the correlation between total mass and luminosity distance for BBH events.
#The color gradient indicates the number of AGNs associated with each event.

#Uses locally accessed gw simulation data: combined_bbh_O5.csv, 

############################################################################### Imports, constants

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

#events = np.loadtxt(r"/Users/mck/Desktop/palmese_research/combined_bbh_O5.csv", skiprows=1)
events = pd.read_csv(r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/combined_bbh_O5.csv", delim_whitespace=True)
print(events)

rho_agn = 10**(-4.75) # 1/MpC^3 # volume density of AGN
rho_event = 1e-4 #proportion of AGN flares associated with BBH mergers

############################################################################### Variables and Plotting

D = events['distance']
M = events['mass1'] + events['mass2']
V_90 = events['vol90']
log_AGN = np.log(V_90 * rho_agn * rho_event)

plt.figure(figsize=(5,5))
plt.scatter(D, M, c=log_AGN)

plt.colorbar(label='log10(N_AGN)')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('dL [MpC]')
plt.ylabel('m1+m2 [M0]')
plt.title('M_total vs dL')
plt.show()
plt.close()

plt.hist(log_AGN, bins=np.linspace(-17,10,30))
plt.show()