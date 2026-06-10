#Scatterplot illustrating the correlation between total mass and luminosity distance for BBH events.
#The color gradient indicates the number of AGNs associated with each event.

#Uses locally accessed gw simulation data: combined_bbh_O5.csv, 

############################################################################### Imports, constants

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math

events = pd.read_csv(r"/path/to/data.csv", delim_whitespace=True)

float rho_agn # 1/MpC^3 # volume density of AGN
float rho_event # proportion of AGN flares associated with BBH mergers

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
plt.xlabel('x axis')
plt.ylabel('y axis')
plt.title('Title')
plt.show()