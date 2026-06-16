#Cumulative Distribution Function plot, shows the distribution of simulated AGN flares
#associated with black hole mergers, for various black hole binary mass bins. 

#Uses locally accessed gw simulation data in csv

############################################################################### Imports & physical parameters

import numpy as np
import matplotlib.pyplot as plt

events = np.loadtxt(r"/path/to/file.csv", skiprows=1)

float rho_agn # 1/MpC^3 density of AGN per unit volume
float rho_event # proportion of AGN flares expected to be associated with BBH mergers?

############################################################################### Mass bins & Plotting
float cutoff1
float cutoff2

group1 = events[events[:,5]<cutoff1, :]
group2 = events[(cutoff1<events[:,5]) & (events[:,5]<cutoff2), :]
group3 = events[cutoff2<events[:,5], :]

g1_AGN90 = group1[:,18] * rho_agn * rho_event #group 1, #AGN based off 90% confidence V
g2_AGN90 = group2[:,18] * rho_agn * rho_event
g3_AGN90 = group3[:,18] * rho_agn * rho_event



N = int(1e7) #plot bin size

hist1, bins1 = np.histogram(g1_AGN90, bins=N, range=(10**-2.5, 10**2.5), density=True)
hist2, bins2 = np.histogram(g2_AGN90, bins=N, range=(10**-2.5, 10**2.5), density=True)
hist3, bins3 = np.histogram(g3_AGN90, bins=N, range=(10**-2.5, 10**2.5), density=True)

cdf1 = np.cumsum(hist1) * (bins1[1] - bins1[0])
cdf2 = np.cumsum(hist2) * (bins2[1] - bins2[0])
cdf3 = np.cumsum(hist3) * (bins3[1] - bins3[0])

plt.figure(figsize=(5,5))
plt.plot(bins1[1:], cdf1, c='red', label=f'm1 < {cutoff1} M$_{\odot}$')
plt.plot(bins2[1:], cdf2, c='blue', label=f'{cutoff1} M$_{\odot}$ < m1 < {cutoff2} M$_{\odot}$')
plt.plot(bins3[1:], cdf3, c='orange', label=f'{cutoff2} M$_{\odot}$ < m1')

plt.legend()
plt.xscale('log')
plt.xlabel('x axis label')
plt.ylabel('y axis label')
plt.title('Title')
plt.show()