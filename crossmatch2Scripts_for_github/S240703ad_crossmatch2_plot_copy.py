# Incomplete skymap plot to serve as a visual aid for S240703ad_crossmatch2.py (in repo)
# was meant to plot the S240703ad skymap with 50% and 90% confidence regions, along with
# reticles for EP240703c and the two nearby AGN objects 
# skymap and EP transient + AGN object positions are saved locally

# at the moment, the coordinate systems for EP240703c and the AGN objects are diffferent
# additionally, the confidence regions are probability level curves instead of cumulative probability regions

# an environment with python=3.10, astropy, ligo.skymap, and pandas installed was used

############################################################################### IMPORTS

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

import ligo.skymap.plot as ligoPlot
from ligo.skymap.io.fits import read_sky_map
from ligo.skymap.plot.allsky import AutoScaledWCSAxes as WCSAxes
from ligo.skymap.plot.marker import reticle

from healpy import read_map

import pandas as pd
from matplotlib import pyplot as plt

############################################################################### EVENTS AND OBJECTS

filename = 'S240703ad_Bilby.fits.gz,0'
filepath = r"/Users/mck/Desktop/palmese_research/crossmatch2_files/" + filename 

gw_flat = {'filename':filename,
           'filepath':filepath,
           'fits':read_sky_map(filepath),
           'healpix':read_map(filepath)}

filename = 'S240703ad_4_Bilby.multiorder.fits,0'
filepath = r"/Users/mck/Desktop/palmese_research/crossmatch_skymaps/" + filename

gw_multi = {'filename':filename,
           'filepath':filepath,
           'fits':read_sky_map(filepath),
           'healpix':None}

csv_path = r"/Users/mck/Desktop/palmese_research/crossmatch2_files/crossmatch2_objects.csv"
obj_data = pd.read_csv(csv_path).dropna(how='any', axis='index') #comment='#' removes the units line of the csv from docs

csv_path = r"/Users/mck/Desktop/palmese_research/crossmatch_analysis_scripts/EP_data1.csv"
ep_data = pd.read_csv(csv_path, comment='#').dropna(how='any', axis='index')

ep_event = ep_data[ep_data['source name']=='EP240703c']

ep_format = ep_event[['source name', 'RA', 'DEC']]
ep_format.columns = ['ID', 'RA', 'DEC']

obj_format = obj_data[['ID 1', 'RA', 'DEC']]
obj_format.columns = ['ID', 'RA', 'DEC']

points = pd.concat([ep_format, obj_format])
points.index = [0,1,2,3]
points['COL'] = pd.Series(pd.Series(['r','b','b','b']))


############################################################################### 

fig = plt.figure(figsize=(6,6), dpi=75)
ax = plt.axes([0.05, 0.05, 0.9, 0.9], #left, bottom, width, height as fractions of figure
              projection='astro hours mollweide')
ax.grid()

hdul = fits.open(gw_flat['filepath'], ignore_missing_simple=True)
hdr = hdul[1].header
    
ax.contour_hpx(data=gw_flat['healpix'], )

ax.imshow_hpx(gw_flat['filepath'], cmap='cylon')

for marker in points.iterrows():
    point = marker[1]
    ax.plot(point['RA'], point['DEC'], transform=ax.get_transform('world'), markersize=10, markeredgewidth=2, marker=True)
