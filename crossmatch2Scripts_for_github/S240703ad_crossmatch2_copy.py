# A secondary crossmatch for GW event S240703ad and AGN objects from the QUAIA and MILLIQUAS catalogs
# AGN objects within 3arcmin, the localization, of S240703ad from both catalogs were saved to a local csv file
# A similar crossmatch procedure for crossmatch_space (in repo) is used,
# except also uses 3D crossmatching with the AGN redshift being used as the distance metric.

# environment with ligo.skymap, astropy, and pandas (python=3.10) was used

import pandas as pd
import numpy as np

from ligo.skymap.postprocess.crossmatch import crossmatch
from ligo.skymap.io.fits import read_sky_map
from astropy.coordinates import SkyCoord
import astropy.units as units
from astropy.cosmology.realizations import Planck15

###############################################################################

gw_filename = 'S240703ad_4_Bilby.multiorder.fits,0'
csv_path = r"/Users/mck/Desktop/palmese_research/crossmatch2_files/crossmatch2_objects.csv"
obj_data = pd.read_csv(csv_path).dropna(how='any', axis='index') #comment='#' removes the units line of the csv from docs

def crossmatch2(gw, obj_data):
    gw_path = r"/Users/mck/Desktop/palmese_research/crossmatch_skymaps/" + gw
    
    ra,dec,dLum = np.array(obj_data['RA'])*units.deg, np.array(obj_data['DEC'])*units.deg, Planck15.luminosity_distance( np.array(obj_data['z']) )

    print(np.array(obj_data['z']))
    
    obj_coords = SkyCoord(ra, dec, distance=dLum, frame='icrs')
            
    results = crossmatch(sky_map=read_sky_map(filename=gw_path, moc=True),
                         coordinates=obj_coords,
                         cosmology=True)
            
    probs = results.searched_prob_vol

    return probs
    
probs = crossmatch2(gw_filename, obj_data)
print(probs)