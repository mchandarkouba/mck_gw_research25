import pandas as pd
import numpy as np

from os import path, listdir
import subprocess
from skymap2_corr import get_skymap

from astropy import units
from astropy.io import fits as fits_astropy
from astropy.table import Table, QTable, join, vstack
from astropy.coordinates import SkyCoord
import astropy_healpix as ah
from astropy.cosmology.realizations import Planck15

from ligo.skymap.io import fits as fits_ligo
from ligo.skymap.postprocess.crossmatch import crossmatch
from ligo.gracedb.rest import GraceDb
import requests




TOP_DIR = f'{path.dirname(__file__)}/'
SKYMAP_DIR = TOP_DIR + 'corr_skymaps/'
OUTPUT_DIR = TOP_DIR + 'corr_output/'

###############################################################################

class GW:
    """
    Stores gravitational wave event files, directories easily
    _flat denotes a flattened healpix skymap, mainly for matplotlib
    
    """
    
    def __init__(self, name:str):
        path_dir = SKYMAP_DIR
        
        region = self.ConfidenceRegion

        self.name = region.name = name
        self.path = region.path = path_dir + self.name
        
        self.name_flat = self.name.replace('.multiorder', '.flattened')
        self.path_flat = self.path.replace('.multiorder', '.flattened')
        subprocess.run(['ligo-skymap-flatten', self.path, self.path_flat]) #uses cmd to flatten the skymap from a filepath
        
        self.fits = fits_ligo.read_sky_map(self.path, moc=True)
        #self.fits_flat = fits_ligo.read_sky_map(self.path_flat)
        
        self.prob = None
        self.img = None

    def __repr__(self):
        return self.name
            
    class ConfidenceRegion():
        """
        Confidence region and area retrieval
        Documentation here: https://emfollow.docs.ligo.org/userguide/tutorial/multiorder_skymaps.html
        
        """
            
        def __init__(self, c):
            self.level = c
                        
            skymap = QTable.read(self.path, format='fits')
            skymap.sort('PROBDENSITY', reverse=True)
            
            level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
            
            pix_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
            
            prob = pix_area * skymap['PROBDENSITY']
            cumprob = np.cumsum(prob)
            
            # find the ith pixel of the sorted pixels where 0th through ith sum cumulatively to prob c
            i = cumprob.searchsorted(c)
            
            #adding 0th through ith pixels' areas to find area of c% confidence region
            area_c = pix_area[:i].sum()
            self.area = area_c.to(units.deg**2)
                
            #Creating, flattening, then saving a skymap
            self.skymap_mod = skymap
            skymap = skymap[:i]
            skymap.sort('UNIQ')
            #skymap = skymap['UNIQ',]
            
            #Making FITS skymap, flattening and saving data to the GW object
            self.name = f"{np.round(100*c)}percent_{self.name}" # <-- technically not yet flattened
            self.name_flat = self.name.replace('.multiorder', '.flattened')
            
            self.path = SKYMAP_DIR + self.name
            self.path_flat = self.path.replace('.multiorder', '.flattened')
            
            skymap.write(self.path, overwrite=True, format='fits')
            
            subprocess.run(['ligo-skymap-flatten', self.path, self.path_flat]) #uses cmd to flatten the skymap from a filepath
            
            #self.fits_flat = fits_ligo.read_sky_map(self.path_flat)
            
        def roundArea(self):
            n = self.area
            g = lambda n: int(-np.round( np.log10(n.value)) +2)
            return np.round(n, g(n))
 

def inRegion(gw, region, objects):
    results = crossmatch(sky_map = gw.fits,
                         coordinates=objects['COORD'], 
                         cosmology=True
                         )
    
    return results.searched_prob < region.level

###############################################################################

def loadCatalog(name:str):
    key = {
            'milliquas' : {'url': 'https://quasars.org/milliquas.fits.zip',
                          'keys': ['NAME', 'RA', 'DEC', 'Z']
                          },
           
           'quaia' : {'url': 'https://zenodo.org/records/8060755/files/quaia_G20.0.fits?download=1',
                      'keys': ['source_id', 'ra', 'dec', 'redshift_quaia']
                     },
          }
    
    catalog_dir = f'{TOP_DIR}{name}.fits'
    address=catalog_dir if path.exists(catalog_dir) else key[name]['url']
    
    try:
        with fits_astropy.open(address, cache=False) as hdul:
            hdul.writeto(catalog_dir, overwrite=True)
            
            data = hdul[1].data
            data = Table( data )
            return data, key[name]['keys']
    
    except Exception as e:
        print(f'Error downloading or retrieving: {e}')

def formatCatalog(c, keys):
    ID, RA, DEC, Z = keys
    c = c[ np.where(c[Z]>=0) ]
    if not isinstance(c[ID].dtype, str): c[ID] = c[ID].astype(str)
    
    return c

def findCandidates(graceID, catalogs, conf):
    gw = GW(get_skymap(graceID)) #'S240807h' for testing
    region = gw.ConfidenceRegion(conf)
    
    colnames = ['ID', 'COORD', 'ORIGIN']
    candidates = Table([ np.array([]) ]*len(colnames),
                       names=colnames
                       )
    
    for name in catalogs:
        catalog, keys = loadCatalog(name)
        ID, RA, DEC, Z = keys
        catalog = formatCatalog(catalog, keys)
        
        ID, RA, DEC, Z = catalog[ID], catalog[RA]*units.deg, catalog[DEC]*units.deg, catalog[Z]
        dL = Planck15.luminosity_distance(Z)
        COORD = SkyCoord(RA, DEC, distance=dL, frame='icrs')
        ORIGIN = Table.Column([name]*len(ID))

        objects = Table([ID, COORD, ORIGIN], names=colnames)
        bools = inRegion(gw, region, objects)
        candidates = vstack([candidates, objects[bools]]) if np.size(candidates)!=0 else objects[bools]
    
    return candidates

def gw_dataFrame():
    
    gw_arr_name = 'gw_dataFrame.csv'
    gw_arr_path = TOP_DIR + gw_arr_name
    
    if path.exists(gw_arr_path):
        print('\tGW File already exists!')
        gw_arr = pd.read_csv(gw_arr_path)
    
    else:
    
        filename = 'event-versions.csv'
        filepath = TOP_DIR + filename
        event_data = pd.read_csv(filepath)
        axes = ['final_mass_source_upper', 'detail_url']#list(set([(axis_name if 'mass' in axis_name else None) for axis_name in event_data.axes]).remove(None))
        data = event_data[axes]
        f_cover = pd.Series([0.9]*len(data['final_mass_source_upper']))
        
        superevents = []
        for i,url in enumerate(event_data['detail_url']):
            graceID = requests.get(url).json()['grace_id']
            superevents += [graceID]
            if i%10==9: print(f'\t Progress: {100*(i+1)/len(event_data["detail_url"])//1}%')
        
        superevents = pd.Series(superevents, name='superevents')
        gw_arr = pd.concat([f_cover, data, superevents], axis=1)
    
        gw_arr.to_csv(gw_arr_path)
        
    return gw_arr

###############################################################################

def main():
    catalogs = {'milliquas', 'quaia'}    
    far = 3.175e-8 #1/yr in Hz
    conf = 0.90
    
    print('Retrieving GW Event data...')
    gw_data = gw_dataFrame()
    superevents = gw_data['superevents']
    print(superevents)
    
    print(f'Retrieved! {len(superevents)} events found.\n')   
    
    continue_retrieval = True

    if continue_retrieval:
        print('Finding Candidates for GW events...')
        
        for n,graceID in enumerate(superevents): 
            
            print(f'\n\tTrying event {graceID} ({n+1}/{len(superevents)})...')
            filename = f"{graceID}_{'&'.join(catalogs)}_candidates_{int(conf*100)}.fits"
            filepath = OUTPUT_DIR + filename
            
            if path.exists(filepath):
                print(f"Path {filepath} already exists. Using this ({n+1}/{len(superevents)}).")
                     
            else:
            
                try:
                    data = findCandidates(graceID, catalogs, conf)
                    print(f"\t{len(data['ID'])} candidates found for event {graceID}!")
                
                except Exception as e:
                    raise
                    print(graceID)
                    print(f"\tException raised retrieving candidates for {graceID}: {e}")
                    continue
                
                print(f"\tSaving {filename} to {OUTPUT_DIR}\n")
                data.write(filepath, format='fits', overwrite=True)
            
        
    
        
    return gw_data
        


main()


