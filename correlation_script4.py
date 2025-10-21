import pandas as pd
import numpy as np
import json

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




TOP_DIR = f'{path.dirname(__file__)}/' #all corr script files
MOD_DIR = TOP_DIR + 'modified_skymaps/' #for flattened and confidence region skymaps
OUTPUT_DIR = TOP_DIR + 'corr_output/' #unused at the moment
LIGO_DIR = TOP_DIR + 'LIGO_runs/' #where to access GW files


###############################################################################

class GW:
    """
    Stores gravitational wave event files, directories easily
    _flat denotes a flattened healpix skymap, mainly for matplotlib
    
    """
    
    def __init__(self, name:str, subpath:str):
        path_dir = LIGO_DIR
        
        region = self.ConfidenceRegion

        self.name = region.name = name
        self.path = region.path = path_dir + subpath
        
        self.path_flat = MOD_DIR + self.name + "_flattened.fits"
        
        if path.exists(self.path_flat):
            print(f"\t\t{self.path_flat} already found, using this.")
        
        else:
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
            
            self.name = f"{np.round(100*c)}percent_{self.name}" # <-- technically not yet flattened
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
            skymap = skymap[:i]
            skymap.sort('UNIQ')
            #skymap = skymap['UNIQ',]
            
            #Making FITS skymap, flattening and saving data to the GW object
            
            self.path = MOD_DIR + self.name
            self.path_flat = MOD_DIR + self.name + "_flattened.fits"
            
            skymap.write(self.path, overwrite=True, format='fits')
            
            if path.exists(self.path_flat):
                print(f'\t\t{self.path_flat} already found, using this.')
            
            else:
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

def findCandidates(gw_name, gw_path, catalogs, conf):
    gw = GW(gw_name, gw_path) #'S240807h' for testing
    region = gw.ConfidenceRegion(conf)
    
    colnames = ['ID', 'COORD', 'CATALOG']
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
        CATALOG = Table.Column([name]*len(ID))

        objects = Table([ID, COORD, CATALOG], names=colnames)
        bools = inRegion(gw, region, objects)
        candidates = vstack([candidates, objects[bools]]) if np.size(candidates)!=0 else objects[bools]
    
    return candidates

def gw_dataFrame(updateDF=False):
    
    gw_arr_name = 'gw_dataFrame.csv'
    gw_arr_path = TOP_DIR + gw_arr_name
    
    if path.exists(gw_arr_path) and not updateDF:
        print('\tGW File already exists!')
        gw_arr = pd.read_csv(gw_arr_path)
    
    else:
        
        #for subdir in LIGO_DIR
        
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

def get_gw_name(s):
    for i in range(len(s)-1):
        c0 = s[i]
        c1 = s[i+1]
        if c0+c1 == 'GW':
            try:
                date = str(int(s[i+2:i+8])) #checking for int following "GW"
                superevent = c0+c1+date
                
                try:
                    time = str(int(s[i+9: i+15])) #checking for int following "GWxxxxxx_"
                    time = '0'*(6-len(time)) + time
                    superevent = superevent + '_' + time
                    return superevent
                    
                except Exception:
                    return superevent
            
            except Exception as e:
                continue
    
    return None       
          
def pick_fits():
        
    ligo_runs = [ s+'/' for s in listdir(LIGO_DIR) ]
    ligo_runs.remove(".DS_Store/")
    eventNames = set()
    
    for subdir in ligo_runs:
        files = listdir(LIGO_DIR+subdir)
        
        for file in files:
            eventNames |= set([get_gw_name(file)])
    
    eventNames.remove(None)
        
    
    
    eventMaster = dict()
    for subdir in ligo_runs:
        sub_filenames = listdir(LIGO_DIR+subdir)
        
        for name in eventNames:
            for filename in sub_filenames:
                
                if name in filename: eventMaster[name] = eventMaster.get(name, set()) | set((subdir+filename,))
    
    
    
    eventPreferred = dict()
    for event in eventMaster:
        for filename in eventMaster[event]:
            
            hierarchy = ['NRSur7dq4', 'Mixed', 'sco']
            fitsFound = False
            for fitsType in hierarchy:
                if fitsType in filename and not fitsFound:
                    eventPreferred[event] = filename
                    fitsFound = True
            
            if not fitsFound:
                eventPreferred[event] = tuple(eventMaster[event])[0]
    
    return eventPreferred

###############################################################################

def main():
    catalogs = {'milliquas', 'quaia'}    
    far = 3.175e-8 #1/yr in Hz
    conf = 0.90
    
    gw_skymaps = pick_fits()
    
    candFilename = f"candidates_{'&'.join(catalogs)}_{np.round(100*conf)}percent.fits"
    eventFilename = f"event_cache_{'&'.join(catalogs)}_{np.round(100*conf)}percent.json"
    
    candFilepath = OUTPUT_DIR + candFilename
    eventFilepath = OUTPUT_DIR + eventFilename
    
    allSavedCands = None
    allUsedEvents = {}
    
    if path.exists(candFilepath):
        print(f"Found candidate table from catalogs {catalogs} at confidence {conf}, loading...")
        with fits_astropy.open(candFilepath, cache=False) as hdul:
            allSavedCands = Table( hdul[1].data )
            negativeDists = np.where(allSavedCands["COORD.distance"] < 0)
            RA, DEC, dL = allSavedCands["COORD.ra"], allSavedCands["COORD.dec"], allSavedCands["COORD.distance"]
            
            COORD = SkyCoord(
                             RA*units.deg,
                             DEC*units.deg,
                             dL*units.Mpc,
                             
                             frame="icrs", 
                             )
            
            allSavedCands.remove_columns(["COORD.ra", "COORD.dec", "COORD.distance"])
            allSavedCands.add_column(COORD, index=1, name="COORD")
            
        print("\tLoaded and formatted!\n")
    
    else:
        print(f"Creating table of candidates for catalogs {catalogs} at confidence {conf} .\n")

            
    if path.exists(eventFilepath):
        print(f"Found JSON of saved GW-AGN pairs crossmatched from catalogs {catalogs} at confidence {conf}, loading...\n")
        with open(eventFilepath, 'r') as file: allUsedEvents = json.load(file)
        
    else:
        print(f"Creating JSON for GW-AGN pairs using catalogs {catalogs} at confidence {conf} .\n")

            
    for tally, event in enumerate(gw_skymaps):
        
        if event in allUsedEvents:
            print(f"\tEvent {event} has already been crossmatched according to JSON {eventFilename}, skipping... ({tally}/{len(gw_skymaps)}, {100*np.round(tally/len(gw_skymaps), 3)}%)\n")
            continue
        
        eventCands = findCandidates(event, gw_skymaps[event], catalogs, conf)

        if allSavedCands is not None:
            newCands = np.setdiff1d( eventCands["ID"], allSavedCands["ID"], assume_unique=True)
            newCands = eventCands[ [x in newCands for x in eventCands["ID"]] ]
            allSavedCands = vstack([allSavedCands, newCands])
            
        else:
            allSavedCands = eventCands
        
        
        allSavedCands.write(candFilepath, format="fits", overwrite=True)
        
        allUsedEvents[event] = str([ str(x) for x in eventCands["ID"] ])
        with open(eventFilepath, 'w') as file: json.dump(allUsedEvents, file)
        print(f"\t\tCrossmatched and cached to JSON! ({tally}/{len(gw_skymaps)}, {100*np.round(tally/len(gw_skymaps), 3)}%)\n")

        
###############################################################################

main()


