###############################################################################
### IMPORTS
###############################################################################

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.cosmology.realizations import Planck15
from astropy.io import fits as fits_astropy
from astropy.table import Table, vstack, unique

from gw_class import GW, recursive_all_skymaps, retrieve_skymaps, get_gw_name, default_hierarchy

import json

import numpy as np

from os import path

import pandas as pd

###############################################################################
### DIRECTORIES
###############################################################################

TOP_DIR = f"{path.dirname(__file__)}/" #all corr script files
MOD_DIR = TOP_DIR + "modified_skymaps/" #for flattened and confidence region skymaps
OUTPUT_DIR = TOP_DIR + "corr_output/" #for final agn and gw dataframes
LIGO_DIR = TOP_DIR + "GWTC_5.0/" #where to access GW files

###############################################################################
### HELPER FUNCTIONS
###############################################################################

def messenger(*keys, fvars:list=[]):
    """
    Parameters
    ----------
    *keys : int > 0
        keys for nested message dictionaries.
    fvars : list, optional
        Required variables for some messages which use f-strings. 
        Args in fvars are assigned in order to inserted variables in messages.
        The default is [].

    Effects
    ---------
    Messenger is toggleable to print stored messages in the interpreter using global verbose.
    All messages are stored into a .txt log file to track progression.
    
    """
    
    logpath = TOP_DIR + "correlation_script6_log.txt"
    hbar = '\n' + '-'*100 + '\n'
    
    master = {1 : {0 : hbar+"Assigning .fits skymaps to GWEvents in GWTC 5.0"+hbar,
                   1 : {0 : lambda L: f"No .fits files found matching {L[0]}; skipping...",
                        },
                   2 : {0 : lambda L: f"Number of GWEvents skipped while loading .fits files: {L[0]}",
                        },
                   },
              
              2 : {0 : hbar+"Creating or loading GW CSV"+hbar,
                   1 : {0 : "GW DataFrame already exists! Using this.",
                        },
                   2 : {0 : "Creating GW DataFrame from GWTC 5.0 and previous step...",
                        },
                   3 : {0 : lambda L: f"Progress: {100*(L[0]+1)/L[1]//1}%",
                        },
                   },
              
              3 : {0 : lambda L: hbar+f"Crossmatching GWTC 5.0 events with {L[0]} to create AGN FITS and JSON"+hbar,
                   1 : {0 : lambda L: f"Found candidate table from catalogs {L[0]} at confidence {L[1]}, loading...",
                        },
                   2 : {0 : "Loaded and formatted!",
                        },
                   3 : {0 : lambda L: f"Creating table of candidates for catalogs {L[0]} at confidence {L[1]} .",
                        },
                   4 : {0 : "Found matching JSON of saved GW-AGN pairs, loading...",
                        },
                   5 : {0 : "Creating JSON for GW-AGN pairs for given catalogs and confidence level",
                        },
                   6 : {0 : lambda L: f"\tEvent {L[0]} has already been crossmatched according to JSON {L[1]}, skipping... {L[2]}",
                        },
                   7 : {0 : lambda L: f"\tBeginning crossmatch: event {L[0]} and input catalogs...",
                        },
                   8 : {0 : lambda L: f"Crossmatched and cached to JSON! {L[0]}",
                        },
                   },
              
              4 : {0 : hbar+"Loading AGN FITS and converting to AGN csv"+hbar,
                    1 : {0 : "Formatted and saved!",
                         },
                    },
              }
    
    d = master
    for k in keys: d = d[k]
    
    message = d[0]
    
    if isinstance(message, type(lambda x: None)): 
        try: 
            message = message(fvars) 
        except IndexError:
            print("ERROR compiling message, check if len(fvars) matches message inputs")
            message = message(fvars+['',]*30)
    
    message = '\t'*(len(keys)-1) + message
    
    if path.exists(logpath): 
        with open(logpath, 'r') as file: txt = file.read()
    else: txt = ''
    
    with open(logpath, 'w') as file: file.write(txt+'\n'+message)
    
    if verbose: print(message)
    return

###############################################################################

def load_catalog(name:str) -> (Table,list):
    """
    Parameters
    ----------
    name : str
        Currently accepts keys "milliquas" or "quaia". 
        Catalogs not on file will be downloaded from stored links.

    Returns
    -------
    (Table,list)
        Returns astropy Table of the catalog and the list keys 
        corresponding to colnames for object id, ra, dec, and redshift.
    """
    
    key = {
            'milliquas' : {'url': 'https://quasars.org/milliquas.fits.zip',
                          'keys': ['NAME', 'RA', 'DEC', 'Z']
                          },
           
           'quaia' : {'url': 'https://zenodo.org/records/8060755/files/quaia_G20.0.fits?download=1',
                      'keys': ['source_id', 'ra', 'dec', 'redshift_quaia']
                     },
           
           'desi_agngal_dr1' : {'url': '',
                                'keys' : ["TARGETID", "TARGET_RA", "TARGET_DEC", "Z"],
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

###############################################################################

def format_catalog(c:np.ndarray, keys:dict) -> np.ndarray:
    ID, RA, DEC, Z = keys
    c = c[ np.where(c[Z]>=0) ]
    if not isinstance(c[ID].dtype, str): c[ID] = c[ID].astype(str)
    
    return c

###############################################################################

def find_cands(gwPath:str, catalogs:set, conf:float,) -> Table:
    """
    Parameters
    ----------
    gwPath : str
        Path to GWEvent skymap.
    catalogs : set
        Object catalogs.
    conf : float
        Crossmatch confidence level.

    Returns
    -------
    Table
        Subset of catalogs crossmatched and found to be within the conf% confidence volume
        of the GWEvent associated with gwPath.

    """
    colnames = ['ID', 'COORD', 'Z', 'CATALOG']
    colnames += ['t_peak', 'f_peak', 't_rise', 't_decay', 'f_base']
    candidates = Table([ np.array([]) ]*len(colnames),
                       names=colnames
                       )
    gw = GW(path=gwPath)
    
    
    for name in catalogs:
        catalog, keys = load_catalog(name)
        ID, RA, DEC, Z = keys
        catalog = format_catalog(catalog, keys)
        
        ID, RA, DEC, Z = catalog[ID], catalog[RA]*units.deg, catalog[DEC]*units.deg, catalog[Z]
        dL = Planck15.luminosity_distance(Z)
        COORD = SkyCoord(RA, DEC, distance=dL, frame='icrs')
        CATALOG = Table.Column([name]*len(ID))
        NA = Table.Column([np.nan]*len(ID))

        objects = Table([ID, COORD, Z, CATALOG]+[NA]*5, names=colnames)
        bools = gw.crossmatch_space(COORD) <= conf
        candidates = vstack([candidates, objects[bools]]) if np.size(candidates)!=0 else objects[bools]
    
    return candidates   

###############################################################################
### MAIN FUNCTIONS
###############################################################################

def pick_fits():
    """
    Returns
    -------
    gwSkymaps : dict
        { gw name : skymap filepath } for all events in GWTC 5.0 csv file.

    """
    gwSkymaps = dict()
    allSkymaps = recursive_all_skymaps(LIGO_DIR)

    filename = 'event-versions_gwtc5.csv'
    filepath = TOP_DIR + filename
    event_data = pd.read_csv(filepath)
    names = np.unique([ get_gw_name(url) for url in event_data["detail_url"] ])
    vera_hierarchy = ["/hildafs/projects/phy220048p/chandark/corr_files/GWTC_5.0/GWTC-Preferred-Skymaps"] + default_hierarchy
    
    skipped = []
    for name in names: 
        try:
            gwSkymaps[name] = retrieve_skymaps(name,
                                               hierarchy=vera_hierarchy,
                                               skymaps=allSkymaps,
                                               )[0]
        except IndexError:
            messenger(1,1, fvars=[name,])
            gwSkymaps[name] = None
            skipped.append(name.item())
        
    messenger(1,2, fvars=[len(skipped),])
    print(sorted(skipped))
    
    return gwSkymaps
    
###############################################################################

def gw_DataFrame(conf, gwSkymaps:dict, updateDF=True):
    """
    Parameters
    ----------
    conf : float
        Confidence level for crossmatching step.
    gwSkymaps : dict
        see pick_fits().
    updateDF : bool, optional
        When False, any existing DataFrame on file is used. 
        The default is True.

    Effects
    ----------
    Reformats GWTC 5.0 event-versions csv for lambda posterior code
    into an original and clean (nan rows dropped) csv files.
    
    """
    
    dfName_origin = 'gw_dataFrame_original.csv'
    dfName_clean = 'gw_dataFrame_clean.csv'
    
    gwDataFrames = () #final output

    
    if all(path.exists(OUTPUT_DIR+dfName) for dfName in [dfName_origin, dfName_clean]) and not updateDF:
        messenger(2,1)
        
        for dfName in [dfName_origin, dfName_clean]:
            dfPath = OUTPUT_DIR + dfName
            gwDataFrames += ( pd.read_csv(dfPath), )
    
    else:
        messenger(2,2)        
        filename = 'event-versions_gwtc5.csv'
        filepath = TOP_DIR + filename
        event_data = pd.read_csv(filepath)
        
        massCols = ['mass_1_source', 'mass_2_source','final_mass_source_upper', ]#list(set([(axis_name if 'mass' in axis_name else None) for axis_name in event_data.axes]).remove(None))
        MASSES = event_data[massCols]
        
        dfPath_origin = OUTPUT_DIR + dfName_origin
                
        FREQ = pd.Series([conf]*len(MASSES))
        ID = []
        PATH = []
        for i,url in enumerate(event_data['detail_url']):
            graceID = get_gw_name(url)
            ID += [graceID]
            
            skymapPath = gwSkymaps.get(graceID, '')
            if skymapPath != None:
                PATH += [skymapPath]
            else:
                PATH += [np.nan]
     
            if i%10==9: messenger(2,3, fvars=[i,len(event_data["detail_url"]),])
        
        print()
        ID = pd.Series(ID, name='grace_id')
        PATH = pd.Series(PATH, name='skymap_path')
        gwDataFrame_origin = pd.concat([FREQ, MASSES, ID, PATH], axis=1)
        gwDataFrame_origin.to_csv(dfPath_origin)
        
        noNulls = np.all(np.array( [ ~gwDataFrame_origin[col].isnull() for col in gwDataFrame_origin ] ),
                         axis=0,
                         )
        gwDataFrame_clean = gwDataFrame_origin[noNulls]
        dfPath_clean = OUTPUT_DIR + dfName_clean
        gwDataFrame_clean.to_csv(dfPath_clean)
        
        #gwDataFrames = (gwDataFrame_origin, gwDataFrame_clean)

    return

###############################################################################

def agn_Table(catalogs:set, conf:float, gwSkymaps:dict) -> Table:
    """
    Parameters
    ----------
    catalogs : set
        Object catalogs to be selected from.
    conf : float
        Crossmatch step confidence level.
    gwSkymaps : dict
        see pick_fits().

    Returns
    -------
    Table
        Subset of combined catalogs, of all events included within the conf% confidence volume 
        of any GWEvent in gwSkymaps.

    """
    
    candFilename = f"candidates_{'&'.join(catalogs)}_{np.round(100*conf)}percent.fits"
    eventFilename = f"event_cache_{'&'.join(catalogs)}_{np.round(100*conf)}percent.json"
    
    candFilepath = OUTPUT_DIR + candFilename
    eventFilepath = OUTPUT_DIR + eventFilename
    
    allSavedCands = None
    allUsedEvents = dict()
    
    ###########################################################################
    
    if path.exists(candFilepath):
        messenger(3,1, fvars=[catalogs,conf])
        with fits_astropy.open(candFilepath, cache=False) as hdul:
            allSavedCands = Table( hdul[1].data )
            RA, DEC, dL = allSavedCands["COORD.ra"], allSavedCands["COORD.dec"], allSavedCands["COORD.distance"]
            
            COORD = SkyCoord(
                             RA*units.deg,
                             DEC*units.deg,
                             dL*units.Mpc,
                             
                             frame="icrs", 
                             )
            
            allSavedCands.remove_columns(["COORD.ra", "COORD.dec", "COORD.distance"])
            allSavedCands.add_column(COORD, index=1, name="COORD")
            
        messenger(3,2)
    
    else:
        messenger(3,3, fvars=[catalogs,conf])
            
    if path.exists(eventFilepath):
        messenger(3,4)
        with open(eventFilepath, 'r') as file: allUsedEvents = json.load(file)
        
    else:
        messenger(3,5)

    ###########################################################################      
      
    for tally, name in enumerate(gwSkymaps):
        progress = f"({tally}/{len(gwSkymaps)}, {100*np.round(tally/len(gwSkymaps), 3)}%)"
        
        if name in allUsedEvents:
            messenger(3,6, fvars=[name, eventFilename, progress])
            continue
        
        messenger(3,7, fvars=[name])
        if gwSkymaps[name]==None: continue
        eventCands = find_cands(gwSkymaps[name], catalogs, conf)
        
        if allSavedCands!=None :
            newCands = np.setdiff1d(eventCands["ID"], allSavedCands["ID"])
            newCands = eventCands[ np.where(np.isin(eventCands["ID"], newCands))[0] ]
            if np.size(newCands["ID"])>0: allSavedCands = vstack([allSavedCands, newCands])
            
        else:
            
            allSavedCands = eventCands  

        allSavedCands.write(candFilepath, format="fits", overwrite=True)
        allUsedEvents[name] = str([ str(x) for x in eventCands["ID"] ])
        with open(eventFilepath, 'w') as file: json.dump(allUsedEvents, file)
        messenger(3,8, fvars=[progress])

    return

###############################################################################

def agnTable_to_csv(catalogs, conf):
    """
    Parameters
    ----------
    catalogs : set
    conf : float

    Returns
    -------
    .csv formatted AGN Table, final step after all crossmatching.

    """
    candFilename = f"candidates_{'&'.join(catalogs)}_{np.round(100*conf)}percent.fits"
    candFilepath = OUTPUT_DIR + candFilename
    
    csvFilename = candFilename[:-4]+"csv"
    csvFilepath = OUTPUT_DIR + csvFilename

    
    with fits_astropy.open(candFilepath, cache=False) as hdul:
        data = Table( hdul[1].data )
        
        data["COORD.ra"].name = 'ra'
        data["COORD.dec"].name = 'dec'
        
        data.write(csvFilepath, format="csv", overwrite=True,)
        
    messenger(4,1)

###############################################################################
### MAIN
###############################################################################

def run_all():
    catalogs = {'desi_agngal_dr1',}    
    conf = 0.95

    messenger(1)
    gwSkymaps = pick_fits()
    
    messenger(2)
    gw_DataFrame(conf, gwSkymaps, updateDF=True)
    
    messenger(3, fvars=[catalogs,])
    agn_Table(catalogs, conf, gwSkymaps)
    
    messenger(4)
    agnTable_to_csv(catalogs, conf)
    
###############################################################################

verbose = True
run_all()
