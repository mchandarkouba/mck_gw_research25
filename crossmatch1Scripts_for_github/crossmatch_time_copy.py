# Crossmatches Einstein Probe events (EPs) with GraceDB gravitational wave events (GWs)
# by comparing a locally saved list of GWs from crossmatch_GraceDB (in repo) with
# another locally-saved dataset of EP events of interest,
# returning GW-EP pairs which were detected within 1 day of one another.

# used an environment with python=3.10, astropy, and pandas installed

import pandas as pd
import datetime
import numpy as np
import os
from astropy.io import fits

###############################################################################

def get_skymapDates(): # returns all observation dates for the skymaps in a locally saved folder of FITS GW skymaps
    
    folder_path =  r"/Users/mck/Desktop/palmese_research/crossmatch_skymaps/"
    dates = []
    
    for filename in os.listdir(folder_path):
        if '.fits' in filename:
            
            file_path = folder_path + filename
            dates.append((filename, retrieveDate(file_path)))
        
    return dates
        
def retrieveDate(path): # reads FITS file for the observation date/time
    with fits.open(path, ignore_missing_simple=True) as hduList:
        metadata = hduList[1].header
        obs_time = metadata["DATE-OBS"]
        return obs_time
    
def crossmatch_time(gws, eps): # crossmatches and returns GW-EP pairs using numpy broadcasting and index preservation
    
    dt = np.timedelta64(1, 'D')
    gwT = gws['trigger time'].values.reshape(-1,1)
    epT = eps['trigger time'].values.reshape(1,-1)        
        
    pairs = pd.DataFrame( (gwT-dt <= epT) & (epT <= gwT+dt),
                           index=list(gws['file name']),
                           columns=list(eps['source name']) )
            
    return pairs

###############################################################################

def getPairs(): # formats datasets and calls crossmatch function

    csv_path = r"/Users/mck/Desktop/palmese_research/crossmatch_analysis_scripts/EP_data1.csv"
    ep_data = pd.read_csv(csv_path, comment='#').dropna(how='any', axis='index') #comment='#' removes the units line of the csv from docs
    
    ep_dates = ep_data[['source name','trigger time']].dropna()
    ep_dates['trigger time'] = pd.to_datetime( ep_dates['trigger time'], format='mixed' )
    
    gw_dates = pd.DataFrame(get_skymapDates(), columns=['file name', 'trigger time'])
    gw_dates['trigger time'] = pd.to_datetime( gw_dates['trigger time'], format='mixed' )
    
    return crossmatch_time(gw_dates, ep_dates)
