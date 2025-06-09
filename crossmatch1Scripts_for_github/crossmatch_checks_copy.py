# Performs checks to crossmatching functions for individual GW/EP pairs for hand-selected events

import pandas as pd
import numpy as np
import datetime
import os

#import crossmatch_space
import crossmatch_time

from astropy.io import fits
from ligo.skymap.postprocess.crossmatch import crossmatch
from ligo.skymap.io.fits import read_sky_map
from astropy.coordinates import SkyCoord
import astropy.units as units

###############################################################################

def get_data(retrieve=None):
    retrieve = retrieve or None
    
    csv_path = r"/Users/mck/Desktop/palmese_research/crossmatch_analysis_scripts/EP_data1.csv"
    ep_data = pd.read_csv(csv_path, comment='#').dropna(how='any', axis='index') #comment='#' removes the units line of the csv from docs

    if retrieve=='time':
        ep_dates = ep_data[['source name','trigger time']].dropna()
        ep_dates['trigger time'] = pd.to_datetime( ep_dates['trigger time'], format='mixed' )
        
        gw_dates = pd.DataFrame(crossmatch_time.get_skymapDates(), columns=['file name', 'trigger time'])
        gw_dates['trigger time'] = pd.to_datetime( gw_dates['trigger time'], format='mixed' )
    
        return ep_dates, gw_dates
    
    if retrieve=='space':
        ep_coords = ep_data[['source name','RA', 'DEC']].dropna()
        return ep_coords
    
def isPair_time(ep, gw):
    ep_dates, gw_dates = get_data(retrieve='time')
    dt = np.timedelta64(1, 'D')
    
    epT = ep_dates[ep_dates['source name'] == ep]['trigger time'].values
    gwT = gw_dates[gw_dates['file name'] == gw]['trigger time'].values
        
    return (gwT-dt <= epT) & (epT <= gwT+dt)

def isPair_space(ep, gw, condition):
    ep_coords = get_data(retrieve='space')
    ep_coords = np.array(ep_coords[['RA', 'DEC']][ep_coords['source name']==ep]) * units.deg
    
    gw_filepath = r"/Users/mck/Desktop/palmese_research/crossmatch_skymaps/" + gw
    results = crossmatch(sky_map=read_sky_map(filename=gw_filepath, moc=True), 
                                              coordinates=ep_coords)
    
    prob = results.searched_prob[0]
    return b(prob, *condition)

def b(p, o, l):
    if o=='==':
        return p==l
    
    elif o=='>':
        return p>l
    
    elif o=='<':
        return p<l

def tests():
    f = isPair_time
    checkFrame = [ [('EP241103a' , 'S241102cy_4_Bilby.offline0.multiorder.fits,0'), True],
                   [('EP240919a' , 'S241102cy_4_Bilby.offline0.multiorder.fits,0'), False],
                 ]
    
    for i in range(len(checkFrame)):
        args, a = checkFrame[i]
        q = f(*args)[0]
        print( f"{f.__name__} #{i+1}: {'passed' if q==a else 'failed' }")
    
    ###############################################################################
    
    f = isPair_space
    checkFrame = [ [('EP241103a' , 'S241102cy_4_Bilby.offline0.multiorder.fits,0'), False, ('<', 0.9)],
                   [('EP240919a' , 'S241102cy_4_Bilby.offline0.multiorder.fits,0'), False, ('<', 0.9)],
                 ]
    
    for i in range(len(checkFrame)):
        args, a, condition = checkFrame[i]
        q = f(*args, condition)
        print( f"{f.__name__} #{i+1}: {'passed' if q==a else 'failed' }")
    
tests()
