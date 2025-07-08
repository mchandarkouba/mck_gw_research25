"""
ENVIRONMENT DETAILS:
    - python=3.11 with conda
    - install***: gcn-kafka, slack_sdk, astropy, astropy-healpix, healpy,*** ligo.skymap (healpy BEFORE ligo.skymap if on macbook), lxml
        FOR MAC: all packages EXCEPT healpy can be installed with pip: Use conda install healpy FOLLOWED by pip install ligo.skymap

This program is intended to be used for crossmatching Einstein Probe WXT alerts with the GraceDB gravitational wave event catalog.
It finds events each that are temporally & spatially coincident, then exports a matplotlib skymap of the gravitational wave PDF and 
Einstein probe location to Slack.
    
"""

import time
import json
from os import path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import subprocess
from datetime import datetime, timedelta
from io import BytesIO

from skymap2_VERA import get_skymap #MODULE FOUND HERE https://github.com/skyportal/gwemopt/blob/2b60b391c84e8f067e50fd7b1c9439393c14652a/gwemopt/io/skymap.py#L51

from gcn_kafka import Consumer
import slack_sdk as slack

from astropy.coordinates import SkyCoord
from astropy import units
from astropy.table import QTable
import astropy_healpix as ah

from ligo.gracedb.rest import GraceDb
from ligo.skymap.io import fits as fits_ligo
from ligo.skymap.postprocess.crossmatch import crossmatch
from ligo.skymap.plot.marker import sun

###############################################################################

TOP_DIR = f'{path.dirname(__file__)}/'
CRED_DIR = TOP_DIR
SKYMAP_DIR = TOP_DIR + 'alert_skymaps/'

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
            
            # find the ith pixel of the sorted pixels where 0th through ith sum cumulatively to provability c
            i = cumprob.searchsorted(c)
            
            #adding 0th through ith pixels' areas to find area of c% confidence region
            area_c = pix_area[:i].sum()
            self.area = area_c.to(units.deg**2)
                
            #Creating, flattening, then saving a skymap
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
        
###############################################################################

def loadCreds(path, keys):
    """
    Takes a .txt path and a list of keys, then matches the lines of the text to the keys by order
    Used for token credentials
    
    """

    creds = {}
    with open(path, 'r') as file: data = file.read().split()
    for n,s in enumerate(data): creds[keys[n]] = s
    
    return creds

def crossmatch_graceDB_time(ep):
    """
    Combines the simplified functions of crossmatch_graceDB and crossmatch_time in crossmatch_analysis_scripts dir
    Looks for GW files with a flase alarm rate less than 1/yr and within 1 day of the EP transient being recorded
    Uses altered skymap2_VERA module, taken from github
    
    """
    
    client = GraceDb()

    t, dt = ep['trigger time'], timedelta(days=1)
    far = 3.175e-8 #1/yr in Hz #False Alarm Rate of a grb event, in this case we expect events like this to not be significant 
    query = f'FAR < {far} created: {t-dt} .. {t+dt}'
    
    event_iter = client.superevents(query=query)
    superevents = [event['superevent_id'] for event in event_iter]
    superevents = [GW(get_skymap(graceID)) for graceID in superevents]
        
    return superevents

def crossmatch_space(ep, gw):
    """
    Uses ligos crossmatch to find the associated confidence area the EP event is found within for the coincident GW event
    
    """
    
    results = crossmatch(sky_map = gw.fits, 
                         coordinates=ep['ICRS'])
        
    return results.searched_prob

def followUp(ep, cands): 
    if cands==[]: 
        messageSlack(ep, None)
    
    else:
        for gw in cands: gw.img = skymapImg(ep, gw)
        messageSlack(ep, cands)
        
def skymapImg(ep, gw):
    """
    Plotting and saving the skymap image as a BytesIO object, which can be read by Slack
    
    """
    
    fig = plt.figure(figsize=(6,6), dpi=150)
    ax = plt.axes([0.05, 0.05, 0.9, 0.9], #left, bottom, width, height as fractions of figure
         projection='astro hours mollweide')
    ax.grid()
    ax.imshow_hpx(gw.path_flat, cmap='cylon')
    
    levels = (0.1, 0.5, 0.9)
    regions = []
    for n, c in enumerate(levels): 
        region = gw.ConfidenceRegion(c)
        regions += [region]
        #ax.imshow_hpx(region.path_flat, cmap='cylon')
        ax.contour_hpx(region.path_flat, levels=0, colors='black', alpha=1-n/len(levels), linewidths=0.5)
        
    for i, region in enumerate(regions): ax.text(0,-30*(i+1.5), f'{np.round(region.level*100)}% confidence region: {region.roundArea()}')
    ax.set_title(f"Skymap PDF for event {gw.name[:9]} \nwith {ep['source name']} superimposed")    
    
    ax.plot(#PLOTS RETICLE
            ep['ICRS'].ra.deg,
            ep['ICRS'].dec.deg,
            transform=ax.get_transform('world'),
            markersize=10,
            markeredgewidth=2,
            marker=sun
            )
    
    buffer = BytesIO()
    plt.savefig(buffer, format='png')
    buffer.seek(0)
    plt.close()
    
    return buffer.getvalue()    

def messageSlack(ep, cands):
    """
    Sends info to Slack app, to be posted into a channel with notice
    
    """
    
    creds = loadCreds(CRED_DIR + 'crossmatchAlert_bot_credentials.txt',
                      ['app id', 'client id', 'client secret', 'signing secret', 'verification token', 'oauth token', '#max-test id'] )
    
    client = slack.WebClient(token=creds['oauth token'])
    
    if cands==[]:
        client.chat_postMessage(channel=creds['#max-test id'],
                                text=f"EP event {ep['source name']} detected on {ep['trigger time']}, but no coincident GW events were found.")
        
    else:    
        for gw in cands:
            client.chat_postMessage(channel=creds['#max-test id'],
                                    text=f"EP event {ep['source name']} detected coincident with GW event {gw.name} within the GWs {np.round(gw.prob*100)}% confidence region...")
            
            if True: # <--for troubleshooting #gw.prob <= 0.9: # <--actual condition #
                client.files_upload_v2(channel=creds['#max-test id'],
                                       content=gw.img
                                      )
            else:
                client.chat_postMessage(channel=creds['#max-test id'],
                                        text='Due to the associated confidence area being too high (>90%), it is unlikely they are from the same event.')
                
            time.sleep(1)

###############################################################################

def main(): 
    """
    Logs into GCN notices, reads in EP events, crossmatches and sends message to Slack
    
    """
    
    creds = loadCreds(CRED_DIR + 'gcn_credentials.txt', #credentials need to be changed
                      ['id', 'secret'] )

    # Connect as a consumer (client "Max")
    consumer = Consumer(client_id=creds['id'],
                        client_secret=creds['secret'])

    # Subscribe to topics and receive alerts
    consumer.subscribe([
                        #'gcn.circulars',
                        #'igwn.gwalert',
                        #'gcn.notices.swift.bat.guano',
                        'gcn.notices.einstein_probe.wxt.alert',
                        'gcn.heartbeat'
                        ])
    
    id_letter = ''
    last_trigger = ''
    
    while True:
        for message in consumer.consume(timeout=1):
            
            if message.error():
                print(message.error())
                continue
    
            # Print the topic and message ID
            print(f'topic={message.topic()}, offset={message.offset()}')
            
            #Information Extract
            value = message.value()
            value = value.decode('utf-8') #bytes --> string
            value = json.loads(value) #string --> json
            
            ep = {#Data needed for crossmatching and documentation
                    'source name':None, # <-- this one will take an f-string of EP + formatted date + index letter
                    'ICRS':SkyCoord(value['ra'], value['dec'], unit=units.degree),
                    'trigger time':pd.to_datetime(value['trigger_time'], format='mixed'),
                    }
            
            ep_date = ep['trigger time'].strftime('%Y:%m:%d')
            id_letter = chr(ord(id_letter)+1) if last_trigger==ep_date else 'a'
            last_trigger = ep_date
                
            ep['source name'] = 'EP' + ep['trigger time'].strftime('%Y%m%s')[2:] + id_letter
            
            cands = crossmatch_graceDB_time(ep)
            for gw in cands: gw.prob = crossmatch_space(ep, gw)

            followUp(ep, cands)
    
def tests(use_preset:bool):
    """
    For latter sections of the main(), using pre-crossmatched EP & GW events
    
    """

    if use_preset:
        ep = {#A previously crossmatched object which is temporally coincident, but not spatially.
                'source name':'EP241126a',
                'ICRS':SkyCoord(300.97,-68.777, unit=units.deg),
                'trigger time':datetime.strptime('2024-08-07 05:42:48', '%Y-%m-%d %H:%M:%S')
                }
    
    else: #for the same as above, use EP241126a,300.97,-68.777,2024-08-07 05:42:48
        keys = ('source name', 'RA', 'DEC', 'trigger time')
        vals = input(f'enter params in this format: {keys}\n').split(',')
        ep = {k:vals[n] for n,k in enumerate(keys)}
        
        ep['ICRS'] = SkyCoord(float(ep['RA']), float(ep['DEC']), unit=units.deg)
        ep['trigger time'] = pd.to_datetime(ep['trigger time'], format='mixed')
                   
    cands = crossmatch_graceDB_time(ep)
    for gw in cands: gw.prob = crossmatch_space(ep, gw)

    followUp(ep, cands)

tests(use_preset=False) #