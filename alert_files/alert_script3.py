"""
ENVIRONMENT DETAILS:
    python=3.11
    install: gcn-kafka, slackclient, slack_sdk astropy, healpy,*** ligo.skymap (healpy BEFORE ligo.skymap if on macbook), lxml
    ***FOR MAC: all packages EXCEPT healpy can be installed with pip: Use conda install healpy FOLLOWED by pip install ligo.skymap
    
"""

import time
import json
import numpy as np
from matplotlib import pyplot as plt
import subprocess
from datetime import datetime, timedelta
from io import BytesIO

from skymap2_VERA import get_skymap #MODULE FOUND HERE https://github.com/skyportal/gwemopt/blob/2b60b391c84e8f067e50fd7b1c9439393c14652a/gwemopt/io/skymap.py#L51

from gcn_kafka import Consumer
import slack_sdk as slack

from astropy.coordinates import SkyCoord
from astropy import units

from ligo.gracedb.rest import GraceDb
from ligo.skymap.io.fits import read_sky_map
from ligo.skymap.postprocess.crossmatch import crossmatch
from ligo.skymap.plot.marker import sun

###############################################################################

SKYMAP_DIR = '/Users/mck/Desktop/palmese_research/alert_skymaps/'
CRED_DIR = '/Users/mck/Desktop/palmese_research/alert_files/'

###############################################################################

class GW:
    """
    Stores gravitational wave event files, directories easily
    _flat denotes a flattened healpix skymap, mainly for matplotlib
    
    """
    
    def __init__(self, name:str):
        path_dir = SKYMAP_DIR
        
        self.name = name
        self.path = path_dir + self.name
        
        self.name_flat = self.name.replace('.multiorder', '.flattened')
        self.path_flat = self.path.replace('.multiorder', '.flattened')
        subprocess.run(['ligo-skymap-flatten', self.path, self.path_flat]) #uses cmd to flatten the skymap from a filepath

        self.fits = read_sky_map(self.path, moc=True)
        self.fits_flat = read_sky_map(self.path_flat)
        
        self.prob = None
        self.img = None
        
    def __repr__(self):
        return self.name

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
            
            if True: # <--for troubleshooting #if gw.prob <= 0.9: # <--actual condition
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
                        'gcn.notices.einstein_probe.wxt.alert'
                        ])
    
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
                    'trigger time':datetime.datetime(value['trigger_time'], format='mixed'),
                    }
            
            cands = crossmatch_graceDB_time(ep)
            for gw in cands: gw.prob = crossmatch_space(ep, gw)

            followUp(ep, cands)
    
def tests():
    """
    For latter sections of the main(), using pre-crossmatched EP & GW events
    
    """

    ep = {#A previously crossmatched object which is temporally coincident, but not spatially.
            'source name':'EP241126a',
            'ICRS':SkyCoord(281.33,-13.167, unit=units.degree),
            'trigger time':datetime.strptime('2024-09-18 18:06:47', '%Y-%m-%d %H:%M:%S')
            }
        
    cands = crossmatch_graceDB_time(ep)
    for gw in cands: gw.prob = crossmatch_space(ep, gw)

    followUp(ep, cands)

tests() #main()