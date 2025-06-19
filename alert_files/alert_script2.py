import time
import json
import numpy as np
from datetime import datetime, timedelta
from gcn_kafka import Consumer
from matplotlib import pyplot as plt
from io import BytesIO
import slack_sdk as slack
import subprocess
from skymap2_VERA import get_skymap
from astropy.coordinates import SkyCoord
from astropy import units
from ligo.gracedb.rest import GraceDb
from ligo.skymap.io.fits import read_sky_map
from ligo.skymap.postprocess.crossmatch import crossmatch
from ligo.skymap.plot.marker import sun

###############################################################################

SKYMAP_DIR = '/Users/mck/Desktop/palmese_research/alert_skymaps/'
IMAGE_DIR = '/Users/mck/Desktop/palmese_research/alert_images/'
CRED_DIR = '/Users/mck/Desktop/palmese_research/alert_files/'

###############################################################################

class GW:
    def __init__(self, name:str):
        path_dir = SKYMAP_DIR
        
        self.name = name
        self.path = path_dir + self.name
        
        self.name_flat = self.name.replace('.multiorder', '.flattened')
        self.path_flat = self.path.replace('.multiorder', '.flattened')
        subprocess.run(['ligo-skymap-flatten', self.path, self.path_flat])

        self.fits = read_sky_map(self.path, moc=True)
        self.fits_flat = read_sky_map(self.path_flat)
        
        self.prob = None
        self.img = None
        
    def __repr__(self):
        return self.name

###############################################################################

def loadCreds(path, keys):
    creds = {}
    with open(path, 'r') as file: data = file.read().split()
    for n,s in enumerate(data): creds[keys[n]] = s
    
    return creds

def crossmatch_graceDB_time(ep):
    client = GraceDb()

    t, dt = ep['trigger time'], timedelta(days=1)
    far = 3.175e-8 #1/yr in Hz #False Alarm Rate of a grb event, in this case we expect events like this to not be significant 
    query = f'FAR < {far} created: {t-dt} .. {t+dt}'
    
    event_iter = client.superevents(query=query)
    superevents = [event['superevent_id'] for event in event_iter]
    superevents = [GW(get_skymap(graceID)) for graceID in superevents]
        
    return superevents

def crossmatch_space(ep, gw):
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
    fig = plt.figure(figsize=(6,6), dpi=75)
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
    
    savepath = f"{IMAGE_DIR}{ep['source name']}_{gw.name_flat}.png"
    plt.savefig(savepath)
    plt.close()
    
    return savepath

def messageSlack(ep, cands):
    
    creds = loadCreds(CRED_DIR + 'crossmatchAlert_bot_credentials.txt', #credentials address need to be changed
                      ['app id', 'client id', 'client secret', 'signing secret', 'verification token', 'oauth token', '#max-test id'] )
    
    client = slack.WebClient(token=creds['oauth token'])
    
    if cands==[]:
        client.chat_postMessage(channel=creds['#max-test id'],
                                text=f"EP event {ep['source name']} detected on {ep['trigger time']}, but no coincident GW events were found.")
        
    else:    
        for gw in cands:
            client.chat_postMessage(channel=creds['#max-test id'],
                                    text=f"EP event {ep['source name']} detected coincident with GW event {gw.name} within the GWs {np.round(gw.prob*100)}% confidence region...")
            
            if True: #prob <= 0.9:
                client.files_upload_v2(channel=creds['#max-test id'],
                                       file=gw.img
                                      )
            else:
                client.chat_postMessage(channel=creds['#max-test id'],
                                        text='Due to the associated confidence area being too high (>90%), it is unlikely they are from the same event.')
                
            time.sleep(1)
###############################################################################

def main(): 
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

    ep = {#Data needed for crossmatching and documentation
            'source name':'EP241101a',
            'ICRS':SkyCoord(37.763, 22.731, unit=units.degree),
            'trigger time':datetime.strptime('2024-11-01 23:54:17', '%Y-%m-%d %H:%M:%S')
            }
        
    cands = crossmatch_graceDB_time(ep)
    for gw in cands: gw.prob = crossmatch_space(ep, gw)

    followUp(ep, cands)

tests() #main()