from gcn_kafka import Consumer
import json
import datetime
from astropy.coordinates import SkyCoord

from ligo.gracedb.rest import GraceDb
from skymap2_NERA import get_skymap

from ligo.skymap.postprocess.crossmatch import crossmatch
from ligo.skymap.io.fits import read_sky_map

filepath = '/Users/mck/Desktop/palmese_research/alert_files/gcn_credentials.txt' #credentials need to be changed
with open(filepath, 'r') as file: creds = file.read().split()
for n,s in enumerate(creds): creds[n] = s[1:-1]

###############################################################################

def crossmatch_graceDB_time(ep):
    client = GraceDb()

    t, dt = ep['trigger time'], datetime.timedelta(days=1)
    far = 3.175e-8 #1/yr in Hz #False Alarm Rate of a grb event, in this case we expect events like this to not be significant 
    query = f'FAR < {far} created: {t-dt} .. {t+dt}'
    
    event_iter = client.superevents(query=query)
    superevents = [event['superevent_id'] for event in event_iter]
    for graceID in superevents: get_skymap(graceID)
    
    return superevents

def crossmatch_space(ep, gw):
    gw_filepath = r'/Users/mck/Desktop/palmese_research/alert_skymaps/' + gw # <-- NEEDS TO BE CHANGED FOR NERA
    
    results = crossmatch(sky_map = read_sky_map(filename=gw_filepath, moc=True), 
                                                  coordinates=ep['ICRS'])
        
    return results.searched_prob

def interpret_probs(ep, probs):
    pass

###############################################################################

# Connect as a consumer (client "Max")
consumer = Consumer(client_id=creds[0],
                    client_secret=creds[1])

# Subscribe to topics and receive alerts
consumer.subscribe([
                    #'gcn.circulars',
                    #'igwn.gwalert',
                    #'gcn.notices.swift.bat.guano',
                    'gcn.notices.einstein_probe.wxt.alert'
                    ])

def main():
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
            
            source = {#Data needed for crossmatching and documentation
                    'source name':None, # <-- this one will take an f-string of EP + formatted date + index letter
                    'ICRS':SkyCoord(value['ra'], value['dec'], unit='deg'),
                    'trigger time':datetime.datetime(value['trigger_time'], format='mixed'),
                    }
            
            candidates = crossmatch_graceDB_time(source)
            results = {gw:crossmatch_space(source, gw) for gw in candidates}
            
            interpret_probs(source, results)
    
        
        