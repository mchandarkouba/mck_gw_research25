from gcn_kafka import Consumer
import json
import datetime
from astropy.coordinates import SkyCoord
import pandas as pd

filepath = 

# Connect as a consumer (client "Max")
consumer = Consumer(client_id=,
                    client_secret=)

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
        
        data = {#Data needed for crossmatching and documentation
                'source name':None, # <-- this one will take an f-string of EP + formatted date + index letter
                'RA':value['ra'], 'DEC':value['dec'],
                'trigger time':datetime.datetime(value['trigger_time']),
                }
        
        