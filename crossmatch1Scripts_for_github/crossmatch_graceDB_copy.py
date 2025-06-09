#Searches through GraceDB events for GW events within the date range between the 
#earliest Einstein Probe event of interest (minus one day), to the latest EP event (plus one day);
#Secondary search criteria is for the False Alarm Rate for GW events to be less than 1/yr
#Selected multiorder skymap files are saved locally for further crossmatching

#note: an environment with pandas, ligo was used for this program.
#Additionally, the module skymap2 was saved/modified locally and can be found here: 
    #https://github.com/skyportal/gwemopt/blob/2b60b391c84e8f067e50fd7b1c9439393c14652a/gwemopt/io/skymap.py#L51

############################################################################### Imports and data

import urllib.request
import pandas as pd
import numpy as np

from ligo.gracedb.rest import GraceDb
from skymap2 import get_skymap
import datetime

csv_path = r"/Users/mck/Desktop/palmese_research/crossmatch_analysis_scripts/EP_data1.csv" #locally saved dataset of EP events that are of interest
ep_data = pd.read_csv(csv_path, comment='#') #comment='#' removes the units line of the csv from docs

############################################################################### GraceDB catalog searching and selecting

client = GraceDb()

dates = pd.to_datetime( ep_data['trigger time'].dropna(), format='mixed' )
dt = datetime.timedelta(days=1)
first, last = min(dates)-dt, max(dates)+dt
far = 3.175e-8 #1/yr in Hz #False Alarm Rate of a grb event, in this case we expect events like this to not be significant 
query = f'FAR < {far} created: {first} .. {last}'
event_iter = client.superevents(query=query)

superevents = [event['superevent_id'] for event in event_iter]

for graceid in superevents: 
    try: 
        get_skymap(graceid)
        print()
    except ValueError:
        print(f'Event {graceid} was skipped for raising an error, and was likely retracted')
