from gw_class import GW, contour_union, contour_intersection
from astropy.coordinates import SkyCoord
from astropy import units
from matplotlib import pyplot as plt
import ligo.skymap.plot

name1 = "GW190424_180648"
name2 = "GW230529_181500"
name3 = "GW231102_071736"

gws = [
        GW(name=name1, find_valid_region=False),
        GW(name=name3, find_valid_region=False),
        GW(name=name2, find_valid_region=False),
       ]

###############################################################################

markers1 = {
            "Milliquas": SkyCoord(ra=[0, 90, 180, 270, 360]*units.deg, 
                                  dec=[0, 0, -30, 0, 90]*units.deg, 
                                  frame='icrs',
                                  ),
            "Quaia": SkyCoord(ra=[242, 145, 83, 97, 5]*units.deg, 
                                  dec=[-70, 55, -32, 30, 0]*units.deg, 
                                  frame='icrs',
                                  ),
          }

markers2 = [
            SkyCoord(ra=[0, 90, 180, 270, 360]*units.deg, 
                     dec=[0, 0, -30, 0, 90]*units.deg, 
                     frame='icrs',
                     ),
            
            SkyCoord(ra=[242, 145, 83, 97, 5]*units.deg, 
                     dec=[-70, 55, -32, 30, 0]*units.deg, 
                     frame='icrs',
                     ),
            
            SkyCoord(ra=[100, 110, 120, 130, 140]*units.deg, 
                     dec=[-50, -40, -30, -20, -10]*units.deg, 
                     frame='icrs',
                     ),
            
            SkyCoord(ra=[190]*units.deg, 
                     dec=[20]*units.deg, 
                     frame='icrs',
                     ),
            ]

###############################################################################

def contour_union_test():
    level = 0.5 
    for gw in gws: gw.ConfidenceRegion(level)
    
    #for gw in gws: gw.matplotlib_plot()
    
    union_num = contour_union(gws,
                          constant_level=level,
                          sum_type="prob",
                          nside_final=2**5,
                          )
    
    union_bin = contour_union(gws,
                          constant_level=level,
                          sum_type="binary",
                          nside_final=2**5,
                          )
    
    fig = plt.figure(figsize=(6,6), dpi=150)
    ax = fig.add_subplot(111, projection="astro degrees mollweide",)
    
    ax.grid()
    ax.contour_hpx(union_bin, linewidths=0.5, colors="black")
    ax.imshow_hpx(union_num, cmap='cylon')

###############################################################################

def ang_to_pixinfo_test():
    gw = gws[0]
    ra1, dec1 = 100*units.degree, -10*units.degree
    coord1 = (ra1, dec1)
    coord2 = SkyCoord(ra1, dec1)
    print(gw.ang_to_pixinfo(coord2))
        
###############################################################################

def corr_candidates_test():
    events = [GW(name="GW150914"),
              GW(path="/Users/mck/Desktop/palmese_research/corr_scripts/GWTC_5.0/GWTC_4.0/O3a_skymaps/IGWN-GWTC2p1-v2-GW151012_095443_PEDataRelease_cosmo_reweight_C01:Mixed.fits")
              ] + gws
    events = events[::-1]
    c = 0.95
    for gw in events: gw.ConfidenceRegion(0.95)
    
    union = contour_union(events,
                              constant_level=c,
                              nside_final=2**5,
                              sum_type="prob",
                              )

    fig = plt.figure(figsize=(6,6), dpi=150)
    ax = fig.add_subplot(111, projection="astro degrees mollweide",)
    
    ax.grid()
    ax.imshow_hpx(union, cmap="cylon")

def nan_events_plot():
    events = [GW(name="GW200322_091133"), 
              GW(name="GW200308_173609"), 
              GW(name="GW231123_135430"),
              ]
    for gw in events: 
        gw.ConfidenceRegion(0.95)
        gw.matplotlib_plot()

nan_events_plot()
