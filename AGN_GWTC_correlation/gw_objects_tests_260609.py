from gw_class import GW, contour_union, contour_intersection
from astropy.coordinates import SkyCoord
from astropy import units
from matplotlib import pyplot as plt
import ligo.skymap.plot

name1 = "GW190424_180648"
name2 = "GW230529_181500"
name3 = "GW231102_071736"

gws = (
       GW(name=name1, find_valid_region=False),
       GW(name=name3, find_valid_region=False),
       GW(name=name2, find_valid_region=False),
       )

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

ang_to_pixinfo_test()