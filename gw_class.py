###############################################################################
### INFO
###############################################################################
"""
- python=3.11 conda environment
- install***: astropy, astropy-healpix, healpy, ligo.skymap, pandas

    *** All packages EXCEPT healpy installed with pip: 
        used conda install healpy FOLLOWED by pip install ligo.skymap
        
- download skymap module (linked below) off Github
_______________________________________________________________________________

This package was made for gravitational wave astronomy research.
The main objects, GW and ConfidenceRegion, load Ligo-Virgo-KAGRA FITS files.
Several functions are defined below object classes.

"""
###############################################################################
### IMPORTS
###############################################################################

from astropy import units
from astropy.coordinates import SkyCoord
from astropy.table import QTable
from astropy.cosmology.realizations import Planck15
import astropy_healpix as ah

from datetime import datetime

import healpy as hp

import ligo.skymap
import ligo.skymap.plot
from ligo.skymap.io import fits as fits_ligo
from ligo.skymap.postprocess.crossmatch import crossmatch
from ligo.gracedb.rest import GraceDb

from matplotlib import pyplot as plt
from matplotlib.patches import Patch

import numpy as np

import os

from scipy.stats import norm as normal_distribution

from skymap import get_skymap
# module link: https://github.com/skyportal/gwemopt/blob/2b60b391c84e8f067e50fd7b1c9439393c14652a/gwemopt/io/skymap.py#L51

import subprocess

###############################################################################
### GLOBALS
###############################################################################

default_hierarchy = ['gstlal', 'NRSur7dq4', 'Mixed', 'sco']
### algorithm vs template hierarchies?
default_cosmo = Planck15

###############################################################################
### DIRECTORIES
###############################################################################

TOP_DIR = f'{os.path.dirname(__file__)}/'
SKYMAP_DIR = TOP_DIR # Assumes that skymaps are in some directory or sub-directory of the current one
MOD_DIR = TOP_DIR + "modified_skymaps/"

###############################################################################
### GW CLASS
###############################################################################

class GW:
    ############################################################
    
    def __init__(self, 
                 path=None, # path or name required
                 name=None,
                 skymaps=None,
                 find_valid_region=True, 
                 preferences=default_hierarchy,
                 verbose=True,
                 ):
        """
        Parameters
        ----------
        # path or name required                 
        path : str, optional
            A path to the GWEvent .fits skymap, accepts multiordered files.
        name : "GWddmmyy" or "GWddmmyy_hhmmss", optional
            The GWEvent name. GW.__init__ will search for .fits filenames containing name in SKYMAP_DIR
        skymaps : iterable, optional
            An iterable containing .fits filepaths. 
            Can serve as a cache so recursive_all_skymaps() is called once
            When kwarg name is provided.
            
        find_valid_region : bool, optional
            When True, __init__ will try to generate test confidence regions for skymaps on file,
            until a valid file is found or an IndexError arises (from lack of valid skymaps on file).
            When False, __init__ will make self.fits the first skymap by preferences, regardless of potential errors.
            The default is True when name is given, and is always False when path is given.
            
        preferences : list, optional
            An ordered list of preferred filetypes. 
            The initialized GW object will search for preferred filetypes first, in order of preference.
            The default is ['NRSur7dq4', 'Mixed', 'sco'].

        Raises
        ------
        TypeError: invalid filepath or name is given (no matching skymaps in SKYMAP_DIR or sub-directories).
        Exception: various errors creating GW.ConfidenceRegion object (corrupt file, missing data)
        
        Returns None
        
        Attributes
        ----------
        self.path: str, filepath to skymap file
        self.name: str, "GWddmmyy" or "GWddmmyy_hhmmss"
        
        list() on self.__init__ :
            self.regions: list, contains all associated ConfidenceRegion objects
        
        None on self.__init__ :
            self.path_flat: str, filepath to flattened skymap file. Call self.flatten()
            self.fits: loaded and read self.path. Call self.read_fits()
            self.fits_flat: loaded and read self.path_flat. Call self.read_fits(flat=True)
            self.plt_tup: tuple, contains (fig,ax). Call self.matplotlib_plot()
        """
        self.verbose = verbose
        cands = []
        if path!=None:
            if ("_flattened" in path) or ("percent_" in path): print("GW path contains '_flattened' or 'percent_', this file may cause errors.")
            cands = [path]
            find_valid_region = False
        
        elif name!=None:
            cands = retrieve_skymaps(name, 
                                     hierarchy=preferences, # optional kwarg, default is None
                                     skymaps=skymaps,
                                     )
        if cands==[]:
            print("No path was given or no skymap matching that name could be found.\nCheck the input name/path, or the SKYMAP_DIR variable in this package.")
            raise TypeError
        
        i=0
        while True:
            self.path = cands[i]
            self.name = name or get_gw_name(path) or self.path
            self.path_flat = None
            self.fits = None
            self.fits_flat = None
            self.regions = []
            
            if find_valid_region:
                try:
                    self.ConfidenceRegion(0.50)
                    break
                    
                except Exception:
                    print("\n\tFailed generating ConfidenceRegion object for this file, trying another...\n")
                    i+=1
                    continue
                    
            else:
                break
                
        self.regions = []
        self.plt_tup = (None, None) # (plt.fig, plt.ax)
        if self.verbose: print(f"\nUsing {self.path} \nFor given value of kwarg find_valid_region, filetype preferences, and input type\n")
        
        return
    
    ############################################################
    
    def __repr__(self):
        return self.name
    
    ############################################################
    
    def ConfidenceRegion(self, C, overwrite=False):
        """
        Documentation here: https://emfollow.docs.ligo.org/userguide/tutorial/multiorder_skymaps.html
        
        Parameters
        ----------
        C : float
            Associated confidence, a value in (0., 1.).
        name : str
            gw.name for associated GW.
        path : str
            gw.path for associated GW.

        Returns None
        
        Attributes
        -------
        self.level: float, associated confidence decimal
        self.area: float, area (in deg^2) of confidence region in the sky
        self.name: str, gw.name for associated GW + self.level
        
        self.path: str, filepath to CR skymap file
        self.path_flat: str, filepath to flattened CR skymap file
        self.fits: loaded and read self.path
        self.fits_flat: loaded and read self.path_flat
        """
        result = tuple()

        C = [C,] if isinstance(C, float) else C
        for region in self.regions: 
            if region.level in C: 
                result.append(self.get_region(level=region.level))
                C.remove(region.level)
                                                
                    
        skymap = QTable.read(self.path, format='fits')
        skymap.sort(['PROBDENSITY'], reverse=True)
        
        level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
        pix_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
        
        prob = pix_area * skymap['PROBDENSITY']
        cumprob = np.cumsum(prob)
        
        # find the ith pixel of the sorted pixels where 0th through ith sum cumulatively to probability c
        I = [cumprob.searchsorted(c) for c in C]
        
        #adding 0th through ith pixels' areas to find area of c% confidence region
        #AREA = [pix_area[:i].sum().to(units.deg**2) for i in I]
            
        #Creating, flattening, then saving a skymap
        SKYMAP = [skymap[:i] for i in I]
        #for skymap in SKYMAP: skymap.sort("UNIQ")
        
        #Making FITS skymap, flattening and saving data to the GW object
        NAME = [f"{self.name}_{np.round(100*c)}percent_" for c in C]
        PATH = [ self.path.replace(self.name, name) for name in NAME]        
        
        for n in range(len(PATH)):
            path, skymap, name, c = PATH[n], SKYMAP[n], NAME[n], C[n]
            
            if not os.path.exists(path) or overwrite:
                skymap.write(path, overwrite=True, format='fits')
                
            else:
                if self.verbose: print(f"kwarg overwrite was set False and the following files were found: \n\t{self.path} \n\t{self.path_flat}")
            
            kwargs = {"c":c,
                      "name":name,
                      "path":path,
                      }
        
            result = result + ( ConfidenceRegion(**kwargs), )
        
        self.regions += list(result)
        
        return result if len(result)>1 else result[0]    
    
    ############################################################
    
    def matplotlib_plot(self, markers=None, plot_regions=True):
        """
        Parameters
        ----------
        markers : dict or list or SkyCoord, optional
            Iterable or single inputs must contain all-SkyCoord values.
            any SkyCoord can contain multiple points,
            
            If an input (list, dict) contains multiple SkyCoords 
            their points will be color-grouped on the skymap.
            
            For dictionaries, SkyCoords will be labeled with their keys in the legend.
            For lists, SkyCoords will be labeled generically in the legend.
            The default is None.
        plot_regions : bool, optional
            When True, contours for confidence regions in self.regions are added to skymap plot. 
            The default is True.

        Returns
        -------
        fig, ax: matplotlib skymap figure and axes objects
        
        Plots a skymap with marker points overlaid 
        and optional confidence regions stored in self.regions
        """
        
        fig = plt.figure(figsize=(6,6), dpi=150)
        ax = fig.add_subplot(111, projection="astro degrees mollweide",)
        
        ax.grid()
        ax.imshow_hpx(self.path_flat, cmap='cylon')
        ax.set_title(self.name)
        
        patches = []
        if plot_regions:
            for n,region in enumerate(self.regions):
                col = (0.8/len(self.regions)) *n
                ax.contour_hpx(region.path_flat, 
                               levels=0, 
                               colors=[col, col, col, 1.0], 
                               linewidths=0.75,
                               )
                patches += [Patch(label=f"{region.level*100 //1}% confidence region",
                                  facecolor=(col, col, col, 1.0),
                                  edgecolor='black',
                                  )]
                
            legend_conf = fig.legend(handles=patches,
                                    loc="upper right",
                                    bbox_to_anchor=(0.275, 0.1, 0.2, 0.2),
                                    )
            ax.add_artist(legend_conf)
        
        mtype = type(markers)
        patches = []
        
        if mtype==SkyCoord:
            markers = [markers]
            mtype = list
        
        if mtype in (list, tuple, dict):
            items = [markers[k] for k in markers] if mtype==dict else [x for x in markers]
            if all(isinstance(x, SkyCoord) for x in items):
                for n,k in enumerate(markers):
                    label = k if mtype==dict else f"group {n+1}"
                    col = 1 - (1/len(self.regions))*n
                    coords = markers[k] if mtype==dict else k
                    
                    ra = coords.ra.wrap_at(180*units.deg) #.radian
                    dec = coords.dec #.radian 
                    ax.scatter(ra, dec, c=[[0., col/2, col, 1.]],
                               transform=ax.get_transform('world'),
                               s=10,)
                    
                    patches += [Patch(label=label,
                                      facecolor=(0., col/2, col, 1.),
                                      edgecolor='black',
                                      )]
                
                legend_markers = fig.legend(handles=patches,
                                        loc="upper left",
                                        bbox_to_anchor=(0.475, 0.1, 0.2, 0.2),
                                        )
                
                ax.add_artist(legend_markers)
                
            else:
                print("kwarg markers is not the proper SkyCoord object(s).")

        elif not mtype==type(None):
            print("kwarg markers is not the proper SkyCoord object(s).")
        
        self.plt_tup = fig,ax
        return (fig, ax,)
    
    ############################################################
    
    def crossmatch_space(self, coords:SkyCoord, cosmology=default_cosmo):
        """
        Parameters
        ----------
        coords : SkyCoord
            A SkyCoord (can contain multiple points) to be crossmatched with self.fits skymap.
            Points can be three-dimensional (contain a redshift coordinate)
        cosmology : from astropy.cosmology, optional
            An object containing cosmological parameters for distance calculations from redshift.
            Necessary for 3D crossmatching, not for 2D.
            The default package is Planck15.

        Returns
        -------
        list
            An ordered list of associated confidence levels for the points in coords.
        """
        if self.fits==None: self.read_fits()
        results = crossmatch(sky_map = self.fits, 
                            coordinates=coords,
                            cosmology=cosmology,
                            )
        return results.searched_prob
      
    def get_region(self, level=None):
        """
        Parameters
        ----------
        level : float, optional
            Confidence level in (0., 1.) to be searched for in self.regions. The default is None.

        Returns
        -------
        region : self.ConfidenceRegion
            region in self.regions meets kwargs
        
        None otherwise
        """
        
        for region in self.regions:
            if region.level==level: return region
    
    ############################################################
    
    def flatten(self):
        """
        Returns
        -------
        self.path_flat : str
            Filepath to a flattened skymap. Sets self.path_flat.
        """
        self.path_flat = self.path.replace('.fits', '_flattened.fits').replace(TOP_DIR, MOD_DIR)
        flat_dir = os.path.dirname(self.path_flat)
        if not os.path.exists( flat_dir ): os.makedirs(flat_dir, exist_ok=True)
        subprocess.run(['ligo-skymap-flatten', self.path, self.path_flat],) #uses cmd to flatten the skymap from a filepath
        
        return self.path_flat
    
    ############################################################
    
    def read_fits(self, flat=False):
        """
        Parameters
        ----------
        flat : bool, optional
            Read the flattened skymap or multiorder. The default is False.

        Returns
        -------
        fits : FITS file
            Skymap file, can be read using astropy QTable or other modules.
            Depending on **flat, saved to self.fits or self.fits_flat .
        """
        path=self.flatten() if flat and not os.path.exists(self.path_flat) else self.path
        fits = fits_ligo.read_sky_map(path, moc=not flat)
        
        if flat: self.fits_flat = fits
        else: self.fits = fits
        return fits
    
    ############################################################
    ### Documentation for the following functions (flat skymaps) found here: 
    ### https://iopscience.iop.org/article/10.3847/0067-0049/226/1/10/pdf
    ############################################################
    
    def ang_to_pixinfo(self, *args, return_kwargs=None,):
        """
        Parameters
        ----------
        *args : (ra, dec) or SkyCoord
            Must be a single coordinate, in degrees. 
        return_kwargs : None or str or iterable, optional
            When None, returns dict with keys "dp_dr", "dp_dA", "dp_dV", "pix_area", "pix_prob".
            When str, returns the value associated with that key of the full dict.
            When iterable of keys, returns the subset of full dict.
            The default is None.

        Returns
        -------
        dict or value.
        
        values :
            dp_dr or dp_dA or dp_dV : function
            pix_area : float, in astropy.units.degrees**2
            pix_prob : float
        """
        if len(args)==2: # ra, dec explicitly provided in degrees
            ra, dec = args
        
        elif len(args)==1 and isinstance(args[0], SkyCoord):
            coord = args[0]
            ra, dec = coord.ra, coord.dec
        
        table = QTable.read(self.path)
        
        LEVEL, IPIX = ah.uniq_to_level_ipix(table["UNIQ"])
        NSIDE = ah.level_to_nside(LEVEL)
        ipix = ah.lonlat_to_healpix(ra, dec, NSIDE, order="nested")
        i = np.where(IPIX==ipix)[0]
        
        result = dict()
        
        result["dp_dA"] = table["PROBDENSITY"][i].to(units.degree**-2)
        result["pix_area"] = ah.nside_to_pixel_area(NSIDE[i]).to(units.degree**2)
        result["pix_prob"] = result["prob_density"] * result["pix_area"]
        
        norm, mu, sigma = table["DISTNORM"][i], table["DISTMU"][i], table["DISTSIGMA"][i]
        result["dp_dr"] = lambda r: r**2 * norm*normal_distribution(mu,sigma).pdf(r)
        result["dp_dV"] = lambda r: table["PROBDENSITY"] * norm*normal_distribution(mu,sigma).pdf(r)
        
        if return_kwargs==None: return result
        elif isinstance(return_kwargs, str): return result[return_kwargs]
        else: return {k:result[k] for k in return_kwargs}
        
    ############################################################
    
    def marginal_dp_dr(self, R=None, suggest_bounds=False):
        """
        Parameters
        ----------
        R : float or np.array containing floats, optional
            input distance(s). Output will then be associated probabilities.
            The default is None.
        suggest_bounds : bool, optional
            If R==None, suggest_bounds=True will return a tuple (dp_dr(), (r0, r1))
            With suggested lower/upper bounds at position 1. 
            suggest_bounds=False will return only the dp_dr() function.
            The default is False.

        Returns
        -------
        (dp_dr(), (r0, r1)) or dp_dr() or dp_dr <-- array
            dp_dr() being the marginal distance probability function,
            averaged over the entire sky.
        """
        table = QTable.read(self.path)
        MU = table["DISTMU"] # in Mpc
        SIG = table["DISTSIGMA"] # in Mpc
        NORM = table["DISTNORM"] # in 1/Mpc^2
        DP_DA = table["PROBDENSITY"] # in 1/steradian
        
        LEVEL, IPIX = ah.uniq_to_level_ipix(table['UNIQ'])
        NSIDE = ah.level_to_nside(LEVEL)
        AREA = ah.nside_to_pixel_area(NSIDE) # in steradian
        
        min_r0, max_r1 = -1, -1
        dp_dr_pix_wlist = []
        
        for ipix in range(len(DP_DA)):
            
            mu, sig, norm = MU[ipix], SIG[ipix], NORM[ipix]
            dp_da, area = DP_DA[ipix], AREA[ipix]
            r0, r1 = mu-2.5*sig, mu+2.5*sig
            
            if r0==np.inf or r1==np.inf: pass
            elif min_r0==-1: min_r0, max_r1 = max(0*units.Mpc, r0), max(0*units.Mpc, r1)
            else: min_r0, max_r1 = max(0*units.Mpc, min(r0, min_r0)), min(np.inf*units.Mpc, max(r1, max_r1))
            
            dp_dr_pix_wlist += [ lambda R: (dp_da * area) * R**2 * norm * normal_distribution(mu, sig).pdf(R) ]
        
        if R!=None: 
            return sum([dp_dr(R) for dp_dr in dp_dr_pix_wlist])
        
        else:
            result = lambda R: sum([dp_dr(R) for dp_dr in dp_dr_pix_wlist])
            return result, (min_r0, max_r1) if suggest_bounds else result
    
###############################################################################
### GW CONFIDENCE REGION CLASS
###############################################################################

class ConfidenceRegion:
    """
    See GW.ConfidenceRegion for generation function
    """
    ############################################################
    
    def __init__(self,
                 c=None,
                 name=None, 
                 path=None,
                 ):
        self.level = c
        self.name = name
        self.path = path
        self.path_flat = None
        self.fits = None
        self.fits_flat = None
        self.area = None
        
    ############################################################
    
    def get_area(self):
        """
        Returns self.area after calculating it.
        """
        table = QTable.read(self.path, format='fits')
        table.sort('PROBDENSITY', reverse=True)
        LEVEL, IPIX = ah.uniq_to_level_ipix(table['UNIQ'])
        AREA = ah.nside_to_pixel_area(ah.level_to_nside(LEVEL))
        PROB = AREA * table['PROBDENSITY']
        cumprob = np.cumsum(PROB)
        i = cumprob.searchsorted(self.level)
        self.area = AREA[:i].sum().to(units.deg**2)
        
        return self.area
    
    ############################################################
    
    def __repr__(self):
        return self.name
    
    ############################################################
    
    def flatten(self):
        """
        See GW.flatten()
        """
        self.path_flat = self.path.replace('.fits', '_flattened.fits').replace(TOP_DIR, MOD_DIR)
        flat_dir = os.path.dirname(self.path_flat)
        if not os.path.exists( flat_dir ): os.makedirs(flat_dir, exist_ok=True)
        subprocess.run(['ligo-skymap-flatten', self.path, self.path_flat]) #uses cmd to flatten the skymap from a filepath
        
        return self.path_flat
    
    ############################################################
    
    def read_fits(self, flat=False):
        """
        See GW.read_fits()
        """
        path=self.flatten() if flat and not os.path.exists(self.path_flat) else self.path
        fits = fits_ligo.read_sky_map(path, moc=not flat)
        
        if flat: self.fits_flat = fits
        else: self.fits = fits
        return fits
            
###############################################################################
### FUNCTIONS
###############################################################################

def get_gw_name(s):
    """
    Parameters
    ----------
    s : str
        A string (typically a path).

    Returns
    -------
    name : str
        GWEvent name, either "GWddmmyy" or "GWddmmyy_xxxxxx".
        (ddmmyy and xxxxxx are ints)
        
    None: filepath format does not contain GWddmmyy format.

    """
    for i in range(len(s)-1):
        c0 = s[i]
        c1 = s[i+1]
        if c0+c1 == 'GW':
            try:
                date = str(int(s[i+2:i+8])) #checking for int following "GW"
                superevent = c0+c1+date
                
                try:
                    time = str(int(s[i+9: i+15])) #checking for int following "GWxxxxxx_"
                    time = '0'*(6-len(time)) + time
                    superevent = superevent + '_' + time
                    return superevent
                    
                except Exception:
                    return superevent
            
            except Exception:
                continue
            
    return None

############################################################

"""
NOTE: gw_search_gracedb() is not tested.
"""
def gw_search_gracedb(query=None, 
                      max_far=None, 
                      creation_dates=None,
                      graceid_list=None,
                      superevent_list=None,
                      ):
    """
    Parameters
    ----------
    query : str, optional
        A complete query formatted properly for the GraceDB client. 
        The default is None.
    max_far : float, optional
        Include a maximum False Alarm Rate (units 1/yr) in the query. 
        The default is None.
    creation_dates : tuple, optional
        A 2-item tuple, being a date range to search the GraceDB catalog between.
        (date1, date2) must be datetime.datetime objects, or None (serving as unbounded on one or both ends).
        ie: (date1, None), (None, date2). (None, None) would return all objects in GraceDB.
        The default is None.
    graceid_list : list, optional
        A list of graceIDs to search, different from superevent IDs ("GWxxxxxx...). 
        The default is None.
    superevent_list : list, optional
        A list of superevents to search, retrieved superevents can contain multiple GWEvents.
        The default is None.

    Returns None <-- should be a list of GW names which can be individually initialized as GW objects

    OBJECTIVE: to save searched GWEvents on GraceDB into SKYMAP_DIR
    """
    print("GraceDB docs for creating a custom query (kwarg): \n\thttps://gracedb.ligo.org/documentation/queries.html")
    print("This function requires the skymap module found here: \n\thttps://github.com/skyportal/gwemopt/blob/2b60b391c84e8f067e50fd7b1c9439393c14652a/gwemopt/io/skymap.py#L51")
    print("\nMake sure the SKYMAP_DIR variable at the beginning of skymap.py matches SKYMAP_DIR in this one.")
    
    client = GraceDb()

    if query==None:
        qstring_list = []
        
        if max_far!=None:
            qstring_list.append(f'FAR < {max_far}')
            
        if creation_dates!=None:
            first, last = creation_dates
            first = first or datetime(2015, 9, 14).strftime("%Y-%m-%d %I:%M:%S")
            last = last or "today"
            qstring_list.append(f"created: {first} .. {last}")
            
        if graceid_list!=None:
            qstring_list.append((" ").join(graceid_list))
            
        if superevent_list!=None:
            qstring_list.append((" superevent_id: ").join(superevent_list))
    
        qstring = (' ').join(qstring_list)
        
    else:
        qstring = query
        
    event_iter = client.superevents(query=qstring)
    superevents = [event['superevent_id'] for event in event_iter]

    for graceid in superevents: 
        try: 
            get_skymap(graceid)
        except ValueError:
            print(f'Event {graceid} was skipped for raising an error, and was likely retracted')

    return
    
############################################################

def retrieve_skymaps(name, hierarchy=None, skymaps=None):
    """
    Searching function with helper recursive_all_skymaps
    Used when for GW.__init__ is passed a value for name
    
    Parameters
    ----------
    name : str
    hierarchy : list, optional
        preferences. The default is None.

    Returns
    -------
    cand_list : list
        list of .fits filepaths.
    """
    hierarchy = hierarchy or []
    all_skymaps = skymaps or recursive_all_skymaps(SKYMAP_DIR)
    cand_dict = {None:[]} | {k:[] for k in hierarchy}
    
    for s in all_skymaps:
        if (name in s) and ("_flattened" not in s) and ("percent_" not in s):
    
            if hierarchy!=[]:
                in_hierarchy = False
                
                for fitsType in hierarchy:
                    if fitsType in s:
                        in_hierarchy = True
                        cand_dict[fitsType].append(s)
                
                if not in_hierarchy: cand_dict[None].append(s)
                
        
            else:
                cand_dict[None].append(s)
    
    cand_list = []
    for key in hierarchy+[None,]:
        cand_list += cand_dict[key]
    
    return cand_list


############################################################
    
def recursive_all_skymaps(ldir):
    result = []
    for s in os.listdir(ldir):
        path = ldir+s
        
        if '.fits' in s:
            result += [path]
            
        else:
            try: 
                result += recursive_all_skymaps(path+'/')
            except NotADirectoryError:
                continue
            
    return result

############################################################

def contour_union(gws,
                  nside_final=2**5,
                  constant_level=None,
                  index_level=None,
                  sum_type="binary", # "prob", "number", "binary"
                  verbose=False,
                  ):
    """
    Parameters
    ----------
    gws : iter(GW) or iter(ConfidenceRegion)
        Iterable containing GW or ConfidenceRegion objects.
        If GW objects, then constant_level or index_level kwarg is required.
    nside_final : int, optional
        Pixel resolution of output skymap array, must be 2**n (n:int < 52) The default is 2**5.
    constant_level : float, optional
        Confidence level in (0., 1.). If specified, will look for a ConfidenceRegion of level confidence_level in gw.regions to use. 
        The default is None.
    index_level : dict, optional
        Maps gw.name for gw in gws to a specific confidence level in gw.regions. 
        The default is None.
    sum_type : str, optional
        When "binary", all included pixel values are 1. , excluded 0. . Recommended when passing into contour_hpx.
        When "freq", included pixels are int>0 for each ConfidenceRegion overlapped, excluded 0 .
        When "prob", pixels are the sum of all probdensity functions, re-normalized.
        The default is "binary".
    verbose : bool, optional
        Prints some information while running. The default is True.

    Raises
    ------
    ValueError
        levels given in constant_level or index_level do not match regions in gw.regions.

    Returns
    -------
    FPIX : np.array
        An array of length 12*nside_final**2, representing healpixels in a flat skymap.
        Values will be 1. or 0., with 1 representing pixels in ANY ConfidenceRegion specified in gws, **kwargs.
        Can be passed into ax.contour_hpx or ax.contourf_hpx for contour superimposition.
    """
        
    if verbose: 
        print(f"Loading contour union for {len(gws)} GWEvents")
        print(f"NSide resolution: {nside_final}")
        if constant_level!=None: print(f"\tlevel: {constant_level}")
    
    if constant_level!=None:
        conf_map = {gw.name:constant_level for gw in gws}
        
    elif index_level!=None:
        conf_map = {gws[i].name:index_level[i] for i in range(len(gws))}
        
    else:
        raise ValueError
    
    fnside = nside_final
    fnpix = hp.nside2npix(fnside)
    FPIX = np.array([0.,]*fnpix)
    FIPIX = np.array([])

    for n,gw in enumerate(gws):
        c = conf_map[gw.name]
        region = None
        
        for r in gw.regions:
            if r.level==c: region=r
        
        if region==None: raise ValueError
        
        #######################################################################
            
        table = QTable.read(region.path)
        
        LEVEL, IPIX = ah.uniq_to_level_ipix(table["UNIQ"])
        NSIDE = ah.level_to_nside(LEVEL)
        
        PHI, THETA = ah.healpix_to_lonlat(IPIX, NSIDE, order='nested')
        THETA = (units.rad*np.pi/2) - THETA
        
        FIPIX = hp.ang2pix(fnside, THETA.value, PHI.value)
        if sum_type=="binary": FPIX[FIPIX] = 1.
        
        elif sum_type=="freq":
            UFIPIX = np.unique(FIPIX)
            FPIX[UFIPIX] = FPIX[UFIPIX] + 1.
            
        elif sum_type=="prob":
            PROB = table["PROBDENSITY"] * ah.nside_to_pixel_area(NSIDE)
            COUNTS = np.array([0,]*np.size(FPIX))
            for i,ipix in enumerate(FIPIX): 
                FPIX[ipix] += PROB[i].value
                COUNTS[ipix] += 1
            
            FPIX[COUNTS!=0] /= COUNTS[COUNTS!=0]
    
    if sum_type=="freq": FPIX *= 1/np.max(FPIX)
    if sum_type=="prob": FPIX *= len(gws)/np.sum(FPIX)
        
    return FPIX

#######################################################################

def contour_intersection(gws,
                         nside_final=2**5,
                         constant_level=None,
                         index_level=None,
                         verbose=True
                         ):
    """
    Parameters
    ----------
    gws : iter(GW) or iter(ConfidenceRegion)
        Iterable containing GW or ConfidenceRegion objects.
        If GW objects, then constant_level or index_level kwarg is required.
    nside_final : int, optional
        Pixel resolution of output skymap array, must be 2**n (n:int < 52) The default is 2**5.
    constant_level : float, optional
        Confidence level in (0., 1.). If specified, will look for a ConfidenceRegion of level confidence_level in gw.regions to use. 
        The default is None.
    index_level : dict, optional
        Maps gw.name for gw in gws to a specific confidence level in gw.regions. 
        The default is None.
    verbose : bool, optional
        Prints some information while running. The default is True.

    Raises
    ------
    ValueError
        levels given in constant_level or index_level do not match regions in gw.regions.

    Returns
    -------
    FPIX : np.array
        An array of length 12*nside_final**2, representing healpixels in a flat skymap.
        Values will be 1. or 0., with 1 representing pixels in ALL ConfidenceRegion specified in gws, **kwargs.
        Can be passed into ax.contour_hpx or ax.contourf_hpx for contour superimposition.
    """
    
    if verbose: 
        print(f"Loading contour intersection for {len(gws)} GWEvents")
        print(f"NSide resolution: {nside_final}")
        if constant_level!=None: print(f"\tlevel: {constant_level}")
    
    if constant_level!=None:
        conf_map = {gw.name:constant_level for gw in gws}
        
    elif index_level!=None:
        conf_map = {gws[i].name:index_level[i] for i in range(len(gws))}
        
    else:
        raise ValueError
    
    fnside = nside_final
    fnpix = hp.nside2npix(fnside)
    FPIX = np.array([1.,]*fnpix)
    FIPIX = np.array([])

    for n,gw in enumerate(gws):
        print(gw.name)
        c = conf_map[gw.name]
        region = None
        
        for r in gw.regions:
            if r.level==c: region=r
        
        if region==None: raise ValueError
        
        #######################################################################
            
        table = QTable.read(region.path)
        
        LEVEL, IPIX = ah.uniq_to_level_ipix(table["UNIQ"])
        NSIDE = ah.level_to_nside(LEVEL)
        
        PHI, THETA = ah.healpix_to_lonlat(IPIX, NSIDE, order='nested')
        THETA = (units.rad*np.pi/2) - THETA
        
        FIPIX = hp.ang2pix(fnside, THETA.value, PHI.value)
        NFIPIX = np.setdiff1d(np.arange(0,fnpix), FIPIX)
        FPIX[NFIPIX] = 0.
        
    return FPIX

#######################################################################