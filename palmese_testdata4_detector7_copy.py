#Plot showing the luminosity distance vs confidence area relationship for simulated BBH events, 
#color-codes based off of categories defined using search parameters (below) 
#includes a searching method for "detectors of interest" given various conditions:
    #any: at least detector within "detectors of interest" must be used
    #all: all detectors within "detectors of interest" must be used
    #none: no detectors within "detectors of interest" must be used
    #unique: exactly one detector from "detectors of interest" must be used
#given the number of detectors used to observe an event is within the "number of interest" set
#includes histograms on both axes to show distribution of each population 
#per luminosity distance and 70% confidence sky area.

#reformatted to make adding code easier,
#plus points in overlapping populations will now be colored bright red

#Uses locally accessed gw simulation data: combined_bbh_O5.csv, 
#and NERSC-retrieved detector data: detector_data6.npy

############################################################################### Imports, data, scatterplot variables

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as matpatch
import math

eventPath = r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/combined_bbh_O5.csv"
detectorPath = r"/Users/mck/Desktop/palmese_research/detector_analysis_scripts/detector_data6.npy"

############################################################################### Search criteria

def colorPops(param, detectorMatch, colors='colors', opaques='opaques'):
    detectorMatchBools = []
    
    for m in detectorMatch: #detectorMatch is strings of detectors used for each event from data
        ofInterestSet = [True]*len(param)
        m = set(m.split(','))

        for i, population in enumerate(param): #checking membership of all populations for each item in detectorMatch, then storing as a nested list of bools
            detectorsInterest, conditionInterest, numberInterest, color, opaque = population
            ofInterest = ofInterestSet[i]

            n = detectorsInterest
     
            if conditionInterest=='any': ofInterest = n&m!=set()
            if conditionInterest=='all': ofInterest = n|m==m
            if conditionInterest=='none': ofInterest = n&m==set()
            if conditionInterest=='uniques': ofInterest = len(n&m)==1    
        
            ofInterest &= len(m) in numberInterest
            ofInterestSet[i] = ofInterest
                        
        detectorMatchBools.append(ofInterestSet)
    
    for i, L in enumerate(detectorMatchBools):
        n = [j for j, b in enumerate(L) if b] #returns the indices of True in list
        
        if len(n)==0:
            pass #no recoloring necessary, is not part of population so we use the initialized color and opacity
        
        elif len(n)==1:
            colors[i] = param[n[0]][3]
            opaques[i] = param[n[0]][4]
            
        elif len(n) not in (0,1): #checking for points in overlapping populations
            colors[i] = 'red'
            opaques[i] = 1
            

    return {'colors':colors, 'opaques':opaques}
    
def scatterHists(x, y, fig, axs, param, mapping):
    ax, xHist, yHist = axs['scatter'], axs['histx'], axs['histy']    
    colors, opaques = mapping['colors'], mapping['opaques']
    
    xHist.tick_params(axis='x')
    yHist.tick_params(axis='y')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.scatter(x,y, c=colors, alpha=opaques)
    
    binwidth = 0.2
    xmin, xmax = math.log10(np.min(x)), math.log10(np.max(x))
    ymin, ymax = math.log10(np.min(y)), math.log10(np.max(y))

    xbins = np.arange(xmin,xmax, binwidth)
    ybins = np.arange(ymin,ymax, binwidth)

    colorsArray = np.array(colors)
    for pop in param:
        c = pop[3]
        xHist.hist(np.log10(x[colorsArray[:]==c]), bins=xbins, histtype='step', color=c, linewidth=1.5)
        yHist.hist(np.log10(y[colorsArray[:]==c]), bins=ybins, histtype='step', color=c, orientation='horizontal', linewidth=1.5)

    patches = [ matpatch.Patch(color=param[n][3], label=f'{param[n][2]} detectors including {param[n][1]} of {param[n][0]}') for n in range(len(param)) ]
    fig.legend(handles=patches, fontsize=12, loc='outside lower left') # fig.legend(loc='outside lower left') or ax.legend(loc='lower right')

############################################################################### population definition and color-coding

def generateMain(eventPath, detectorPath):
    events = np.loadtxt(eventPath, skiprows=1)
    detectorData = np.loadtxt(detectorPath, dtype=str)

    D = events[:,1] #x axis
    A_70 = events[:,11] #y axis

    simID = detectorData[:, 2].astype(int)
    sharedID = np.intersect1d(simID, events[:,0], return_indices=True)

    #colorMatch = detectorData[sharedID[1], 1].astype(int)
    detectorMatch = detectorData[sharedID[1], 0].astype(str)

    ################################################################################ Initial state and Parameters

    mapping = {'colors': ['gray' for _ in events], #default color is gray for non-populations
              'opaques': [1 for _ in events] } #default opacity is 1 for non-populations
    
    param = [ # set(*detectors), condition, set(*number), dotColor, opacity
                [set(['H1','L1','K1','V1']), 'any', set([1]), 'gold', 0.5],
                [set(['H1','L1']), 'uniques', set([2,3]), 'green', 0.5],
                [set(['H1','L1']), 'all', set([2,3,4]), 'purple', 0.5],
                [set(['K1','V1']), 'all', set([2]), 'blue', 0.5]  
            ]
    
    ################################################################################ Scatterplot and histograms 
    
    fig, axs = plt.subplot_mosaic([['histx', '.'],
                                   ['scatter', 'histy']],
                                   figsize=(10,10),
                                   width_ratios=(4,1), height_ratios=(1,4),
                                   layout='constrained')
    
    mapping = colorPops(param, detectorMatch, **mapping) # *args, **kwargs

    args = (D, A_70, fig, axs, param, mapping)
    scatterHists(*args)
    
    fig.suptitle(f"""Area (70% CL) vs dL by Population""", fontsize=24, fontweight='semibold')
    axs['scatter'].set_xlabel('dL [Mpc]', fontsize=20)
    axs['scatter'].set_ylabel('Area (70% CL) [deg^2]', fontsize=20)
    
    
    plt.show()
    
generateMain(eventPath, detectorPath)