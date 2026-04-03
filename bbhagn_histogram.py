import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde as kde

###############################################################################
    
def get_intervals(f, bins, conf_list):
    """
    Parameters
    ----------
    f : fitting function
        Returns numpy arrays/singletons for inputted arrays/values,
        this project uses a gaussian kde fit.
    bins : list
        A list of x-bins for the fitting function 
        (used mostly for graphing and establishing domain bounds).
    conf_list : float[]
        a list of confidence levels in (0,1) to be used to make intervals.

    Returns
    -------
    dict
        Keys are the same values as in conf_lists. 
        Values are dictionaries containing data useful for plotting.

    """
    result = dict()

    ###########################################################################

    def integrater(f, bounds):
        """
        Parameters
        ----------
        f : fitting function
        bounds : x-axis bounds (float:a, float:b)

        Returns
        -------
        float
            A trapezoid-integral of only the POSITIVE regions of f within (a,b).
        """
        
        a,b = bounds
        X = np.linspace( a, b, num=int(1E3) )
        Y = f(X)
        YP = Y
        YP[YP<0] = 0
        return np.trapezoid(YP, X)
    
    ###########################################################################
    
    def binary_est(target_val, interval, err_bound):
        """
        Parameters
        ----------
        target_val : float
            The target value to optimize.
        interval : (float:lo, float:hi)
            The lower and upper bounds of binary search (continuous distribution).
        err_bound : float
            Once f(median(lo,hi)) is within err_bound of target_val, 
            the function returns median(lo,hi).

        Returns
        -------
        float
            The result of binary searching the domain of a continuous function
            for a target value in its codomain.

        """
        lo,hi = interval
        med = (lo+hi)/2
        val = integrater(lambda x: A0*(f(x)-med), bounds)
                        
        if abs(val - target_val) < err_bound:
            return med
         
        elif val<target_val: 
            """
            For probability distributions, this condition checks overestimation
            val is the probability of the interval under med
            if val is less than target_val, then med is too high.
            
            Bottom of range is closer to minimizing error
            Range becomes (lo, med) and val_lo stays the same
            """
            return binary_est(target_val=target_val, 
                              interval=(lo,med), 
                              err_bound=err_bound, 
                              )

        else:
            """
            Top of range is closer to minimizing error
            Range becomes (med, hi) and val_lo becomes val_med
            """
            return binary_est(target_val=target_val, 
                              interval=(med,hi), 
                              err_bound=err_bound, 
                              )
        
    ###########################################################################

    def get_intercepts(f, y0, bounds):
        """
        Parameters
        ----------
        f : fit function
        y0 : float
            y=y0 line of intercept
        bounds : (float:a, float:b)
            Bounds on domain of f.

        Returns
        -------
        zeros : np.array()
            Array of intercepts with f(x) and y=y0.
        """
        X = np.linspace(bounds[0], bounds[1], num=int(1E3))
        Y = f(X)-y0
        
        zeros = []
        w = 1
        for i in range(len(X[::w])-1):
            x0, x1 = X[ i*w ], X[ (i+1)*w ]
            y0, y1 = Y[ i*w ], Y[ (i+1)*w ]
            if y0*y1<0: zeros.append((x0+x1)/2)
        
        return zeros
            
    ###########################################################################
    
    bounds = (bins[0], bins[-1])
    A0 = integrater(f, bounds) **-1 #normalization coefficient
    y_vals = f(bins)
    y0_max, y0_min = np.max( y_vals ), np.min( y_vals ) #bounds on the threshold line
    err_bound = 1E-3
    
    for conf in conf_list:
        y0 = binary_est(target_val=conf, 
                        interval=(y0_min,y0_max),
                        err_bound=err_bound,
                        )
                     
        result[conf] = dict()
        result[conf]["threshold"] = y0
        
        zeros = get_intercepts(f, y0, bounds)
        positive_interval = f(bounds[0])[0] > y0
        int_bools = np.array([positive_interval]*np.size(bins))    

        for zero in zeros:
            int_bools[bins>zero] = ~int_bools[bins>zero]
        
        result[conf]["bools"] = int_bools
    
    return result
        
###############################################################################

def main():
    
    path = "/Users/mck/Desktop/palmese_research/corr_scripts/O4_samples_graham23.dat"
    data = np.loadtxt(path)
    bins = np.linspace(0,0.24,num=24)
    fit_bins = np.linspace(0,0.24,num=1000)
    fit_bins = fit_bins[fit_bins <= bins[-2]]
    
    #reflecting the distribution across x-axis  for cleaner probdens fit
    symm_bins = -fit_bins[::-1]
    symm_bins = np.append(symm_bins, fit_bins)
    reflect_data = -data
    symm_data = np.append(reflect_data, data)
    
    ###########################################################################
    # Fitting functions for linear and logarithmic scale graphs
    
    fit_linax = kde(symm_data,
              bw_method = 5E-2
              )

    fit_logax = kde(symm_data,
              bw_method = 2E-1
              )

    ###########################################################################
    # Generating confidence intervals
    
    conf_intervals = {
                      0.50:{},
                      0.75:{},
                      0.90:{},
                      0.95:{},
                     }
    
    colors = ["#cc"+f"{(153//len(conf_intervals) * n)//16}{(153//len(conf_intervals) * n) %16//9}"*2 for n in range(len(conf_intervals))]
    
    conf_intervals = get_intervals(fit_logax, fit_bins, list(conf_intervals))
    for c,conf in enumerate(conf_intervals):
        conf_intervals[conf]["color"] = colors[c]
    
    
    ###########################################################################
    # Plotting
    
    fig,ax = plt.subplots(2, figsize=(5, 5), dpi=100)

    subplots_config = {
                        'linear':   {
                                    'fit':fit_linax,
                                    'log':False,
                                    'ax':0,
                                    },
                
                       'logarithmic':   {
                                         'fit':fit_logax,
                                         'log':True,
                                         'ax':1,
                                         },
                       }

    for fit in subplots_config:
        params = subplots_config[fit]
        
        ax[params["ax"]].hist(data,
                 bins=bins,
                 log=params["log"],
                 alpha=0.5,
                 density=True
                 )
        
        fit = ax[params["ax"]].plot(fit_bins,
                 params["fit"](fit_bins),
                 label = "kde fit line",
                 )
        
        filled_intervals = []
        for level in conf_intervals:
            interval_data = conf_intervals[level]
            interval = ax[params["ax"]].fill_between(fit_bins,
                                                     params["fit"](fit_bins),
                                                     where=interval_data["bools"],
                                                     color=interval_data["color"],
                                                     alpha=0.5,
                                                     label=f"{round(level*100)}% confidence region",
                                                     )
            
            if params["ax"]==0: filled_intervals += [interval]
            
        for level in conf_intervals:
            interval_data = conf_intervals[level]
            ax[params["ax"]].plot(fit_bins,
                                  [interval_data["threshold"]]*len(fit_bins),
                                  color=interval_data["color"],
                                  alpha=0.5,
                                  dashes=(level*10, (1-level)*10),
                                  )
            
        if params["ax"]==0: 
            ax[params["ax"]].legend(loc=0,)
    
    plt.xlabel("$\lambda$")
    fig.suptitle("MCMC Posterior Distribution")
    
    plt.show()

###############################################################################
    
main()