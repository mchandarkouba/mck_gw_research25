import numpy as np
import math

def binary_est(target_val, interval, err_bound, prev_val=None):  
    lo,hi = interval
    med = (lo+hi)/2
    
    val_lo = prev_val or integrater(lambda x: A0*(f(x)-lo), bounds)
    val_med = integrater(lambda x: A0*(f(x)-med), bounds)
    
    if abs(val_med - target_val) < err_bound:
        return med
     
    elif abs(val_lo - target_val) <= abs(val_med - target_val): 
        """
        Bottom of range is closer to minimizing error
        Range becomes (lo, med) and val_lo stays the same
        """
        return binary_est(target_val=target_val, 
                          interval=(lo,med), 
                          err_bound=err_bound, 
                          prev_val=val_lo,
                          )
     
    else:
        """
        Top of range is closer to minimizing error
        Range becomes (med, hi) and val_lo becomes val_med
        """
        return binary_est(target_val=target_val, 
                          interval=(med,hi), 
                          err_bound=err_bound, 
                          prev_val=val_med,
                          )
    
def integrater(f, bounds):
    a,b = bounds
    X = np.linspace( a, b, num=int(1E3) )
    Y = f(X)
    YP = Y
    YP[YP<0] = 0
    return np.trapezoid(YP, X)

###########################################################################

def binaryAdjust(lo=-math.pi/2, 
                 hi=math.pi/2, 
                 prevErrLo=None,
                 prevErrHi=None,
                 n=0, N=100
                 ):
    
    med = (lo+hi)/2
    q1 = (lo+med)/2
    q3 = (med+hi)/2
    errLo = prevErrLo or error(q1)
    errHi = prevErrHi or error(q3)
    errMed = error(med)
    
    if n>=N:
        return med
     
    elif errLo <= errHi: 
        """
        Bottom of range is closer to minimizing error
        Range becomes (lo, med) and lo stays the same
        """
        return binaryAdjust(lo, med, 
                            prevErrLo=errLo, 
                            prevErrHi=errMed,
                            n=n+1, N=N
                            )
     
    else:
        """
        Top of range is closer to minimizing error
        Range becomes (med, hi) and lo becomes med
        """
        return binaryAdjust(med, hi,
                            prevErrLo=errMed,
                            prevErrHi=errHi,
                            n=n+1, N=N
                            ) 
    
def f(x):
    pass

def g(a,x):
    pass

def error(a):
    pass
    

 ###########################################################################

f = lambda x: (x-1)**3 -2*(x-1)**2 -(x-1) +2
bounds = (0,5)
A0 = integrater(f, bounds) **-1 #normalization coefficient
y0_max, y0_min = (0,2.12) #starting tolerance  appx. the abs. max of the distribution
err_bound = 1E-2
conf = 0.95
 
y0 = binary_est(target_val=conf,
                interval=(y0_min,y0_max),
                err_bound=err_bound
                )

print(y0)