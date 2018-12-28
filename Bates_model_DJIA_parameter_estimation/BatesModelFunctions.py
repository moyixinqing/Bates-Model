# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 10:07:50 2018

@author: renxi
"""
# In[GenerateGaussLaguerre]
def GenerateGaussLaguerre(n):
    # =============================================================================
    #     # Generate abscissas (x) and weights (w) for Gauss Laguerre integration
    #     # The Laguerre polynomial
    #     L=np.zeros(n)
    #     for k in range(n):
    #     	L[k]=(((-1)**k)/math.factorial(k)*nchoosek(n,k))
    #     
    #     L.ndim
    #     # Need to flip the vector to get the roots in the correct order
    #     L = np.flip(L);
    # 
    #     # Find the roots.  
    #     x = np.flipud(np.roots(L));
    #     
    # =============================================================================
    import math
    from math import exp
    import numpy as np
    
    def nchoosek(n,r):
        f = math.factorial
        return f(n) / f(r) / f(n-r)
        
    [x,w]=np.polynomial.laguerre.laggauss(n) 
    
    # Find the weights
    w = np.zeros(n)
    dL = np.zeros((n,len(x)));    
    for j in range(len(x)):
    	# The coefficients of the derivative of the Laguerre polynomial
    	for k in range(n):
    		dL[k,j] = (-1)**(k+1)/math.factorial(k)*nchoosek(n,k+1)*x[j]**(k);
    	# The weight w[j]
    	w[j] = 1/x[j]/sum(dL[:,j])**2;
    	# The total weight w[j]exp(x(j))
    	w[j] = w[j]*exp(x[j]);
  
    return [x, w]

# In[]
def BisecBSIV(PutCall,S,K,rf,q,T,a,b,MktPrice,Tol,MaxIter):

        from scipy.stats import norm
        from math import exp,sqrt,log
        # Function for the Black Scholes call and put
        def BSC(s,K,rf,q,v,T):
            BSC=(s*exp(-q*T)*norm.cdf((log(s/K) + (rf-q+v**2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*norm.cdf((log(s/K) + (rf-q+v**2/2)*T)/v/sqrt(T) - v*sqrt(T)));
            return BSC
        
        def BSP(s,K,rf,q,v,T): 
            BSP= (K*exp(-rf*T)*norm.cdf(-(log(s/K) + (rf-q+v**2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*norm.cdf(-(log(s/K) + (rf-q+v**2/2)*T)/v/sqrt(T)));
            return BSP
        
        if 'C' in PutCall:
            lowCdif  = MktPrice - BSC(S,K,rf,q,a,T);
            highCdif = MktPrice - BSC(S,K,rf,q,b,T);
        else:
            lowCdif  = MktPrice - BSP(S,K,rf,q,a,T);
            highCdif = MktPrice - BSP(S,K,rf,q,b,T);
        
        if lowCdif*highCdif > 0:
            y = -1;
        else:
            for x in range(MaxIter):
                midP = (a+b)/2;
        		
                if 'C' in PutCall:
                    midCdif = MktPrice - BSC(S,K,rf,q,midP,T);
                else:
                    midCdif = MktPrice - BSP(S,K,rf,q,midP,T);
        		
                if abs(midCdif)<Tol:
                    break
                else:
                    if midCdif>0:
                        a = midP;
                    else:
                        b = midP;
            
            y = midP;
        
        return y


# In[time]
import time

def TicTocGenerator():
    # Generator that returns time differences
    ti = 0           # initial time
    tf = time.time() # final time
    while True:
        ti = tf
        tf = time.time()
        yield tf-ti # returns the time difference

TicToc = TicTocGenerator() # create an instance of the TicTocGen generator

# This will be the main function through which we define both tic() and toc()
def toc(tempBool=True):
    tempTimeInterval = next(TicToc)
    return tempTimeInterval

def tic():
    # Records a time in TicToc, marks the beginning of a time interval
    toc(False)