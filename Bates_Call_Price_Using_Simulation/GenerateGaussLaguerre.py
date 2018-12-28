# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 23:44:02 2018

@author: renxi
"""

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
