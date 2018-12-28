# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 09:54:46 2018

@author: renxi
"""

def BisecBSIV(PutCall,S,K,rf,q,T,a,b,MktPrice,Tol,MaxIter):
        # =============================================================================
        # The earlier code invokes theMatlab function BisecBSIV.m that uses the bisection
        # algorithm to find the implied volatilities once the prices are generated. The first part
        # of the function defines a function handle for calculating the price of calls under the
        # Black-Scholes model. The main part of the function calculates the Black-Scholes call
        # price at the midpoint of the two endpoints, and compares the difference between this
        # call price and themarket price. If the difference is greater than the user-specified tolerance,
        # the function updates one of the endpoints with the midpoint and passes through
        # the iteration again. To conserve space, parts of the function have been omitted.
        # =============================================================================
        # Bisection algorithm
        # PutCall = "P"ut or "C"all
        # S = Spot price
        # K = Strike price
        # rf = Risk free rate
        # T = Maturity
        # a = Starting point lower value of vol
        # b = Starting point upper value of vol
        # MktPrice = Market price
        # Tol = Tolerance
        # MaxIter = maximum number of iterations
    
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

