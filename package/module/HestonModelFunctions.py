# -*- coding: utf-8 -*-
"""
Created on Sat Dec 15 15:03:44 2018

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

# In[HestonProb]    
def HestonProb(phi,kappa,theta,lambd,rho,sigma,tau,K,S,r,q,v,Pnum,Trap):
    # Returns the integrand for the risk neutral probabilities P1 and P2.
    # phi = integration variable
    # Pnum = 1 or 2 (for the probabilities)
    # Heston parameters:
    #    kappa  = volatility mean reversion speed parameter
    #    theta  = volatility mean reversion level parameter
    #    lambda = risk parameter
    #    rho    = correlation between two Brownian motions
    #    sigma  = volatility of variance
    #    v      = initial variance
    # Option features.
    #    PutCall = 'C'all or 'P'ut
    #    K = strike price
    #    S = spot price
    #    r = risk free rate
    #    q = dividend yield
    #    Trap = 1 "Little Trap" formulation 
    #           0  Original Heston formulation
    
    from math import pi
    from cmath import exp, sqrt,log
    import cmath
    import numpy as np
    
    i=1j
    # Log of the stock price.
    x = log(S);
    
    # Parameter "a" is the same for P1 and P2.
    a = kappa*theta;
    
    # Parameters "u" and "b" are different for P1 and P2.
    if Pnum==1:
    	u = 0.5;
    	b = kappa + lambd - rho*sigma;
    else:
    	u = -0.5;
    	b = kappa + lambd;
    d = sqrt((rho*sigma*i*phi - b)**2 - (sigma**2)*(2*u*i*phi - phi**2));
    g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
    if Trap==1:
    	# "Little Heston Trap" formulation
    	c = 1/g;
    	D = (b - rho*sigma*i*phi - d)/(sigma**2)*((1-np.exp(-d*tau))/(1-c*np.exp(-d*tau)));
    	G = (1 - c*np.exp(-d*tau))/(1-c);
    	C = (r-q)*i*phi*tau + a/(sigma**2)*((b - rho*sigma*i*phi - d)*tau - 2*np.log(G));
    elif Trap==0:
    	# Original Heston formulation.
    	G = (1 - g*np.exp(d*tau))/(1-g);
    	C = (r-q)*i*phi*tau + a/sigma**2*((b - rho*sigma*i*phi + d)*tau - 2*np.log(G));
    	D = (b - rho*sigma*i*phi + d)/sigma**2*((1-np.exp(d*tau))/(1-g*np.exp(d*tau)));
    # The characteristic function.
    f =np.exp(C + D*v + i*phi*x);
    # Return the real part of the integrand.
    y = (np.exp(-i*phi*np.log(K))*f/i/phi)
    y=y.real
    return y

# In[HestonPriceGaussLaguerre]    
def HestonPriceGaussLaguerre(PutCall,S,K,T,r,q,kappa,theta,sigma,lambd,V0,rho,trap,x,w):
    

    # Heston (1993) call or put price by Gauss-Laguerre Quadrature
    # Uses the original Heston formulation of the characteristic function,
    # or the "Little Heston Trap" formulation of Albrecher et al.
    # INPUTS -------------------------------------------------------
    #   PutCall = 'C' Call or 'P' Put
    #   S = Spot price.
    #   K = Strike
    #   T = Time to maturity.
    #   r = Risk free rate.
    #   kappa  = Heston parameter: mean reversion speed.
    #   theta  = Heston parameter: mean reversion level.
    #   sigma  = Heston parameter: volatility of vol
    #   lambda = Heston parameter: risk.
    #   v0     = Heston parameter: initial variance.
    #   rho    = Heston parameter: correlation
    #   trap:  1 = "Little Trap" formulation
    #          0 = Original Heston formulation
    #   x = Gauss Laguerre abscissas
    #   w = Gauss Laguerre weights
    # OUTPUT -------------------------------------------------------
    #   The Heston call or put price
# =============================================================================
# T=tau
# S=S0
# trap=1
#  
# =============================================================================
# =============================================================================
#     kappa  = param[0]; 
#     theta  = param[1]; 
#     sigma  = param[2]; 
#     v0     = param[3]; 
#     rho    = param[4];
#     lambd = 0;
# =============================================================================

    from math import exp,pi
    import numpy as np    
    # Numerical integration
    int1=np.zeros(len(x))
    int2=np.zeros(len(x))
    for k in range(len(x)):
    	int1[k] = w[k]*HestonProb(x[k],kappa,theta,lambd,rho,sigma,T,K,S,r,q,v0,1,trap);
    	int2[k] = w[k]*HestonProb(x[k],kappa,theta,lambd,rho,sigma,T,K,S,r,q,v0,2,trap);
    # Define P1 and P2
    P1 = 1/2 + 1/pi*sum(int1);
    P2 = 1/2 + 1/pi*sum(int2);
    
    # The call price
    HestonC = S*exp(-q*T)*P1 - K*exp(-r*T)*P2;
    
    # The put price by put-call parity
    HestonP = HestonC - S*exp(-q*T) + K*exp(-r*T);
    
    # Output the option price
    if 'C' in PutCall:
    	y = HestonC;
    else:
    	y = HestonP;
    return y

# In[HestonCF]    
def HestonCF(phi,kappa,theta,lambd,rho,sigma,tau,S,r,q,v0,Trap):
    # Returns the Heston characteristic function,f2 (the second CF)
    # Uses original Heston formulation 1,or formulation 2 from "The Little Heston Trap"
    # phi = integration variable
    # Heston parameters:
    #    kappa  = volatility mean reversion speed parameter
    #    theta  = volatility mean reversion level parameter
    #    lambd = risk parameter
    #    rho    = correlation between two Brownian motions
    #    sigma  = volatility of variance
    #    v0     = initial variance
    # Option features.
    #    tau = maturity
    #    S = spot price
    #    r = risk free rate
    #    q = dividend yield
    # Trap = 0 --> Original Heston Formulation
    # Trap = 1 --> Little Trap formulation
    
    from math import pi
    from cmath import exp, sqrt,log
    import cmath
    import numpy as np
    #phi = 0.0444894
    #tau =1.5
    #v=0.05
    
    i=1j
    # Log of the stock price.
    x = log(S);
    
    # Parameter a
    a = kappa*theta;
    
    # Parameters u,b,d,and b
    u = -0.5;
    b = kappa + lambd;
    d = np.sqrt((rho*sigma*i*phi - b)**2 - (sigma**2)*(2*u*i*phi - phi**2));
    g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
    if Trap==1:
    	# "Little Heston Trap" formulation
    	c = 1/g;
    	G = (1 - c*np.exp(-d*tau))/(1-c);
    	C = (r-q)*i*phi*tau + a/sigma**2*((b - rho*sigma*i*phi - d)*tau - 2*np.log(G));
    	D = (b - rho*sigma*i*phi - d)/sigma**2*((1-np.exp(-d*tau))/(1-c*np.exp(-d*tau)));
    elif Trap==0:
    	# Original Heston formulation.
    	G = (1 - g*np.exp(d*tau))/(1-g);
    	C = (r-q)*i*phi*tau + a/sigma**2*((b - rho*sigma*i*phi + d)*tau - 2*np.log(G));
    	D = (b - rho*sigma*i*phi + d)/sigma**2*((1-np.exp(d*tau))/(1-g*np.exp(d*tau)));
    
    # The characteristic function.
    y = np.exp(C + D*v0 + i*phi*x);
    return y

# In[BisecBSIV]
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
