def  EulerMilsteinPrice(scheme,negvar,params,PutCall,S0,K,Mat,r,q,T,N,alpha):

    # Heston Call or Put price using Euler, Milstein, 
    # and Implicit Milstein discretization of the Heston model.
    # INPUTS
    #   scheme = 'E' Euler, 'M' Milstein, or 'IM' Implicit Milstein
    #   negvar = 'R' Reflection or 'T' Truncation of negative variances
    #   params = Heston parameters
    #   PutCall = 'C' for Call option, 'P' for put
    #   S0  = Spot Price
    #   K = Strike price
    #   Mat = Maturity
    #   r = risk free rate
    #   q = dividend yield
    #   T   = Number of time steps
    #   N   = Number of stock price paths
    #   alpha = weight for implicit-explicit scheme
    # OUTPUTS
    #   S = Vector of simulated stock prices
    #   V = Vector of simulated variances
    #   F = Matrix of overriden negative variances
    #   SimPrice = option price by simulation
    import numpy as np
    from math import exp, log, sqrt
    
    # Obtain the simulated stock price and simulated variance
    [S, V, F] = EulerMilsteinSim(scheme,negvar,params,S0,Mat,r,q,T,N,alpha);
    
    # Terminal stock prices
    ST = S[-1,:];
    
    # Payoff vectors
    if 'C' in PutCall:
        Payoff = (ST - K)*(ST - K > 0);
    else: 
        Payoff = (K - ST)*(K - ST > 0);
        
	
    # Simulated price
    SimPrice = exp(-r*Mat)*np.mean(Payoff);
    
    return [S, V, F, SimPrice]