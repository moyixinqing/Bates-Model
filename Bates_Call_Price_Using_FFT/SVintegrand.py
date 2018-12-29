# =============================================================================
# # Logarithmic normal jump characteristic function
# def JumpCF(phi,param,T):
#     from cmath import exp
#     i=1j
#     lambdJ  = param[5];        # Annual jump frequency 
#     muJ      = param[6];        # Random percentage jump
#     sigmaJ   = param[7];        # Jump volatility
#     jcf = exp(-lambdJ*muJ*i*phi*T + lambdJ*T*((1+muJ)**(i*phi)*exp(0.5*sigmaJ**2*i*phi*(i*phi-1))-1));
#     return jcf
# 
# # Heston (1993) second characteristic function (f2) using the "Little Trap" formulation
# def HestonCF(phi,param,T,S,r,q):
#     from math import pi
#     from cmath import exp,sqrt ,log
#     import math
#     import cmath
#     import numpy as np
# 
#     i=1j
#     kappa = param[0];
#     theta = param[1];
#     sigma = param[2];
#     v0    = param[3];
#     rho   = param[4];
#     lambd = 0;
#     
#     # Log of the stock price.
#     x = log(S);
#     
#     # Required parameters.
#     a = kappa*theta;
#     u = -0.5;
#     b = kappa + lambd;
#     
#     d = np.sqrt((rho*sigma*i*phi - b)**2 - sigma**2*(2*u*i*phi - phi**2));
#     g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);
# 
#     # "Little Heston Trap" formulation
#     if g.real == np.inf:
#         g = np.inf
#     else:
#         g
#     c = 1/g;
#     D = (b - rho*sigma*i*phi - d)/sigma**2*((1-np.exp(-d*T))/(1-c*np.exp(-d*T)));
#     G = (1 - c*np.exp(-d*T))/(1-c);
#     C = (r-q)*i*phi*T + a/sigma**2*((b - rho*sigma*i*phi - d)*T - 2*np.log(G));
#     
#     # The characteristic function.
#     f2 = np.exp(C + D*v0 + i*phi*x);
#     return f2
# 
# =============================================================================
# Stochastic volatility integrand
def SVintegrand(phi,Model,param,S,K,r,q,T,Pnum):
        
    from math import pi
    from cmath import exp, sqrt,log
    import cmath
    import numpy as np
    i=1j
    
    if Pnum==2:
        if Model =='Heston':
            f = HestonCF(phi,param,T,S,r,q);
        elif Model =='Bates':
            f = HestonCF(phi,param,T,S,r,q) * JumpCF(phi,param,T);
        
        svi = (exp(-i*phi*log(K))*f/i/phi).real;
    elif Pnum==1:
        if Model =='Heston':
            fnum = HestonCF(phi-i,param,T,S,r,q);
            fden = HestonCF(   -i,param,T,S,r,q);
        elif Model =='Bates':
            fnum = HestonCF(phi-i,param,T,S,r,q) * JumpCF(phi-i,param,T);
            fden = HestonCF(   -i,param,T,S,r,q) * JumpCF(   -i,param,T);
        
        svi = (exp(-i*phi*log(K))*fnum/i/phi/fden).real;
  
    return svi
    