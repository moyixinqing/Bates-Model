def JumpCF(phi,param,T):
    from cmath import exp
    import numpy as np
    i=1j
    lambdJ  = param[5];        # Annual jump frequency 
    muJ      = param[6];        # Random percentage jump
    sigmaJ   = param[7];        # Jump volatility
    jcf = np.exp(-lambdJ*muJ*i*phi*T + lambdJ*T*((1+muJ)**(i*phi)*np.exp(0.5*sigmaJ**2*i*phi*(i*phi-1))-1));
    return jcf

# Heston (1993) second characteristic function (f2) using the "Little Trap" formulation
def HestonCF(phi,param,T,S,r,q):
    from math import pi
    from cmath import exp,sqrt ,log
    import math
    import cmath
    import numpy as np

    i=1j
    kappa = param[0];
    theta = param[1];
    sigma = param[2];
    v0    = param[3];
    rho   = param[4];
    lambd = 0;
    
    # Log of the stock price.
    x = log(S);
    
    # Required parameters.
    a = kappa*theta;
    u = -0.5;
    b = kappa + lambd;
    
    d = np.sqrt((rho*sigma*i*phi - b)**2 - sigma**2*(2*u*i*phi - phi**2));
    g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);

    # "Little Heston Trap" formulation
    if np.any(g.real == math.inf):
        g = np.inf
    else:
        g
    c = 1/g;
    D = (b - rho*sigma*i*phi - d)/sigma**2*((1-np.exp(-d*T))/(1-c*np.exp(-d*T)));
    G = (1 - c*np.exp(-d*T))/(1-c);
    C = (r-q)*i*phi*T + a/sigma**2*((b - rho*sigma*i*phi - d)*T - 2*np.log(G));
    
    # The characteristic function.
    f2 = np.exp(C + D*v0 + i*phi*x);
    return f2

# Bates characteristic function
def BatesCF(phi,param,T,S,rf,q):
    bcf = HestonCF(phi,param,T,S,rf,q) * JumpCF(phi,param,T);
    return bcf 
 
def BatesObjFunSVC(param,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,ObjFun,a,b,Tol,MaxIter):

    import cmath
    from cmath import exp, log
    from math import sqrt,pi
    import numpy as np
    from scipy.stats import norm
    i =1j    
    [NK,NT] = MktPrice.shape;
    
    int1=np.zeros(len(x))
    int2=np.zeros(len(x))
    f=np.zeros(len(x),dtype=complex)
    fnum=np.zeros(len(x),dtype=complex)
    fden=np.zeros(len(x),dtype=complex)
    ModelPrice = np.zeros((NK,NT));
    error  = np.zeros((NK,NT));
    Vega   = np.zeros((NK,NT));
    
    for t in range(NT):
        for j in range(len(x)):    
            phi = x[j];
            f[j]    = BatesCF(phi  ,param,T[t],S,rf,q);
            fnum[j] = BatesCF(phi-i,param,T[t],S,rf,q);
            fden[j] = BatesCF(   -i,param,T[t],S,rf,q);
        
        for k in range(NK):
            for j in range(len(x)):
                phi = x[j];
                int2[j] = w[j] * (exp(-i*phi*log(K[k]))*f[j]/i/phi).real;
                int1[j] = w[j] * (exp(-i*phi*log(K[k]))*fnum[j]/i/phi/fden[j]).real;
            
            # The probabilities P1 and P2
            P1 = 1/2 + 1/pi*sum(int1);
            P2 = 1/2 + 1/pi*sum(int2);
            # The call price
            CallPrice = S*exp(-q*T[t])*P1 - K[k]*exp(-rf*T[t])*P2;
            # Output the option price (put price by put call parity)
            if 'C' in PutCall:
                ModelPrice[k,t] = CallPrice;
            else:
                ModelPrice[k,t] = CallPrice - S*exp(-q*T[t]) + K[k]*exp(-rf*T[t]);
            
            # Select the objective function
            if ObjFun == 1:
                # MSE
                error[k,t] = (MktPrice[k,t] - ModelPrice[k,t])**2;
            elif ObjFun == 2:
                # RMSE
                error[k,t] = (MktPrice[k,t] - ModelPrice[k,t])**2 / MktPrice[k,t];
            elif ObjFun == 3:
                # IVMSE
                ModelIV = BisecBSIV(PutCall[k,t],S,K[k],rf,q,T[t],a,b,ModelPrice[k,t],Tol,MaxIter);
                error[k,t] = (ModelIV - MktIV[k,t])**2;
            elif ObjFun == 4:
                # IVRMSE Christoffersen, Heston, Jacobs proxy
                d = (log(S/K[k]) + (rf-q+MktIV[k,t]**2/2)*T[t])/MktIV[k,t]/sqrt(T[t]);
                Vega[k,t] = S*norm.pdf(d)*sqrt(T[t]);
                error[k,t] = (ModelPrice[k,t] - MktPrice[k,t])**2 / Vega[k,t]**2;
            
    y = sum(sum(error));
    
    
    return y