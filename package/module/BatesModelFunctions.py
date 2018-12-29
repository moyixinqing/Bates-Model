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
# In[BatesPrice]
def BatesPrice(PutCall,param,T,K,S,r,q,x,w):
    from math import exp,pi
    import numpy as np
    
    int1 = np.zeros(len(x));
    int2 = np.zeros(len(x));

    # Build the integrands for P1 and P2;
    for k in range(len(x)):
        int1[k] = w[k] * BatesIntegrand(x[k],param,S,K,r,q,T,1);
        int2[k] = w[k] * BatesIntegrand(x[k],param,S,K,r,q,T,2);

    # The integrals
    I1 = sum(int1);
    I2 = sum(int2);
    
    # The probabilities P1 and P2
    P1 = 1/2 + 1/pi*I1;
    P2 = 1/2 + 1/pi*I2;
    
    # The call price
    CallPrice = S*exp(-q*T)*P1 - K*exp(-r*T)*P2;
    
    # The put price by put-call parity
    PutPrice = CallPrice - S*exp(-q*T) + K*exp(-r*T);
    
    # Output the option price
    if 'C' in PutCall:
    	y = CallPrice;
    else:
    	y = PutPrice;

    return y
# In[BatesCF]
# Logarithmic normal jump characteristic function
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
    if g.real == math.inf:
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
# In[BatesIntegrand]
# Bates integrand
def BatesIntegrand(phi,param,S,K,r,q,T,Pnum):
    from math import pi
    from cmath import exp, sqrt,log
    import cmath
    import numpy as np
    i=1j
    
    if Pnum==2:
        f = HestonCF(phi,param,T,S,r,q) * JumpCF(phi,param,T);
        svi = (exp(-i*phi*log(K))*f/i/phi).real;
    elif Pnum==1:
        fnum = HestonCF(phi-i,param,T,S,r,q) * JumpCF(phi-i,param,T);
        fden = HestonCF(   -i,param,T,S,r,q) * JumpCF(   -i,param,T);
        svi = (exp(-i*phi*log(K))*fnum/i/phi/fden).real;
    return svi
# In[BatesCF]
# Bates characteristic function
def BatesCF(phi,param,T,S,rf,q):
    bcf = HestonCF(phi,param,T,S,rf,q) * JumpCF(phi,param,T);
    return bcf 
# In[SVintegrand]
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
# In[SVprice]    
def SVprice(Model,PutCall,param,T,K,S,r,q,x,w):
    from math import exp,pi
    import numpy as np
    
    int1 = np.zeros(len(x));
    int2 = np.zeros(len(x));
    # Build the integrands for P1 and P2;
    for k in range(len(x)):
        int1[k] = w[k] * SVintegrand(x[k],Model,param,S,K,r,q,T,1);
        int2[k] = w[k] * SVintegrand(x[k],Model,param,S,K,r,q,T,2);
    
    # The integrals
    I1 = sum(int1);
    I2 = sum(int2);
    
    # The probabilities P1 and P2
    P1 = 1/2 + 1/pi*I1;
    P2 = 1/2 + 1/pi*I2;
    
    # The call price
    CallPrice = S*exp(-q*T)*P1 - K*exp(-r*T)*P2;
    
    # The put price by put-call parity
    PutPrice = CallPrice - S*exp(-q*T) + K*exp(-r*T);
    
    # Output the option price
    if 'C' in PutCall:
    	y = CallPrice;
    else:
    	y = PutPrice;

    return y 
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
# In[BatesCallFFT]
def BatesCallFFT(N,uplimit,S0,r,q,tau,param,alpha,fast,rule):
    
    # Fast Fourier Transform of the Bates Model
    # INPUTS
    #   N  = number of discretization points
    #   uplimit = Upper limit of integration
    #   S0 = spot price
    #   r = risk free rate
    #   q = divid yield
    #   tau = maturity
    #   sigma = volatility
    #   alpha = dampening factor
    #   fast = fast versus slow algorithm.
    #     fast = 1 fast version.  Uses vectorization.
    #     fast = 0 slow version.  Uses loops
    #   rule = integration rule
    #     rule = 1 --> Trapezoidal
    #     rule = 2 --> Simpson's Rule
    # --------------------------------------
    # Outputs
    #   CallFFT = Black Scholes call prices using FFT
    #   CallBS  = Black Scholes call prices using closed form
    #         K = Strike prices
    #       eta = increment for integration range
    # lambdainc = increment for log-strike range
    # --------------------------------------
    from math import log, pi
    import numpy as np
    from cmath import exp
    i =1j
    # log spot price
    s0 = log(S0);
    
    # Specify the increments
    eta = uplimit/N;
    lambdainc = 2*pi/N/eta;
    
    # Initialize and specify the weights
    w = np.ones((N));
    if 'T' in rule:            # Trapezoidal rule
        w[0] = 1/2;
        w[N-1] = 1/2;
    elif 'S' in rule:         # Simpson's rule
        w[0] = 1/3;
        w[N-1] = 1/3;
        for k in range(1,N-1):
            if k%2==0:
                w[k] = 2/3;
            else:
                w[k] = 4/3;
    
    # Specify the b parameter
    b = N*lambdainc/2;
    
    # Create the grid for the integration
    v = eta*np.arange(N);
    
    # Create the grid for the log-strikes
    
    k = -b + lambdainc*np.arange(N) + s0;
    k = k.reshape((N,1))
    # Create the strikes and identify ATM 
    K = np.exp(k);
    
    # Initialize the price vector;
    CallFFT = np.zeros((N,1));
    
    if fast==1:
        # Implement the FFT - fast algorithm
        U = np.arange(N).reshape((N, 1));
        J = np.arange(N).reshape((N, 1));
        #U = np.matrix(U)
        #J = np.matrix(J)
        
        psi = BatesFTCF(v-(alpha+1)*i,param,S0,r,q,tau);
        #psi = np.conj(psi);
        phi = exp(-r*tau)*psi / (alpha**2 + alpha - v**2 + i*v*(2*alpha+1));
        x = np.exp(i*(b-s0)*v)*phi*w;
        #e = np.exp(-i*2*pi/N*(U.T*J))*np.matrix(x).T;
        x = x.reshape((N, 1))
        e = np.exp(-i*2*pi/N*(U.T*J))@x;
        CallFFT = eta*np.exp(-alpha*k)/pi *np.array(e.real);
    
    elif fast==0:
        # Implement the FFT - slow algorithm
        for u in range(N):
            for j in range(N):
                psi[j] = BatesFTCF(v[j]-(alpha+1)*i,param,S0,r,q,tau);
                phi[j] = exp(-r*tau)*psi[j]/(alpha^2 + alpha - v[j]^2 + i*v[j]*(2*alpha+1));
                x[j] = exp(i*(b-s0)*v[j])*phi[j]*w[j];
                e[j] = exp(-i*2*pi/N*(j-1)*(u-1))*x[j];
            
            CallFFT[u] = eta*exp(-alpha*k[u])/pi * e.sum().real;
        
    return [CallFFT, K, lambdainc, eta] 

def JumpFTCF(phi,param,T):
    from cmath import exp
    import numpy as np
    i=1j
    lambdJ  = param[5];        # Annual jump frequency 
    muJ      = param[6];        # Random percentage jump
    sigmaJ   = param[7];        # Jump volatility
    jcf = np.exp(-lambdJ*muJ*i*phi*T + lambdJ*T*((1+muJ)**(i*phi)*np.exp(0.5*sigmaJ**2*i*phi*(i*phi-1))-1));
    return jcf

# Heston (1993) second characteristic function (f2) using the "Little Trap" formulation
def HestonFTCF(phi,param,T,S,r,q):
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
def BatesFTCF(phi,param,S,r,q,T):
    bcf = HestonFTCF(phi,param,T,S,r,q) * JumpFTCF(phi,param,T);
    return bcf 

# In[ExtractRND]
def ExtractRND(K,CallPrice):
    # Function to extract the Risk Neutral Density from Call prices.
    # Inputs
    #   K = vector of strikes (size N)
    #   CallPrice = vector of call prices (size N).
    # Outputs
    #   K2 = vector of strikes (size N-4).
    #   RND = vector of RND (size N-4).
    
    # Calculate the first derivatives of calls w.r.t. strike K. 
    # Use central differences
    import numpy as np
    
    dCdK = np.zeros(len(K)-2)
    for i in range(1,len(K)-1):
    	    dK = K[i+1] - K[i-1];
    	    dC = CallPrice[i+1] - CallPrice[i-1];
    	    dCdK[i-1] = dC/dK;   
    # Calculate the risk neutral density by central finite differences.
    RND = np.zeros(len(dCdK)-2)
    for i in range(1,len(dCdK)-1):
    	dK  = K[i+1] - K[i-1];
    	dC2 = dCdK[i+1] - dCdK[i-1];
    	RND[i-1] = dC2/dK;   
    K2 = K[2:-2];
    return [RND, K2]
# In[BatesObjFunSVC] 
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
# In[EulerMilsteinSim]
def EulerMilsteinSim(scheme,negvar,param,S0,Mat,r,q,NT,NS,alpha):
    
    # Euler, Milstein, Implicit Milstein, and Weighted Implicit Milstein
    # discretization of the Heston model.
    # INPUTS
    #   scheme = Scheme for the variance and stock price
    #            'E' Euler, 'M' Milstein, 'IM' Implicit Milstein, or
    #            'WM' Weighted Explicit-Implicit Scheme
    #   negvar = 'R' Reflection or 'T' Truncation of negative variances
    #   params = Heston parameters
    #   S0  = Spot Price
    #   Mat = Maturity
    #   r = riskf ree rate
    #   q = dividend yield
    #   NT   = Number of time steps
    #   NS   = Number of stock price paths
    #   alpha = Weight for explicit-implicit scheme
    # OUTPUTS
    #   S = Vector of simulated stock prices
    #   V = Vector of simulated variances
    #   F = Matrix of overriden negative variances
    import numpy as np
    from math import exp, sqrt
    import random
    
    # Heston parameters
    kappa   = param[0];
    theta   = param[1];
    sigma   = param[2];
    v0      = param[3];
    rho     = param[4];
    lambdaJ = param[5];
    muJ     = param[6];
    sigmaJ  = param[7];
    
    # Time increment
    dt = Mat/NT;
    
    # Expected value of k, and drift term
    kappa2 = exp(muJ) - 1;
    drift = r - q - lambdaJ*kappa2;
    
    # Initialize the variance and stock processes
    V = np.zeros((NT,NS));
    S = np.zeros((NT,NS));
    
    # Flags for negative variances
    F = 0;
    
    # Starting values for the variance and stock processes
    S[0,:] = S0;       # Spot price
    V[0,:] = v0;       # Heston v0 initial variance
    
    # Generate the stock and volatility paths
    for i in range(NS):
        for t in range(1,NT):
            # Generate two dependent N(0,1) variables with correlation rho
            Zv = random.normalvariate(0,1);
            Zs = rho*Zv + sqrt(1-rho**2)*random.normalvariate(0,1);
    
            if scheme =='E':
                # Euler discretization for the variance
                V[t,i] = V[t-1,i] + kappa*(theta-V[t-1,i])*dt \
                    + sigma*sqrt(V[t-1,i]*dt)*Zv;
            elif scheme =='M':
                # Milstein discretization for the variance.
                V[t,i] = V[t-1,i] + kappa*(theta-V[t-1,i])*dt \
                    + sigma*sqrt(V[t-1,i]*dt)*Zv \
                    + (1/4)*sigma**2*dt*(Zv**2-1);
            elif scheme =='IM':
                # Implicit Milstein for the variance.
                V[t,i] = (V[t-1,i] + kappa*theta*dt \
                    + sigma*sqrt(V[t-1,i]*dt)*Zv \
                    + sigma**2*dt*(Zv**2-1)/4) / (1+kappa*dt);
            elif scheme =='WM':
                # Weighted Explicit-Implicit Milstein Scheme
                V[t,i] = (V[t-1,i] + kappa*(theta-alpha*V[t-1,i])*dt \
                    + sigma*sqrt(V[t-1,i]*dt)*Zv\
                    + sigma**2*dt*(Zv**2-1)/4) / (1+(1-alpha)*kappa*dt);
    
            # Apply the full truncation or reflection scheme to the variance
            if V[t,i] <= 0:
                F = F+1;
                if 'R' == negvar:          # Reflection: take -V
                    V[t,i] = abs(V[t,i]);
                elif 'T' in negvar:
                    V[t,i] = max(0, V[t,i]);   # Truncation: take max(0,V)
    
            # Simulate the lognormal jumps
            J = 0;
            if lambdaJ != 0:
                Nt = np.random.poisson(lambdaJ*dt);
                if Nt > 0:
                    for x in range(Nt):
                        J = J + random.normalvariate(muJ - sigmaJ**2/2,sigmaJ);

            # Discretize the log stock price
            S[t,i] = S[t-1,i]*exp((drift-V[t-1,i]/2)*dt + J + sqrt(V[t-1,i]*dt)*Zs);

    return [S, V, F]    

# In[EulerMilsteinPrice]
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
    
