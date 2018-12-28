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