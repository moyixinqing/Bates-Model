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
        
        psi = BatesCF(v-(alpha+1)*i,param,S0,r,q,tau);
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
                psi[j] = BatesCF(v[j]-(alpha+1)*i,param,S0,r,q,tau);
                phi[j] = exp(-r*tau)*psi[j]/(alpha^2 + alpha - v[j]^2 + i*v[j]*(2*alpha+1));
                x[j] = exp(i*(b-s0)*v[j])*phi[j]*w[j];
                e[j] = exp(-i*2*pi/N*(j-1)*(u-1))*x[j];
            
            CallFFT[u] = eta*exp(-alpha*k[u])/pi * e.sum().real;
        
    return [CallFFT, K, lambdainc, eta] 