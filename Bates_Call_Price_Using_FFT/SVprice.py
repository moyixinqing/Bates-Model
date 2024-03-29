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