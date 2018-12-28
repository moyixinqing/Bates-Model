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