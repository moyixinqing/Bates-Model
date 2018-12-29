# Reproduction of Table 1 of Bates (1996)

# Fabrice Douglas Rouah, FRouah.com and Volopta.com
from module import BatesModelFunctions as bmf
import numpy as np

# Integration range settings for numerical integration
Lphi = 1e-50;  # Lower limit
Uphi = 100;    # Upper limit
N = 1000;      # Number of integration points
dphi = (Uphi-Lphi)/(N-1);

# Construct the weights and abscissas for the trapezoidal rule
#x = [Lphi:dphi:Uphi];
#w = [0.5 ones(1,length(x)-2) 0.5] .* dphi;
x = np.arange(Lphi,Uphi+dphi,dphi)
w = (np.ones(N-2)).tolist()
w = [1/2] + w +[1/2]
w = dphi*np.array(w)
    
# Global settings for the first entries of the five rows of Table 1
PutCall = 'P';
S = 40;
T = 0.25;
r = 0.08;
b = 0.02;
q = r - b;
kappa = 4;

# Parameters for the first entries of the five rows of Table 1
# Heston parameters
theta   = [0.0225, 0.0225, 0.0225, 0.0225, 0.0125];
sigma   = [0.15, 0.15, 0.30, 0.15, 0.20];
v0      = [0.0225, 0.04, 0.0225, 0.0225, 0.0125];
rho     = [0, 0, 0, 0.1, 0];
# Jump parameters
lambdaJ = [0, 0, 0, 0, 2];
muJ = 0;
sigmaJ  = [0, 0, 0, 0, 0.07];

# Strike range
#K = 38:42;
K = np.arange(38,43)

## Generate the prices for the first entries in the rows of Table 1
price = np.zeros((5,5))
for j in range(5):
    param = [kappa, theta[j], sigma[j], v0[j], rho[j], lambdaJ[j], muJ, sigmaJ[j]];
    for k in range(5):
        price[j,k] = bmf.SVprice('Bates',PutCall,param,T,K[k],S,r,q,x,w);

## Printthe results
print('Table 1 of Bates (1996), first entries in each row only \n')
print('European prices \n')
print('--------------------------------------------------------\n')
print('Set    K=38      K=39       K=40       K=41      K=42   \n')
print('--------------------------------------------------------\n')
import pandas as pd
result = pd.DataFrame(price, columns=['K=38','K=39', 'K=40', 'K=41', 'K=42'])
print(result.round(4))

