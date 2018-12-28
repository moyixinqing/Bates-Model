# Effect of jump parameters on Bates implied volatility

# Fabrice Rouah, FRouah.com and Volopta.com
import numpy as np

# Gauss Laguerre abscissas and weights
[x, w] = GenerateGaussLaguerre(32);

# Global settings for the call prices
PutCall = 'C';
S = 40;
T = 0.25;
r = 0.08;
q = 0.02;
kappa = 4;

# Heston parameters
theta   =  0.0225;
sigma   =  0.20;
v0      =  0.04;
rho     = -0.7;
# Jump parameters
lambdaJ = [0, 1, 2]; # Jump frequency parameter
muJ     = 0.1; # Jump mean
sigmaJ  = 0.17; # Jump volatility

# Strike range
K = np.arange(25,55+.25,.25)

## Increasing lambda increases volatility
BatesPrice1 =  np.zeros(len(K));
BatesPrice2 =  np.zeros(len(K));
BatesPrice3 =  np.zeros(len(K));
for k in range(len(K)):
    param1 = [kappa, theta, sigma, v0, rho, lambdaJ[0], muJ, sigmaJ];
    BatesPrice1[k] = SVprice('Bates', PutCall,param1,T,K[k],S,r,q,x,w);
    param2 = [kappa, theta, sigma, v0, rho, lambdaJ[1], muJ, sigmaJ];
    BatesPrice2[k] = SVprice('Bates', PutCall,param2,T,K[k],S,r,q,x,w);
    param3 = [kappa, theta, sigma, v0, rho, lambdaJ[2], muJ, sigmaJ];
    BatesPrice3[k] = SVprice('Bates', PutCall,param3,T,K[k],S,r,q,x,w);

# Extract the risk neutral densities
[RND1, K2] = ExtractRND(K,BatesPrice1);
[RND2, K2] = ExtractRND(K,BatesPrice2);
[RND3, K2] = ExtractRND(K,BatesPrice3);

import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
png_loc=[]
name='Effect of Jump Parameters on Risk Netural Densities'
path= 'C:/temp/Jump-Process-master/Bates_Effect_of_jump_parameters_on_RND/'
fig = plt.figure(figsize = (10,8))
ax = fig.add_subplot(2,2,1)
ax.plot(K2,RND1,'kx-',label = 'Bates,lambda=0')
ax.plot(K2,RND2,'rx-',label = 'Bates,lambda=1')
ax.plot(K2,RND3,'bx-',label = 'Bates,lambda=2')
ax.set_ylim([0,0.12])
ax.set_xlim([20,60])
ax.grid()
ax.legend()

#%%
## Mean jump size affects skewness
lambdaJ = 1;
muJ = [-0.1, 0, 0.1];

BatesPrice1 =  np.zeros(len(K));
BatesPrice2 =  np.zeros(len(K));
BatesPrice3 =  np.zeros(len(K));
for k in range(len(K)):
    param1 = [kappa, theta, sigma, v0, rho, lambdaJ, muJ[0], sigmaJ];
    BatesPrice1[k] = SVprice('Bates', PutCall,param1,T,K[k],S,r,q,x,w);
    param2 = [kappa, theta, sigma, v0, rho, lambdaJ, muJ[1], sigmaJ];
    BatesPrice2[k] = SVprice('Bates', PutCall,param2,T,K[k],S,r,q,x,w);
    param3 = [kappa, theta, sigma, v0, rho, lambdaJ, muJ[2], sigmaJ];
    BatesPrice3[k] = SVprice('Bates', PutCall,param3,T,K[k],S,r,q,x,w);

# Extract the risk neutral densities
[RND1, K2] = ExtractRND(K,BatesPrice1);
[RND2, K2] = ExtractRND(K,BatesPrice2);
[RND3, K2] = ExtractRND(K,BatesPrice3);

ax = fig.add_subplot(2,2,2)
ax.plot(K2,RND1,'kx-',label = 'Bates,mu = -0.1')
ax.plot(K2,RND2,'rx-',label = 'Bates,mu = 0')
ax.plot(K2,RND3,'bx-',label = 'Bates,mu = +0.1')
ax.set_ylim([0,0.12])
ax.set_xlim([20,60])
ax.grid()
ax.legend()

#%%
## Jump size variance affects kurtosis
lambdaJ = 1;
muJ = 0;
sigmaJ = [0.05, 0.10, 0.20];

BatesPrice1 =  np.zeros(len(K));
BatesPrice2 =  np.zeros(len(K));
BatesPrice3 =  np.zeros(len(K));
for k in range(len(K)):
    param1 = [kappa, theta, sigma, v0, rho, lambdaJ, muJ, sigmaJ[0]];
    BatesPrice1[k] = SVprice('Bates', PutCall,param1,T,K[k],S,r,q,x,w);
    param2 = [kappa, theta, sigma, v0, rho, lambdaJ, muJ, sigmaJ[1]];
    BatesPrice2[k] = SVprice('Bates', PutCall,param2,T,K[k],S,r,q,x,w);
    param3 = [kappa, theta, sigma, v0, rho, lambdaJ, muJ, sigmaJ[2]];
    BatesPrice3[k] = SVprice('Bates', PutCall,param3,T,K[k],S,r,q,x,w);
    
# Extract the risk neutral densities
[RND1, K2] = ExtractRND(K,BatesPrice1);
[RND2, K2] = ExtractRND(K,BatesPrice2);
[RND3, K2] = ExtractRND(K,BatesPrice3);

ax = fig.add_subplot(2,2,3)
ax.plot(K2,RND1,'kx-',label = 'Bates,sigma = 0.05')
ax.plot(K2,RND2,'rx-',label = 'Bates,sigma = 0.10')
ax.plot(K2,RND3,'bx-',label = 'Bates,sigma = 0.20')
ax.set_ylim([0,0.12])
ax.set_xlim([20,60])
ax.grid()
ax.legend(loc='upper left')
fig.savefig(path +'{}.jpg'.format(name))
png_loc.append(path + '{}.jpg'.format(name))