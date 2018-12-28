# Effect of mean jump parameter muJ on Bates implied volatility

# Fabrice Rouah, FRouah.com and Volopta.com
import numpy as np
# Gauss Laguerre abscissas and weights
[x, w] = GenerateGaussLaguerre(32);

# Settings for bisection algorithm
a = 0.01;
b = 3.0;
Tol = 1e-5;
MaxIter = 1000;

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
rho     =  0;
# Jump parameters
lambdaJ = 1;
muJ     = [-0.20, -0.05, 0.05, 0.20];
sigmaJ  = 0.17;

# Strike range
K = np.arange(25,51)
T = np.arange(.18,1+.1,.1)

IV = np.zeros((len(K),len(T),4));
## Mean jump size affects skewness
for i in range(4):
    for t in range(len(T)):
        for k in range(len(K)):
            param = [kappa, theta, sigma, v0, rho, lambdaJ, muJ[i], sigmaJ];
            BatesPrice1 = SVprice('Bates', PutCall,param,T[t],K[k],S,r,q,x,w);
            IV[k,t,i] = BisecBSIV(PutCall,S,K[k],r,q,T[t],a,b,BatesPrice1,Tol,MaxIter);
            if IV[k,t,i] < 0:
                IV[k,t,i] = np.nan;
#%%
## Plot the implied vols
X = T
Y = K
X, Y = np.meshgrid(X, Y)
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
path= 'C:/temp/Jump-Process-master/Bates_Effect_of_jump_parameters_on_RND/'
name='Implied and Local Volatility'
fig = plt.figure(figsize=(12,10))
plt.title('Bates_Effect_of_jump_parameters_on_RND')
plt.axis('off')
for j in range(4):
    ax = fig.add_subplot(2,2,j+1, projection='3d')
    #ax.plot_wireframe(X, Y, IV[:,:,j], rstride=1, cstride=1)
    ax.plot_surface(X, Y, IV[:,:,j], rstride=2, cstride=2, alpha=0.25,cmap=cm.coolwarm, linewidth=0, antialiased=False)
    #ax.set_xlim3d(T[0], T[-1]);
    #ax.set_ylim3d(K[0], K[-1]);
    ax.set_zlim3d(0, 0.4);
    ax.set_xlabel('Maturity');
    ax.set_ylabel('Strike')
    ax.set_zlabel('Volatility')
    ax.set_title('Mu = {}'.format(muJ[j]))
# =============================================================================
#     for angle in range(0, 360):
#         ax.view_init(30, angle)
#         plt.draw()
#         plt.pause(.01)
# =============================================================================
plt.show()
fig.savefig(path +'{}.jpg'.format(name))
