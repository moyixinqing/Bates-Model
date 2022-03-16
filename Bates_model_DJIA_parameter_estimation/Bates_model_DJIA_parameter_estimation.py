import numpy as np
from scipy.stats import norm
from math import exp,sqrt,log

# Black Scholes call and put
#BSC = @(s,K,rf,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));
#BSP = @(s,K,rf,q,v,T) (K*exp(-rf*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)));
def BSC(s,K,rf,q,v,T):
    BSC=(s*exp(-q*T)*norm.cdf((log(s/K) + (rf-q+v**2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*norm.cdf((log(s/K) + (rf-q+v**2/2)*T)/v/sqrt(T) - v*sqrt(T)));
    return BSC

def BSP(s,K,rf,q,v,T): 
    BSP= (K*exp(-rf*T)*norm.cdf(-(log(s/K) + (rf-q+v**2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*norm.cdf(-(log(s/K) + (rf-q+v**2/2)*T)/v/sqrt(T)));
    return BSP

## Load the DIA data
MktIV = np.array([[
    19.62,    19.47,    20.19,    21.15],
    [19.10,    19.05,    19.80,    20.82],
    [18.60,    18.61,    19.43,    20.57],
    [18.10,    18.12,    19.07,    20.21],
    [17.61,    17.74,    18.71,    20.00],
    [17.18,    17.43,    18.42,    19.74],
    [16.71,    17.06,    18.13,    19.50],
    [16.44,    16.71,    17.83,    19.27],
    [16.45,    16.41,    17.60,    18.99],
    [16.61,    16.25,    17.43,    18.84],
    [17.01,    16.02,    17.26,    18.62],
    [17.55,    16.10,    17.16,    18.46],
    [17.96,    16.57,    17.24,    18.42]])/100;        
         
K = np.arange(124,137);
Days = [37, 72, 135, 226, 278]
T = list(np.array(Days)/365);
[NK, NT] = MktIV.shape;
S = 129.14;
rf = 0.0010;
q  = 0.0068;
PutCall = np.tile('P',(NK,NT));
         
## Find the market prices
MktPrice=np.zeros((NK,NT))
for k in range(NK):
    for t in range(NT):
        if PutCall[k,t] == 'C':
            MktPrice[k,t] = BSC(S,K[k],rf,q,MktIV[k,t],T[t]);
        else: 
            MktPrice[k,t] = BSP(S,K[k],rf,q,MktIV[k,t],T[t]);

## Parameter Estimation
# Starting values
#     kappa,theta,sigma,v0,rho,lambdaJ,muJ,sigmaJ
#start = [5, 0.05, 2, 0.04, -0.7, 0.1, 0.05, 0.10];
start = [3, 0.25, 4, 0.06,  0.13,-0.02, 3, 3] # <-- local minimization is sensitive to initial start point. 

# Specify the objective function
# 1 = MSE | 2 = RMSE | 3 = IVMSE | 4 = Christoffersen, Heston, Jacobs
ObjFun = 4;

## Find the parameter estimates using Strike Vector Computation
# Settings for the Bisection method to find the Model IV
a = .001;
b = 3;
Tol = 1e-7;
MaxIter = 1000;
e = 1e-2;
inf=np.inf
#    kappa theta sigma v0  rho   lambdaJ  muJ  sigmaJ
lb = [e,   e,    e,    e, -.999,  -inf,   -inf, e ];  # Lower bound on the estimates
ub = [20,  5,    5,    5,  .999,   inf,    inf, 10];  # Upper bound on the estimates
#         kappa theta sigma   v0   rho           lambdaJ     muJ     sigmaJ
bounds= ((e,20),(e,5),(e,5),(e,5),(-0.999,0.999),(-inf,inf),(-inf,inf),(e,10))
# Gauss-Laguerre weights and abscissas
[x, w] = GenerateGaussLaguerre(32);

#%%
## Find the parameter estimates
tic()
import scipy.optimize
def obj(p):
    return BatesObjFunSVC(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,ObjFun,a,b,Tol,MaxIter)
OPresult = scipy.optimize.minimize(fun=obj, x0=start, method='SLSQP', bounds=bounds, options={'gtol': 1e-6, 'disp': True})
param = OPresult.x
t1 = toc()

## Print out the estimates

print('Bates parameter estimates for the DJIA \n')
print('Parameter   Estimate \n')
print('-------------------- \n')
print('kappa      {:.4f} \n'.format(param[0]))
print('theta      {:.4f} \n'.format(param[1]))
print('sigma      {:.4f} \n'.format(param[2]))
print('v0         {:.4f} \n'.format(param[3]))
print('rho        {:.4f} \n'.format(param[4]))
print('lambdaJ    {:.4f} \n'.format(param[5]))
print('numJ       {:.4f} \n'.format(param[6]))
print('sigmaJ     {:.4f} \n'.format(param[7]))
print('-------------------- \n')
print('Estimation time {:.3f} seconds \n'.format(t1));
#%%
#param =    [2.7997,    0.2757,    4.8352,    0.0646,    0.1248,   -0.0234,    3.1484,    3.1797]

## Find the Bates prices and the implied volatilities
BPrice = np.zeros((NK,NT));
BatesIV = np.zeros((NK,NT));
for k in range(NK):
    for t in range(NT):
        # Bates Call price with Gauss Laguerre integration
        BPrice[k,t] = BatesPrice(PutCall[k,t],param,T[t],K[k],S,rf,q,x,w);
        BatesIV[k,t] = BisecBSIV(PutCall[k,t],S,K[k],rf,q,T[t],a,b,BPrice[k,t],Tol,MaxIter);

## Plot the implied volatilities
import matplotlib.pyplot as plt
png_loc=[]
name='Market and Bates Implied Volatilities from Puts on DIA'
path= ''
fig = plt.figure(figsize = (10,8))
plt.text(0.5, 1.08, name,
         horizontalalignment='center',
         fontsize=20)
plt.axis('off')
for t  in range(NT):
    ax = fig.add_subplot(2,2,t+1)
    ax.plot(K,MktIV[:,t],'kx-',label = 'Market IV')
    ax.plot(K,BatesIV[:,t],'ro-',label = 'Bates IV')
    ax.set_title('Maturity {} days'.format(int(T[t]*365)))
    ax.grid()
    ax.legend()
plt.show()
fig.savefig(path +'{}.jpg'.format(name))
png_loc.append(path + '{}.jpg'.format(name))
