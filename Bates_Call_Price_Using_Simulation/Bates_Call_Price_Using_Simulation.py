# Euler or Milstein simulation for the Bates Model
# Fabrice Douglas Rouah, FRouah.com and Volopta.com 

# Spot price, strike, risk free rate, dividend yield, maturity, flavor
S0 = 50;
K  = 50;
r = 0.03;
q = 0.01;
tau = .25;
PutCall = 'C';

# Heston parameters
kappa  =  0.2;
theta  =  0.05;
sigma  =  0.3;
rho    = -0.7;
v0     =  0.05;

# Jump parameters
lambdaJ =  1;
muJ     = 0.1;
sigmaJ  = 0.60;
params = [kappa, theta, sigma, v0, rho, lambdaJ, muJ, sigmaJ];

# Number of time steps, number of stock price paths
NT = 90;
NS = 20000;

# Simulation scheme, negative variance overrides, alpha for weighted scheme
scheme = 'M';
negvar = 'R';
alpha = 0;

# Obtain the simulated price
[S, V, F, SimPrice] = EulerMilsteinPrice(scheme,negvar,params,PutCall,S0,K,tau,r,q,NT,NS,alpha);

# Gauss-Laguerre weights and abscissas
[x, w] = GenerateGaussLaguerre(32);

# Obtain the closed form price
ClosedPrice = BatesPrice(PutCall,params,tau,K,S0,r,q,x,w);

## Output everything
print('Bates price by Euler or Milstein Simulation simulation \n')
print('-------------------------------------\n')
print('Method         Price     DollarError \n')
print('-------------------------------------\n')
print('Closed Form   {:.4f}      n/a   \n'.format(ClosedPrice))
print('Simulation    {:.4f}   {:.4f}  \n'.format(SimPrice,ClosedPrice-SimPrice))
print('-------------------------------------\n')

