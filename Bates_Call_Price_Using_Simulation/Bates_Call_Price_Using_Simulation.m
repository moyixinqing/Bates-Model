% Euler or Milstein simulation for the Bates Model
% Fabrice Douglas Rouah, FRouah.com and Volopta.com 

clc; clear;

% Spot price, strike, risk free rate, dividend yield, maturity, flavor
S0 = 50;
K  = 50;
r = 0.03;
q = 0.01;
tau = .25;
PutCall = 'C';

% Heston parameters
kappa  =  0.2;
theta  =  0.05;
sigma  =  0.3;
rho    = -0.7;
v0     =  0.05;

% Jump parameters
lambdaJ =  1;
muJ     = 0.1;
sigmaJ  = 0.60;
params = [kappa theta sigma v0 rho lambdaJ muJ sigmaJ];

% Number of time steps, number of stock price paths
NT = 90;
NS = 20000;

% Simulation scheme, negative variance overrides, alpha for weighted scheme
scheme = 'M';
negvar = 'R';
alpha = 0;

% Obtain the simulated price
[S V F SimPrice] = EulerMilsteinPrice(scheme,negvar,params,PutCall,S0,K,tau,r,q,NT,NS,alpha);

% Gauss-Laguerre weights and abscissas
xw = [...
    0.0445    0.1142;
    0.2345    0.2661;
    0.5769    0.4188;
    1.0724    0.5725;
    1.7224    0.7276;
    2.5283    0.8845;
    3.4922    1.0436;
    4.6165    1.2053;
    5.9040    1.3702;
    7.3581    1.5388;
    8.9829    1.7116;
   10.7830    1.8895;
   12.7637    2.0730;
   14.9314    2.2633;
   17.2914    2.4622;
   19.8586    2.6600;
   22.6270    2.8909;
   25.6295    3.1235;
   28.8726    3.3645;
   32.3193    3.6697;
   36.1371    4.1866;
   40.1256    4.1695;
   44.4829    4.4523;
   49.3076    4.9941;
   54.2188    5.8800;
   59.9992    5.3379;
   65.9033    6.8224;
   72.7216    6.8573;
   80.1762    8.0512;
   88.7377    9.1847;
   98.8293   11.1658;
  111.7514   15.3900];
x = xw(:,1);
w = xw(:,2);

% Obtain the closed form price
ClosedPrice = BatesPrice(PutCall,params,tau,K,S0,r,q,x,w);

%% Output everything
fprintf('Bates price by Euler or Milstein Simulation simulation \n')
fprintf('-------------------------------------\n')
fprintf('Method         Price     DollarError \n')
fprintf('-------------------------------------\n')
fprintf('Closed Form   %8.4f      n/a   \n',ClosedPrice)
fprintf('Simulation    %8.4f   %8.4f  \n',SimPrice,ClosedPrice-SimPrice)
fprintf('-------------------------------------\n')

