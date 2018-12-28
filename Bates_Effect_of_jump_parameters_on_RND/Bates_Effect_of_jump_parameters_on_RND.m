% Effect of jump parameters on Bates implied volatility

% Fabrice Rouah, FRouah.com and Volopta.com

clc; clear;
% Gauss Laguerre abscissas and weights
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


% Global settings for the call prices
PutCall = 'C';
S = 40;
T = 0.25;
r = 0.08;
q = 0.02;
kappa = 4;

% Heston parameters
theta   =  0.0225;
sigma   =  0.20;
v0      =  0.04;
rho     = -0.7;
% Jump parameters
lambdaJ = [0 1 2];
muJ     = 0.1;
sigmaJ  = 0.17;

% Strike range
K = 25:.25:55;

%% Increasing lambda increases volatility
for k=1:length(K)
    param1 = [kappa theta sigma v0 rho lambdaJ(1) muJ sigmaJ];
    BatesPrice1(k) = SVprice('Bates', PutCall,param1,T,K(k),S,r,q,x,w);
    param2 = [kappa theta sigma v0 rho lambdaJ(2) muJ sigmaJ];
    BatesPrice2(k) = SVprice('Bates', PutCall,param2,T,K(k),S,r,q,x,w);
    param3 = [kappa theta sigma v0 rho lambdaJ(3) muJ sigmaJ];
    BatesPrice3(k) = SVprice('Bates', PutCall,param3,T,K(k),S,r,q,x,w);
end

% Extract the risk neutral densities
[RND1, K2] = ExtractRND(K,BatesPrice1);
[RND2, K2] = ExtractRND(K,BatesPrice2);
[RND3, K2] = ExtractRND(K,BatesPrice3);

subplot(2,2,1)
plot(K2,RND1,'kx-',K2,RND2,'rx-',K2,RND3,'bx-')
legend('Bates,lambda=0','Bates,lambda=1','Bates,lambda=2','Location','NorthWest')
ylim([0 0.12])

%% Mean jump size affects skewness
lambdaJ = 1;
muJ = [-0.1 0 0.1];

for k=1:length(K)
    param1 = [kappa theta sigma v0 rho lambdaJ muJ(1) sigmaJ];
    BatesPrice1(k) = SVprice('Bates', PutCall,param1,T,K(k),S,r,q,x,w);
    param2 = [kappa theta sigma v0 rho lambdaJ muJ(2) sigmaJ];
    BatesPrice2(k) = SVprice('Bates', PutCall,param2,T,K(k),S,r,q,x,w);
    param3 = [kappa theta sigma v0 rho lambdaJ muJ(3) sigmaJ];
    BatesPrice3(k) = SVprice('Bates', PutCall,param3,T,K(k),S,r,q,x,w);
end

% Extract the risk neutral densities
[RND1, K2] = ExtractRND(K,BatesPrice1);
[RND2, K2] = ExtractRND(K,BatesPrice2);
[RND3, K2] = ExtractRND(K,BatesPrice3);

subplot(2,2,2)
plot(K2,RND1,'kx-',K2,RND2,'rx-',K2,RND3,'bx-')
legend('Bates,mu = -0.1','Bates,mu = 0','Bates,mu = +0.1','Location','NorthWest')
ylim([0 0.12])


%% Jump size variance affects kurtosis
lambdaJ = 1;
muJ = 0;
sigmaJ = [0.05 0.10 0.20];

for k=1:length(K)
    param1 = [kappa theta sigma v0 rho lambdaJ muJ sigmaJ(1)];
    BatesPrice1(k) = SVprice('Bates', PutCall,param1,T,K(k),S,r,q,x,w);
    param2 = [kappa theta sigma v0 rho lambdaJ muJ sigmaJ(2)];
    BatesPrice2(k) = SVprice('Bates', PutCall,param2,T,K(k),S,r,q,x,w);
    param3 = [kappa theta sigma v0 rho lambdaJ muJ sigmaJ(3)];
    BatesPrice3(k) = SVprice('Bates', PutCall,param3,T,K(k),S,r,q,x,w);
end

% Extract the risk neutral densities
[RND1, K2] = ExtractRND(K,BatesPrice1);
[RND2, K2] = ExtractRND(K,BatesPrice2);
[RND3, K2] = ExtractRND(K,BatesPrice3);

subplot(2,2,3)
plot(K2,RND1,'kx-',K2,RND2,'rx-',K2,RND3,'bx-')
legend('Bates,sigma = 0.05','Bates,sigma = 0.10','Bates,sigma = 0.20','Location','NorthWest')
ylim([0 0.12])
