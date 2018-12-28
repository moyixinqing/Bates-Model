% Effect of mean jump parameter muJ on Bates implied volatility

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

% Settings for bisection algorithm
a = 0.01;
b = 3.0;
Tol = 1e-5;
MaxIter = 1000;

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
rho     =  0;
% Jump parameters
lambdaJ = 1;
muJ     = [-0.20 -0.05 0.05 0.20];
sigmaJ  = 0.17;

% Strike range
K = 25:50;
T = .18:.1:1;

%% Mean jump size affects skewness
for i=1:4
    for t=1:length(T)
        for k=1:length(K)
            param = [kappa theta sigma v0 rho lambdaJ muJ(i) sigmaJ];
            BatesPrice1 = SVprice('Bates', PutCall,param,T(t),K(k),S,r,q,x,w);
            IV(k,t,i) = BisecBSIV(PutCall,S,K(k),r,q,T(t),a,b,BatesPrice1,Tol,MaxIter);
            if IV(k,t,i) < 0
                IV(k,t,i) = nan;
            end
        end
    end
end

%% Plot the implied vols
for j=1:4
    subplot(2,2,j)
    surf(IV(:,:,j))
    zlim([0 .4]);
    title(['Mu = ' num2str(muJ(j))])
    ylabel('Strike')
    xlabel('Maturity')
end
