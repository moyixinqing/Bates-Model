clc; clear;

% Black Scholes call and put
BSC = @(s,K,rf,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));
BSP = @(s,K,rf,q,v,T) (K*exp(-rf*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)));
 

%% Load the DIA data
MktIV = [...
    19.62	19.47	20.19	21.15;
    19.10	19.05	19.80	20.82;
    18.60	18.61	19.43	20.57;
    18.10	18.12	19.07	20.21;
    17.61	17.74	18.71	20.00;
    17.18	17.43	18.42	19.74;
    16.71	17.06	18.13	19.50;
    16.44	16.71	17.83	19.27;
    16.45	16.41	17.60	18.99;
    16.61	16.25	17.43	18.84;
    17.01	16.02	17.26	18.62;
    17.55	16.10	17.16	18.46;
    17.96	16.57	17.24	18.42]./100;
K = (124:136);
T = [37 72 135 226 278]./365;
S = 129.14;
rf = 0.0010;
q  = 0.0068;

[NK NT] = size(MktIV);
PutCall = repmat('P',NK,NT);

%% Find the market prices
for k=1:NK
    for t=1:NT
        if PutCall(k,t) == 'C';
            MktPrice(k,t) = BSC(S,K(k),rf,q,MktIV(k,t),T(t));
        else 
            MktPrice(k,t) = BSP(S,K(k),rf,q,MktIV(k,t),T(t));
        end
    end
end


%% Parameter Estimation
% Starting values
% kappa,theta,sigma,v0,rho,lambdaJ,muJ,sigmaJ
start = [5, 0.05, 2, 0.04, -0.7, 0.1, 0.05, 0.10];

% Specify the objective function
% 1 = MSE | 2 = RMSE | 3 = IVMSE | 4 = Christoffersen, Heston, Jacobs
ObjFun = 4;


%% Find the parameter estimates using Strike Vector Computation
% Settings for the Bisection method to find the Model IV
a = .001;
b = 3;
Tol = 1e-7;
MaxIter = 1000;
e = 1e-2;
%    kappa theta sigma v0  rho   lambdaJ  muJ  sigmaJ
lb = [e    e     e     e  -.999  -inf    -inf  e ];  % Lower bound on the estimates
ub = [20   5     5     5   .999   inf     inf  10];  % Upper bound on the estimates

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

%% Find the parameter estimates
tic
[param feval] = fmincon(@(p) BatesObjFunSVC(p,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,ObjFun,a,b,Tol,MaxIter),start,[],[],[],[],lb,ub);
t1 = toc;

%% Print out the estimates
clc;
fprintf('Bates parameter estimates for the DJIA \n')
fprintf('Parameter   Estimate \n')
fprintf('-------------------- \n')
fprintf('kappa      %8.4f \n',param(1))
fprintf('theta      %8.4f \n',param(2))
fprintf('sigma      %8.4f \n',param(3))
fprintf('v0         %8.4f \n',param(4))
fprintf('rho        %8.4f \n',param(5))
fprintf('lambdaJ    %8.4f \n',param(6))
fprintf('numJ       %8.4f \n',param(7))
fprintf('sigmaJ     %8.4f \n',param(8))
fprintf('-------------------- \n')
fprintf('Estimation time %6.3f seconds \n',t1);
fprintf('Value of loss function %10.8f \n',feval);


%% Find the Bates prices and the implied volatilities
for k=1:NK
	for t=1:NT
        % Bates Call price with Gauss Laguerre integration
        BPrice(k,t) = BatesPrice(PutCall(k,t),param,T(t),K(k),S,rf,q,x,w);
        BatesIV(k,t) = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,BPrice(k,t),Tol,MaxIter);
	end
end

%% Plot the implied volatilities
for t=1:NT
	subplot(2,2,t)
	plot(K,MktIV(:,t),'kx-',K,BatesIV(:,t),'ro-')
	title(['Maturity ' num2str(T(t)*365) ' days'])
 	legend('Market IV', 'Bates IV')
end


