function y = BatesObjFunSVC(param,S,rf,q,MktPrice,K,T,PutCall,MktIV,x,w,ObjFun,a,b,Tol,MaxIter)

% S  = Spot price
% rf = risk free rate
% q  = dividend yield
% MktPrice = vector of market prices.
% K = vector of strikes.
% T = vector of maturities.
% PutCall = vector of put/call indicator 'C' or 'P'
% MktIV = vector of market implied volatilities
% x = abscissas for Gauss Laguerre integration
% w = weights for Gauss Laguerre integration
% ObjFun = Type of Objective Function
%    1 = MSE
%    2 = RMSE
%    3 = IVMSE
%    4 = Christoffersen, Heston, Jacobs
% a = Bisection algorithm, small initial estimate
% b = Bisection algorithm, large initial estimate
% Tol = Bisection algorithm, tolerance
% MaxIter = Bisection algorithm, maximum iterations

[NK,NT] = size(MktPrice);

for t=1:NT
    for j=1:length(x);
        phi = x(j);
        f(j)    = BatesCF(phi  ,param,T(t),S,rf,q);
        fnum(j) = BatesCF(phi-i,param,T(t),S,rf,q);
        fden(j) = BatesCF(   -i,param,T(t),S,rf,q);
    end
    for k=1:NK
        for j=1:length(x);
            phi = x(j);
            int2(j) = w(j) * real(exp(-i*phi*log(K(k)))*f(j)/i/phi);
            int1(j) = w(j) * real(exp(-i*phi*log(K(k)))*fnum(j)/i/phi/fden(j));
        end
        % The probabilities P1 and P2
        P1 = 1/2 + 1/pi*sum(int1);
        P2 = 1/2 + 1/pi*sum(int2);
        % The call price
        CallPrice = S*exp(-q*T(t))*P1 - K(k)*exp(-rf*T(t))*P2;
        % Output the option price (put price by put call parity)
        if strcmp(PutCall,'C')
            ModelPrice(k,t) = CallPrice;
        else
            ModelPrice(k,t) = CallPrice - S*exp(-q*T(t)) + K(k)*exp(-rf*T(t));
        end
        % Select the objective function
        if ObjFun == 1
            % MSE
            error(k,t) = (MktPrice(k,t) - ModelPrice(k,t))^2;
        elseif ObjFun == 2
            % RMSE
            error(k,t) = (MktPrice(k,t) - ModelPrice(k,t))^2 / MktPrice(k,t);
        elseif ObjFun == 3
            % IVMSE
            ModelIV = BisecBSIV(PutCall(k,t),S,K(k),rf,q,T(t),a,b,ModelPrice(k,t),Tol,MaxIter);
            error(k,t) = (ModelIV - MktIV(k,t))^2;
        elseif ObjFun == 4
            % IVRMSE Christoffersen, Heston, Jacobs proxy
            d = (log(S/K(k)) + (rf-q+MktIV(k,t)^2/2)*T(t))/MktIV(k,t)/sqrt(T(t));
            Vega(k,t) = S*normpdf(d)*sqrt(T(t));
            error(k,t) = (ModelPrice(k,t) - MktPrice(k,t))^2 / Vega(k,t)^2;
        end
    end
end
y = sum(sum(error));


% Bates characteristic function
function y = BatesCF(phi,param,T,S,rf,q)

y = HestonCF(phi,param,T,S,rf,q) * JumpCF(phi,param,T);


% Logarithmic normal jump characteristic function
function jcf = JumpCF(phi,param,T)
lambdaJ  = param(6);        % Annual jump frequency
muJ      = param(7);        % Random percentage jump
sigmaJ   = param(8);        % Jump volatility
jcf = exp(-lambdaJ*muJ*i*phi*T + lambdaJ*T*((1+muJ)^(i*phi)*exp(0.5*sigmaJ^2*i*phi*(i*phi-1))-1));


% Heston (1993) second characteristic function (f2) using the "Little Trap" formulation
function f2 = HestonCF(phi,param,T,S,rf,q)
kappa = param(1);
theta = param(2);
sigma = param(3);
v0    = param(4);
rho   = param(5);
lambda = 0;

% Log of the stock price.
x = log(S);

% Required parameters.
a = kappa*theta;
u = -0.5;
b = kappa + lambda;

d = sqrt((rho*sigma*i*phi - b)^2 - sigma^2*(2*u*i*phi - phi^2));
g = (b - rho*sigma*i*phi + d) / (b - rho*sigma*i*phi - d);

% "Little Heston Trap" formulation
c = 1/g;
D = (b - rho*sigma*i*phi - d)/sigma^2*((1-exp(-d*T))/(1-c*exp(-d*T)));
G = (1 - c*exp(-d*T))/(1-c);
C = (rf-q)*i*phi*T + a/sigma^2*((b - rho*sigma*i*phi - d)*T - 2*log(G));

% The characteristic function.
f2 = exp(C + D*v0 + i*phi*x);

