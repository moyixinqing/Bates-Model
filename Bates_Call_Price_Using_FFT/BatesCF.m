% Bates characteristic function
function bcf = BatesCF(phi,param,S,r,q,T)
bcf = HestonCF(phi,param,T,S,r,q) .* JumpCF(phi,param,T);
 
% Logarithmic normal jump characteristic function
function jcf = JumpCF(phi,param,T)
lambdaJ  = param(6);        % Annual jump frequency 
muJ      = param(7);        % Random percentage jump
sigmaJ   = param(8);        % Jump volatility
jcf = exp(-lambdaJ.*muJ.*i.*phi.*T + lambdaJ.*T.*((1+muJ).^(i.*phi).*exp(0.5.*sigmaJ.^2.*i.*phi.*(i.*phi-1))-1));

% Heston (1993) second characteristic function (f2) using the "Little Trap" formulation
function f2 = HestonCF(phi,param,T,S,r,q)
kappa = param(1);
theta = param(2);
sigma = param(3);
v0    = param(4);
rho   = param(5);
lambda = 0;

% Log of the stock price.
x = log(S);

% Required parameters.
a = kappa.*theta;
u = -0.5;
b = kappa + lambda;

d = sqrt((rho.*sigma.*i.*phi - b).^2 - sigma.^2.*(2.*u.*i.*phi - phi.^2));
g = (b - rho.*sigma.*i.*phi + d) ./ (b - rho.*sigma.*i.*phi - d);

% "Little Heston Trap" formulation
c = 1./g;
D = (b - rho.*sigma.*i.*phi - d)./sigma.^2.*((1-exp(-d.*T))./(1-c.*exp(-d.*T)));
G = (1 - c.*exp(-d.*T))./(1-c);
C = (r-q).*i.*phi.*T + a./sigma.^2.*((b - rho.*sigma.*i.*phi - d).*T - 2.*log(G));

% The characteristic function.
f2 = exp(C + D.*v0 + i.*phi.*x);

