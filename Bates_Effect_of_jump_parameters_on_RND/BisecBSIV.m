function y = BisecBSIV(PutCall,S,K,rf,q,T,a,b,MktPrice,Tol,MaxIter)

% Function for the Black Scholes call and put
BSC = @(s,K,rf,q,v,T) (s*exp(-q*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)) - K*exp(-rf*T)*normcdf((log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) - v*sqrt(T)));
BSP = @(s,K,rf,q,v,T) (K*exp(-rf*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T) + v*sqrt(T)) - s*exp(-q*T)*normcdf(-(log(s/K) + (rf-q+v^2/2)*T)/v/sqrt(T)));
  
if strcmp(PutCall,'C')
	lowCdif  = MktPrice - BSC(S,K,rf,q,a,T);
	highCdif = MktPrice - BSC(S,K,rf,q,b,T);
else
	lowCdif  = MktPrice - BSP(S,K,rf,q,a,T);
	highCdif = MktPrice - BSP(S,K,rf,q,b,T);
end

if lowCdif*highCdif > 0
	y = -1;
else
	for x=1:MaxIter
		midP = (a+b)/2;
		if strcmp(PutCall,'C')
			midCdif = MktPrice - BSC(S,K,rf,q,midP,T);
		else
			midCdif = MktPrice - BSP(S,K,rf,q,midP,T);
		end
		if abs(midCdif)<Tol
			break
		else
			if midCdif>0
				a = midP;
			else
				b = midP;
			end
		end
	end
	y = midP;
end

