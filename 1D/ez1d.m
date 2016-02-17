function [out] = ez1d(hy,ez,x,t,deltax,deltat,eps)
    out = ez(x,t) + deltat*(hy(x+1,t+1)-hy(x-1,t+1))/(eps*deltax);
end

