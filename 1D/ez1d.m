function [out] = ez1d(hy,ez,x,t,deltax,deltat,eps)
    out = ez(x,t) + deltat*(hy(x+0.5,t+0.5)-hy(x-0.5,t+0.5))/(eps*deltax);
end

