function [out] = hy1d(hy,ez,x,t,deltax,deltat,u)
    out = hy(x+1,t-1) + deltat*(ez(x+2,t)-ez(x,t))/(u*deltax);
end

