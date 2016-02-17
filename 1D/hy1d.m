function [out] = hy1d(hy,ez,x,t,deltax,deltat,u)
    out = hy(x+0.5,t-0.5) + deltat*(ez(x+1,t)-ez(x,t))/(u*deltax);
end

