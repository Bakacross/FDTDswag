clc
clear

% Constant (void)
u = 4*pi*1e-12;
eps = 8.8541878176e-12;

% Option
deltax = 0.001;
deltat = 0.001;
size = 100;

% Initial condition

hy = zeros(100,100);
ez = zeros(100,100);
ez(2,2) = 1;

% hy(3,1) = 0;
% hy(1,3) = 42;
% ez(2,2) = 1;
% ez(4,2) = 2;
%%
for t=2:2:4
    x=t
    hy(x+1,t+1)= hy(x+1,t-1) + deltat*(ez(x+2,t)-ez(x,t))/(u*deltax);
    ez(x,t+2)= ez(x,t) + deltat*(hy(x+1,t+1)-hy(x-1,t+1))/(eps*deltax);
    
end
