clc
clear

% Constant (void)
u = 1.25663706*1e-6;
eps = 8.8541878176e-12;
c = 299792458;
f = 2.4*1e9;
lambda = c/f;

% Option
deltax = lambda/10;
deltat = deltax/c/10;
size = 2000;

% Initial condition

hy = zeros(size,size);
ez = zeros(size,size);

% hy(3,1) = 0;
% hy(1,3) = 42;
% ez(2,2) = 1;
% ez(4,2) = 2;
%%
for t=2:2:size
    t
    if t < 500
        ez(2,t) = cos(2*pi*f*deltat*(t-2));
    end
    for x=2:2:size-2
        hy(x+1,t+1)= hy(x+1,t-1) + deltat*(ez(x+2,t)-ez(x,t))/(u*deltax);
    end
    for x=4:2:size-2
        ez(x,t+2)= ez(x,t) + deltat*(hy(x+1,t+1)-hy(x-1,t+1))/(eps*deltax);
    end
end
