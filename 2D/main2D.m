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
deltay = lambda/10;
deltat = deltax/c/10;
size = 2000;

% Initial condition
hx = zeros(size,size);
hy = zeros(size,size);
ez = zeros(size,size,size);

% hy(3,1) = 0;
% hy(1,3) = 42;
% ez(2,2) = 1;
% ez(4,2) = 2;
%%
for t=2:2:size
    t
    ez(2,2,t) = cos(2*pi*f*deltat*(t-2));
    for x=2:2:size-2
        hy(x+1,t+1)= hy(x+1,t-1) + deltat*(ez(x+2,t)-ez(x,t))/(u*deltax);
    end
    
    clc
clear

% Constant (void)
u = 1.25663706*1e-6;
eps = 8.8541878176e-12;
c = 299792458;
f = 2.4*1e9;
lambda = c/f;
om = 1; %%% A definir becuz bs 
o = 1; %%% corriger bs fouad becuz it is bs hihihi
% Option
deltax = lambda/10;
deltat = deltax/c/10;
size = 2000;
c1 = (1-(om*deltat)/(2*u))/(1+(om*deltat)/(2*u));
c2 = deltat/(u*deltay)/(1+(om*deltat)/(2*u));
c3 = deltat/(u*deltax)/(1+(om*deltat)/(2*u));
c4 = (1-(om*deltat)/(2*eps))/(1+(om*deltat)/(2*eps));
c5 = 1/(1+(o*deltat/(2*eps)));
% Initial condition
hx = zeros(size,1);
hy = zeros(size,1);
ez = zeros(size,size,1);

% hy(3,1) = 0;
% hy(1,3) = 42;
% ez(2,2) = 1;
% ez(4,2) = 2;
%%
for t=2:2:size
    t
    ez(2,2,t) = cos(2*pi*f*deltat*(t-2));
    for y=2:2:size-2
        hx(x)= c1*hx(x) + c2*(ez(x)-ez(x,t))/(u*deltax);
        for x=2:2:size-2
            ez(x,t+2)= ez(x,t) + deltat*(hy(x+1,t+1)-hy(x-1,t+1))/(eps*deltax);
        end
    end
    
end
    for x=4:2:size-2
        ez(x,t+2)= ez(x,t) + deltat*(hy(x+1,t+1)-hy(x-1,t+1))/(eps*deltax);
    end
end