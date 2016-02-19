clc
clear
close all
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

hy = zeros(size,1);
ez = zeros(size,1);

% hy(3,1) = 0;
% hy(1,3) = 42;
% ez(2,2) = 1;
% ez(4,2) = 2;
%%

figure
ylim([-1.5,1.5])
hold on
for t=1:size
    t
    if t < 500
        ez(1) = cos(2*pi*f*deltat*(t-1));
    else
        ez(1) = 0;
    end
    for x=1:size-1
        hy(x)= hy(x) + deltat*(ez(x+1)-ez(x))/(u*deltax);
    end
    for x=2:size
        ez(x)= ez(x) + deltat*(hy(x)-hy(x-1))/(eps*deltax);
    end
    plot(t,ez(50),'.')
    drawnow
end
