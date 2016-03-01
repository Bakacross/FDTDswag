
    
clc
clear
close all

% Constant (void)
u = 1.25663706*1e-6;
eps = 8.8541878176e-12;
c = 299792458;
f = 2.4*1e9;
lambda = c/f;
om = 0; %%% A definir becuz bs 
o = 0; %%% corriger bs fouad becuz it is bs hihihi
% Option
deltax = lambda/10;
deltay = lambda/10;
deltat = deltax/c/10;
size = 50;
x=size;
y=x;
c1 = (1-(om*deltat)/(2*u))/(1+(om*deltat)/(2*u));
c2 = deltat/(u*deltay)/(1+(om*deltat)/(2*u));
c3 = deltat/(u*deltax)/(1+(om*deltat)/(2*u));
c4 = (1-(om*deltat)/(2*eps))/(1+(om*deltat)/(2*eps));
c5 = 1/(1+(o*deltat/(2*eps)));
% Initial condition
hx = zeros(x,y);
hy = zeros(x,y);
ez = zeros(x,y);

% hy(3,1) = 0;
% hy(1,3) = 42;
% ez(2,2) = 1;
% ez(4,2) = 2;
%%

for t=1:10000
    t
    if(t<100)
        ez(size/2,size/2) = cos(2*pi*f*deltat*(t-1));
%     else
%         ez(size/2,size/2) = 0;
    end
    for x=1:size-2
        for y=1:size-2
            hx(x,y)=c1*hx(x,y)-c2*(ez(x,y+1)-ez(x,y));
            hy(x,y)=c1*hy(x,y)+c3*(ez(x+1,y)-ez(x,y));
        end
    end
    for x=2:size-1
        for y=2:size-1
%             if(x~=size/2 && y~=size/2)
                ez(x,y)=c4*ez(x,y)+c5*((deltat/(eps*deltax))*(hy(x,y)-hy(x-1,y)))-c5*((deltat/(eps*deltay))*(hx(x,y)-hx(x,y-1)));
%             end
        end
    end
%     imagesc([1:size]*deltax,[1:size]*deltay,ez)
    surf([1:size]*deltax,[1:size]*deltay,ez)
%     view(2);
    drawnow
%     hold on;
end





% 
%     ez(2,2,t) = cos(2*pi*f*deltat*(t-2));
%     for y=2:2:size-2
%         hx(y)= c1*hx(y) + c2*(ez(x)-ez(x,t))/(u*deltax);
%         for x=2:2:size-2
%             ez(x,t+2)= ez(x,t) + deltat*(hy(x+1,t+1)-hy(x-1,t+1))/(eps*deltax);
%         end
%     end
%     
% end
%     for x=4:2:size-2
%         ez(x,t+2)= ez(x,t) + deltat*(hy(x+1,t+1)-hy(x-1,t+1))/(eps*deltax);
%     end