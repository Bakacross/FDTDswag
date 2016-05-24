

clc
clear
close all

% Constant (void)
u_0 = 1.25663706*1e-6;
eps_0 = 8.8541878176e-12;
c = 299792458;
f = 2.4*1e9;
lambda = c/f;

% Option
deltax = lambda/10;
deltay = lambda/10;
deltat = deltax/c/10;
size = 50;
loss_layer = 2;
size_x = size+loss_layer;
size_y = size+loss_layer;

% Loss layer
loss = 0.01;
eps_r = 4;
eps = eps_0*eps_r;
o = 0.1;
u_r = 1;
u = u_0*u_r;
om = 0;
c1l = (1-(om*deltat)/(2*u))/(1+(om*deltat)/(2*u));
c4l = (1-(om*deltat)/(2*eps))/(1+(om*deltat)/(2*eps));
cezh_l = deltat/(eps*deltax);
c2l = deltat/(u*deltay)/(1+(om*deltat)/(2*u));
c5l = cezh_l/(1+(o*deltat/(2*eps)));

% Head
% brain : eps_r = 43 sigma = 1,3
eps_r = 43;
eps = eps_0*eps_r;
o = 1.3; %%% corriger bs fouad becuz it is bs hihihi

u_r = 1;
u = u_0*u_r;
om = 0;

pos_x = size_x/2 - 10; % real x is pos_x*deltax
pos_y = size_y/2;
radius = 0.08/deltax; % 16 cm diameter
c1 = (1-(om*deltat)/(2*u))/(1+(om*deltat)/(2*u));

% c3 = deltat/(u*deltax)/(1+(om*deltat)/(2*u));
c4 = (1-(om*deltat)/(2*eps))/(1+(om*deltat)/(2*eps));
cezh_ = deltat/(eps*deltax);
c2 = deltat/(u*deltay)/(1+(om*deltat)/(2*u));
c5 = cezh_/(1+(o*deltat/(2*eps)));



% draw head
c_t = linspace(0,2*pi);
% figure
% plot(radius*cos(c_t)+pos_x,radius*sin(c_t)+pos_y,'k*')

% Void om/o = 0
c0 = deltat/(u*deltax);
cezh = deltat/(eps_0*deltax);
% cez2 = deltat/(eps_0*deltay);

% Initial condition
hx = zeros(size_x,size_y);
hy = zeros(size_x,size_y);
ez = zeros(size_x,size_y);

% Second order ABC
ezLeft = zeros(3,2,size_y); % space - time - all the left column
ezRight = zeros(3,2,size_y);
ezTop = zeros(3,2,size_x);
ezBottom = zeros(3,2,size_x);
temp1 = c*deltat/deltax; % void
temp2 = 1/temp1 + 2 + temp1;
coef0 = -(1 / temp1 - 2 + temp1) / temp2; 
coef1 = -2 * (temp1 - 1 / temp1) / temp2; 
coef2 = 4 * (temp1 + 1 / temp1) / temp2;
ezold1 = zeros(size_y,1);
ezold2 = zeros(size_y,2);

ezold1Right = zeros(size_y,1);
ezold2Right = zeros(size_y,2);

ezold1Bottom = zeros(size_x,1);
ezold2Bottom = zeros(size_x,2);

ezold1Top = zeros(size_x,1);
ezold2Top= zeros(size_x,2);

for t=1:10000
    
    t;
    %if(t<100)
    % additive source needed
%     ez(size/2,size/2) = cos(2*pi*f*deltat*(t-1));    
   % end
   
   % ABC left
   for i = 1:size_y
       ez(1,i) = 2*ezold1(i)-ezold2(i,2);
       ezold1(i,1) = ez(2,i);
       ezold2(i,2) = ezold2(i,1);
       ezold2(i,1) = ez(3,i);
       
       ez(size_x,i) = 2*ezold1Right(i)-ezold2Right(i,2);
       ezold1Right(i) = ez(size_x-1,i);
       ezold2Right(i,2) = ezold2Right(i,1);
       ezold2Right(i,1) = ez(size_x-2,i);
%        ez(1,i) = coef0 * ez(3,i) + ezLeft(1,2,i)... 
%        + coef1 * (ezLeft(1,1,i) + ezLeft(3,1,i) - ez(2,i) - ezLeft(2,2,i))...
%        + coef2 * ezLeft(2,1,i) - ezLeft(3,2,i);
%        for j = 1:3
%            ezLeft(j,2,i) = ezLeft(j,1,i);
%            ezLeft(j,1,i) = ez(j,i);
%        end
   end
   
   for i = 1:size_x
       ez(i,1) = 2*ezold1Bottom(i)-ezold2Bottom(i,2);
       ezold1Bottom(i) = ez(i,2);
       ezold2Bottom(i,2) = ezold2Bottom(i,1);
       ezold2Bottom(i,1) = ez(i,3);
              
       ez(i,size_x) = 2*ezold1Top(i)-ezold2Top(i,2);
       ezold1Top(i) = ez(i,size_x-1);
       ezold2Top(i,2) = ezold2Top(i,1);
       ezold2Top(i,1) = ez(i,size_x-2);
   end
   
   
    for x=1:size_x
        for y=1:size_y
            if((x-pos_x)^2 + (y-pos_y)^2<radius^2) 
                if(y ~= size_y)
                    hx(x,y)=c1*hx(x,y)-c2*(ez(x,y+1)-ez(x,y));
                end
                if( x~= size_x)
                    hy(x,y)=c1*hy(x,y)+c2*(ez(x+1,y)-ez(x,y));
                end
            elseif( x > size || y > size || x < loss_layer || y < loss_layer) 
                if(y ~= size_y)
                    hx(x,y)=c1l*hx(x,y)-c2l*(ez(x,y+1)-ez(x,y));
                end
                if( x~= size_x)
                    hy(x,y)=c1l*hy(x,y)+c2l*(ez(x+1,y)-ez(x,y));
                end
            else
                if(y ~= size_y)
                    hx(x,y)=hx(x,y)-c0*(ez(x,y+1)-ez(x,y));
                end
                if( x~= size_x)
                    hy(x,y)=hy(x,y)+c0*(ez(x+1,y)-ez(x,y));
                end
            end
        end
    end
    for x=2:size_x-1
        for y=2:size_y-1
            if((x-pos_x)^2 + (y-pos_y)^2<radius^2)
                ez(x,y)=c4*ez(x,y)+c5*((hy(x,y)-hy(x-1,y))-(hx(x,y)-hx(x,y-1)));
            elseif ( x > size || y > size || x < loss_layer || y < loss_layer) 
                ez(x,y)=c4l*ez(x,y)+c5l*((hy(x,y)-hy(x-1,y))-(hx(x,y)-hx(x,y-1)));
            else %if (x ~= size/2 && y ~= size/2)
                ez(x,y)=ez(x,y)+cezh*((hy(x,y)-hy(x-1,y))-(hx(x,y)-hx(x,y-1)));  
            end
        end
    end
    ez(size_x/2,size_y/2) = ez(size_x/2,size_y/2)+ cos(2*pi*f*deltat*(t-1));
%     1ez(10:30,30)=0;
%     ez(30,10:30)=0;
%     ez(10:30,20)=0;
%     ez(20,10:30)=0;
    
    imagesc([1:size_x]*deltax,[1:size_y]*deltay,ez')
    title(['Time elapsed :', num2str((t-1)*deltat)])
%     axis equal
%     colorbar
    set(gca,'XDir','normal')
    set(gca,'YDir','normal')
    xlabel('x');
    ylabel('y');
    caxis([-0.01, 0.01])
    hold on
    plot((radius*cos(c_t)+pos_x)*deltax,(radius*sin(c_t)+pos_y)*deltay,'k','LineWidth',6)
    hold off

%     surf([1:size]*deltax,[1:size]*deltay,ez)
%     zlim([-1 1]);
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