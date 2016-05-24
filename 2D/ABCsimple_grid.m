clc
clear
close all
%------------------------
% Constant void
%------------------------
u_0 = 1.25663706*1e-6;
eps_0 = 8.8541878176e-12;
c = 299792458;
f = 2.4*1e9;
lambda = c/f;

%------------------------
% Parameters
%------------------------
% Simulation
deltax = lambda/10;
deltay = lambda/10;
deltat = deltax/c/10;
size = 50;
size_x = size;
size_y = size;

% Brain
sigma = 1.3;
sigma_m = 0;
eps_r = 43;
u_r = 1;

%------------------------
% Initial condition
%------------------------
hx = zeros(size_x,size_y);
hy = zeros(size_x,size_y);
ez = zeros(size_x,size_y);

%------------------------
% Grid
%------------------------
Chxh = zeros(size_x,size_y);
Chxe = zeros(size_x,size_y);
Chyh = zeros(size_x,size_y);
Chye = zeros(size_x,size_y);
Ceze = zeros(size_x,size_y);
Cezh = zeros(size_x,size_y);

% void
temp1 = deltat/(u_0*deltax); 
temp2 = deltat/(eps_0*deltax);
% if not void, divide temp1/2 by u_r/eps_r and use temp3 and temp4
eps = eps_0*eps_r;
u = u_0*u_r;
temp3 = sigma_m*deltat/(2*u);
temp4 = sigma*deltat/(2*eps);

% test
eps_r = 9;

for x=1:size_x
    for y=1:size_y
        if(y ~= size_y)
            Chxh(x,y) = 1;
            Chxe(x,y) = temp1;
        end
        if( x~= size_x)
            Chyh(x,y) = 1;
            Chye(x,y) = temp1;
        end
    end
end

for x=2:size_x-1
    for y=2:size_y-1
%         if(x > 10 && x < 40 && y > 10 && y < 40)
            Ceze(x,y) = 1;
            Cezh(x,y) = temp2;
%         else
%             Ceze(x,y) = 1;
%             Cezh(x,y) = temp2/eps_r;
%         end
        
    end
end
  

%------------------------
% ABC (same because void)
%------------------------
tempABC = sqrt(temp1*temp2);
ABCcoeff = (tempABC-1)/(tempABC+1);
ezOldLeft = zeros(size_x,1);
ezOldRight = zeros(size_x,1);
ezOldTop = zeros(1,size_y);
ezOldBot = zeros(1,size_x);



for t=1:10000
    
    for x=1:size_x
        for y=1:size_y
            if(y ~= size_y)
                hx(x,y)=Chxh(x,y)*hx(x,y)-Chxe(x,y)*(ez(x,y+1)-ez(x,y));
            end
            if( x~= size_x)
                hy(x,y)=Chyh(x,y)*hy(x,y)+Chye(x,y)*(ez(x+1,y)-ez(x,y));
            end
        end
    end
    for x=2:size_x-1
        for y=2:size_y-1
            ez(x,y)=Ceze(x,y)*ez(x,y)+Cezh(x,y)*((hy(x,y)-hy(x-1,y))-(hx(x,y)-hx(x,y-1)));  
        end
    end
    ez(size_x/2,size_y/2) = 3*cos(2*pi*f*deltat*(t-1));

    ez(:,1) = ezOldLeft + ABCcoeff*(ez(:,2)-ez(:,1));
    ezOldLeft = ez(:,2);
    ez(:,size_y) = ezOldRight + ABCcoeff*(ez(:,size_y-1)-ez(:,size_y));
    ezOldRight = ez(:,size_y-1);
    ez(1,:) = ezOldTop + ABCcoeff*(ez(2,:)-ez(1,:));
    ezOldTop = ez(2,:);
    ez(size_x,:) = ezOldBot + ABCcoeff*(ez(size_x-1,:)-ez(size_x,:));
    ezOldBot = ez(size_x-1,:);
    

    imagesc([1:size_x]*deltax,[1:size_y]*deltay,ez')
    title(['Time elapsed :', num2str((t-1)*deltat)])
    axis equal
%     colorbar
    colormap(jet);
    set(gca,'XDir','normal')
    set(gca,'YDir','normal')
    xlabel('x');
    ylabel('y');
    caxis([-1, 1])
    drawnow

end


