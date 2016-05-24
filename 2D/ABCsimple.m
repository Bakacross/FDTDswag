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
size_x = size;
size_y = size;

% Void om/o = 0
c0 = deltat/(u_0*deltax);
cezh = deltat/(eps_0*deltax);

% ABC (same because void)
temp = sqrt(c0*cezh);
ABCcoeff = (temp-1)/(temp+1);
ezOldLeft = zeros(size_x,1);
ezOldRight = zeros(size_x,1);
ezOldTop = zeros(1,size_y);
ezOldBot = zeros(1,size_x);

% Initial condition
hx = zeros(size_x,size_y);
hy = zeros(size_x,size_y);
ez = zeros(size_x,size_y);

for t=1:10000
    
    t;
    
    for x=1:size_x
        for y=1:size_y
            if(y ~= size_y)
                hx(x,y)=hx(x,y)-c0*(ez(x,y+1)-ez(x,y));
            end
            if( x~= size_x)
                hy(x,y)=hy(x,y)+c0*(ez(x+1,y)-ez(x,y));
            end
        end
    end
    for x=2:size_x-1
        for y=2:size_y-1
            ez(x,y)=ez(x,y)+cezh*((hy(x,y)-hy(x-1,y))-(hx(x,y)-hx(x,y-1)));  
        end
    end
    ez(size_x/2,size_y/2) = ez(size_x/2,size_y/2)+ cos(2*pi*f*deltat*(t-1));

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
    set(gca,'XDir','normal')
    set(gca,'YDir','normal')
    xlabel('x');
    ylabel('y');
    caxis([-0.1, 0.1])
    drawnow

end


