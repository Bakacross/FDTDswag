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
deltax = lambda/50;
deltay = lambda/50;
deltat = deltax/c/2;
% size = 100;
size_x =200;
size_y = 100;

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

% Young
dist = 20;
line1 = 1:size_x/2-dist;
line2 = size_x/2-dist+2:size_x/2+dist-2;
line3 = size_x/2+dist:size_x;
y_young = size_y/2-10;
Cezh([line1 line2 line3],y_young) = 0;
  

%------------------------
% ABC (same because void)
%------------------------
temp1 = sqrt(Cezh(2,2)*Chye(1,1));
temp2 = 1/temp1 + 2 + temp1;
coef0 = -(1 / temp1 - 2 + temp1) / temp2; 
coef1 = -2 * (temp1 - 1 / temp1) / temp2; 
coef2 = 4 * (temp1 + 1 / temp1) / temp2;
ezOldLeft = zeros(size_x,3,2); %1) all column 2) which column 3) n-1 or n-2
ezOldRight = zeros(size_x,3,2);
ezOldTop = zeros(3,size_y,2);
ezOldBot = zeros(3,size_y,2);

for t=1:10000
    
    % Update Hx, Hy and Ez
    hx(:,1:end-1) = Chxh(:,1:end-1).*hx(:,1:end-1)-Chxe(:,1:end-1).*(ez(:,2:end)-ez(:,1:end-1));
    hy(1:end-1,:) = Chyh(1:end-1,:).*hy(1:end-1,:)+Chye(1:end-1,:).*(ez(2:end,:)-ez(1:end-1,:));
    ez(2:end-1,2:end-1) = Ceze(2:end-1,2:end-1).*ez(2:end-1,2:end-1)+Cezh(2:end-1,2:end-1).*((hy(2:end-1,2:end-1)-hy(1:end-2,2:end-1))-(hx(2:end-1,2:end-1)-hx(2:end-1,1:end-2)));
    
    % Hardwire source
    ez(size_x/2,2) = 5*cos(2*pi*f*deltat*(t-1));
    
    % ABC Left
    ez(:,1) = coef0*(ez(:,3) + ezOldLeft(:,1,2)) + coef1*(ezOldLeft(:,1,1) + ezOldLeft(:,3,1) - ez(:,2) - ezOldLeft(:,2,2)) + coef2*ezOldLeft(:,2,1) - ezOldLeft(:,3,2);
    ezOldLeft(:,:,2) = ezOldLeft(:,:,1);
    ezOldLeft(:,:,1) = ez(:,1:3);
    
    % ABC Right
    ez(:,size_y) = coef0*(ez(:,size_y-2) + ezOldRight(:,1,2)) + coef1*(ezOldRight(:,1,1) + ezOldRight(:,3,1) - ez(:,size_y-1) - ezOldRight(:,2,2)) + coef2*ezOldRight(:,2,1) - ezOldRight(:,3,2);
    ezOldRight(:,:,2) = ezOldRight(:,:,1);
    ezOldRight(:,:,1) = ez(:,size_y:-1:size_y-2);
    
    % ABC Top
    ez(1,:) = coef0*(ez(3,:) + ezOldTop(1,:,2)) + coef1*(ezOldTop(1,:,1) + ezOldTop(3,:,1) - ez(2,:) - ezOldTop(2,:,2)) + coef2*ezOldTop(2,:,1) - ezOldTop(3,:,2);
    ezOldTop(:,:,2) = ezOldTop(:,:,1);
    ezOldTop(:,:,1) = ez(1:3,:);
    
    % ABC Bot
    ez(size_x,:) = coef0*(ez(size_x-2,:) + ezOldBot(1,:,2)) + coef1*(ezOldBot(1,:,1) + ezOldBot(3,:,1) - ez(size_x-1,:) - ezOldBot(2,:,2)) + coef2*ezOldBot(2,:,1) - ezOldBot(3,:,2);
    ezOldBot(:,:,2) = ezOldBot(:,:,1);
    ezOldBot(:,:,1) = ez(size_x:-1:size_x-2,:);
  
    imagesc([1:size_x]*deltax,[1:size_y]*deltay,ez')
    title(['Time elapsed :', num2str((t-1)*deltat)])
    axis equal
%     colorbar
    colormap(jet);
    set(gca,'XDir','normal')
    set(gca,'YDir','normal')
    xlabel('x');
    ylabel('y');
%     xlim([size_x/2-15 size_x/2+15]*deltax)
    caxis([-0.01, 0.01])
    hold on
    plot(line1*deltax,y_young*ones(1,length(line1))*deltay,'k','LineWidth',2)
    plot(line2*deltax,y_young*ones(1,length(line2))*deltay,'k','LineWidth',2)
    plot(line3*deltax,y_young*ones(1,length(line3))*deltay,'k','LineWidth',2)
    hold off
    
    drawnow

end


