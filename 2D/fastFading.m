clc
clear
close all
%------------------------
% Constant void
%------------------------
u_0 = 1.25663706*1e-6;
eps_0 = 8.8541878176e-12;
c = 299792458;
f = 2.45*1e9;
lambda = c/f;

%------------------------
% Parameters
%------------------------
% Simulation
deltax = lambda/10;
deltay = lambda/10;
deltat = deltax/c/2;
% size = 100;
size_x = 502;
size_y = 502;

% Copper
eps_r = 1;
u_r = 1;
sigma = 5.8e7;
sigma_m = 0;

%------------------------
% Initial condition
%------------------------
hx = zeros(size_x,size_y);
hy = zeros(size_x,size_y);
ez = zeros(size_x,size_y);

%------------------------
% Grid
%------------------------
Chxh = ones(size_x,size_y);
Chxe = zeros(size_x,size_y);
Chyh = ones(size_x,size_y);
Chye = zeros(size_x,size_y);
Ceze = ones(size_x,size_y);
Cezh = zeros(size_x,size_y);

% void
temp1 = deltat/(u_0*deltax); 
temp2 = deltat/(eps_0*deltax);
% if not void, divide temp1/2 by u_r/eps_r and use temp3 and temp4
eps = eps_0*eps_r;
u = u_0*u_r;
temp3 = sigma_m*deltat/(2*u);
temp4 = sigma*deltat/(2*eps);

Chxe(:,:) = temp1;
Chye(:,:) = temp1;
Cezh(:,:) = temp2;

% Walls in copper
sigma_copper = 5.8e7;
temp_copper = sigma_copper*deltat/(2*eps_0);
Ceze(2,:) = (1-temp_copper)/(1+temp_copper);
Ceze(501,:) = (1-temp_copper)/(1+temp_copper);
Ceze(:,2) = (1-temp_copper)/(1+temp_copper);
Ceze(:,501) = (1-temp_copper)/(1+temp_copper);
Cezh(2,:) = temp2/temp_copper;
Cezh(501,:) = temp2/temp_copper;
Cezh(:,2) = temp2/temp_copper;
Cezh(:,501) = temp2/temp_copper;


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

%-----------------------
% Zone
%-----------------------

zone1 = ones(1,100);
z1 = 1;
zone2 = ones(1,100);
z2 = 1;
zone3 = ones(1,100);
z3 = 1;
zone4 = ones(1,1500);
z4 = 1;

for t=1:5000
    
    % Update Hx, Hy and Ez
    hx(:,1:end-1) = Chxh(:,1:end-1).*hx(:,1:end-1)-Chxe(:,1:end-1).*(ez(:,2:end)-ez(:,1:end-1));
    hy(1:end-1,:) = Chyh(1:end-1,:).*hy(1:end-1,:)+Chye(1:end-1,:).*(ez(2:end,:)-ez(1:end-1,:));
    ez(2:end-1,2:end-1) = Ceze(2:end-1,2:end-1).*ez(2:end-1,2:end-1)+Cezh(2:end-1,2:end-1).*((hy(2:end-1,2:end-1)-hy(1:end-2,2:end-1))-(hx(2:end-1,2:end-1)-hx(2:end-1,1:end-2)));

    % Hardwire source
    ez(size_x/2,size_y/2) = 3*cos(2*pi*f*deltat*(t-1));
    
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
    
    % Plot
    imagesc([1:size_x]*deltax,[1:size_y]*deltay,ez')
    title(['Time step :', num2str(t)])
    axis equal
    colormap(jet(50));
    set(gca,'XDir','normal')
    set(gca,'YDir','normal')
    xlabel('x');
    ylabel('y');
    caxis([-1, 1])
    
    % Study Zone
    if(t>=150 && t < 250)
        zone1(z1) = ez(250,200);
        hold on
        plot(250*deltax,200*deltay,'ko','LineWidth',3,'MarkerSize',5);
        hold off
        z1 = z1+1;
    elseif(t>=450 && t < 550)
        zone2(z2) = ez(250,50);
        hold on
        plot(250*deltax,50*deltay,'ko','MarkerSize',10);
        hold off
        z2 = z2+1;
    elseif(t>=650 && t < 750)
        zone3(z3) = ez(50,50);
        hold on
        plot(50*deltax,50*deltay,'ko');
        hold off
        z3 = z3+1;
    elseif(t>= 3000 && t < 4500)
        zone4(z4) = ez(350,50);
        hold on
        plot(350*deltax,50*deltay,'ko');
        hold off
        z4 = z4+1;
    end
    
    drawnow

end


