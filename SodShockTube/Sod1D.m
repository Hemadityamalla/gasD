clear;clc;format long;

%Initial data (used SI units- kg, m, sec)
pL = 1.0; % Pa
rhoL = 1.0; %kg/m3
uL = 0.0; %m/sec


pR = 0.1;
rhoR = 0.125;
uR = 0.0;

gamma = 1.4;
aL = sqrt(gamma*(pL/rhoL));
aR = sqrt(gamma*(pR/rhoR));

xL = -10.0; xR = 10.0;
x = linspace(xL,xR,500);
x0 = 0.0;
for t = 0:0.1:3

MachShock = @(x) x - (1.0/x) - aL*((gamma+1)/(gamma-1))*(1.0 - ((pR/pL)*((2*gamma*x^2)/(gamma+1) - (gamma-1)/(gamma+1)))^((gamma-1)/(2*gamma)));


%Solving for the Mach shock number

Ms = fzero(MachShock, 10);

%Computing the data in region 1

p1 = pR*((2.0*gamma*Ms^2)/(gamma+1) - (gamma-1)/(gamma+1));
u1 = (2/(gamma+1))*(Ms - 1.0/Ms);
rho1 = rhoR/((2.0)/((gamma+1)*Ms^2) + (gamma-1)/(gamma+1));

%Data in region 2

p2 = p1;
u2 = u1;
rho2 = rhoL*(p2/pL)^(1.0/gamma);
a2 = sqrt(gamma*(p2/rho2));

%Data in region E
x1 = x0 - aL*t;
x2 = x0 + (u2 - a2)*t;

xELog = ((x1 < x) & (x < x2));
uE = (2.0/(gamma+1))*(aL*xELog + (x.*(xELog) - x0*xELog)/t);
aE = aL*xELog - (gamma-1)*uE*0.5;
pE = pL*(aE/aL).^((2.0*gamma)/(gamma-1));
rhoE = (gamma*pE)./aE.^2;
rhoE(isnan(rhoE)) = 0.0; %problems due to vectorization

%Other positiions
x3 = x0 + u2*t;
x4 = x0 + Ms*t; %Position of shock(!)


Pressure = (x < x1)*pL + ((x1 < x) & (x < x2)).*pE + ((x2 < x) & (x < x4))*p2 + (x > x4)*pR; 
density = (x < x1)*rhoL + ((x1 < x) & (x < x2)).*rhoE + ((x2 < x) & (x < x3))*rho2 + ((x3 < x) & (x < x4))*rho1 +  (x > x4)*rhoR;
velocity = (x < x1)*uL + ((x1 < x) & (x < x2)).*uE + ((x2 < x) & (x < x4))*u2 + (x > x4)*uR;
energy = (x < x1)*uL + ((x1 < x) & (x < x2)).*uE + ((x2 < x) & (x < x4))*u2 + (x > x4)*uR;

figure(1)
subplot(2,2,1),
plot(x,density,'-b','LineWidth',2);
xlabel('x (m)');
ylabel('Density (kg/m^3)');
title('Plot of Density vs Position');
grid on;
subplot(2,2,2),
plot(x,Pressure,'-g','LineWidth',2);
xlabel('x (m)');
ylabel('Pressure (Pa)');
title('Plot of Pressure vs Position');
grid on;
subplot(2,2,3),
plot(x,velocity,'-r','LineWidth',2);
xlabel('x (m)');
ylabel('Velocity (m/s)');
title('Plot of Velocity vs Position');
grid on;
subplot(2,2,4),
plot(x,energy,'-k','LineWidth',2);
xlabel('x (m)');
ylabel('Specific Internal Energy (J/kg)');
title('Plot of Internal Energy vs Position');
grid on;
pause(0.1);
end

