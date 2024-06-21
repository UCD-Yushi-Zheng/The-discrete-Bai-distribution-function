%------------------------------%
close all;
clear;
clc;
%---------Initialize the test signal------------%
N=101; % Number of samples
x0=-(N-1)/2:1:(N-1)/2; % initialize an array for the corrdinate
L=10; % size of the coordinate 
T=L/N; % Sampling frequency 
x=x0*T; % coordinate in the space domain 
f=(x0/(N*T*2)); % coordinate in the frequency domain 

y1=exp(-(x + 2).^2/sqrt(2));
y2=exp(-(x - 2).^2/sqrt(2));
y=y1+y2; % Test signal 

figure(1); % Plot the test signal 
plot(x,abs(y),'linewidth',2);
colorbar
set(gca,'FontSize',15)
set(gcf,'unit','centimeters','position',[3 5 12 10])
xlabel('x');
ylabel('A');
title('Gaussian signal');


%---------Calculate the test signal's BDF-------------%
a=2.5; % Initialize the LCT parameter 
b=10;
d=1;

[BDFA,aA]=BDF(y,T,a,b,d); % Calculate the test signal's BDF
BDFA=BDFA*T*2;

figure % Plot the instantaneous auto-correlation function
subplot(2,2,1)
contourf(x,2*x,abs(aA));
colorbar
set(gca,'FontSize',15)
xlabel('x');
ylabel('tau');
title('IAF');

subplot(2,2,2) % Plot the BDF
contourf(x,f*b,abs(BDFA));
colorbar
set(gca,'FontSize',15)
xlabel('x');
ylabel('u');
title('Simulated BDF');


%--------Compare the discrete BDF and its analytical solution------------%
[t,~]=meshgrid(x); % initialize the coordinate 
[~,u]=meshgrid(f*b);

GT=exp (1).^((-0.565685E1).*t+(-0.141421E1).*t.^2+((-0.470369E-1)+ ... % Analytical solution for the BDF when, a=2.5, b=10 and d=1 calculated using Mathematica.
    sqrt(-1)*(-0.10449E0)).*u.^2).*((0.560305E-2+sqrt(-1)*0.362233E-2) ...
    .*exp (1).^((0.E-323+sqrt(-1)*0.314159E0).*u.^2)+((-0.155815E-1)+ ...
    sqrt(-1)*0.753217E-2).*exp (1).^(0.565685E1.*t+(((-0.940738E0)+ ...
    sqrt(-1)*0.423481E0)+(0.E-323+sqrt(-1)*0.314159E0).*u).*u)+(( ...
    -0.155815E-1)+sqrt(-1)*0.753217E-2).*exp (1).^(0.565685E1.*t+(( ...
    0.940738E0+sqrt(-1)*(-0.423481E0))+(0.E-323+sqrt(-1)*0.314159E0).* ...
    u).*u)+(0.560305E-2+sqrt(-1)*0.362233E-2).*exp (1).^(0.113137E2.*t+ ...
    (0.E-323+sqrt(-1)*0.314159E0).*u.^2));

subplot(2,2,3) % Plot the analytical BDF
contourf(x,f*b,abs(GT));
colorbar
set(gca,'FontSize',15)
xlabel('x');
ylabel('u');
title('Analytical BDF');

subplot(2,2,4) % Plot the error 
contourf(x,f*b,abs(BDFA)-abs(GT));
colorbar
set(gca,'FontSize',15)
xlabel('x');
ylabel('u');
title('Error');

mean_squared_error=immse(BDFA,GT) % Calculate the mean-squared error
