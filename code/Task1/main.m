% ===================================
% HOMEWORK 3B IN COMP.PHYS. - TASK 1
% ===================================
% By Victor Nilsson, Simon Nilsson

clear all, clc, close all

% ------ SIMULATION PARAMETERS ---------
hbar        = 1;
d           = 0.5;
m           = 1;
p_0         = sqrt(0.1*2*m);
x_0         = 0;
dx          = 1;
n_points    = 1e3;


% ----------- VARIABLES ------------
x = dx*(0:n_points-1);
p = 1/(n_points*dx)*(0:n_points-1);
% Functions handles
Gaussian_Wave_Packet = @(x)1/(pi*d^2)^(1/4)*exp(-(x-x_0).^2/(2*d^2)).*exp(1i*p_0*(x-x_0)/hbar);
Gaussian_Fourier_Packet = @(p)(exp(1i*p*x_0 - (d^2*(p_0 + p*hbar).^2)./(2*hbar^2))./((d^(-2))^(1/4)*pi^(1/4)));
% ----

gauss = Gaussian_Wave_Packet(x);
fourier = abs(Gaussian_Fourier_Packet(p)).^2;
prob = abs(gauss).^2;
momentum = abs(fft(prob)).^2;
momentum = momentum(1:floor(length(x)/2))*1/length(x);


figure(1); clf;
%plot(prob)
plot(p(1:floor(length(x)/2)), momentum)
hold on
plot(p, fourier)
hold off