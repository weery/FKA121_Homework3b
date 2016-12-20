% ===================================
% HOMEWORK 3B IN COMP.PHYS. - TASK 1
% ===================================
% By Victor Nilsson, Simon Nilsson
% 2015
%
% Length scale: 1 Å
% Time scale:   1 fs = 1e-15 s
% Energy scale: 1 eV

clear all, clc, close all

% ------ SIMULATION PARAMETERS ---------
hbar        = 1.054/1.602; % JS -> f eV s
d           = 0.5;
m           = 1.66/1.6*1e2;
p_0         = sqrt(0.1*2*m);
x_0         = 0;
dx          = 0.01;
n_points    = 1024;


% ----------- VARIABLES ------------
x = dx*(0:n_points-1);
p = 1/(n_points*dx)*(0:n_points/2);
% Functions handles
Gaussian_Wave_Packet = @(x)1/(pi*d^2)^(1/4)*exp(-(x-x_0).^2/(2*d^2)).*exp(1i*p_0*(x-x_0)/hbar);
Gaussian_Fourier_Packet = @(p)(exp(1i*p*x_0 - (d^2*(p_0 + p*hbar).^2)./(2*hbar^2))./((d^(-2))^(1/4)*pi^(1/4)));
% ----

gauss = Gaussian_Wave_Packet(x);
fourier = abs(Gaussian_Fourier_Packet(p)).^2;
prob = abs(gauss).^2;
momentum = abs(fft(gauss)).^2;
momentum = momentum(1:floor(length(x)/2))*1/length(x);

temp = abs(fft(prob))/n_points;
momentum = temp(1:n_points/2+1);
momentum(2:end-1) = 2*momentum(2:end-1);
momentum = momentum.^2;

figure(1); clf;
%plot(prob)
plot(p, momentum)
hold on
%figure(2); clf;
plot(p, fourier.*1e4)
hold off