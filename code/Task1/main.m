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


% ----------- VARIABLES ------------
x=0:0.01:5;
% Functions handles
Gaussian_Wave_Packet = @(x)1/(pi*d^2)^(1/4)*exp(-(x-x_0).^2/(2*d^2)).*exp(1i*p_0*(x-x_0)/hbar);
% ----

gauss = Gaussian_Wave_Packet(x);
prob = abs(gauss).^2;
momentum = abs(fft(prob)).^2;
momentum = momentum(1:floor(length(x)/2))*1/length(x);

plot(prob)
hold on
plot(momentum)