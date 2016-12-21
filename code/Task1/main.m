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
% space samples
x = dx*(0:n_points-1);
% and the corresponding samples in momentum space
p = 2*pi/(n_points*dx)*(0:n_points/2);
% ---- Functions handles ----
Gaussian_Wave_Packet = @(x)1/(pi*d^2)^(1/4)*exp(-(x-x_0).^2/(2*d^2)).*exp(1i*p_0*(x-x_0)/hbar);
% Fourier transform obtained via Mathematica
Gaussian_Packet_Fourier = @(p)(exp(1i*p*x_0 - (d^2*(p_0 - p*hbar).^2)./(2*hbar^2))./((d^(-2))^(1/4)*pi^(1/4)));
% ---------------------------------

% Sample-discretize the wave packet function
wave_packet = Gaussian_Wave_Packet(x)*dx;
fourier_prob = abs(Gaussian_Packet_Fourier(p)).^2;

%prob = abs(wave_packet).^2;
%momentum_prob = abs(fft(wave_packet)).^2;
%momentum_prob = momentum_prob(1:floor(length(x)/2))*1/length(x);

temp = abs(fft(wave_packet))/n_points;
momentum_prob = temp(1:length(p));
% Compute double-sided contribution
momentum_prob(2:end-1) = 2*momentum_prob(2:end-1);
momentum_prob = momentum_prob.^2;

figure(1); clf;
plot(p, momentum_prob*1e5)
hold on
plot(p, fourier_prob)
hold off
xlabel('Momentum')
legend('fourier prob', 'theoretic prob')