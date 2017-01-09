% ===================================
% HOMEWORK 3B IN COMP.PHYS. - TASK 2
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
x_0         = 0;
dx          = 0.1;
n_points    = 2^10;
dp          = 2*pi/(n_points*dx);
dt          = 1;
v_0         = 0.08;
p_0         = sqrt(v_0*2*m);
alpha       = 2.0;


% ----------- VARIABLES ------------
% space samples
x = x_0+dx*((0:n_points-1)-n_points/2);
% and the corresponding samples in momentum space
p = ((0:n_points-1)-n_points/2)*dp;
% Functions handles
Gaussian_Wave_Packet = @(x)1/(pi*d^2)^(1/4)*exp(-(x-x_0).^2/(2*d^2)).*exp(1i*p_0*(x-x_0)/hbar);
Potential_Function = @(x) v_0*cosh(x/alpha).^(-2);
% ----
step_three=Gaussian_Wave_Packet(x);

potential = Potential_Function(x);
exp_potential = exp(-1i/hbar.*potential*dt);
inv_pot = exp(-1i/hbar*(hbar^2*p.^2./(2*m))*dt);

for j=1:n_points/4
    step_one = step_three;
    step_two = fftshift(fft(step_one.*exp_potential));
    step_three = ifft(ifftshift(inv_pot.*step_two));
    plot(x,abs(step_three).^2)
    pause(0.01)
end