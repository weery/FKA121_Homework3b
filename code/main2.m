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
p_0         = sqrt(0.1*2*m);
x_0         = 0;
dx          = 0.1;
n_points    = 2^10;
dp          = 2*pi/(n_points*dx);
dt          = 1;


% ----------- VARIABLES ------------
% space samples
x = x_0+dx*(0:n_points-1);
% and the corresponding samples in momentum space
p = ((0:n_points-1)-n_points/2)*dp;
% Functions handles
Gaussian_Wave_Packet = @(x)1/(pi*d^2)^(1/4)*exp(-(x-x_0).^2/(2*d^2)).*exp(1i*p_0*(x-x_0)/hbar);
Potential_Function = @(x) 0;
% ----
step_three=Gaussian_Wave_Packet(x);

potential = Potential_Function(x);
exp_potential = exp(-1i/hbar.*potential*dt);
inv_pot = exp(-1i/hbar*(hbar^2*p.^2./(2*m))*dt);


% Plot initial prob.distr.
figure(1); clf;
plot(x(1:n_points/2), abs(step_three(1:n_points/2).^2))
xlim([0 2])
xlabel('Position / [\AA]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Probability distribution', 'fontsize', 14)
title('$\left| \psi (0) \right|^2$', 'interpreter', 'latex', 'fontsize', 18)

% Plot the rest in a separate figure
figure(2); clf;
step_one = step_three;
step_two = fftshift(fft(step_one.*exp_potential));
step_three = ifft(ifftshift(inv_pot.*step_two));

plotHandle = plot(x(1:n_points/2), abs(step_three(1:n_points/2).^2));
xlabel('Position / [\AA]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Probability distribution', 'fontsize', 14)
title('$\left| \psi (t = 256 \mathrm{fs}) \right|^2$', 'interpreter', 'latex', 'fontsize', 18)

for j=1:n_points/4-1
    step_one = step_three;
    
    step_two = fftshift(fft(step_one.*exp_potential));
    step_three = ifft(ifftshift(inv_pot.*step_two));
    set(plotHandle, 'YData', abs(step_three(1:n_points/2).^2))
    pause(0.01)
end

