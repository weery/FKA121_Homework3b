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
p_0         = sqrt(0.1*2*m);
dx          = 0.1;
n_points    = 2^10;
dp          = 2*pi/(n_points*dx);
dt          = 1;
v_0         = 0.10;
alpha       = 0.5;


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
bar = waitbar(0,'Calculating...');
for j=1:n_points/2
    step_one = step_three;
    step_two = fftshift(fft(step_one.*exp_potential));
    step_three = ifft(ifftshift(inv_pot.*step_two));
    waitbar(j/(n_points/2))
end

final_prob = abs(step_three).^2;

plot(x, final_prob)
close(bar)
pause(0.01)

% -- Compute the integrals --
left_integral = dx*sum(final_prob(1:n_points/2-1))
right_integral = dx*sum(final_prob(n_points/2+1:end))


xlabel('Position / [\AA]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Probability distribution', 'fontsize', 14)
s=sprintf('$\\left| \\psi (t = %i \\; \\mathrm{fs}) \\right|^2$',j*dt);
title(s, 'interpreter', 'latex', 'fontsize', 18)