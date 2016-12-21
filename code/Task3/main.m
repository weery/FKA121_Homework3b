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
dx          = 0.01;
n_points    = 2^12;
x_0         = n_points/2*dx;
dt          = 0.0001;
dk          = 2*pi/(n_points*dx);
v_0         = 0.1;
a           = 0.5;

% ----------- VARIABLES ------------
x = dx*(1:n_points);
k = dk*((1:n_points)-n_points/2);
% Functions handles
Gaussian_Wave_Packet = @(x)1/(pi*d^2)^(1/4)*exp(-(x-x_0).^2/(2*d^2)).*exp(1i*p_0*(x-x_0)/hbar);
Potential_Function = @(x) v_0*cosh(x/a).^(-2);
% ----
step_three=Gaussian_Wave_Packet(x);
plot(abs(step_three).^2)

for j=1:100
    step_one = step_three;
    potential = Potential_Function(x);
    exp_potential = exp(-1i/hbar.*potential*dt);
    
    step_two=fft(step_one.*exp_potential);
    inv_pot = exp(-i/hbar*(hbar^2*k.^2./(2*(1:n_points)))*dt);

    step_three = ifft(inv_pot.*step_two);
    momentum = fft(step_three);
    figure(1)
    plot(x,abs(step_three).^2)
    pause(0.01)
    figure(2)
    plot(k,abs(momentum).^2)
    pause(0.01)
    
   
end

