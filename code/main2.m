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
x = x_0+dx*((0:n_points-1)-n_points/2);
% and the corresponding samples in momentum space
p = ((0:n_points-1)-n_points/2)*dp;
% Functions handles
Gaussian_Wave_Packet = @(x)1/(pi*d^2)^(1/4)*exp(-(x-x_0).^2/(2*d^2)).*exp(1i*p_0*(x-x_0)/hbar);
Potential_Function = @(x) x*0;
Gaussian_Wave_Time_Evolution=@(x,t) (pi^(1/2)*(d+1i*hbar*t/(m*d))).^(-1/2).*...
    exp(-(x-p_0*t/m).^2/(2*d^2*(1+1i*hbar*t/(m*d^2)))).*...
    exp(i*p_0/hbar*(x-p_0*t/(2*m)));
Gaussian_Time_Probability=@(x,t) 1/(pi^(1/2)*(d^2+hbar^2*t^2/(m^2*d^2))^(1/2))*exp(-(x-p_0/m*t).^2/(d^2+hbar^2*t^2/(m^2*d^2))); 
% ----
step_three=Gaussian_Wave_Packet(x);

potential = Potential_Function(x);
exp_potential = exp(-1i/hbar.*potential*dt);
inv_pot = exp(-1i/hbar*(hbar^2*p.^2./(2*m))*dt);


% Plot initial prob.distr.
figure(1); clf;
plotX=1:n_points;
plot(x(plotX), abs(step_three(plotX).^2))
xlim([0 2])
xlabel('Position / [\AA]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Probability distribution', 'fontsize', 14)
title('$\left| \psi (0) \right|^2$', 'interpreter', 'latex', 'fontsize', 18)

% Plot the rest in a separate figure

figure(2); clf;
plotHandle_1 = plot(x(plotX), abs(step_three(plotX).^2));
hold on
plotHandle_2 = plot(x(plotX), abs(Gaussian_Wave_Time_Evolution(plotX,0)).^2,'--');
hold off

xlabel('Position / [\AA]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Probability distribution', 'fontsize', 14)
s=sprintf('$\\left| \\psi (t = %i \\; \\mathrm{fs}) \\right|^2$',j*dt);
title(s, 'interpreter', 'latex', 'fontsize', 18)

bar= waitbar(0,'Calculating...')

for j=1:n_points/2
    step_one = step_three;    
    step_two = fftshift(fft(step_one.*exp_potential));
    step_three = ifft(ifftshift(inv_pot.*step_two));
    set(plotHandle_1, 'YData', abs(step_three(plotX)).^2)
    set(plotHandle_2, 'YData', Gaussian_Time_Probability(x(plotX),j*dt))
    waitbar(j/(n_points/2))
    pause(0.01)
end

close(bar)

s=sprintf('$\\left| \\psi (t = %i \\; \\mathrm{fs}) \\right|^2$',j*dt);
title(s, 'interpreter', 'latex', 'fontsize', 18)


figure(3); clf;
plotHandle_1 = plot(p(plotX)/p_0, abs(fftshift(fft(step_three(plotX)))).^2);
hold on
plotHandle_2 = plot(p(plotX)/p_0, abs(fftshift(fft(Gaussian_Wave_Time_Evolution(x(plotX),j*dt)))).^2,'--');
hold off
xlabel('P/p_0 / [\AA]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Probability distribution', 'fontsize', 14)
s=sprintf('$\\left| \\psi (t = %i \\; \\mathrm{fs}) \\right|^2$',j*dt);
title(s, 'interpreter', 'latex', 'fontsize', 18)
