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
n_points    = 2^12;
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
Theoretical_Width_Position=@(t)d/2^(1/2)*(1+hbar^2*t.^2/(m^2*d^4)).^(1/2);
Theoretical_Width_Momentum=@(t)1./(2*Theoretical_Width_Position(t));
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
plotHandle_1 = plot(x(plotX), abs(step_three(plotX)).^2);
hold on
plotHandle_Mean = plot([p_0*0*dt/m,p_0*0*dt/m],[0, max(abs(step_three(plotX)).^2)]);
plotHandle_2 = plot(x(plotX), abs(Gaussian_Wave_Time_Evolution(plotX,0)).^2,'--');
hold off

xlabel('Position / [\AA]', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Probability distribution', 'fontsize', 14)
s=sprintf('$\\left| \\psi (t = %i \\; \\mathrm{fs}) \\right|^2$',j*dt);
title(s, 'interpreter', 'latex', 'fontsize', 18)


gaussian_momentum = abs(fftshift(fft(step_three(plotX)))).^2;
gaussian_position = abs(step_three(plotX)).^2;

width_position= calculate_width(gaussian_position,dx);
width_momentum= calculate_width(gaussian_momentum,dp);

gaussian_width_position=width_position;
gaussian_width_momentum=width_momentum;

for j=1:n_points/2
    step_one = step_three;    
    step_two = fftshift(fft(step_one.*exp_potential));
    step_three = ifft(ifftshift(inv_pot.*step_two));    
    gaussian_momentum = fftshift(fft(abs(step_three(plotX)).^2));
    gaussian_position = abs(step_three(plotX)).^2;
    width_position= calculate_width(gaussian_position,dx);
    width_momentum= calculate_width(gaussian_momentum,dp);
    gaussian_width_momentum=[gaussian_width_momentum width_momentum];
    gaussian_width_position=[gaussian_width_position width_position];
end

set(plotHandle_1, 'YData', abs(step_three(plotX)).^2)
set(plotHandle_Mean, 'YData',[0 max(Gaussian_Time_Probability(x(plotX),j*dt))])
set(plotHandle_Mean, 'XData', [p_0*j*dt/m,p_0*j*dt/m])
set(plotHandle_2, 'YData', Gaussian_Time_Probability(x(plotX),j*dt))

s=sprintf('$\\left| \\psi (t = %i \\; \\mathrm{fs}) \\right|^2$',j*dt);
title(s, 'interpreter', 'latex', 'fontsize', 18)



figure(3); clf;
plotHandle_1 = plot(p(plotX)/p_0, abs(step_three(plotX)).^2);
hold on
plotHandle_2 = plot(p(plotX)/p_0, abs(fftshift(fft(Gaussian_Wave_Time_Evolution(x(plotX),j*dt)))).^2,'--');
hold off
xlabel('$P/p_0$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('Probability distribution', 'fontsize', 14)
s=sprintf('$\\left| \\psi (t = %i \\; \\mathrm{fs}) \\right|^2$',j*dt);
title(s, 'interpreter', 'latex', 'fontsize', 18)

fig = figure(4);
left_color = [0 0 0];
right_color = [0 0 0];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
yyaxis left
plotHandle_wPosition = plot(gaussian_width_position,'b');
hold on
plotHandle_wPosition_Theoretical = plot(Theoretical_Width_Position(((1:length(gaussian_width_position))-1)*dt),'--r');
hold off
xlabel('$itererings variabel$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('\Delta X(t)', 'fontsize', 14)
yyaxis right 
plotHandle_wMomentum = plot(gaussian_width_momentum,'r');
hold on
plotHandle_wMomentum_Theoretical = plot(Theoretical_Width_Momentum(((1:length(gaussian_width_position))-1)*dt),'--g');
hold off
ylabel('\Delta P(t)', 'fontsize', 14)
s=sprintf('$Uncertainty of Wave Packets$',j*dt);
title(s, 'interpreter', 'latex', 'fontsize', 18)