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
p_0         = sqrt(0.12*2*m);
x_0         = 0;
dx          = 0.1;
n_points    = 2^10;
dp          = 2*pi/(n_points*dx);
dt          = 1;
a           = 0.3;
b           = 0.4;
d           = 0.7;
c           = 0.05; % 0.1;

V_11 = @(x) 
if (x > 0) 
    a*(2-exp(-x/b)) 
else
    a*(2*exp(-x/b)) 
end;

V_12 = @(x) ;
V_21 = @(x) a*(2*exp(-x/b));
V_22 = @(x) a*(2*exp(-x/b));
V_d =[];
