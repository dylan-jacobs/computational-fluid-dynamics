% u_t + (u^2/2)_x = 0, x in [-pi, pi]
% IC = u0(x)
% Use Lax-Friedrichs Scheme to approximate nonlinear advection equation

clear variables; clc; close all;

% ---- Riemann Problems --------
Nx = 160; % N+1 points, N intervals
I = [-pi, pi]; % x interval
tf = 1; % final time
alpha = 1; % alpha = max_u(|f'(u)|) --> f'(u) = u0(x), which has a max val of 1

% 1
u0 = @(x) (x <= 0) + (x > 0)*(-0.5); % initial condition
u_exact = @(x) (x<=0.25) + (x>0.25)*(-0.5);
burgers_equation(Nx, I, u0, tf, alpha, u_exact);

% 2
u0 = @(x) (x <= 0)*(-0.5) + (x > 0)*(1); % initial condition
u_exact = @(x) (x<=-0.5)*(-0.5) + (x>-0.5 & x <= 1).*(x) + (x>1);
burgers_equation(Nx, I, u0, tf, alpha, u_exact);
