%clear variables; close all; clc;
% ---- MUSCL Problems --------
Nx = 100; % N+1 points, N intervals
I = [-pi, pi]; % x interval
tf = 2; % final time
xvals = linspace(I(1), I(2), Nx + 1)';
dx = xvals(2) - xvals(1);
xmid = xvals(1:end-1) + dx/2;
f = @(u) 0.5*u.^2;

% 1
u0 = @(x) sin(x-pi) + 1; % initial condition
u_exact_eqn = @(y) y-(sin(xmid-(tf*y)-pi)+1);
u_exact = fsolve(u_exact_eqn, xmid, optimset('Display', 'off'));
MUSCL(xvals, u0, tf, u_exact, f);

% 2 - Riemann Problem
u0 = @(x) (x <= 0)*(2); % initial condition
u_exact = (xmid<=(tf))*(2);
MUSCL(xvals, u0, tf, u_exact, f);

% 3 - Riemann Problem
u0 = @(x) (x > 0)*(2); % initial condition
u_exact = (xmid>0 & xmid <= tf).*(xmid/tf) + (xmid>tf)*(2);
MUSCL(xvals, u0, tf, u_exact, f);








