clear variables; close all; clc;

% approximate solution curves
% u_t + u_x = 0
% u0(x) = {1, x<= 0.5; 0, x>0.5}
% tf = 1, 201 gridpoints, 500 timesteps

Nx = 200; % N intervals, N+1 points
Nt = 500;
t0 = 0; tf = 1;
x0 = 0; xf = 1;

u0 = @(x) (x-0.5)<=0; % Heaviside step function about x=0.5

xvals = linspace(x0, xf, Nx + 1)'; 
tvals = linspace(t0, tf, Nt)'; % assume stability conditions are met
dx = xvals(2) - xvals(1);
dt = tvals(2) - tvals(1);

%% UPWIND SCHEME
Dx = gallery('tridiag', Nx + 1, -1, 1, 0);
Dx(1, end) = -1; % period BC

u = u0(xvals); % initial condition

for n = 2:numel(tvals)
    u = u - (dt/dx)*(Dx*u);
end

figure(1); clf;
plot(xvals, u0(xvals), 'black--', 'LineWidth', 1); hold on;
plot(xvals, u, 'b-', 'LineWidth', 1.5);
ylim([-0.2, 1.2]);
title('Upwind scheme at t=1'); xlabel('x'); ylabel('u');
legend('Heaviside Step Function around x=0.5', 'Upwind Approximation');

%% LAX-WENDROFF
Dx = gallery('tridiag', Nx, -1, 0, 1);
Dx(1, end) = -1; Dx(end, 1) = 1; % periodic BCs

Dxx = gallery('tridiag', Nx, 1, -2, 1);
Dxx(1, end) = 1; Dxx(end, 1) = 1; % periodic BCs

xvals = xvals(1:end-1);
u = u0(xvals); % initial condition

for n = 2:numel(tvals)
    u = u - (dt/(2*dx))*(Dx*u) + (dt^2/(2*dx^2))*(Dxx*u); % Lax-Wendroff
end

figure(2); clf;
plot(xvals, u0(xvals), 'black--', 'LineWidth', 1); hold on;
plot(xvals, u, 'b-', 'LineWidth', 1.5);
ylim([-0.3, 1.3]);
title('Lax-Wendroff scheme at t=1'); xlabel('x'); ylabel('u');
legend('Heaviside Step Function around x=0.5', 'Lax-Wendroff Approximation');








