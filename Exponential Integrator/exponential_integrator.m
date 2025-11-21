% Time exponential integrator

clc; clear variables; close all;

Nx = 100;
interval = [0, pi];
xvals = linspace(interval(1), interval(2), Nx+1)';
dx = xvals(2) - xvals(1);

Tf = 15;
dt = dx;
tvals = 0:dt:Tf;

% U_t = k*U_xx ===> F(t, u) = k*U_xx
k = 1;
u0 = @(x) sin(2*x);
u_exact = @(x, t) exp(-4*k^2*t).*sin(2*k*x);

u = u0(xvals);
phi = @(z) (expm(z) - 1) ./ z;
g = @(u) 0;
figure(1); clf;
plot(xvals, u);
Dxx = (1/dx^2)*gallery("tridiag", Nx+1, 1, -2, 1);
A = k*Dxx;
for n = 2:numel(tvals)
    t = tvals(n);
    dt = t - tvals(n-1);

    u = expm(dt*A)*u + dt*phi(-dt*A).*g(u);
    u(1) = 0; u(end) = 0; % homo. Dirichlet BCs

    figure(2); clf;
    plot(u);
end







