% u_t + (u^2/2)_x = 0, x in [-pi, pi]
% IC = u0(x)
% Use Lax-Friedrichs Scheme to approximate nonlinear advection equation

clear variables; clc; close all;

% ----------- SPATIAL CONVERGENCE -------------
u0 = @(x) sin(x);
a = -1;
I = [-pi, pi];
f = @(u) a*(u.^2)/2;
tf = 1; % final time
alpha = 1; % alpha = max_u(|f'(u)|) --> f'(u) = u0(x), which has a max val of 1
Nxvals = [10, 20, 40, 80]; % N+1 points, N intervals
errors = zeros(numel(Nxvals), 2);

for k = 1:size(errors, 1)
    Nx = Nxvals(k);
    xvals = linspace(I(1), I(2), Nx + 1)';
    dx = xvals(2) - xvals(1);
    dt = 0.5*dx; % CFL condition: dt = 0.5dx
    tvals = (0:dt:tf)';
    if tvals(end)~=tf
        tvals = [tvals; tf];
    end
    
    xmid = xvals(1:end-1) + dx/2;
    u_bar = (1/dx)*(-cos(xvals(2:end)) + cos(xvals(1:end-1))); % initial condition
    
    u_exact_eqn = @(y) y-sin(xmid-sign(a)*(tf*y));
    u_exact = fsolve(u_exact_eqn, xmid, optimset('Display', 'off'));
    
    for n = 2:numel(tvals)
        u_mid = u_bar;
        u_pos = [u_bar(2:end) ; u_bar(1)];
        f_pos = 0.5*(f(u_bar) + f(u_pos) - (alpha)*(u_pos - u_bar));
        f_neg = [f_pos(end); f_pos(1:end-1)];

        u_bar = u_bar - (dt/dx)*(f_pos - f_neg);
    end

    errors(k, 1) = dx*(sum(abs(u_bar - u_exact))); % L1 error
    errors(k, 2) = sqrt(dx*sum((u_bar - u_exact).^2)); % L2 error
end

L1_error = errors(:, 1);
L2_error = errors(:, 2);
L1_order = [0; log2(errors(1:end-1, 1) ./ errors(2:end, 1))];
L2_order = [0; log2(errors(1:end-1, 2) ./ errors(2:end, 2))];

table(L1_error, L1_order, L2_error, L2_order)


figure(1); clf;
plot(xmid, u_exact); hold on;
plot(xmid, u_bar);
legend('Exact', 'Approx');













