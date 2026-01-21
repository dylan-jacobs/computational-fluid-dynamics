%% PART A
% Forward Euler + Upwind Scheme
% u_t - u_x = 0
% u0(x) = sin(2*pi*x)
% period BCs
% u^(n+1)_j = u^(n)_j + (∆t/∆x)(u^n_j − u^n_(j−1)))

clc; clear variables; close all;


% ---------- Temporal convergence --------------
tf = 3; % final time
interval = [0, 1];
a = -1; % advection speed
Nx = 100;
xvals = linspace(interval(1), interval(2), Nx + 1)';
dx = xvals(2) - xvals(1);

u_exact = sin(2*pi*(xvals - (a*tf))); % exact soln

lambdavals = 0.01:0.02:1; % "convergent" for dt <= dx
errors = zeros(numel(lambdavals), 1);
for k = 1:numel(lambdavals)
    dt = lambdavals(k)*dx;

    tvals = (0:dt:tf)';
    if (tvals(end) ~= tf)
        tvals = [tvals; tf];
    end

    Dx = (gallery('tridiag',Nx+1,-1,1,0)); Dx(1, end) = -1; % periodic BCs

    u = sin(2*pi.*xvals); % u0

    for t = 2:numel(tvals)
        u = u - (a*dt/(dx))*(Dx*u);
    end

    errors(k) = dx*sum(abs(u - u_exact));
end


% Temporal error plot
figure(1);
loglog(lambdavals, errors, 'r-'); hold on;
loglog(lambdavals, lambdavals, 'r--');
title('Part A error plot');
legend('Approximation Error', 'Order 1');

% Solution curves
figure(2);
plot(xvals, u); hold on;
plot(xvals, u_exact);
title('Part A solution curves');
xlabel('x'); ylabel('u');
legend('approx', 'exact');




% ----------- Spatial convergence ----------
Nvals = [10, 20, 40, 80]';
errors = zeros(numel(Nvals), 1);

for k = 1:numel(Nvals)
    Nx = Nvals(k);
    xvals = linspace(interval(1), interval(2), Nx + 1)'; % n+1 points, n intervals
    u_exact = sin(2*pi*(xvals - (a*tf)));
    dx = xvals(2) - xvals(1);
    dt = dx;

    tvals = (0:dt:tf)';
    if (tvals(end) ~= tf)
        tvals = [tvals; tf];
    end

    Dx = (gallery('tridiag',Nx+1,-1,1,0)); Dx(1, end) = -1; % periodic BCs

    u = sin(2*pi.*xvals); % u0

    for t = 2:numel(tvals)
        u = u - (a*dt/(dx))*(Dx*u);
    end

    errors(k) = dx*sum(abs(u - u_exact));

end

disp('Order = ')
disp(log2(errors(1:end-1)./errors(2:end)));


%% PART B
% Forward Euler + Downwind Scheme
% u_t - u_x = 0
% u0(x) = sin(2*pi*x)
% period BCs
% u^(n+1)_j = u^(n)_j + (∆t/∆x)(u^n_j − u^n_(j−1)))

clc; clear variables; close all;

% -------------- Temporal convergence -----------------
tf = 3; % final time
interval = [0, 1];
a = -1; % advection speed
Nx = 500;
xvals = linspace(interval(1), interval(2), Nx + 1)';
dx = xvals(2) - xvals(1);

u_exact = sin(2*pi*(xvals - (a*tf)));

lambdavals = 0.01:0.02:1; % convergent for dt <= dx
errors = zeros(numel(lambdavals), 1);
for k = 1:numel(lambdavals)
    dt = lambdavals(k)*dx;

    tvals = (0:dt:tf)';
    if (tvals(end) ~= tf)
        tvals = [tvals; tf];
    end

    Dx = (gallery('tridiag',Nx+1, 0, -1, 1)); Dx(end, 1) = 1; % periodic BC

    u = sin(2*pi.*xvals);

    for t = 2:numel(tvals)
        u = u - (a*dt/(dx))*(Dx*u);
    end

    errors(k) = dx*sum(abs(u - u_exact));
end

% Temporal error plot
figure(1);
loglog(lambdavals, errors, 'r-'); hold on;
loglog(lambdavals, lambdavals, 'r--');
title('Part B error plot');
legend('Approximation Error', 'Order 1');

% Solution curves
figure(2);
plot(xvals, u); hold on;
plot(xvals, u_exact);
title('Part B solution curves');
xlabel('x'); ylabel('u');
legend('approx', 'exact');



% ------------ Spatial convergence ----------------------
Nvals = [10, 20, 40, 80]';
errors = zeros(numel(Nvals), 1);

for k = 1:numel(Nvals)
    Nx = Nvals(k);
    xvals = linspace(interval(1), interval(2), Nx + 1)'; % n+1 points, n intervals
    u_exact = sin(2*pi*(xvals - (a*tf)));
    dx = xvals(2) - xvals(1);
    dt = dx;

    tvals = (0:dt:tf)';
    if (tvals(end) ~= tf)
        tvals = [tvals; tf];
    end

    Dx = (gallery('tridiag',Nx+1, 0, -1, 1)); Dx(end, 1) = 1; % periodic BC

    u = sin(2*pi.*xvals); % u0

    for t = 2:numel(tvals)
        u = u - (a*dt/(dx))*(Dx*u);
    end

    errors(k) = dx*sum(abs(u - u_exact));

end

disp('Order = ')
disp(log2(errors(1:end-1)./errors(2:end)));






