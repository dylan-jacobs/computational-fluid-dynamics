% solves advection-diffusion eqn Ut + a1Ux + a2Uy = d1Uxx + d2Uyy + phi(t)

%% TEST 1 -- Accuracy
clc; clear variables; close all;

d = 1/5;
A = cell(2, 3);
    A{1, 1} = @(x) x.^0;
    A{1, 2} = @(t) -1;
    A{1, 3} = @(y) y;
    A{2, 1} = @(x) x;
    A{2, 2} = @(t) 1;
    A{2, 3} = @(y) y.^0;
tf = 0.5;
interval = [-2*pi, 2*pi, -2*pi, 2*pi];
Nx = 200; Ny = 200;
tolerance = 1e-8;
phi = cell(1, 3);
    phi{1, 1} = @(x, t) [exp(-(x.^2)), x.*exp(-(x.^2)), (x.^2).*exp(-(x.^2)), exp(-(x.^2))];
    phi{1, 2} = @(t) diag([6*d*exp(-2*d*t), -4*exp(-2*d*t), -4*d*exp(-2*d*t), -36*d*exp(-2*d*t)]);
    phi{1, 3} = @(y, t) [exp(-3*(y.^2)), y.*exp(-3*(y.^2)), exp(-3*(y.^2)), (y.^2).*exp(-3*(y.^2))];
u_exact = @(x, y, t) exp(-((x.^2) + (3*y.^2) + (2*d*t)));
u0 = @(x, y) u_exact(x, y, 0);
[X, Y, dx, dy] = GetXY(Nx, Ny, interval);

lambdavals = (0.1:0.02:2)';
errors = zeros(numel(lambdavals), 3); % L1 norms for IMEX(1, 1, 1), (2, 2, 2), (4, 4, 3)

for k = 1:numel(lambdavals)
    dt = lambdavals(k)*dx;
    disp([num2str(k), '/', num2str(numel(lambdavals))]);

    % [U1, ~] = IMEX('111', u0(X, Y), dt, Nx, Ny, tf, interval, A, d, d, phi, tolerance);
    % U_exact = u_exact(X, Y, tf);
    % errors(k, 1) = dx*dy*(sum(sum(abs(U1 - U_exact))));
    % 
    [U2, ~] = IMEX('222', u0(X, Y), dt, Nx, Ny, tf, interval, A, d, d, phi, tolerance);
    U_exact = u_exact(X, Y, tf);
    errors(k, 2) = dx*dy*(sum(sum(abs(U2 - U_exact))));

    % [U3, ~] = DIRK3(u0(X, Y), dt, Nx, Ny, tf, interval, d1, d2, tolerance);
    % U_exact = u_exact(X, Y, tf);
    % errors(k, 3) = dx*dy*(sum(sum(abs(U3 - U_exact))));
end

figure(1); clf; surf(X, Y, U2);
colorbar;
shading flat; % removes gridlines
legend(sprintf('N_x = %s, N_y = %s', num2str(Nx, 3), num2str(Ny, 3)), 'Location','northwest');
xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('RAIL approximation at time %s', num2str(tf, 4))]);

figure(2); clf; surf(X, Y, U_exact);
colorbar;
shading flat; % removes gridlines
xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('U_{Exact} at time %s', num2str(tf, 4))]);

figure(3); clf;
loglog(lambdavals, errors(:, 1), 'black-', 'LineWidth', 1.5); hold on; % IMEX111
loglog(lambdavals, errors(:, 2), 'blue-', 'LineWidth', 1.5); % IMEX222
loglog(lambdavals, errors(:, 3), 'green-', 'LineWidth', 1.5); % DIRK3
loglog(lambdavals, 0.08*lambdavals, 'black--', 'LineWidth', 1); % Order 1
loglog(lambdavals, 0.0003*lambdavals.^2, 'blue--', 'LineWidth', 1); % Order 2
loglog(lambdavals, 0.000008*lambdavals.^3, 'green--', 'LineWidth', 1); % Order 3
title('RAIL Temporal Convergence at tf=0.5, Nx = Ny = 200'); xlabel('\lambda'); ylabel('L1 Error');
legend('IMEX111', 'IMEX222', 'IMEX443', 'Order 1', 'Order 2', 'Order 3');

%% TEST 2 -- Rank Test
clc; clear variables; close all;

d1 = 1/4; d2 = 1/9;
tf = 5;
interval = [0, 14, 0, 14];
tolerance = 1e-6;
Nx = 200; Ny = 200;
u0 = @(x, y) 0.8*exp(-15 * ((x - 6.5).^2 + (y - 6.5).^2)) + (0.5*exp(-15 * ((x - 7.5).^2 + (y - 7).^2)));
u_exact = @(x, y, t) (0.8 / (4 * pi * sqrt(d1 * d2 * t))) .* ...
    exp(-((x - 6.5).^2 / (4 * d1 * t) + (y - 6.5).^2 / (4 * d2 * t))) + ...
    (0.5 / (4 * pi * sqrt(d1 * d2 * t))) .* ...
    exp(-((x - 7.5).^2 / (4 * d1 * t) + (y - 7).^2 / (4 * d2 * t)));
[X, Y, dx, dy] = GetXY(Nx, Ny, interval);

lambdavals = (0.1:0.02:1)';
errors = zeros(numel(lambdavals), 3); % L1 norm
U_exact = u_exact(X, Y, tf);

for k = 1:1%numel(lambdavals)
    dt = lambdavals(k)*dx;
    disp([num2str(k), '/', num2str(numel(lambdavals))]);

    [U1, ranks1] = Backward_Euler(u0(X, Y), dt, Nx, Ny, tf, interval, d1, d2, tolerance);

    [U2, ranks2] = DIRK2(u0(X, Y), dt, Nx, Ny, tf, interval, d1, d2, tolerance);

    [U3, ranks3] = DIRK3(u0(X, Y), dt, Nx, Ny, tf, interval, d1, d2, tolerance);
end

figure(1); clf; surf(X, Y, U3);
colorbar;
shading flat; % removes gridlines
legend(sprintf('N_x = %s, N_y = %s', num2str(Nx, 3), num2str(Ny, 3)), 'Location','northwest');
xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('RAIL approximation at time %s', num2str(tf, 4))]);

figure(2); clf; surf(X, Y, U_exact);
colorbar;
shading flat; % removes gridlines
xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('U_{Exact} at time %s', num2str(tf, 4))]);

figure(3); clf;
loglog(lambdavals, errors(:, 1), 'black--', 'LineWidth', 1); hold on;
loglog(lambdavals, 0.05*lambdavals, 'b-', 'LineWidth', 1.5);
title('RAIL Temporal Convergence at tf=0.5'); xlabel('\lambda'); ylabel('L1 Error');
legend('Nx = Ny = 200', 'Order 1');

figure(4); clf;
plot(ranks1(:, 1), ranks1(:, 2), 'black-', 'LineWidth', 1.5); hold on;
plot(ranks2(:, 1), ranks2(:, 2), 'blue-', 'LineWidth', 1.5);
plot(ranks3(:, 1), ranks3(:, 2), 'green-', 'LineWidth', 1.5);
xlabel('time'); ylabel('rank'); title('Rank plot over time');
legend('Backward Euler', 'DIRK2', 'DIRK3');




