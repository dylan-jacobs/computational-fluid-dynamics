% solves advection-diffusion eqn Ut + a1Ux + a2Uy = d1Uxx + d2Uyy + phi(t)

%% TEST 1 -- Rigidbody Rotation Accuracy
clc; clear variables; close all;

d = 1/5;
A = cell(2, 3);
    A{1, 1} = @(x) x.^0;
    A{1, 2} = @(t) -1;
    A{1, 3} = @(y) y;
    A{2, 1} = @(x) x;
    A{2, 2} = @(t) 1;
    A{2, 3} = @(y) y.^0;
tf = 2.5;
interval = [-2*pi, 2*pi, -2*pi, 2*pi];
Nx = 100; Ny = 100;
tolerance = 1e-8;
r0 = 15;
phi = cell(1, 3);
    phi{1, 1} = @(x, t) [exp(-(x.^2)), x.*exp(-(x.^2)), (x.^2).*exp(-(x.^2)), exp(-(x.^2))];
    phi{1, 2} = @(t) diag([6*d*exp(-2*d*t), -4*exp(-2*d*t), -4*d*exp(-2*d*t), -36*d*exp(-2*d*t)]);
    phi{1, 3} = @(y, t) [exp(-3*(y.^2)), y.*exp(-3*(y.^2)), exp(-3*(y.^2)), (y.^2).*exp(-3*(y.^2))];
% turn source term off
% phi = cell(1, 3);
%     phi{1, 1} = @(x, t) 0;
%     phi{1, 2} = @(t) 0;
%     phi{1, 3} = @(y, t) 0;
u_exact = @(x, y, t) exp(-((x.^2) + (3*y.^2) + (2*d*t)));
u0 = @(x, y) u_exact(x, y, 0);
[X, Y, dx, dy] = GetXY(Nx, Ny, interval);

lambdavals = 0.05;%(0.1:0.02:2)';
errors = zeros(numel(lambdavals), 3); % L1 norms for IMEX(1, 1, 1), (2, 2, 2), (4, 4, 3)

for k = 1:numel(lambdavals)
    dt = lambdavals(k)*dx;
    disp([num2str(k), '/', num2str(numel(lambdavals))]);

    [U1, ranks1] = IMEX('111', u0(X, Y), dt, Nx, Ny, tf, interval, A, d, d, phi, tolerance, r0);
    U_exact = u_exact(X, Y, tf);
    errors(k, 1) = dx*dy*(sum(sum(abs(U1 - U_exact))));

    [U2, ranks2] = IMEX('222', u0(X, Y), dt, Nx, Ny, tf, interval, A, d, d, phi, tolerance, r0);
    U_exact = u_exact(X, Y, tf);
    errors(k, 2) = dx*dy*(sum(sum(abs(U2 - U_exact))));

    [U3, ranks3] = IMEX('443', u0(X, Y), dt, Nx, Ny, tf, interval, A, d, d, phi, tolerance, r0);
    U_exact = u_exact(X, Y, tf);
    errors(k, 3) = dx*dy*(sum(sum(abs(U3 - U_exact))));
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
loglog(lambdavals, errors(:, 1), 'black-', 'LineWidth', 1.5); hold on; % IMEX111
loglog(lambdavals, errors(:, 2), 'blue-', 'LineWidth', 1.5); % IMEX222
loglog(lambdavals, errors(:, 3), 'green-', 'LineWidth', 1.5); % DIRK3
loglog(lambdavals, 0.08*lambdavals, 'black--', 'LineWidth', 1); % Order 1
loglog(lambdavals, 0.003*(lambdavals.^2), 'blue--', 'LineWidth', 1); % Order 2
loglog(lambdavals, 0.00008*(lambdavals.^3), 'green--', 'LineWidth', 1); % Order 3
title('RAIL Temporal Convergence at tf=0.5, Nx = Ny = 100'); xlabel('\lambda'); ylabel('L1 Error');
legend('IMEX111', 'IMEX222', 'IMEX443', 'Order 1', 'Order 2', 'Order 3');

figure(4); clf;
plot(ranks1(:, 1), ranks1(:, 2), 'black-', 'LineWidth', 1.5); hold on;
plot(ranks2(:, 1), ranks2(:, 2), 'blue-', 'LineWidth', 1.5);
plot(ranks3(:, 1), ranks3(:, 2), 'green-', 'LineWidth', 1.5);
xlabel('time'); ylabel('rank'); title('Rank plot over time');
legend('IMEX111', 'IMEX222', 'IMEX443');


%% TEST 2: Swirling Deformation with Diffusion
clc;  close all;

tf = 2.5;
interval = [-pi, pi, -pi, pi];
Nx = 100; Ny = 100;
tolerance = 1e-8;
r0 = 30;

d = 1;
gt = @(t) cos(pi*t/tf)*pi;
A = cell(2, 3);
    A{1, 1} = @(x) (cos(x./2).^2);
    A{1, 2} = @(t) -gt(t);
    A{1, 3} = @(y) sin(y);
    A{2, 1} = @(x) sin(x);
    A{2, 2} = @(t) gt(t);
    A{2, 3} = @(y) (cos(y./2).^2);
% f = @(u, x, y, t) -(cos(x./2).^2).*(sin(y).*gt(t).*u);
% g = @(u, x, y, t) (sin(x).*(cos(y./2).^2).*gt(t).*u);
u0 = @(x, y) cos_bell(x, y);
phi = cell(1, 3);
    phi{1, 1} = @(x, t) 0;
    phi{1, 2} = @(t) 0;
    phi{1, 3} = @(y, t) 0;

% U_exact(X_exact_soln, Y_exact_soln, dx_exact_soln) = IMEX('443', u0(X_exact_soln, Y_exact_soln), 0.05*dx_exact_soln, 300, 300, tf, interval, A, d, d, phi, tolerance);

[X, Y, dx, dy] = GetXY(Nx, Ny, interval);

lambdavals = (0.15:0.02:2)';
errors = zeros(numel(lambdavals), 3); % L1 norms for IMEX(1, 1, 1), (2, 2, 2), (4, 4, 3)

for k = 1:1%numel(lambdavals)
    dt = lambdavals(k)*dx;
    disp([num2str(k), '/', num2str(numel(lambdavals))]);

    [U1, ranks1] = IMEX('111', u0(X, Y), dt, Nx, Ny, tf, interval, A, d, d, phi, tolerance, r0);
    % errors(k, 1) = dx*dy*(sum(sum(abs(U1 - U_exact))));

    [U2, ranks2] = IMEX('222', u0(X, Y), dt, Nx, Ny, tf, interval, A, d, d, phi, tolerance, r0);
    % errors(k, 2) = dx*dy*(sum(sum(abs(U2 - U_exact))));

    [U3, ranks3] = IMEX('443', u0(X, Y), dt, Nx, Ny, tf, interval, A, d, d, phi, tolerance, r0);
    % errors(k, 3) = dx*dy*(sum(sum(abs(U3 - U_exact))));
end

figure(1); clf; surf(X, Y, U3);
colorbar;
shading flat; % removes gridlines
legend(sprintf('N_x = %s, N_y = %s', num2str(Nx, 3), num2str(Ny, 3)), 'Location','northwest');
xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('RAIL approximation at time %s', num2str(tf, 4))]);

% figure(2); clf; surf(X, Y, U_exact);
% colorbar;
% shading flat; % removes gridlines
% xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('U_{Exact} at time %s', num2str(tf, 4))]);

% figure(3); clf;
% loglog(lambdavals, errors(:, 1), 'black-', 'LineWidth', 1.5); hold on; % IMEX111
% loglog(lambdavals, errors(:, 2), 'blue-', 'LineWidth', 1.5); % IMEX222
% loglog(lambdavals, errors(:, 3), 'green-', 'LineWidth', 1.5); % DIRK3
% loglog(lambdavals, 0.08*lambdavals, 'black--', 'LineWidth', 1); % Order 1
% loglog(lambdavals, 0.003*(lambdavals.^2), 'blue--', 'LineWidth', 1); % Order 2
% loglog(lambdavals, 0.00008*(lambdavals.^3), 'green--', 'LineWidth', 1); % Order 3
% title('RAIL Temporal Convergence at tf=0.5, Nx = Ny = 100'); xlabel('\lambda'); ylabel('L1 Error');
% legend('IMEX111', 'IMEX222', 'IMEX443', 'Order 1', 'Order 2', 'Order 3');

figure(2); clf;
plot(ranks1(:, 1), ranks1(:, 2), 'black-', 'LineWidth', 1.5); hold on;
plot(ranks2(:, 1), ranks2(:, 2), 'blue-', 'LineWidth', 1.5);
plot(ranks3(:, 1), ranks3(:, 2), 'green-', 'LineWidth', 1.5);
xlabel('time'); ylabel('rank'); title('Rank plot over time');
legend('IMEX111', 'IMEX222', 'IMEX443');



function [output] = cos_bell(x, y)
    rb0 = 0.3*pi;
    rb = sqrt((x-(0.3*pi)).^2 + ((y).^2));
    output = (rb < rb0).*(rb.*(cos(rb.*pi./(2.*rb0))).^6);
end

function [output] = swirl(x, y)    
    rb0 = 8*pi/5;
    rb = sqrt((x-(pi)).^2 + ((y-pi).^2));
    output = (rb < rb0).*1;
end



%% TEST 3 - Lenard-Bernstien-Fokker-Planck Equation

clc; close all;

tf = 15;
interval = [-8, 8, -8, 8];
Nx = 300; Ny = 300;
tolerance = 1e-6;
r0 = 30;

% model-specific parameters
R = 1/6; T = 3; 
n = pi; vx_bar = 0; vy_bar = 0;
fM = @(vx, vy, n, vx_bar, vy_bar, T) (n/(2*pi*R*T))*exp(-( ((vx-vx_bar).^2) + ((vy - vy_bar).^2) ) ./ (2*R*T));

d = 0.5;
A = cell(2, 3);
    A{1, 1} = @(vx) (vx - vx_bar);
    A{1, 2} = @(t) -1;
    A{1, 3} = @(vy) (1);
    A{2, 1} = @(vx) (1);
    A{2, 2} = @(t) -1;
    A{2, 3} = @(vy) (vy - vy_bar);
    
u0 = @(Vx, Vy) fM(Vx, Vy, 1.990964530353041, 0.4979792385268875, 0, 2.46518981703837) + fM(Vx, Vy, 1.150628123236752, -0.8616676237412346, 0, 0.4107062104302872);
phi = cell(1, 3);
    phi{1, 1} = @(vx, t) 0;
    phi{1, 2} = @(t) 0;
    phi{1, 3} = @(vy, t) 0;

[Vx, Vy, dx, dy] = GetXY(Nx, Ny, interval);
U_exact = fM(Vx, Vy, n, vx_bar, vy_bar, T);

lambdavals = 0.15;
errors = zeros(3, 1877); % L1 norms for IMEX(1, 1, 1), (2, 2, 2), (4, 4, 3)

for k = 1:numel(lambdavals)
    dt = lambdavals(k)*dx;
    disp([num2str(k), '/', num2str(numel(lambdavals))]);

    [U1, ranks1] = IMEX('111', u0(Vx, Vy), dt, Nx, Ny, tf, interval, A, d, d, phi, tolerance, r0);
    errors(1, :) = dx*dy*(sum(sum(abs(U1 - U_exact))));

    [U2, ranks2] = IMEX('222', u0(Vx, Vy), dt, Nx, Ny, tf, interval, A, d, d, phi, tolerance, r0);
    errors(2, :) = dx*dy*(sum(sum(abs(U2 - U_exact))));

    [U3, ranks3] = IMEX('443', u0(Vx, Vy), dt, Nx, Ny, tf, interval, A, d, d, phi, tolerance, r0);
    errors(3, :) = dx*dy*(sum(sum(abs(U3 - U_exact))));
end

figure(1); clf; surf(Vx, Vy, U1(:, :, end));
colorbar;
shading flat; % removes gridlines
legend(sprintf('N_x = %s, N_y = %s', num2str(Nx, 3), num2str(Ny, 3)), 'Location','northwest');
xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('RAIL approximation at time %s', num2str(tf, 4))]);

figure(2); clf; surf(Vx, Vy, U_exact);
colorbar;
shading flat; % removes gridlines
xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('U_{Exact} at time %s', num2str(tf, 4))]);

figure(3); clf;
semilogy(errors(1, :), 'black-', 'LineWidth', 1.5); hold on; % IMEX111
semilogy(errors(2, :), 'blue-', 'LineWidth', 1.5); % IMEX222
semilogy(errors(3, :), 'green-', 'LineWidth', 1.5); % IMEX443
title('IMEX Temporal Convergence at tf=15, Nx = Ny = 100'); xlabel('\lambda'); ylabel('L1 Error');
legend('IMEX111', 'IMEX222', 'IMEX443');

figure(2); clf;
plot(ranks1(:, 1), ranks1(:, 2), 'black-', 'LineWidth', 1.5); hold on;
plot(ranks2(:, 1), ranks2(:, 2), 'blue-', 'LineWidth', 1.5);
plot(ranks3(:, 1), ranks3(:, 2), 'green-', 'LineWidth', 1.5);
xlabel('time'); ylabel('rank'); title('Rank plot over time');
legend('IMEX111', 'IMEX222', 'IMEX443');













