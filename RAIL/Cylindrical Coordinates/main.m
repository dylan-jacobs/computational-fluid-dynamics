clear variables; close all; clc;

tf = 0.5;
Rf = 1; Zf = 1;
interval = [0, Rf, 0, 1];
Nr = 100; Nz = 100;
r0 = 100;
tolerance = 1e-12;

j01 = 2.40482555769577; % first root of bessel function
u0 = @(r, z) (besselj(0, (j01/Rf).*r)) + (sin((pi/Rf).*r));
u_exact = @(r, z, t) (exp(-t*(((j01/Rf).*r) + (pi/Rf).*r)).*u0(r, z));
[R, Z, dr, dz] = GetRZ(Nr, Nz, interval);

lambdavals = (0.1:0.02:1)';
errors = zeros(numel(lambdavals), 3); % L1 norm for B. Euler, DIRK2, and DIRK3

for k = 1:numel(lambdavals)
    dt = lambdavals(k)*dr;
    disp([num2str(k), '/', num2str(numel(lambdavals))]);

    [U1, ranks1] = HeatEqnSolver('1', u0(R, Z), dt, Nr, Nz, tf, interval, 1, 1, tolerance, r0);
    U_exact = u_exact(R, Z, tf);
    errors(k, 1) = dr*dz*(sum(sum(abs(U1 - U_exact))));
end

figure(1); clf; surf(R, Z, U1);
colorbar;
shading flat; % removes gridlines
legend(sprintf('N_x = %s, N_y = %s', num2str(Nr, 3), num2str(Nz, 3)), 'Location','northwest');
xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('RAIL approximation at time %s', num2str(tf, 4))]);

figure(2); clf; surf(R, Z, U_exact);
colorbar;
shading flat; % removes gridlines
xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('U_{Exact} at time %s', num2str(tf, 4))]);

figure(3); clf;
loglog(lambdavals, errors(:, 1), 'black-', 'LineWidth', 1.5); hold on; % B. Euler
% loglog(lambdavals, errors(:, 2), 'blue-', 'LineWidth', 1.5); % DIRK2
% loglog(lambdavals, errors(:, 3), 'green-', 'LineWidth', 1.5); % DIRK3
loglog(lambdavals, 0.08*lambdavals, 'black--', 'LineWidth', 1); % Order 1
loglog(lambdavals, 0.0003*lambdavals.^2, 'blue--', 'LineWidth', 1); % Order 2
loglog(lambdavals, 0.000008*lambdavals.^3, 'green--', 'LineWidth', 1); % Order 3
title('RAIL Temporal Convergence at tf=0.5, Nx = Ny = 200'); xlabel('\lambda'); ylabel('L1 Error');
legend('Backward Euler', 'DIRK2', 'DIRK3', 'Order 1', 'Order 2', 'Order 3');