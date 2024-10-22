clear variables; close all; clc;

tf = 1;
Rf = 1; L = 1;
interval = [0, Rf, 0, L];
Nr = 40; Nz = 80;
r0 = 1;
tolerance = 1e-6;

j01 = 2.40482555769577; % first root of bessel function
u0 = @(r, z) (besselj(0, (j01/Rf)*r)) .* (sin((2*pi/L)*z));
u_exact = @(r, z, t) (exp(-t*(((j01/Rf)^2) + (2*pi/L)^2))*u0(r, z));
[R, Z, dr, dz] = GetRZ(Nr, Nz, interval);

lambdavals = 0.05;%(0.02:0.02:10)';
errors = zeros(numel(lambdavals), 3); % L1 norm for B. Euler, DIRK2, and DIRK3
U_exact = u_exact(R, Z, tf);

for k = 1:numel(lambdavals)
    dt = lambdavals(k)/(1/dr + 1/dz);
    disp([num2str(k), '/', num2str(numel(lambdavals))]);

    [U1, ranks1] = HeatEqnSolver('1', u0(R, Z), dt, Nr, Nz, tf, interval, tolerance, r0);
    errors(k, 1) = dr*dz*(sum(sum(abs(U1 - U_exact))));

    [U2, ranks2] = HeatEqnSolver('2', u0(R, Z), dt, Nr, Nz, tf, interval, tolerance, r0);
    errors(k, 2) = dr*dz*(sum(sum(abs(U2 - U_exact))));

    [U3, ranks3] = HeatEqnSolver('3', u0(R, Z), dt, Nr, Nz, tf, interval, tolerance, r0);
    errors(k, 3) = dr*dz*(sum(sum(abs(U3 - U_exact))));
end

figure(1); clf; surf(R, Z, U1);
colorbar;
shading flat; % removes gridlines
legend(sprintf('N_r = %s, N_z = %s', num2str(Nr, 3), num2str(Nz, 3)), 'Location','northwest');
xlabel('R'); ylabel('Z'); zlabel('U(R, Z)'); title([sprintf('RAIL approximation at time %s', num2str(tf, 4))]);

figure(2); clf; surf(R, Z, U_exact);
colorbar;
shading flat; % removes gridlines
xlabel('R'); ylabel('Z'); zlabel('U(R, Z)'); title([sprintf('U_{Exact} at time %s', num2str(tf, 4))]);

figure(3); clf; cutoff = 0.2;
loglog(lambdavals, errors(:, 1), 'black-', 'LineWidth', 1.5); hold on; % B. Euler
loglog(lambdavals, errors(:, 2), 'blue-', 'LineWidth', 1.5); % DIRK2
loglog(lambdavals, errors(:, 3), 'green-', 'LineWidth', 1.5); % DIRK3
loglog(lambdavals(ceil(cutoff*end):end), 0.01*lambdavals(ceil(cutoff*end):end), 'black--', 'LineWidth', 1); % Order 1
loglog(lambdavals(ceil(cutoff*end):end), 0.0003*lambdavals(ceil(cutoff*end):end).^2, 'blue--', 'LineWidth', 1); % Order 2
loglog(lambdavals(ceil(cutoff*end):end), 0.000005*lambdavals(ceil(cutoff*end):end).^3, 'green--', 'LineWidth', 1); % Order 3
title(sprintf('RAIL Temporal Convergence at tf=%s, Nr = %s, Nz = %s', num2str(tf), num2str(Nr), num2str(Nz))); xlabel('\lambda'); ylabel('L1 Error');
legend('Backward Euler', 'DIRK2', 'DIRK3', 'Order 1', 'Order 2', 'Order 3', 'Location','northwest');

figure(4); clf;
plot(ranks1(:, 1), ranks1(:, 2), 'black-', 'LineWidth', 1.5); hold on;
plot(ranks2(:, 1), ranks2(:, 2), 'blue-', 'LineWidth', 1.5);
plot(ranks3(:, 1), ranks3(:, 2), 'green-', 'LineWidth', 1.5);
xlabel('time'); ylabel('rank'); title('Rank plot over time');
legend('Backward Euler', 'DIRK2', 'DIRK3');


%% Maxwellian diffusion rank test

clear variables; close all; clc;

tf = 5;
Rf = 14; L = 14;
interval = [0, Rf, 0, L];
Nr = 80; Nz = 80;
r0 = 10;
tolerance = 1e-6;

u0 = @(r, z) 0.8*(exp(-15*((r - 7).^2 + (z - 4).^2))) + 0.5*(exp(-15*((r - 7).^2 + (z - 12).^2)));
% u_exact = @(r, z, t) (exp(-t*(((j01/Rf)^2) + (2*pi/L)^2))*u0(r, z));
[R, Z, dr, dz] = GetRZ(Nr, Nz, interval);

lambdavals = 0.05;%(0.02:0.02:10)';
% errors = zeros(numel(lambdavals), 3); % L1 norm for B. Euler, DIRK2, and DIRK3
% U_exact = u_exact(R, Z, tf);

for k = 1:numel(lambdavals)
    dt = lambdavals(k)/(1/dr + 1/dz);
    disp([num2str(k), '/', num2str(numel(lambdavals))]);

    [U1, ranks1] = HeatEqnSolver('1', u0(R, Z), dt, Nr, Nz, tf, interval, tolerance, r0);
    % errors(k, 1) = dr*dz*(sum(sum(abs(U1 - U_exact))));

    [U2, ranks2] = HeatEqnSolver('2', u0(R, Z), dt, Nr, Nz, tf, interval, tolerance, r0);
    % errors(k, 2) = dr*dz*(sum(sum(abs(U2 - U_exact))));

    [U3, ranks3] = HeatEqnSolver('3', u0(R, Z), dt, Nr, Nz, tf, interval, tolerance, r0);
    % errors(k, 3) = dr*dz*(sum(sum(abs(U3 - U_exact))));
end

figure(1); clf; surf(R, Z, U1);
colorbar;
shading flat; % removes gridlines
legend(sprintf('N_r = %s, N_z = %s', num2str(Nr, 3), num2str(Nz, 3)), 'Location','northwest');
xlabel('R'); ylabel('Z'); zlabel('U(R, Z)'); title([sprintf('RAIL approximation at time %s', num2str(tf, 4))]);

figure(2); clf;
plot(ranks1(:, 1), ranks1(:, 2), 'black-', 'LineWidth', 1.5); hold on;
plot(ranks2(:, 1), ranks2(:, 2), 'blue-', 'LineWidth', 1.5);
plot(ranks3(:, 1), ranks3(:, 2), 'green-', 'LineWidth', 1.5);
xlabel('time'); ylabel('rank'); title('Rank plot over time');
legend('Backward Euler', 'DIRK2', 'DIRK3');

















