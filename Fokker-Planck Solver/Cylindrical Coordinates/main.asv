% solves 1D Fokker-Planck system
%% TEST 1 -- Accuracy
clc; clear variables; close all;

% Fokker-Planck parameters
u = 0; % advection velocity
B = @(vr, t) (vr - u); 
D = @(vr) 0.5*(vr.^0);
tf = 5;
interval = [0, 14, -16, 16];
Nr = 200; Nz = 200;
[rvals, zvals, dr, dz] = GetRZ(Nr, Nz, interval);
B_max = max(0.5*(B([rvals(1:end, 1); rvals(end, 1)+dr], tf).^2) - 0.5*(B([rvals(1, 1)-dr; rvals(1:end, 1)], tf).^2));
D_max = max(D(rvals(:, 1)));

tolerance = 1e-6;

% Maxwellian parameters
R = 1;
n1 = 1.902813990281176;
T1 = 1.10608227396894;
u_vec1 = [0, -2.113532196926305];

n2 = 1.238778663308618;
T2 = 0.1087698976066122;
u_vec2 = [0, 3.246470699196378];

f_M = @(n, R, T, u_vec, vr, vz) (n/((2*pi*R*T)^(3/2))).*exp(-((vr - u_vec(1)).^2 + ((vz - u_vec(2)).^2))./(2*R*T));
f_M1 = f_M(n1, R, T1, u_vec1, rvals, zvals);
f_M2 = f_M(n2, R, T2, u_vec2, rvals, zvals);

f0 = f_M1 + f_M2;
f_exact = @(vr, vz) (pi/(2*pi*R*3)^(3/2)).*exp(-((vr.^2) + (vz.^2))./(2*R*3));
f_exact = f_exact(rvals, zvals);

lambdavals = 0.5;%(0.5:0.1:1)';

% figure(1); clf; surf(rvals, zvals, f0);shading interp;
% legend(sprintf('N_r = %s', num2str(Nr, 3)), 'Location','northwest');
% xlabel('v_r'); ylabel('v_z'); title([sprintf('Forward Euler approximation of 0D2V Fokker-Planck system at time %s', num2str(tf, 4))]);
% return

for k = 1:numel(lambdavals)
    dt = lambdavals(k)*(dr^2)/(2*((B_max*dr) + D_max));
    disp([num2str(k), '/', num2str(numel(lambdavals))]);

    [f, data] = FokkerPlanckSolver('1', f0, dt, Nr, Nz, tf, interval, B, D, true, f_exact, tolerance);
    % errors(k, 1) = dx*(sum(abs(f - f_exact))); % L1 error
end

l1 = data(:, 1);
positivity = data(:, 2);
relative_entropy = data(:, 3);
mass = data(:, 4);
tvals = data(:, 5);

figure(1); clf; surf(rvals, zvals, f, 'blue-', 'LineWidth', 1.5);
legend(sprintf('N_r = %s', num2str(Nr, 3)), 'Location','northwest');
xlabel('v_r'); ylabel('f(v_z)'); title([sprintf('Forward Euler approximation of 0D2V Fokker-Planck system at time %s', num2str(tf, 4))]);

figure(2); clf; surf(rvals, zvals, f_exact, 'blue-', 'LineWidth', 1.5);
xlabel('v_r'); ylabel('f(v_z)'); title([sprintf('f_{exact} at time t=%s', num2str(tf, 4))]);

% l1 decay
figure(3); clf; semilogy(tvals, l1, 'black-', 'LineWidth', 1.5);
xlabel('t'); ylabel('L_1(f(V_r, V_z))'); title('L_1 decay of numerical solution over time');

% Positivity
figure(4); clf; plot(tvals, positivity, 'green-', 'LineWidth', 1.5);
xlabel('t'); ylabel('min(f(V_r, V_z))'); title('Minimum values of numerical solution over time');

% Relative entropy
figure(5); clf; plot(tvals, relative_entropy, 'magenta-', 'LineWidth', 1.5);
xlabel('t'); ylabel('Relative entropy'); title('Relative entropy of numerical solution over time');

% Mass
figure(6); clf; plot(tvals, abs(mass-mass(1))/mass(1), 'red-', 'LineWidth', 1.5);
xlabel('t'); ylabel('relative mass'); title('Relative mass of numerical solution over time');
