% Compute 2D Finite Differences Problems for Conservation Laws
clear variables; close all; clc;

tf = 60; % final time
lambda = 0.3; % dt = lambda*dx
vmax = 5; % max velocity
interval = [0, 4*pi, -vmax, vmax]; % xmin, xmax, vmin, vmax
Nx = 64; Nv = 128;
discretizationType = 'RK3'; % time discretization type (RK1, RK2, RK3, RK4)

test = input('Enter test: \n 1) Weak Landau Damping \n 2) Strong Landau Damping \n 3) 2-Stream Instability \n');

switch test
    case 1 % Weak Landau Damping Test
        alpha = 0.01;
        k = 0.5;
        f0 = @(x, v) (1/sqrt(2*pi)).*(1+(alpha.*cos(k.*x))).*exp(-0.5*(v.^2));

    case 2 % Strong Landau Damping Test
        alpha = 0.5;
        k = 0.5;
        f0 = @(x, v) (1/sqrt(2*pi)).*(1+(alpha.*cos(k.*x))).*exp(-0.5*(v.^2));

    case 3 % 2-Stream Instability Test
        alpha = 0.05;
        k = 0.5;
        f0 = @(x, v) (1/sqrt(2*pi)).*(1+(alpha.*cos(k.*x))).*exp(-0.5*(v.^2)).*(v.^2);
end

[X, V, ~, ~] = GetXY(Nx, Nv, interval);
[f_matrix, EF, mass, L1, L2, energy, entropy, tvals] = VPSolver(discretizationType, Nx, Nv, lambda, interval, tf, f0);
%%
figure(1);
semilogy(tvals, EF, 'LineWidth', 1.5);
xlabel('time'); ylabel('E_2'); title('L_2 norm of Electric Field over Time');
%%
figure(2);
plot(tvals, mass, 'LineWidth', 1.5);
xlabel('time'); ylabel('E_2'); title('Mass over Time');
%%
figure(3);
plot(tvals, L1, 'LineWidth', 1.5);
xlabel('time'); ylabel('L_1'); title('L_1 norm of Probability Distribution f over Time');

figure(4);
plot(tvals, L2, 'LineWidth', 1.5);
xlabel('time'); ylabel('L_2'); title('L_2 norm of Probability Distribution f over Time');

figure(5);
plot(tvals, energy, 'LineWidth', 1.5);
xlabel('time'); ylabel('Energy'); title('Energy over Time');

figure(6);
plot(tvals, entropy, 'LineWidth', 1.5);
xlabel('time'); ylabel('Entropy'); title('Entropy over Time');

figure(7); clf; surf(X, V, f_matrix(:, :, end));
colorbar;
shading flat; % removes gridlines
legend(sprintf('N_x = %s, N_v = %s', num2str(Nx, 3), num2str(Nv, 3)), 'Location','northwest');
xlabel('X'); ylabel('V'); zlabel('F(X, V)'); title([sprintf('2D WENO+%s', discretizationType), sprintf(' approximation at time %s', num2str(tf, 4))]);
view(2); % bird's eye view
xlim([interval(1), interval(2)]);

















