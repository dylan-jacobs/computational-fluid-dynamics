% solves 1D Fokker-Planck system
%% TEST 1 -- Accuracy
clc; clear variables; close all;

% Fokker-Planck parameters
u = 0; % advection velocity
B = @(vx, t) (vx - u); 
D = @(vx) vx.^0;
tf = 10;
interval = [-pi, pi];
Nx = 100;
[xvals, dx] = GetXVals(Nx, interval);
B_max = max(0.5*(B([xvals(2:end); xvals(end)+dx], tf).^2) - 0.5*(B([xvals(1)-dx; xvals(1:end-1)], tf).^2));
D_max = max(D(xvals));

% Maxwellian parameters
R = 1/6;
n1 = 2.427991645004178;
T1 = 1.142108881605655;
vxbar1 = -0.3092290436287916;

n2 = 0.7136010085856154;
T2 = 0.727329792278731;
vxbar2 = 1.052136313276047;

f_M1 = @(vx) (n1/sqrt(2*pi*R*T1)).*exp(-((vx-vxbar1).^2)./(2*R*T1));
f_M2 = @(vx) (n2/sqrt(2*pi*R*T2)).*exp(-((vx-vxbar2).^2)./(2*R*T2));

f0 = @(vx) f_M1(vx) + f_M2(vx);
f_exact = @(vx) (pi/sqrt(2*pi*R*3)).*exp(-((vx).^2)./(2*R*3));
f_exact = f_exact(xvals);

lambdavals = 0.5;%(0.5:0.1:1)';

for k = 1:numel(lambdavals)
    dt = lambdavals(k)*(dx^2)/(2*((B_max*dx) + D_max));
    disp([num2str(k), '/', num2str(numel(lambdavals))]);

    [f, data] = FokkerPlanckSolver('1', f0(xvals), dt, Nx, tf, interval, B, D, true, f_exact);
    % errors(k, 1) = dx*(sum(abs(f - f_exact))); % L1 error
end

l1 = data(:, 1);
positivity = data(:, 2);
relative_entropy = data(:, 3);
mass = data(:, 4);
tvals = data(:, 5);

figure(1); clf; plot(xvals, f, 'blue-', 'LineWidth', 1.5);
legend(sprintf('N_x = %s', num2str(Nx, 3)), 'Location','northwest');
xlabel('v_x'); ylabel('f(v_x)'); title([sprintf('Forward Euler approximation of 1D Fokker-Planck system at time %s', num2str(tf, 4))]);

figure(2); clf; plot(xvals, f_exact, 'blue-', 'LineWidth', 1.5);
xlabel('v_x'); ylabel('f(v_x)'); title([sprintf('f_{exact} at time t=%s', num2str(tf, 4))]);

% l1 decay
figure(3); clf; semilogy(tvals, l1, 'black-', 'LineWidth', 1.5);
xlabel('t'); ylabel('L_1(f(V_x))'); title('L_1 decay of numerical solution over time');

% Positivity
figure(4); clf; plot(tvals, positivity, 'green-', 'LineWidth', 1.5);
xlabel('t'); ylabel('min(f(V_x))'); title('Minimum values of numerical solution over time');

% Relative entropy
figure(5); clf; plot(tvals, relative_entropy, 'magenta-', 'LineWidth', 1.5);
xlabel('t'); ylabel('Relative entropy'); title('Relative entropy of numerical solution over time');

% Mass
figure(6); clf; plot(tvals, abs(mass-mass(1))/mass(1), 'red-', 'LineWidth', 1.5);
xlabel('t'); ylabel('relative mass'); title('Relative mass of numerical solution over time');
