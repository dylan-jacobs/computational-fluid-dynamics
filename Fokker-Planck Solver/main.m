% solves 1D Fokker-Planck system
%% TEST 1 -- Accuracy
clc; clear variables; close all;

% Fokker-Planck parameters
u = 1; % advection velocity
C = @(vx, t) (vx - u); 
D = @(vx) vx.^0;
tf = 3;
interval = [-pi, 2*pi];
Nx = 100;
X = linspace(interval(1), interval(2), Nx + 1)';
dx = X(2) - X(1);
C_max = max(0.5*(C(X).^2));
D_max = max(D(X));

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
f_exact = exp(-X.^2);

lambdavals = 1;%(0.1:0.1:1)';
errors = zeros(numel(lambdavals), 1); % L1 norm for F. Euler

for k = 1:numel(lambdavals)
    dt = lambdavals(k)*(dx^2)/(2*((C_max*dx) + D_max));
    disp([num2str(k), '/', num2str(numel(lambdavals))]);

    [f] = FokkerPlanckSolver('1', f0(X), dt, Nx, tf, interval, C, D);
    errors(k, 1) = dx*(sum(abs(f - f_exact))); % L1 error
end

figure(1); clf; plot(X, f);
legend(sprintf('N_x = %s', num2str(Nx, 3)), 'Location','northwest');
xlabel('V_x'); ylabel('f(V_x)'); title([sprintf('Forward Euler approximation of 1D Fokker-Planck system at time %s', num2str(tf, 4))]);

figure(2); clf; plot(X, f_exact);
xlabel('V_x'); ylabel('f(V_x)'); title([sprintf('f_{exact} at time t=%s', num2str(tf, 4))]);

figure(3); clf;
loglog(lambdavals, errors(:, 1), 'black-', 'LineWidth', 1.5); hold on; % B. Euler
loglog(lambdavals, 0.08*lambdavals, 'black--', 'LineWidth', 1); % Order 1
title(sprintf('Forward Euler Temporal Convergence at tf=%s, Nx=%s', num2str(tf), num2str(Nx))); xlabel('\lambda'); ylabel('L1 Error');
legend('Forward Euler', 'Order 1');
