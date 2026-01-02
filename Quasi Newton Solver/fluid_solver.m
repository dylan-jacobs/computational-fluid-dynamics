% Tests the fluid solver for the FP plasma system
% clc; clear variables; close all;

% Simulation parameters
Nx = 80;
Vr = 50;
Vz = 100;
x_min = 0;
x_max = 200;
interval = [x_min, x_max, 0, 8, -8, 8]; % for standing shock
% interval = [x_min, x_max, 0, 12, -12, 12]; % for smooth IC % 1D in x, 2D in v (first: r interval, second: z interval)
tf = 1;

% Constants
ma = 1; % ion mass
me = 1/1836; % electron mass
qa = 1; % ion charge
qe = -1; % electron charge
R_const = 1/ma; % gas constant


%%% COMPUTE REFERENCE SOLUTION %%%
disp('Computing reference solution...');
dt_ref = 0.001;
tvals_ref = 0:dt_ref:tf;
if tvals_ref(end) ~= tf
    tvals_ref = [tvals_ref, tf];
end
[n, u_para, T_ae] = get_boundaries(); % electron and ion temps equal at t=0 and at boundaries --> T_ae = T_a = T_e
[xvals, v_perp, v_para, dx, dv_perp, dv_para] = GetRZ(Nx, Vr, Vz, interval);

n = n(xvals);
u_para = u_para(xvals);
T_ae = T_ae(xvals); Ta = T_ae; Te = T_ae;

f_ref = maxwellian(n, v_para, v_perp, u_para, T_ae, R_const);

% store macroscopic parameters at every timestep for all xvals
nvals_ref = zeros(numel(tvals_ref), Nx);
uvals_ref = zeros(numel(tvals_ref), Nx);
Te_vals_ref = zeros(numel(tvals_ref), Nx);
Ta_vals_ref = zeros(numel(tvals_ref), Nx);

% ---- set initial vals ----
nvals_ref(1, :) = n;
uvals_ref(1, :) = u_para;
Te_vals_ref(1, :) = Te;
Ta_vals_ref(1, :) = Ta;

for tn = 2:numel(tvals_ref)
    disp(tvals_ref(tn));
    dt = tvals_ref(tn) - tvals_ref(tn-1);
    [n, u_para, Ta, Te] = newton_solver_FP_order_2(f_ref, n, u_para, Ta, Te, dt, dx, dv_para, dv_perp, v_para, v_perp, qa, qe, ma, me, R_const, x_min, x_max);
   
    nvals_ref(tn, :) = n;
    uvals_ref(tn, :) = u_para;
    Te_vals_ref(tn, :) = Te;
    Ta_vals_ref(tn, :) = Ta;

    % reconstruct f
    f_ref = maxwellian(n, v_para, v_perp, u_para, Ta, R_const);  
end

%%

%%% RUN CONVERGENCE TEST %%%
dtvals = [0.2, 0.1, 0.05, 0.025];%, 0.0125];
errors = zeros(numel(dtvals), 6); % compute errors for f, n, u, Ta, Te
for i = 1:numel(dtvals)

dt = dtvals(i)
tvals = [0, 5e-3:dt:tf];
if tvals(end) ~= tf
    tvals = [tvals, tf];
end

[n, u_para, T_ae] = get_boundaries(); % electron and ion temps equal at t=0 and at boundaries --> T_ae = T_a = T_e
[xvals, v_perp, v_para, dx, dv_perp, dv_para] = GetRZ(Nx, Vr, Vz, interval);

n = n(xvals);
u_para = u_para(xvals);
T_ae = T_ae(xvals); Ta = T_ae; Te = T_ae;

f = maxwellian(n, v_para, v_perp, u_para, T_ae, R_const);

% store macroscopic parameters at every timestep for all xvals
nvals = zeros(numel(tvals), Nx);
uvals = zeros(numel(tvals), Nx);
Te_vals = zeros(numel(tvals), Nx);
Ta_vals = zeros(numel(tvals), Nx);

% ---- set initial vals ----
nvals(1, :) = n;
uvals(1, :) = u_para;
Te_vals(1, :) = Te;
Ta_vals(1, :) = Ta;

% figure; clf;
for tn = 2:numel(tvals)
    disp(['t = ', num2str(tvals(tn))])
    dt = tvals(tn) - tvals(tn-1);
    [n, u_para, Ta, Te] = newton_solver_FP(f, n, u_para, Ta, Te, dt, dx, dv_para, dv_perp, v_para, v_perp, qa, qe, ma, me, R_const, x_min, x_max);
   
    nvals(tn, :) = n;
    uvals(tn, :) = u_para;
    Te_vals(tn, :) = Te;
    Ta_vals(tn, :) = Ta;

    % reconstruct f
    f = maxwellian(n, v_para, v_perp, u_para, Ta, R_const);

    plot_freq = 50;
    num_subplots = (tf / plot_freq) + 1;
    num_rows = ceil(sqrt(num_subplots));
    num_cols = ceil(num_subplots / num_rows);
    if mod(tvals(tn), plot_freq) < 0.25
        subplot(num_rows, num_cols, floor(tvals(tn) / plot_freq)+1);
        plot(xvals, nvals(tn, :), "LineWidth", 1.5); hold on;
        plot(xvals, uvals(tn, :)./uvals(tn, 1), "LineWidth",1.5);
        plot(xvals, Te_vals(tn, :), "LineWidth", 1.5);
        plot(xvals, Ta_vals(tn, :), "LineWidth", 1.5);
        title(['tn=', num2str(tvals(tn))]);
        ylim([0, 1.5]);
        pause(0.05);
    end

      
end

% legend('n', 'u/u0', 'T_e', 'T_\alpha')
% nvals = 2*pi*dv_perp*dv_para*sum(sum(f .* v_perp, 1), 2);
% nvals_ref = 2*pi*dv_perp*dv_para*sum(sum(f_ref .* v_perp, 1), 2);
errors(i, 1) = dx*2*pi*dv_perp*dv_para*sum(sum(sum(abs(f - f_ref) .* v_perp)))/(200*18*8);
errors(i, 2) = dx*sum(abs(nvals(end, :) - nvals_ref(end, :)))/200;
errors(i, 3) = dx*sum(abs(uvals(end, :) - uvals_ref(end, :)))/200;
errors(i, 4) = dx*sum(abs(Ta_vals(end, :) - Ta_vals_ref(end, :)))/200;
errors(i, 5) = dx*sum(abs(Te_vals(end, :) - Te_vals_ref(end, :)))/200;
errors(i, 6) = max(max(max(abs(f - f_ref))));
end
%%
figure; clf;
plot(xvals, nvals(tn, :), "LineWidth", 1.5); hold on;
plot(xvals, uvals(tn, :)./uvals(tn, 1), "LineWidth",1.5);
plot(xvals, Te_vals(tn, :), "LineWidth", 1.5);
plot(xvals, Ta_vals(tn, :), "LineWidth", 1.5);
title(['tn=', num2str(tvals(tn))]);
ylim([0, 1.5]);
legend('n', 'u/u0', 'T_e', 'T_\alpha')

% for i = 1:21
%     idx = i*500;
%     figure; clf;
%     plot(xvals, nvals(idx, :), "LineWidth",1.5); hold on;
%     plot(xvals, uvals(idx, :)./uvals(idx, 1), "LineWidth",1.5);
%     plot(xvals, Te_vals(idx, :), "LineWidth",1.5);
%     plot(xvals, Ta_vals(idx, :), "LineWidth",1.5);
%     title(['tn=', num2str(tvals(idx))]);
%     ylim([0, 1.8]);
%     legend('n', 'u/u0', 'T_e', 'T_\alpha')
% end


%%% CONVERGENCE %%%
disp('Errors:');
disp('L1 (global)');
disp(errors(:, 1));
disp(log2(errors(1:end-1, 1)./errors(2:end, 1)));
disp('Mass');
disp(errors(:, 2));
disp(log2(errors(1:end-1, 2)./errors(2:end, 2)));
disp('Momentum');
disp(errors(:, 3));
disp(log2(errors(1:end-1, 3)./errors(2:end, 3)));
disp('Ion Temperature');
disp(errors(:, 4));
disp(log2(errors(1:end-1, 4)./errors(2:end, 4)));
disp('Electron Temperature');
disp(errors(:, 5));
disp(log2(errors(1:end-1, 5)./errors(2:end, 5)));
disp('Linf (global)');
disp(errors(:, 6));
disp(log2(errors(1:end-1, 6)./errors(2:end, 6)));