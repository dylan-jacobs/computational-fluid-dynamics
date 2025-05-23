% Tests the fluid solver for the FP plasma system
clc; clear variables; close all;

Nx = 100;
Vr = 50;
Vz = 100;
x_min = 0;
x_max = 200;
interval = [x_min, x_max, 0, 8, -8, 10]; % 1D in x, 2D in v (first: r interval, second: z interval)

dt = 0.3;
tf = 300;
tvals = [0, 5e-3:dt:tf];
if tvals(end) ~= tf
    tvals = [tvals, tf];
end

ma = 1; % ion mass
me = 1/1836; % electron mass
qa = 1; % ion charge
qe = -1; % electron charge
R_const = 1/ma; % gas constant

[n, u_para, T_ae] = get_boundaries(); % electron and ion temps equal at t=0 and at boundaries --> T_ae = T_a = T_e
[xvals, v_perp, v_para, dx, dv_perp, dv_para] = GetRZ(Nx, Vr, Vz, interval);

n = n(xvals);
u_para = u_para(xvals);
T_ae = T_ae(xvals); Ta = T_ae; Te = T_ae;
U = 0.5*((3*Ta/ma) + (u_para.^2));

% ---- plot ICs ----
% figure; clf;
% plot(xvals, n); hold on;
% plot(xvals, u_para0./u_para0(1));
% plot(xvals, T_ae);
% legend('n', 'u/u0', 'T_e', 'T_\alpha')
% title('tn=0');
% xlabel('x');
% drawnow;
% pause(0.05);

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
tic
figure; clf;
for tn = 2:numel(tvals)
    disp(['t = ', num2str(tvals(tn))])
    dt = tvals(tn) - tvals(tn-1);
    [n, nu_para, nU, Te] = newton_solver_FP(f, n, u_para, U, Te, dt, dx, dv_para, dv_perp, v_para, v_perp, qa, qe, ma, me, R_const, x_min, x_max);

    u_para = nu_para./n;
    U = nU./n;

    % get T_a from macroscopic parameters
    Ta = ((2*U) - (u_para.^2)).*ma/3;
   
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
legend('n', 'u/u0', 'T_e', 'T_\alpha')
toc

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

