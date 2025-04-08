% Tests the fluid solver for the FP plasma system
clc; clear variables; close all;

Nx = 80;
Vx = 120;
x_min = 0;
x_max = 200;
interval = [x_min, x_max, -8, 10, 0, 8]; % 1D in x, 2D in v

dt = 5e-3;
tf = 150;
tvals = [0, dt, dt+0.3:0.3:tf];
if tvals(end) ~= tf
    tvals = [tvals, tf];
end

R_const = 1/6; % gas constant
ma = 1; % ion mass
me = 1/1836; % electron mass
qa = 1; % ion charge
qe = -1; % electron charge

[n, u_para, T_ae] = get_boundaries(); % electron and ion temps equal at t=0 and at boundaries --> T_ae = T_a = T_e
[xvals, v_para, v_perp, dx, dv_para, dv_perp] = GetXV(Nx, Vx, Vx, interval);

n = n(xvals);
u_para = u_para(xvals);
T_ae = T_ae(xvals); T_a = T_ae; Te = T_ae;
U = 0.5*((3*T_a/ma) + (u_para.^2));

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

for tn = 1:numel(tvals)
    disp(['t = ', num2str(tvals(tn))])
    [n, nu_para, nU, Te] = newton_solver_FP(f, n, u_para, U, Te, dt, dx, dv_para, dv_perp, v_para, v_perp, qa, qe, ma, me, R_const, x_min, x_max);

    u_para = nu_para./n;
    U = nU./n;

    % get T_a from macroscopic parameters
    T_a = (2*U - u_para.^2).*ma/3;
   
    nvals(tn, :) = n;
    uvals(tn, :) = u_para;
    Te_vals(tn, :) = Te;
    Ta_vals(tn, :) = T_a;

    % reconstruct f
    f = maxwellian(n, v_para, v_perp, u_para, T_a, R_const);

    if mod(tvals(tn), 25) < 0.25
        figure; clf;
        plot(xvals, nvals(tn, :), "LineWidth",1.5); hold on;
        plot(xvals, uvals(tn, :)./uvals(tn, 1), "LineWidth",1.5);
        plot(xvals, Te_vals(tn, :), "LineWidth",1.5);
        plot(xvals, Ta_vals(tn, :), "LineWidth",1.5);
        legend('n', 'u/u0', 'T_e', 'T_\alpha')
        title(['tn=', num2str(tvals(tn))]);
        xlabel('x');
        drawnow;
        pause(0.05);

        % figure;
        % surf(f(:, :, 1));
        % drawnow;
        % pause(0.05);
    end
end
