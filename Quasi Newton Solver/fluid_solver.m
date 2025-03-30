% Tests the fluid solver for the FP plasma system
clc; clear variables; close all;

Nx = 80;
Vx = 120;
x_min = 0;
x_max = 200;
interval = [x_min, x_max, -8, 10, 0, 8]; % 1D in x, 2D in v

dt = 5e-3;
tf = 25;
tvals = 0:dt:tf;

R = 1/6;
ma = 1;
me = 1/8136;
qa = 1;
qe = -1;

[n, u_para, T_ae] = get_boundaries(); % electron and ion temps equal at t=0 and at boundaries --> T_ae = T_a = T_e
[xvals, v_para, v_perp, dx, dv_para, dv_perp] = GetXV(Nx, Vx, Vx, interval);

n = n(xvals);
u_para = u_para(xvals);
T_ae = T_ae(xvals); Ta = T_ae; Te = T_ae;
U = 0.5*(3*Ta./ma) + (u_para.^2);

f = maxwell(n, v_para, v_perp, u_para, T_ae, R);


for tn = 1:numel(tvals)

    [n, nu_para, nU, T_e] = newton_solver_FP(f, n, u_para, U, Te, dt, dx, dv_para, dv_perp, v_para, v_perp, qa, qe, ma, me, R, x_min, x_max);

    u_para = nu_para./n;
    U = nU./n;
    Ta = (2*U - u_para.^2).*ma/3;

    % reconstruct f
    f = maxwell(n, v_para, v_perp, u_para, T_a, R);

end
