clc; clear variables; close all;

% load reference solution 
soln = load('vdfp_refsoln_imex111.mat');
soln = soln.f_master;

DTvals = [0.4,0.2,0.1,0.05,0.025,0.0125,0.00625];
DTvals = [0.4,0.2,0.1,0.05,0.025,0.0125];
% DTvals = [0.4, 0.2, 0.1, 0.05];
% DTvals = [0.4, 0.2, 0.1];
% DTvals = [0.025];
Nx = 80;
moments_ref_reimann = zeros(numel(DTvals), Nx, 3);
ERRORS = zeros(5 + Nx, numel(DTvals));

for dtIndex = 1:numel(DTvals)
dt = DTvals(dtIndex);

% initial rank
r0 = 30;
% time-stepping method: 1=B.Euler, 2=DIRK2
method = '1';
tolerance = 1e-6;

% mesh parameters
% Nx = 80;
Nr = 50;
Nz = 100;
x_min = 0;
x_max = 200;
interval = [x_min, x_max, 0, 8, -8, 8]; % 1D in x, 2D in v (x interval, r interval, z interval)
tf = 10;

ma = 1; % ion mass
me = 1/1836; % electron mass
qa = 1; % ion charge
qe = -1; % electron charge
R_const = 1; % gas constant

[n_vals, u_para, Tae_vals] = get_boundaries(x_min, x_max); % electron and ion temps equal at t=0 and at boundaries --> T_ae = T_a = T_e
n_IC = n_vals; u_IC = u_para; Te_IC = Tae_vals;
[xvals, Rmat, Zmat, dx, dr, dz] = GetRZ(Nx, Nr, Nz, interval);
leftBC = maxwellian(n_vals(x_min - dx/2), Rmat, Zmat, u_para(x_min - dx/2), Tae_vals(x_min - dx/2), R_const);
rightBC = maxwellian(n_vals(x_max + dx/2), Rmat, Zmat, u_para(x_max + dx/2), Tae_vals(x_max + dx/2), R_const);
rvals = Rmat(:, 1);
zvals = Zmat(1, :)';
[Xmat, Zmat2] = meshgrid(xvals, zvals); Xmat=Xmat'; Zmat2=Zmat2'; % for plotting only

n_vals = n_vals(xvals);
u_perp = zeros(Nx, 1);
u_para = u_para(xvals);
Tae_vals = Tae_vals(xvals); Ta_vals = Tae_vals; Te_vals = Tae_vals;
U = 0.5*((3*Ta_vals/ma) + (u_para.^2));

% initialize f_inf using QCM at each spatial node
f_vals = maxwellian(n_vals, Rmat, Zmat, u_para, Tae_vals, R_const);
f_inf_vals = zeros(Nr, Nz, Nx);
for spatialIndex = 1:Nx
    % f0 = @(vr,vz) f_M(vr,vz,n1,u1r,u1z,T1,R_const) + f_M(vr,vz,n2,u2r,u2z,T2,R_const); % IC
    % f = f0(Rmat,Zmat);
    f0 = f_vals(:, :, spatialIndex);
    
    % discrete moments of f0
    rho0 = sum(sum(f0.*Rmat))*2*pi*dr*dz;
    Jz0  = sum(sum(f0.*Rmat.*Zmat))*2*pi*dr*dz;
    kappa0   = sum(sum(f0.*Rmat.*((Rmat.^2 + Zmat.^2)/2)))*2*pi*dr*dz;
    
    [f_inf, n_inf, uz_inf, T_inf] = QCM(rho0, Jz0, kappa0, R_const, Rmat, Zmat);
    f_inf_vals(:, :, spatialIndex) = f_inf;
        
    % moments at equilibrium
    rhoM = rho0;
    JzM = Jz0;
    kappaM = kappa0;
    
    ur_inf = 0; % no drift in r
    ur = ur_inf;
    uz = uz_inf;
end

% dt = 3/((1/dr) + (1/dz));
% dt = 0.3;
tvals = [0, 5e-3:dt:tf]';
if tvals(end) ~= tf
    tvals = [tvals; tf];
end
Nt = numel(tvals);

Ar = @(u, w, t) w; %cell centers
Br = @(u, w, t) w.*(w - u); %evaluated on cell boundaries
Cr = @(u, w, t, T, m) (T/m)*w; %evaluated cell boundaries
Az = @(u, w, t) w.^0;
Bz = @(u, w, t) w - u;
Cz = @(u, w, t, T, m) (T/m)*w.^0;

% init bases
f_vals_low_rank = cell(Nx, 3); % store Vr, S, Vz at each spatial node
for spatialIndex = 1:Nx
    f = f_vals(:, :, spatialIndex);
    [Vr, S, Vz] = svd2(f, rvals);
    r0 = min(r0, size(Vr, 2));
    Vr = Vr(:, 1:r0); S = S(1:r0, 1:r0); Vz = Vz(:, 1:r0);
    f_vals_low_rank{spatialIndex, 1} = Vr;
    f_vals_low_rank{spatialIndex, 2} = S;
    f_vals_low_rank{spatialIndex, 3} = Vz;  
end

% compute LoMaC parameters once
wr = exp(-(rvals.^2));
wz = exp(-(zvals.^2));
c = (dr*sum(rvals.^2.*wr.*rvals))/(dr*sum(wr.*rvals)) + (dz*sum(zvals.^2.*wz))/(dz*sum(wz));
w_norm_1_squared = 2*pi*dr*dz*sum(rvals .* wr)*sum(wz);
w_norm_v_squared = 2*pi*dr*dz*sum(rvals .* wr)*sum(zvals.^2 .* wz);
w_norm_v2_squared = 2*pi*dr*dz*sum(sum((Rmat.^2 + Zmat.^2 - c).^2 .* exp(-Rmat.^2 - Zmat.^2) .* Rmat));

% store rank, mass, momentum, energy, l1 decay, etc...
mass = zeros(Nt, 1);
momentum = zeros(Nt, 1);
energy = zeros(Nt, 1);
min_vals = zeros(Nt, 1);
ranks = zeros(Nt, 1);

mass(1) = dx*ma*sum(n_vals);
momentum(1) = dx*ma*sum(n_vals .* u_para);
energy(1) = dx*ma*sum(0.5*((3/ma)*n_vals.*Ta_vals + n_vals.*u_para.^2) + (3/2)*n_vals.*Te_vals);
min_vals(1) = min(min(min(f_vals)));
ranks(1) = r0;

flux_diff_n = 0;
flux_diff_nu = 0;
flux_diff_T = 0;

% time-stepping loop
for t_index = 2:Nt
    tval = tvals(t_index);
    disp(tval);
    dt = tval - tvals(t_index-1);
    [n_vals, u_para, Ta_vals, Te_vals, u_hat, nu_hat, S_hat, Q_hat, nTe_hat, kappaTx] = fluid_solver_IMEX111(f_vals_low_rank, n_vals, u_para, Ta_vals, Te_vals, dt, dx, dr, dz, Rmat, Zmat, qa, qe, ma, me, R_const, x_min, x_max);
    % [nvals,upara,Tavals,Tevals,nuhat,Shat,Qhat,uhat,nTehat,KTxhat] = QN(f_vals_low_rank,n_vals,u_para,Ta_vals,Te_vals,qa,qe,ma,me,Nx,dx,dt,leftBC,rightBC,dr,dz,zvals,Rmat,Zmat,n_IC,u_IC,Te_IC);
    % u_para = nu ./ n_vals;
    % U = nU ./ n_vals;
    % Ta_vals = ((2*U) - (u_para.^2)).*ma/3;
    nu = n_vals .* u_para;
    nU = 0.5*((3/ma)*n_vals.*Ta_vals + n_vals.*u_para.^2);

    f_vals_low_rank_temp = cell(Nx, 3);
    % nvals_diff = zeros(Nx,1);
    % uvals_diff = zeros(Nx,1);
    % Tivals_diff = zeros(Nx,1);

    if abs(tval - 1e-5) < dt || abs(tval - 200) < dt
        figure(1); clf;
        plot(xvals, n_vals, "LineWidth", 1.5); hold on;
        plot(xvals, u_para./u_para(1), "LineWidth",1.5);
        plot(xvals, Te_vals, "LineWidth", 1.5);
        plot(xvals, Ta_vals, "LineWidth", 1.5);
        legend('n', 'u_{||}', 'T_e', 'T_a');
        title(sprintf('Mass, momentum, ion temperature, and electron temperature at t=%s', num2str(tval)));
        ylim([0, 1.2]);
        saveas(gcf, sprintf('Plots/moments_time_%s.m', num2str(round(tval))));
        pause(0.05);
    end

    for spatialIndex = 1:Nx
        Vr = f_vals_low_rank{spatialIndex, 1};
        S = f_vals_low_rank{spatialIndex, 2};
        Vz = f_vals_low_rank{spatialIndex, 3};

        rhoM = n_vals(spatialIndex);
        JzM  = nu(spatialIndex);
        kappaM = nU(spatialIndex);
        % i = spatialIndex;
        % JzM = n_vals(i)*u_para(i);
        % kappaM = 0.5*(3*n_vals(i)*Ta_vals(i)/ma + n_vals(i)*u_para(i)^2);

        switch(method)
            case '1'
                [Vr, S, Vz, rank] = IMEX111(Vr, S, Vz, u_perp, u_para, dt, dx, tval, rvals, zvals, x_min, x_max, f_vals_low_rank, Rmat, Zmat, ma, me, qa, qe, Ar, Az, Br, Bz, Cr, Cz, tolerance, n_vals, Ta_vals, Te_vals, rhoM, JzM, kappaM, spatialIndex, R_const, leftBC, rightBC, wr, wz, c, w_norm_1_squared, w_norm_v_squared, w_norm_v2_squared);
        end
        f_vals_low_rank_temp{spatialIndex, 1} = Vr;
        f_vals_low_rank_temp{spatialIndex, 2} = S;
        f_vals_low_rank_temp{spatialIndex, 3} = Vz;  

        % nvals_diff(spatialIndex) = 2*pi*dr*dz*sum(sum((Vr*S*Vz').*Rmat));
        % uvals_diff(spatialIndex) = (1/nvals_diff(spatialIndex))*2*pi*dr*dz*sum(sum(Zmat.*(Vr*S*Vz').*Rmat));
        % Tivals_diff(spatialIndex) = (2*2*pi*dr*dz*sum(sum(((Rmat.^2+Zmat.^2)/2).*(Vr*S*Vz').*Rmat)) - nvals_diff(spatialIndex)*uvals_diff(spatialIndex)^2)/(3*nvals_diff(spatialIndex)/ma);
    end

    % update f_vals simultaneously
    f_vals_low_rank = f_vals_low_rank_temp;

    % update mass, momentum, energy
    flux_diff_n = flux_diff_n + dt*(nu_hat(end) - nu_hat(1));
    flux_diff_nu = flux_diff_nu + dt*((S_hat(end) - S_hat(1)) - (qa/(2*qe*ma))*(n_IC(x_max + dx/2)*Te_IC(x_max + dx/2) + n_vals(end)*Te_vals(end) - n_vals(1)*Te_vals(1) - n_IC(x_min - dx/2)*Te_IC(x_min - dx/2)));
    flux_diff_T = flux_diff_T + dt*((Q_hat(end) - Q_hat(1)) + 2.5*(u_hat(end)*nTe_hat(end) - u_hat(1)*nTe_hat(1)) - (kappaTx(end) - kappaTx(1)));
                                
    mass(t_index) = dx*ma*sum(n_vals) + flux_diff_n;
    momentum(t_index) = dx*ma*sum(u_para .* n_vals) + flux_diff_nu;
    energy(t_index) = dx*ma*sum(0.5*((3/ma)*n_vals.*Ta_vals + n_vals.*u_para.^2) + (3/2)*n_vals.*Te_vals) + flux_diff_T;
    % min_vals(t_index) = min(min(min(f)));
    ranks(t_index) = rank;
 
    % n_vals = nvals_diff;
    % u_para = uvals_diff;
    % Ta_vals = Tivals_diff;

    % figure(2); clf; surf(Xmat, Zmat2, squeeze(f_vals(1, :, :))');
    % colorbar; shading interp;
    % drawnow;
    % PlotF(f_vals_low_rank, Xmat, Rmat, Zmat2, Nz);
end

% - compute errors
error = zeros(Nx, 1);
soln_n = zeros(Nx, 1);
soln_u = zeros(Nx, 1);
soln_Ta = zeros(Nx, 1);

for i = 1:Nx
    f = f_vals_low_rank{i, 1}*f_vals_low_rank{i, 2}*f_vals_low_rank{i, 3}';
    f_exact = soln{i, 1}*soln{i, 2}*soln{i, 3}';
    % error = error + 2*pi*dr*dz*sum(sum(abs(f - f_exact) .* Rmat));
    error(i) = 2*pi*dr*dz*sum(sum(abs(f - f_exact) .* Rmat));

    % % error = max(error, max(max(abs(f-f_exact))));
    % moments_ref_reimann(dtIndex, i, 1) = 2*pi*dr*dz*sum(sum(f.*Rmat));
    % moments_ref_reimann(dtIndex, i, 2) = (1/moments_ref_reimann(dtIndex, i, 1))*2*pi*dr*dz*sum(sum(Zmat.*f.*Rmat));
    % moments_ref_reimann(dtIndex, i, 3) = (2*pi*dr*dz*sum(sum((Rmat.^2 + Zmat.^2).*f.*Rmat)) - moments_ref_reimann(dtIndex, i, 1)*moments_ref_reimann(dtIndex, i, 2)^2)/(3*moments_ref_reimann(dtIndex, i, 1)/ma);

    soln_n(i) = 2*pi*dr*dz*sum(sum(f_exact.*Rmat));
    soln_u(i) = (1/soln_n(i))*2*pi*dr*dz*sum(sum(Zmat.*f_exact.*Rmat));
    soln_Ta(i) = (2*pi*dr*dz*sum(sum((Rmat.^2 + Zmat.^2).*f_exact.*Rmat)) - soln_n(i)*soln_u(i)^2)/(3*soln_n(i)/ma);
end
%ERRORS(dtIndex, 1) = dx*error;
ERRORS(1, dtIndex) = dx*sum(error);
ERRORS(2, dtIndex) = dx*sum(abs(n_vals - soln_n));
ERRORS(3, dtIndex) = dx*sum(abs(u_para - soln_u));
ERRORS(4, dtIndex) = dx*sum(abs(Ta_vals - soln_Ta));
ERRORS(5, dtIndex) = dx*sum(error(51:80));

% accuracy at individual points
for j = 1:Nx
    ERRORS(5 + j, dtIndex) = 2*pi*dr*dz*sum(sum(abs(f_vals_low_rank{j, 1}*f_vals_low_rank{j, 2}*f_vals_low_rank{j, 3}' - soln{j, 1}*soln{j, 2}*soln{j, 3}')));
end
end
%%

figure(1); clf;
plot(xvals, n_vals, "LineWidth", 1.5); hold on;
plot(xvals, u_para./u_para(1), "LineWidth",1.5);
plot(xvals, Te_vals, "LineWidth", 1.5);
plot(xvals, Ta_vals, "LineWidth", 1.5);
legend('n', 'nu_{||}', 'T_e', 'T_a');
title(sprintf('Mass, momentum, ion temperature, and electron temperature at t=%s', num2str(tval)));
ylim([0, 1.2]);
pause(0.05);

figure(4); clf; plot(tvals, abs(mass-mass(1))/mass(1), 'k-', 'LineWidth', 1.5);
xlabel('t'); ylabel('Relative error (mass)'); title('Relative mass of numerical solution over time');

figure(5); clf; plot(tvals, abs(momentum-momentum(1))/momentum(1), 'k-', 'LineWidth', 1.5);
xlabel('t'); ylabel('Absolute error (Uz)'); title('Absolute error of bulk velocity over time');

figure(6); clf; plot(tvals, abs(energy-energy(1))/energy(1), 'k-', 'LineWidth', 1.5);
xlabel('t'); ylabel('Relative error (Energy)'); title('Relative energy of numerical solution over time');

PlotF(f_vals_low_rank, Xmat, Rmat, Zmat2, Nz);
title(sprintf('Numerical solution at time t=%s', num2str(tf)));

%% Errors
disp('Errors:');
disp('L1 (global)');
disp(ERRORS(1, :)');
disp(log2(ERRORS(1, 1:end-1)./ERRORS(1, 2:end))');
disp('Mass');
disp(ERRORS(2, :)');
disp(log2(ERRORS(2, 1:end-1)./ERRORS(2, 2:end))');
disp('Momentum');
disp(ERRORS(3, :)');
disp(log2(ERRORS(3, 1:end-1)./ERRORS(3, 2:end))');
disp('Energy');
disp(ERRORS(4, :)');
disp(log2(ERRORS(4, 1:end-1)./ERRORS(4, 2:end))');
disp('L1 (nodes 51 through 80)');
disp(ERRORS(5, :)');
disp(log2(ERRORS(5, 1:end-1)./ERRORS(5, 2:end))');

for j = 1:Nx
    fprintf('L1 (local, i=%d)\n', j);
    disp(ERRORS(5 + j, :)');
    disp(log2(ERRORS(5 + j, 1:end-1)./ERRORS(5 + j, 2:end))');
end


% figure(1); clf; surf(Rmat, Zmat, f);
% colorbar; shading interp;
% legend(sprintf('N_r = %s', num2str(Nr, 3)), 'Location','northwest');
% xlabel('V_r'); ylabel('V_z'); zlabel('U'); title([sprintf('Backward Euler approximation of 0D2V Fokker-Planck system at time %s', num2str(tf, 4))]);
% 
% figure(2); clf; surf(Rmat, Zmat, f_inf);
% colorbar; shading interp;
% xlabel('V_r'); ylabel('V_z'); zlabel('f(V_r, V_z, t)'); title([sprintf('f_{exact} at time t=%s', num2str(tf, 4))]);
% 
% l1 decay
% figure(3); clf; semilogy(tvals, l1, 'black-', 'LineWidth', 1.5);
% xlabel('t'); ylabel('L_1(f(V_r, V_z))'); title('L_1 decay of numerical solution over time');

% Temporal errors!
return
%%
% Mass
figure(4); clf; hold on;

plot(tvals, abs(mass-mass(1))/mass(1), 'b-', 'LineWidth', 1.5);
xlabel('t'); 

plot(tvals, abs(momentum-momentum(1)), 'k-', 'LineWidth', 1.5);
xlabel('t'); 

plot(tvals, abs(energy-energy(1))/energy(1), 'g-', 'LineWidth', 1.5);
xlabel('t');

legend('Mass', 'Momentum', 'Energy');
title('Relative errors in mass, momentum, energy over time');

return

% Positivity
figure(4); clf; plot(tvals, min_vals, 'green-', 'LineWidth', 1.5);
xlabel('t'); ylabel('min(f(V_r, V_z))'); title('Minimum values of numerical solution over time');

% Relative entropy
figure(5); clf; semilogy(tvals, relative_entropy, 'magenta-', 'LineWidth', 1.5);
xlabel('t'); ylabel('Relative entropy'); title('Relative entropy of numerical solution over time');

% Rank plot
figure(7); clf;
plot(tvals, ranks, 'black-', 'LineWidth', 1.5); hold on;
% plot(ranks2(:, 1), ranks2(:, 2), 'blue-', 'LineWidth', 1.5);
% plot(ranks3(:, 1), ranks3(:, 2), 'green-', 'LineWidth', 1.5);
xlabel('time'); ylabel('rank'); title('Rank plot over time');
legend('Backward Euler', 'RK2', 'RK3');







%%%%%% FUNCTIONS %%%%%%
function [Dx] = GetVlasov(dt, dx, Zmat, f_vals_low_rank, leftBC, rightBC, i)
    f_i = f_vals_low_rank{i, 1}*f_vals_low_rank{i, 2}*f_vals_low_rank{i, 3}';
    if i == 1
        f_i_pos = f_vals_low_rank{i+1, 1}*f_vals_low_rank{i+1, 2}*f_vals_low_rank{i+1, 3}';
        f_i_neg = leftBC;
    elseif i == size(f_vals_low_rank, 1)
        f_i_pos = rightBC;
        f_i_neg = f_vals_low_rank{i-1, 1}*f_vals_low_rank{i-1, 2}*f_vals_low_rank{i-1, 3}';
    else
        f_i_pos = f_vals_low_rank{i+1, 1}*f_vals_low_rank{i+1, 2}*f_vals_low_rank{i+1, 3}';
        f_i_neg = f_vals_low_rank{i-1, 1}*f_vals_low_rank{i-1, 2}*f_vals_low_rank{i-1, 3}';
    end
    Dx = (dt/dx)*max(Zmat, 0).*(f_i - f_i_neg) + (dt/dx)*min(Zmat, 0).*(f_i_pos - f_i);
end

function [Lorentz] = GetLorentz(dx, x_min, x_max, zvals, ma, qe, qa, ne_vals, Te_vals, i)
    Nx = numel(ne_vals);
    Nz = numel(zvals);
    dz = zvals(2) - zvals(1);

    % first compute E_para using Generalized Ohm's Law
    [n0, ~, Te0] = get_boundaries(x_min, x_max);
    n_leftBC = n0(x_min - dx/2);
    n_rightBC = n0(x_max + dx/2);
    Te_leftBC = Te0(x_min - dx/2);
    Te_rightBC = Te0(x_max + dx/2);
    if i == 1
        E = (1/(2*dx*ne_vals(i)*qe))*(ne_vals(i+1)*Te_vals(i+1) - n_leftBC*Te_leftBC); % centered difference
    elseif i == Nx
        E = (1/(2*dx*ne_vals(i)*qe))*(n_rightBC*Te_rightBC - ne_vals(i-1)*Te_vals(i-1)); % centered difference
    else
        E = (1/(2*dx*ne_vals(i)*qe))*(ne_vals(i+1)*Te_vals(i+1) - ne_vals(i-1)*Te_vals(i-1)); % centered difference
    end

    diff_mtx_pos = (1/dz)*(spdiags(ones(Nz, 1), 0, Nz, Nz) + spdiags(-1*ones(Nz, 1), -1, Nz, Nz));
    diff_mtx_neg = -1*diff_mtx_pos';
    Lorentz = max(E, 0)*(qa/ma)*diff_mtx_pos + min(E, 0)*(qa/ma)*diff_mtx_neg; 
end

function [Flux] = GetFokkerPlanckFlux(A, B, C, u, T, m, xvals, t, dx)
    N = numel(xvals);
    A = A(u, xvals, t);
    B = B(u, xvals(1:N-1) + dx/2, t);
    C = C(u, xvals(1:N-1) + dx/2, t, T, m);
    
    w = dx*B./C + 1e-14; % prevent zero division
    delta = (1 ./ w) - (1 ./ (exp(w) - 1));

    F1 = -((1/dx)*C - delta.*B);
    F2 = (1 - delta).*B + (1/dx)*C;

    F_pos = spdiags([F1;0], 0, N, N) + spdiags([0;F2], 1, N, N);
    F_neg = spdiags([0; F2], 0, N, N) + spdiags(F1, -1, N, N);

    Flux = diag(1./A)*(1/dx)*(F_pos - F_neg);
end

function [Vr, S, Vz, rank] = IMEX111(Vr0, S0, Vz0, ur, uz, dt, dx, tval, rvals, zvals, x_min, x_max, f_vals_low_rank, Rmat, Zmat, ma, me, qa, qe, Ar, Az, Br, Bz, Cr, Cz, tolerance, n_vals, Ta_vals, Te_vals, rhoM, JzM, kappaM, spatialIndex, R_const, leftBC, rightBC, wr, wz, c, w_norm_1_squared, w_norm_v_squared, w_norm_v2_squared)
    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);
    Nr = numel(rvals);
    Nz = numel(zvals);

    % discretize Fokker-Planck operators
    nu_aa = n_vals(spatialIndex)/(Ta_vals(spatialIndex)^(3/2));
    nu_ae = sqrt(2*me)*n_vals(spatialIndex)/(Te_vals(spatialIndex)^(3/2));
    C_aa_r1 = nu_aa*GetFokkerPlanckFlux(Ar, Br, Cr, ur(spatialIndex), Ta_vals(spatialIndex), ma, rvals, tval, dr);
    C_ae_r1 = nu_ae*GetFokkerPlanckFlux(Ar, Br, Cr, ur(spatialIndex), Te_vals(spatialIndex), ma, rvals, tval, dr);
    C_aa_z1 = nu_aa*GetFokkerPlanckFlux(Az, Bz, Cz, uz(spatialIndex), Ta_vals(spatialIndex), ma, zvals, tval, dz);
    C_ae_z1 = nu_ae*GetFokkerPlanckFlux(Az, Bz, Cz, uz(spatialIndex), Te_vals(spatialIndex), ma, zvals, tval, dz);
    Cr1 = C_aa_r1 + C_ae_r1;
    Cz1 = C_aa_z1 + C_ae_z1;

    % discretize Lorentz force
    Dz1 = GetLorentz(dx, x_min, x_max, zvals, ma, qe, qa, n_vals, Te_vals, spatialIndex);

    % discretize velocity explicitly
    Dx1 = GetVlasov(dt, dx, Zmat, f_vals_low_rank, leftBC, rightBC, spatialIndex);
    W0 = (Vr0*S0*Vz0') - Dx1;

    Vr0_star = Vr0;
    Vz0_star = Vz0;

    K1 = sylvester(eye(Nr) - (dt*Cr1), dt*((Dz1 - Cz1)*Vz0_star)'*Vz0_star, W0*Vz0_star);
    L1 = sylvester(eye(Nz) + dt*(Dz1 - Cz1), -dt*(Cr1*Vr0_star)'*(rvals.*Vr0_star), W0'*(rvals .* Vr0_star));
    
    [Vr1_ddagger, ~] = qr2(K1, rvals);
    [Vz1_ddagger, ~] = qr(L1, 0);

    [Vr1_hat, Vz1_hat] = reduced_augmentation([Vr1_ddagger, Vr0], [Vz1_ddagger, Vz0], rvals);

    S1_hat = sylvester((eye(size(Vr1_hat, 2)) - (dt*((rvals .* Vr1_hat)')*(Cr1*Vr1_hat))), dt*((Dz1 - Cz1)*Vz1_hat)'*Vz1_hat, ((rvals .* Vr1_hat)'*W0*Vz1_hat));
    [Vr, S, Vz, rank] = LoMaC(Vr1_hat, S1_hat, Vz1_hat, Rmat, Zmat, rvals, zvals, tolerance, rhoM, JzM, kappaM, wr, wz, c, w_norm_1_squared, w_norm_v_squared, w_norm_v2_squared);
    % [Vr, S, Vz, rank] = truncate_svd(Vr1_hat, S1_hat, Vz1_hat, tolerance);
end

function [Vr, S, Vz, rank] = DIRK2Timestep(Vr0, S0, Vz0, ur, uz, dt, tval, rvals, zvals, Rmat, Zmat, Ar, Az, Br, Bz, Cr, Cz, tolerance, rhoM, JzM, kappaM)
    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);
    Nr = numel(rvals);
    Nz = numel(zvals);

    Fr1 = GetFokkerPlanckFlux(Ar, Br, Cr, ur, rvals, tval, dr);
    Fz1 = GetFokkerPlanckFlux(Az, Bz, Cz, uz, zvals, tval, dz);

    nu = 1-(sqrt(2)/2);

    % Stage 1: Backward Euler
    [Vr1, S1, Vz1, ~] = IMEX111(Vr0, S0, Vz0, ur, uz, nu*dt, tval, rvals, zvals, Rmat, Zmat, Ar, Az, Br, Bz, Cr, Cz, tolerance, rhoM, JzM, kappaM);

    W0 = (Vr0*S0*(Vz0')) + ((1-nu)*dt*(((Fr1*(Vr1)*S1*(Vz1')) + ((Vr1)*S1*((Fz1*Vz1)')))));

    % Reduced Augmentation
    % Predict V_dagger using B. Euler for second stage
    [Vr1_dagger, ~, Vz1_dagger, ~] = IMEX111(Vr0, S0, Vz0, ur, uz, dt, tval, rvals, zvals, Rmat, Zmat, Ar, Az, Br, Bz, Cr, Cz, tolerance, rhoM, JzM, kappaM);
    [Vr_star, Vz_star] = reduced_augmentation([Vr1_dagger, Vr1, Vr0], [Vz1_dagger, Vz1, Vz0], rvals);

    % Stage 2: KLS Steps
    Fr2 = GetFokkerPlanckFlux(Ar, Br, Cr, ur, rvals, tval+dt, dr);
    Fz2 = GetFokkerPlanckFlux(Az, Bz, Cz, uz, zvals, tval+dt, dz);
   
    % K/L-Step
    K1 = sylvester(eye(Nr) - (nu*dt*Fr2), -nu*dt*(Fz2*Vz_star)'*Vz_star, W0*Vz_star);
    L1 = sylvester(eye(Nz) - (nu*dt*Fz2), -nu*dt*(Fr2*Vr_star)'*(rvals .* Vr_star), W0'*(rvals .* Vr_star));

    % Get bases
    [Vr_ddagger, ~] = qr2(K1, rvals); [Vz_ddagger, ~] = qr(L1, 0);

    % Reduced Augmentation
    [Vr1_hat, Vz1_hat] = reduced_augmentation([Vr_ddagger, Vr1, Vr0], [Vz_ddagger, Vz1, Vz0], rvals);

    % S-Step
    S1_hat = sylvester(eye(size(Vr1_hat, 2)) - (nu*dt*((rvals .* Vr1_hat)')*Fr2*Vr1_hat), -nu*dt*(Fz2*Vz1_hat)'*Vz1_hat, ((rvals .* Vr1_hat)')*W0*Vz1_hat);
    [Vr, S, Vz, rank] = LoMaC(Vr1_hat, S1_hat, Vz1_hat, Rmat, Zmat, rvals, zvals, tolerance, rhoM, JzM, kappaM);
end

function [Vr, Vz] = reduced_augmentation(Vr_aug, Vz_aug, rvals)
    tolerance = 1e-14;
    [Qr, Rr] = qr2(Vr_aug, rvals);
    [Qz, Rz] = qr(Vz_aug, 0);
    [Ur, Sr, ~] = svd(Rr, 0);
    [Uz, Sz, ~] = svd(Rz, 0);
    rr = find(diag(Sr) > tolerance, 1, 'last');
    rz = find(diag(Sz) > tolerance, 1, 'last');
    rank = max(rr, rz);
    rank = min(rank, min(size(Ur, 2), size(Uz, 2)));
    Vr = Qr*Ur(:, 1:rank);
    Vz = Qz*Uz(:, 1:rank);
end

function [Q, R] = qr2(X, rvals)
    [Q, R] = qr(sqrt(rvals) .* X, 0);
    Q = Q ./ sqrt(rvals);
end

function [U, S, V] = svd2(X, rvals)
    [U, S, V] = svd(sqrt(rvals) .* X, 0);
    U = U./sqrt(rvals);
end

function [X, R, Z, dx, dr, dz] = GetRZ(Nx, Nr, Nz, interval)
    xvals = linspace(interval(1), interval(2), Nx+1)';
    rvals = linspace(interval(3), interval(4), Nr+1)';
    zvals = linspace(interval(5), interval(6), Nz+1)';
    dx = xvals(2) - xvals(1);
    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);
    xmid = xvals(1:end-1) + dx/2; % centered mesh
    rmid = rvals(1:end-1) + dr/2; % centered mesh
    zmid = zvals(1:end-1) + dz/2; % centered mesh
    [R, Z] = meshgrid(rmid, zmid); 
    X = xmid; R = R'; Z = Z';
end



% ------- FLUID SOLVER STUFF -------
function [n1, u_para1, T_a1, T_e1, u_para0_half_nodes, nu_hat1, S_hat1, Q_hat1, nTe_hat1, kappaTx] = fluid_solver_IMEX111(f_vals_low_rank, n0, u_para0, T_a0, T_e0, dt, dx, dv_perp, dv_para, v_perp, v_para, qa, qe, ma, me, R_const, x_min, x_max)
    Nx = numel(n0); % max val of i, j
    x_ghost_left = x_min - dx/2; % left ghost cell
    x_ghost_right = x_max + dx/2; % right ghost cell

    % get boundary conditions
    [n0_boundary, u_para0_boundary, T_ae_boundary] = get_boundaries(x_min, x_max); % Ta = Te at boundary --> T_ae
    n_BC_left = n0_boundary(x_ghost_left); n_BC_right = n0_boundary(x_ghost_right);
    u_para_BC_left = u_para0_boundary(x_ghost_left); u_para_BC_right = u_para0_boundary(x_ghost_right);
    Tae_BC_left = T_ae_boundary(x_ghost_left); Tae_BC_right = T_ae_boundary(x_ghost_right);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ---------- STAGE ONE ----------- %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n0_pos = [n0; n_BC_right]; n0_neg = [n_BC_left; n0];
    u_para0_half_nodes = ([u_para_BC_left; u_para0] + [u_para0; u_para_BC_right])/2; %u_para0_half_nodes = [u_para0_min; u_para0_half_nodes; u_para0_max];
    U0 = 0.5*((3/ma)*T_a0 + u_para0.^2);
    Te0_pos = [T_e0(2:end); Tae_BC_right]; Te0_neg = [Tae_BC_left; T_e0(1:end-1)];

    % first, compute fluxes via summation
    [nu_hat1, S_hat1, Q_hat1] = get_fluxes(f_vals_low_rank, v_perp, v_para, R_const, x_min, x_max, dx, dv_perp, dv_para);
    nTe_hat1 = ((n0_neg.*[Tae_BC_left; T_e0]).*(u_para0_half_nodes > 0) + (n0_pos.*[T_e0; Tae_BC_right]).*(u_para0_half_nodes <= 0)); % upwinding

    % shift bounds to get flux pos/neg
    nu_hat1_pos = nu_hat1(2:end); nu_hat1_neg = nu_hat1(1:end-1);

    S_hat1_pos = S_hat1(2:end); S_hat1_neg = S_hat1(1:end-1);
    Q_hat1_pos = Q_hat1(2:end); Q_hat1_neg = Q_hat1(1:end-1);
    nTe_hat1_pos = nTe_hat1(2:end); nTe_hat1_neg = nTe_hat1(1:end-1);
  
    % explicitly find n via Forward Euler
    n1 = n0 - (dt/dx)*(nu_hat1_pos - nu_hat1_neg); % now find n_i+1, n_i-1
    n_pos = [n1(2:end); n_BC_right]; n_neg = [n_BC_left; n1(1:end-1)];

    % ---- init y_vec, R_norm ----  
    % y = [n0.*u_para0; n0.*U0; T_e0];
    y = [n1.*u_para0; n1.*U0; T_e0];

    nu = y(1:Nx); 
    nU = y(Nx+1:2*Nx); 
    T_e = y(2*Nx+1:end);
    Te_pos = [T_e(2:end); Tae_BC_right]; Te_neg = [Tae_BC_left; T_e(1:end-1)];
    nTe_pos = n_pos.*Te_pos; nTe_neg = n_neg.*Te_neg;

    kappa_pos = (3.2/(2*sqrt(2*me)))*((Te_pos.^(5/2) + T_e.^(5/2)));
    kappa_neg = (3.2/(2*sqrt(2*me)))*((T_e.^(5/2) + Te_neg.^(5/2))); 

    R1 = nu - (n0.*u_para0) + (dt/dx).*(S_hat1_pos - S_hat1_neg) - ((dt*qa)/(2*dx*qe*ma)).*((n_pos.*Te_pos) - (n_neg.*Te_neg));
    R2 = nU - (n0.*U0) + (dt/dx).*(Q_hat1_pos - Q_hat1_neg) - (((dt*qa*nu)./(2*dx*qe*n1)).*((n_pos.*Te_pos) - (n_neg.*Te_neg))) - (((dt.*3.*sqrt(2*me))./((ma.^2).*(T_e.^(3/2)))) .* (((n1.^2).*T_e) - ((ma/3).*((2.*n1.*nU) - (nu.^2)))));
    R3 = (n1.*T_e) - (n0.*T_e0) + ((5*dt)/(3*dx)).*((u_para0_half_nodes(2:end).*nTe_hat1_pos) - (u_para0_half_nodes(1:end-1).*nTe_hat1_neg)) - (((dt*(nu./n1))./(3*dx)).*(nTe_pos - nTe_neg)) - (((2*dt)/(3*dx.^2)) .* ((kappa_pos.*(Te_pos - T_e)) - (kappa_neg.*(T_e - Te_neg)))) - (((dt.*2.*sqrt(2*me))./(ma.*(T_e.^(3/2)))) .* (((ma/3).*((2*n1.*nU) - (nu.^2))) - ((n1.^2).*T_e)));
    R = [R1; R2; R3];

    tol = min(5e-12, max(abs(R))*5e-10); % ensure we don't get worse!
   
    while max(abs(R)) > tol
        % define partial derivatives of residual
        nu_nu = spdiags(ones(Nx), 0, Nx, Nx);
        nu_nU = spdiags(zeros(Nx),0, Nx, Nx);
        nu_Te = gallery('tridiag', (dt*qa*n1(1:end-1))/(2*dx*qe*ma), zeros(Nx,1), -(dt*qa*n1(2:end))/(2*dx*qe*ma));
    
        nU_nu = diag( ((-dt*qa)./(2*dx.*qe.*n1)).*((n_pos.*Te_pos) - (n_neg.*Te_neg)) );
        nU_nU = spdiags(ones(Nx), 0, Nx, Nx);
        nU_Te_mid = ((3*dt*sqrt(2*me))./(2*(ma.^2))).*(((n1.^2)./(T_e.^(3/2))  + (ma./(T_e.^(5/2))).*((2*n1.*nU) - (nu.^2)) ));
        nU_Te = gallery('tridiag', (dt.*qa.*nu(2:end).*(n1(1:end-1)./n1(2:end)))./(2*dx*qe), nU_Te_mid, -(dt.*qa.*nu(1:end-1).*(n1(2:end)./n1(1:end-1)))/(2*dx*qe) );
    
        Te_nu = diag(-(dt.*((n_pos.*Te_pos) - (n_neg.*Te_neg)))./(3*dx*n1));
        Te_nU = spdiags(zeros(Nx), 0, Nx, Nx);
        Te_Te_left = ((dt.*nu.*n_neg)./(3*dx.*n1)) + (((dt*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(Te_neg.^(3/2))/2).*(T_e - Te_neg)) - (T_e.^(5/2) + Te_neg.^(5/2)))); Te_Te_left = Te_Te_left(2:end);
        Te_Te_mid = n1 - (((dt*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(T_e.^(3/2))./2).*(Te_pos - T_e)) - (Te_pos.^(5/2) + T_e.^(5/2)) - ((5.*(T_e.^(3/2))/2).*(T_e - Te_neg)) - (T_e.^(5/2) + Te_neg.^(5/2)))) - ((2*dt*sqrt(2*me)./(2*ma)).*(((-ma./(T_e.^(5/2))).*(2*n1.*nU - (nu.^2)) + ((n1.^2)./(T_e.^(3/2))))));
        Te_Te_right = ((-dt.*nu.*n_pos)./(3*dx.*n1)) - (((dt*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(Te_pos.^(3/2))/2).*(Te_pos - T_e)) + (Te_pos.^(5/2) + T_e.^(5/2)))); Te_Te_right = Te_Te_right(1:end-1);
        Te_Te = gallery('tridiag', Te_Te_left, Te_Te_mid, Te_Te_right);

        P = [nu_nu, nu_nU, nu_Te;
             nU_nu, nU_nU, nU_Te;
             Te_nu, Te_nU, Te_Te]; % holy block matrix
     
        dy = -P\R; % solve for delta y at time l+1
        y = y + dy;

        % update y_vec stuff
        nu = y(1:Nx); 
        nU = y(Nx+1:2*Nx); 
        T_e = y(2*Nx+1:end);
        Te_pos = [T_e(2:end); Tae_BC_right]; Te_neg = [Tae_BC_left; T_e(1:end-1)];
        nTe_pos = n_pos.*Te_pos; nTe_neg = n_neg.*Te_neg;

        kappa_pos = (3.2/(2*sqrt(2*me)))*((Te_pos.^(5/2) + T_e.^(5/2)));
        kappa_neg = (3.2/(2*sqrt(2*me)))*((T_e.^(5/2) + Te_neg.^(5/2)));
        
        R1 = nu - (n0.*u_para0) + (dt/dx).*(S_hat1_pos - S_hat1_neg) - ((dt*qa)/(2*dx*qe*ma)).*((n_pos.*Te_pos) - (n_neg.*Te_neg));
        R2 = nU - (n0.*U0) + (dt/dx).*(Q_hat1_pos - Q_hat1_neg) - (((dt*qa*nu)./(2*dx*qe*n1)).*((n_pos.*Te_pos) - (n_neg.*Te_neg))) - (((dt.*3.*sqrt(2*me))./((ma.^2).*(T_e.^(3/2)))) .* (((n1.^2).*T_e) - ((ma/3).*((2.*n1.*nU) - (nu.^2)))));
        R3 = (n1.*T_e) - (n0.*T_e0) + ((5*dt)/(3*dx)).*((u_para0_half_nodes(2:end).*nTe_hat1_pos) - (u_para0_half_nodes(1:end-1).*nTe_hat1_neg)) - (((dt*(nu./n1))./(3*dx)).*(nTe_pos - nTe_neg)) - (((2*dt)/(3*dx.^2)) .* ((kappa_pos.*(Te_pos - T_e)) - (kappa_neg.*(T_e - Te_neg)))) - (((dt.*2.*sqrt(2*me))./(ma.*(T_e.^(3/2)))) .* (((ma/3).*((2*n1.*nU) - (nu.^2))) - ((n1.^2).*T_e)));
        R = [R1; R2; R3];
    end
  
    % parameters to return
    nu1 = (y(1:Nx));
    nU1 = y(Nx+1:2*Nx);
    T_e1 = y(2*Nx+1:end);

    u_para1 = nu1 ./ n1;
    T_a1 = (ma/3)*(2*nU1./n1 - (nu1.^2)./(n1.^2));

    Te_all = [Tae_BC_left; T_e; Tae_BC_right];
    kappaTx = (3.2/(2*sqrt(2*me)*dx))*(Te_all(1:end-1).^(2.5)+Te_all(2:end).^(2.5)).*(Te_all(2:end)-Te_all(1:end-1));
end

function [n2, u_para2, T_a2, T_e2, u_para0_half_nodes, nu_hat2, S_hat2, Q_hat2, nTe_hat2, kappaTx] = fluid_solver_IMEX222(f_vals_low_rank, n0, u_para0, T_a0, T_e0, dt, dx, dv_perp, dv_para, v_perp, v_para, qa, qe, ma, me, R_const, x_min, x_max)

    gamma = 1 - (sqrt(2)/2);
    delta = 1 - (1/(2*gamma));
    dt1 = gamma*dt;
    dt2 = (1-gamma)*dt;

    Nx = numel(n0); % max val of i, j
    x_ghost_left = x_min - dx/2; % left ghost cell
    x_ghost_right = x_max + dx/2; % right ghost cell

    % get boundary conditions
    [n0_boundary, u_para0_boundary, T_ae_boundary] = get_boundaries(); % Ta = Te at boundary --> T_ae
    n_BC_left = n0_boundary(x_ghost_left); n_BC_right = n0_boundary(x_ghost_right);
    u_para_BC_left = u_para0_boundary(x_ghost_left); u_para_BC_right = u_para0_boundary(x_ghost_right);
    Tae_BC_left = T_ae_boundary(x_ghost_left); Tae_BC_right = T_ae_boundary(x_ghost_right);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ---------- STAGE ONE ----------- %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n0_pos = [n0; n_BC_right]; n0_neg = [n_BC_left; n0];
    u_para0_half_nodes = ([u_para_BC_left; u_para0] + [u_para0; u_para_BC_right])/2; %u_para0_half_nodes = [u_para0_min; u_para0_half_nodes; u_para0_max];
    U0 = 0.5*((3/ma)*T_a0 + u_para0.^2);
    Te0_pos = [T_e0(2:end); Tae_BC_right]; Te0_neg = [Tae_BC_left; T_e0(1:end-1)];

    % first, compute fluxes via summation
    [nu_hat1, S_hat1, Q_hat1] = get_fluxes(f_vals_low_rank, v_perp, v_para, R_const, x_min, x_max, dx, dv_perp, dv_para);
    nTe_hat1 = ((n0_neg.*[Tae_BC_left; T_e0]).*(u_para0_half_nodes > 0) + (n0_pos.*[T_e0; Tae_BC_right]).*(u_para0_half_nodes <= 0)); % upwinding

    % shift bounds to get flux pos/neg
    nu_hat1_pos = nu_hat1(2:end); nu_hat1_neg = nu_hat1(1:end-1);

    S_hat1_pos = S_hat1(2:end); S_hat1_neg = S_hat1(1:end-1);
    Q_hat1_pos = Q_hat1(2:end); Q_hat1_neg = Q_hat1(1:end-1);
    nTe_hat1_pos = nTe_hat1(2:end); nTe_hat1_neg = nTe_hat1(1:end-1);

    % define A0 for later
    An_0 = - (1/dx)*(nu_hat1_pos - nu_hat1_neg);
  
    % explicitly find n via Forward Euler
    n1 = n0 + dt1*An_0; % now find n_i+1, n_i-1
    n_pos = [n1(2:end); n_BC_right]; n_neg = [n_BC_left; n1(1:end-1)];

    % define the rest of the A0s
    Anu_0 = - (1/dx)*(S_hat1_pos - S_hat1_neg);
    AnU_0 = - (1/dx).*(Q_hat1_pos - Q_hat1_neg);
    AnTe_0 = - ((5)/(3*dx)).*((u_para0_half_nodes(2:end).*nTe_hat1_pos) - (u_para0_half_nodes(1:end-1).*nTe_hat1_neg));

    % ---- init y_vec, R_norm ----  
    % y = [n0.*u_para0; n0.*U0; T_e0];
    y = [n1.*u_para0; n1.*U0; T_e0];

    nu = y(1:Nx); 
    nU = y(Nx+1:2*Nx); 
    T_e = y(2*Nx+1:end);
    Te_pos = [T_e(2:end); Tae_BC_right]; Te_neg = [Tae_BC_left; T_e(1:end-1)];
    nTe_pos = n_pos.*Te_pos; nTe_neg = n_neg.*Te_neg;

    kappa_pos = (3.2/(2*sqrt(2*me)))*((Te_pos.^(5/2) + T_e.^(5/2)));
    kappa_neg = (3.2/(2*sqrt(2*me)))*((T_e.^(5/2) + Te_neg.^(5/2))); 

    R1 = nu - (n0.*u_para0) + (dt1/dx).*(S_hat1_pos - S_hat1_neg) - ((dt1*qa)/(2*dx*qe*ma)).*((n_pos.*Te_pos) - (n_neg.*Te_neg));
    R2 = nU - (n0.*U0) + (dt1/dx).*(Q_hat1_pos - Q_hat1_neg) - (((dt1*qa*nu)./(2*dx*qe*n1)).*((n_pos.*Te_pos) - (n_neg.*Te_neg))) - (((dt1.*3.*sqrt(2*me))./((ma.^2).*(T_e.^(3/2)))) .* (((n1.^2).*T_e) - ((ma/3).*((2.*n1.*nU) - (nu.^2)))));
    R3 = (n1.*T_e) - (n0.*T_e0) + ((5*dt1)/(3*dx)).*((u_para0_half_nodes(2:end).*nTe_hat1_pos) - (u_para0_half_nodes(1:end-1).*nTe_hat1_neg)) - (((dt1*(nu./n1))./(3*dx)).*(nTe_pos - nTe_neg)) - (((2*dt1)/(3*dx.^2)) .* ((kappa_pos.*(Te_pos - T_e)) - (kappa_neg.*(T_e - Te_neg)))) - (((dt1.*2.*sqrt(2*me))./(ma.*(T_e.^(3/2)))) .* (((ma/3).*((2*n1.*nU) - (nu.^2))) - ((n1.^2).*T_e)));
    R = [R1; R2; R3];

    err = 1;
    tol = min(5e-12, max(abs(R))*5e-10); % ensure we don't get worse!
   
    while err > tol
        % define partial derivatives of residual
        nu_nu = spdiags(ones(Nx), 0, Nx, Nx);
        nu_nU = spdiags(zeros(Nx),0, Nx, Nx);
        nu_Te = gallery('tridiag', (dt1*qa*n1(1:end-1))/(2*dx*qe*ma), zeros(Nx,1), -(dt1*qa*n1(2:end))/(2*dx*qe*ma));
    
        nU_nu = diag( ((-dt1*qa)./(2*dx.*qe.*n1)).*((n_pos.*Te_pos) - (n_neg.*Te_neg)) );
        nU_nU = spdiags(ones(Nx), 0, Nx, Nx);
        nU_Te_mid = ((3*dt1*sqrt(2*me))./(2*(ma.^2))).*(((n1.^2)./(T_e.^(3/2))  + (ma./(T_e.^(5/2))).*((2*n1.*nU) - (nu.^2)) ));
        nU_Te = gallery('tridiag', (dt1.*qa.*nu(2:end).*(n1(1:end-1)./n1(2:end)))./(2*dx*qe), nU_Te_mid, -(dt1.*qa.*nu(1:end-1).*(n1(2:end)./n1(1:end-1)))/(2*dx*qe) );
    
        Te_nu = diag(-(dt1.*((n_pos.*Te_pos) - (n_neg.*Te_neg)))./(3*dx*n1));
        Te_nU = spdiags(zeros(Nx), 0, Nx, Nx);
        Te_Te_left = ((dt1.*nu.*n_neg)./(3*dx.*n1)) + (((dt1*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(Te_neg.^(3/2))/2).*(T_e - Te_neg)) - (T_e.^(5/2) + Te_neg.^(5/2)))); Te_Te_left = Te_Te_left(2:end);
        Te_Te_mid = n1 - (((dt1*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(T_e.^(3/2))./2).*(Te_pos - T_e)) - (Te_pos.^(5/2) + T_e.^(5/2)) - ((5.*(T_e.^(3/2))/2).*(T_e - Te_neg)) - (T_e.^(5/2) + Te_neg.^(5/2)))) - ((2*dt1*sqrt(2*me)./(2*ma)).*(((-ma./(T_e.^(5/2))).*(2*n1.*nU - (nu.^2)) + ((n1.^2)./(T_e.^(3/2))))));
        Te_Te_right = ((-dt1.*nu.*n_pos)./(3*dx.*n1)) - (((dt1*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(Te_pos.^(3/2))/2).*(Te_pos - T_e)) + (Te_pos.^(5/2) + T_e.^(5/2)))); Te_Te_right = Te_Te_right(1:end-1);
        Te_Te = gallery('tridiag', Te_Te_left, Te_Te_mid, Te_Te_right);

        P = [nu_nu, nu_nU, nu_Te;
             nU_nu, nU_nU, nU_Te;
             Te_nu, Te_nU, Te_Te]; % holy block matrix
     
        dy = -P\R; % solve for delta y at time l+1
        y = y + dy;

        % update y_vec stuff
        nu = y(1:Nx); u = nu ./ n1;
        nU = y(Nx+1:2*Nx); U = nU ./ n1;
        T_e = y(2*Nx+1:end);
        Te_pos = [T_e(2:end); Tae_BC_right]; Te_neg = [Tae_BC_left; T_e(1:end-1)];
        nTe_pos = n_pos.*Te_pos; nTe_neg = n_neg.*Te_neg;

        kappa_pos = (3.2/(2*sqrt(2*me)))*((Te_pos.^(5/2) + T_e.^(5/2)));
        kappa_neg = (3.2/(2*sqrt(2*me)))*((T_e.^(5/2) + Te_neg.^(5/2)));
        
        R1 = nu - (n0.*u_para0) + (dt1/dx).*(S_hat1_pos - S_hat1_neg) - ((dt1*qa)/(2*dx*qe*ma)).*((n_pos.*Te_pos) - (n_neg.*Te_neg));
        R2 = nU - (n0.*U0) + (dt1/dx).*(Q_hat1_pos - Q_hat1_neg) - (((dt1*qa*nu)./(2*dx*qe*n1)).*((n_pos.*Te_pos) - (n_neg.*Te_neg))) - (((dt1.*3.*sqrt(2*me))./((ma.^2).*(T_e.^(3/2)))) .* (((n1.^2).*T_e) - ((ma/3).*((2.*n1.*nU) - (nu.^2)))));
        R3 = (n1.*T_e) - (n0.*T_e0) + ((5*dt1)/(3*dx)).*((u_para0_half_nodes(2:end).*nTe_hat1_pos) - (u_para0_half_nodes(1:end-1).*nTe_hat1_neg)) - (((dt1*(nu./n1))./(3*dx)).*(nTe_pos - nTe_neg)) - (((2*dt1)/(3*dx.^2)) .* ((kappa_pos.*(Te_pos - T_e)) - (kappa_neg.*(T_e - Te_neg)))) - (((dt1.*2.*sqrt(2*me))./(ma.*(T_e.^(3/2)))) .* (((ma/3).*((2*n1.*nU) - (nu.^2))) - ((n1.^2).*T_e)));
        R = [R1; R2; R3];
        err = max(abs(R));
    end
  
    % parameters to return
    nu1 = (y(1:Nx));
    nU1 = y(Nx+1:2*Nx);
    T_e1 = y(2*Nx+1:end);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ---------- STAGE TWO ----------- %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u_para1 = nu1 ./ n1;
    U1 = nU1 ./ n1;
    Ta1 = (ma/3)*(2*U1 - (u_para1.^2));

    n1_pos = [n1(2:end); n_BC_right]; n1_neg = [n_BC_left; n1(1:end-1)];
    u_para1_half_nodes = ([u_para_BC_left; u_para1] + [u_para1; u_para_BC_right])/2;
    Te1_pos = [T_e1(2:end); Tae_BC_right]; Te1_neg = [Tae_BC_left; T_e1(1:end-1)];
    nTe1_pos = n1_pos.*Te1_pos; nTe1_neg = n1_neg.*Te1_neg;

    kappa1_pos = (3.2/(2*sqrt(2*me)))*((Te1_pos.^(5/2) + T_e1.^(5/2)));
    kappa1_neg = (3.2/(2*sqrt(2*me)))*((T_e1.^(5/2) + Te1_neg.^(5/2)));

    % recompute fluxes using updated moments
    %%% SHOULDN'T WE JUST COMPUTE THESE VIA BASIS SEPARATION %%%
    f1_vals = maxwellian(n1, v_para, v_perp, u_para1, Ta1, R_const);
    f1_vals_low_rank = cell(Nx, 3);
    for spatialIndex = 1:Nx
        f = f1_vals(:, :, spatialIndex);
        [Vr, S, Vz] = svd2(f, rvals);
        r0 = min(r0, size(Vr, 2));
        Vr = Vr(:, 1:r0); S = S(1:r0, 1:r0); Vz = Vz(:, 1:r0);
        f1_vals_low_rank{spatialIndex, 1} = Vr;
        f1_vals_low_rank{spatialIndex, 2} = S;
        f1_vals_low_rank{spatialIndex, 3} = Vz;  
    end

    [nu_hat2, S_hat2, Q_hat2] = get_fluxes(f1_vals_low_rank, v_perp, v_para, R_const, x_min, x_max, dx, dv_perp, dv_para);
    nTe_hat2 = (([n_BC_left; n1].*[Tae_BC_left; T_e1]).*(u_para1_half_nodes > 0) + ([n1; n_BC_right].*[T_e1; Tae_BC_right]).*(u_para1_half_nodes <= 0)); % upwinding

    % shift bounds to get flux pos/neg
    nu_hat2_pos = nu_hat2(2:end); nu_hat2_neg = nu_hat2(1:end-1);

    S_hat2_pos = S_hat2(2:end); S_hat2_neg = S_hat2(1:end-1);
    Q_hat2_pos = Q_hat2(2:end); Q_hat2_neg = Q_hat2(1:end-1);
    nTe_hat2_pos = nTe_hat2(2:end); nTe_hat2_neg = nTe_hat2(1:end-1);

    % explicitly find n via Forward Euler
    n2 = n0 + (delta*dt*An_0) - (1-delta)*(dt/dx)*(nu_hat2_pos - nu_hat2_neg);
    n2_pos = [n2(2:end); n_BC_right]; n2_neg = [n_BC_left; n2(1:end-1)];

     % ---- init y_vec, R_norm ----
    % y = [n2.*u_para0; n2.*U0; T_e0];
    y = [nu1; nU1; T_e1];

    nu = y(1:Nx); 
    nU = y(Nx+1:2*Nx);
    T_e = y(2*Nx+1:end);
    Te_pos = [T_e(2:end); Tae_BC_right]; Te_neg = [Tae_BC_left; T_e(1:end-1)];
    nTe_pos = n2_pos.*Te_pos; nTe_neg = n2_neg.*Te_neg;

    n_pos = n2_pos; n_neg = n2_neg;
    kappa_pos = (3.2/(2*sqrt(2*me)))*((Te_pos.^(5/2) + T_e.^(5/2)));
    kappa_neg = (3.2/(2*sqrt(2*me)))*((T_e.^(5/2) + Te_neg.^(5/2))); 

    R1 = nu - (n0.*u_para0)...
            - (delta*dt*Anu_0)...
            + ((1-delta)*dt/dx)*(S_hat2_pos - S_hat2_neg) ... 
            - ((dt2*qa)/(2*dx*qe*ma)).*((n1_pos.*Te1_pos) - (n1_neg.*Te1_neg)) ...
            - ((dt1*qa)/(2*dx*qe*ma)).*((n2_pos.*Te_pos) - (n2_neg.*Te_neg));
    R2 = nU - (n0.*U0)...
            - (delta*dt*AnU_0) ...
            + ((1-delta)*dt/dx).*(Q_hat2_pos - Q_hat2_neg) ...
            - (((dt2*qa*nu1)./(2*dx*qe*n1)).*((n1_pos.*Te1_pos) - (n1_neg.*Te1_neg)))...
            - (((dt2.*3.*sqrt(2*me))./((ma.^2).*(T_e1.^(3/2)))) .* (((n1.^2).*T_e1) - ((ma/3).*((2.*n1.*nU1) - (nu1.^2))))) ...
            - (((dt1*qa*nu)./(2*dx*qe*n2)).*((n2_pos.*Te_pos) - (n2_neg.*Te_neg)))...
            - (((dt1.*3.*sqrt(2*me))./((ma.^2).*(T_e.^(3/2)))) .* (((n2.^2).*T_e) - ((ma/3).*((2.*n2.*nU) - (nu.^2)))));
    R3 = (n2.*T_e) - (n0.*T_e0) - (delta*dt*AnTe_0) + (1-delta)*dt*(5/(3*dx)).*((u_para1_half_nodes(2:end).*nTe_hat2_pos) - (u_para1_half_nodes(1:end-1).*nTe_hat2_neg)) ...
            - (((dt2*u_para1)./(3*dx)).*(nTe1_pos - nTe1_neg)) - (((2*dt2)/(3*dx.^2)) .* ((kappa1_pos.*(Te1_pos - T_e1)) - (kappa1_neg.*(T_e1 - Te1_neg)))) - (((dt2.*2.*sqrt(2*me))./(ma.*(T_e1.^(3/2)))) .* (((ma/3).*((2*n1.*nU1) - (nu1.^2))) - ((n1.^2).*T_e1))) ...
            - (((dt1*(nu./n2))./(3*dx)).*(nTe_pos - nTe_neg)) - (((2*dt1)/(3*dx.^2)) .* ((kappa_pos.*(Te_pos - T_e)) - (kappa_neg.*(T_e - Te_neg)))) - (((dt1.*2.*sqrt(2*me))./(ma.*(T_e.^(3/2)))) .* (((ma/3).*((2*n2.*nU) - (nu.^2))) - ((n2.^2).*T_e)));
    R = [R1; R2; R3];

    while max(abs(R)) > tol

        % define partial derivatives of residual
        nu_nu = spdiags(ones(Nx), 0, Nx, Nx);
        nu_nU = spdiags(zeros(Nx),0, Nx, Nx);
        nu_Te = gallery('tridiag', (dt1*qa*n2(1:end-1))/(2*dx*qe*ma), zeros(Nx,1), -(dt1*qa*n2(2:end))/(2*dx*qe*ma));
    
        nU_nu = diag( ((-dt1*qa)./(2*dx.*qe.*n2)).*((n_pos.*Te_pos) - (n_neg.*Te_neg)) );
        nU_nU = spdiags(ones(Nx), 0, Nx, Nx);
        nU_Te_mid = ((3*dt1*sqrt(2*me))./(2*(ma.^2))).*(((n2.^2)./(T_e.^(3/2))  + (ma./(T_e.^(5/2))).*((2*n2.*nU) - (nu.^2)) ));
        nU_Te = gallery('tridiag', (dt1.*qa.*nu(2:end).*(n2(1:end-1)./n2(2:end)))./(2*dx*qe), nU_Te_mid, -(dt1.*qa.*nu(1:end-1).*(n2(2:end)./n2(1:end-1)))/(2*dx*qe) );
    
        Te_nu = diag(-(dt1.*((n_pos.*Te_pos) - (n_neg.*Te_neg)))./(3*dx*n2));
        Te_nU = spdiags(zeros(Nx), 0, Nx, Nx);
        Te_Te_left = ((dt1.*nu.*n_neg)./(3*dx.*n2)) + (((dt1*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(Te_neg.^(3/2))/2).*(T_e - Te_neg)) - (T_e.^(5/2) + Te_neg.^(5/2)))); Te_Te_left = Te_Te_left(2:end);
        Te_Te_mid = n2 - (((dt1*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(T_e.^(3/2))./2).*(Te_pos - T_e)) - (Te_pos.^(5/2) + T_e.^(5/2)) - ((5.*(T_e.^(3/2))/2).*(T_e - Te_neg)) - (T_e.^(5/2) + Te_neg.^(5/2)))) - ((2*dt1*sqrt(2*me)./(2*ma)).*(((-ma./(T_e.^(5/2))).*(2*n2.*nU - (nu.^2)) + ((n2.^2)./(T_e.^(3/2))))));
        Te_Te_right = ((-dt1.*nu.*n_pos)./(3*dx.*n2)) - (((dt1*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(Te_pos.^(3/2))/2).*(Te_pos - T_e)) + (Te_pos.^(5/2) + T_e.^(5/2)))); Te_Te_right = Te_Te_right(1:end-1);
        Te_Te = gallery('tridiag', Te_Te_left, Te_Te_mid, Te_Te_right);

        P = [nu_nu, nu_nU, nu_Te;
             nU_nu, nU_nU, nU_Te;
             Te_nu, Te_nU, Te_Te]; % holy block matrix
     
        dy = -P\R; % solve for delta y at time l+1
        y = y + dy;

        % update y_vec stuff
        nu = y(1:Nx); 
        nU = y(Nx+1:2*Nx); 
        T_e = y(2*Nx+1:end);
        Te_pos = [T_e(2:end); Tae_BC_right]; Te_neg = [Tae_BC_left; T_e(1:end-1)];
        nTe_pos = n_pos.*Te_pos; nTe_neg = n_neg.*Te_neg;

        kappa_pos = (3.2/(2*sqrt(2*me)))*((Te_pos.^(5/2) + T_e.^(5/2)));
        kappa_neg = (3.2/(2*sqrt(2*me)))*((T_e.^(5/2) + Te_neg.^(5/2)));

        R1 = nu - (n0.*u_para0) - (delta*dt*Anu_0) + ((1-delta)*dt/dx)*(S_hat2_pos - S_hat2_neg) ... 
                - ((dt2*qa)/(2*dx*qe*ma)).*((n1_pos.*Te1_pos) - (n1_neg.*Te1_neg)) ...
                - ((dt1*qa)/(2*dx*qe*ma)).*((n2_pos.*Te_pos) - (n2_neg.*Te_neg));
        R2 = nU - (n0.*U0)      - (delta*dt*AnU_0) + ((1-delta)*dt/dx).*(Q_hat2_pos - Q_hat2_neg) ...
                - (((dt2*qa*nu1)./(2*dx*qe*n1)).*((n1_pos.*Te1_pos) - (n1_neg.*Te1_neg))) - (((dt2.*3.*sqrt(2*me))./((ma.^2).*(T_e1.^(3/2)))) .* (((n1.^2).*T_e1) - ((ma/3).*((2.*n1.*nU1) - (nu1.^2))))) ...
                - (((dt1*qa*nu)./(2*dx*qe*n2)).*((n2_pos.*Te_pos) - (n2_neg.*Te_neg))) - (((dt1.*3.*sqrt(2*me))./((ma.^2).*(T_e.^(3/2)))) .* (((n2.^2).*T_e) - ((ma/3).*((2.*n2.*nU) - (nu.^2)))));
        R3 = (n2.*T_e) - (n0.*T_e0) ...
                - (delta*dt*AnTe_0) ...
                + (1-delta)*dt*(5/(3*dx)).*((u_para1_half_nodes(2:end).*nTe_hat2_pos) - (u_para1_half_nodes(1:end-1).*nTe_hat2_neg)) ...
                - (((dt2*nu1./n1)./(3*dx)).*(nTe1_pos - nTe1_neg))...
                - (((2*dt2)/(3*dx.^2)) .* ((kappa1_pos.*(Te1_pos - T_e1)) - (kappa1_neg.*(T_e1 - Te1_neg))))...
                - (((2*dt2*sqrt(2*me))./(ma.*(T_e1.^(3/2)))) .* (((ma/3).*((2*n1.*nU1) - (nu1.^2))) - ((n1.^2).*T_e1))) ...
                - (((dt1*(nu./n2))./(3*dx)).*(nTe_pos - nTe_neg))...
                - (((2*dt1)/(3*dx.^2)) .* ((kappa_pos.*(Te_pos - T_e)) - (kappa_neg.*(T_e - Te_neg)))) ...
                - (((2*dt1*sqrt(2*me))./(ma.*(T_e.^(3/2)))) .* (((ma/3).*((2*n2.*nU) - (nu.^2))) - ((n2.^2).*T_e)));
        R = [R1; R2; R3];
    end

    % return moments (yay)
    nu2 = y(1:Nx);
    nU2 = y(Nx+1:2*Nx);
    T_e2 = y(2*Nx+1:end);

    u_para2 = nu2 ./ n2;
    T_a2 = (ma/3)*(2*nU2./n2 - (nu2.^2)./(n2.^2));
    Te_all = [T_ae_min; T_e; T_ae_max];
    kappaTx = (3.2/(2*sqrt(2*me)*dx))*(Te_all(1:end-1).^(2.5)+Te_all(2:end).^(2.5)).*(Te_all(2:end)-Te_all(1:end-1));
   
end

function [nu_hat, S_hat, Q_hat] = get_fluxes(f_vals_low_rank, V_perp, V_para, R_const, x_min, x_max, dx, dv_perp, dv_para)
    Nx = size(f_vals_low_rank, 1);
    [n0, u_para0, T0] = get_boundaries(x_min, x_max);

    leftBC = maxwellian(n0(x_min - dx/2), V_perp, V_para, u_para0(x_min - dx/2), T0(x_min - dx/2), R_const);
    rightBC = maxwellian(n0(x_max + dx/2), V_perp, V_para, u_para0(x_max + dx/2), T0(x_max + dx/2), R_const);

    nu_hat = zeros(Nx+1, 1);
    S_hat = zeros(Nx+1, 1);
    Q_hat = zeros(Nx+1, 1);

    v_para_split_idx = find(V_para(1, :) > 0, 1); % upwinding!
    for i = 0:Nx % loop through spatial nodes    
        f_hat = zeros(size(V_perp));
        if i == 0 % use left BC
            f_left = leftBC;
            f_right = f_vals_low_rank{i+1, 1}*f_vals_low_rank{i+1, 2}*f_vals_low_rank{i+1, 3}';
        elseif i == Nx % use right BC
            f_left = f_vals_low_rank{i, 1}*f_vals_low_rank{i, 2}*f_vals_low_rank{i, 3}';
            f_right = rightBC;
        else
            f_left = f_vals_low_rank{i, 1}*f_vals_low_rank{i, 2}*f_vals_low_rank{i, 3}';
            f_right = f_vals_low_rank{i+1, 1}*f_vals_low_rank{i+1, 2}*f_vals_low_rank{i+1, 3}';
        end
        f_hat(:, 1:v_para_split_idx-1) = f_right(:, 1:v_para_split_idx-1);
        f_hat(:, v_para_split_idx:end) = f_left(:, v_para_split_idx:end);

        nu_hat(i+1) = 2*pi*(sum(sum(V_para .* f_hat .* (V_perp.*dv_para.*dv_perp))));
        S_hat(i+1) = 2*pi*(sum(sum(V_para.^2 .* f_hat .* (V_perp.*dv_para.*dv_perp))));
        Q_hat(i+1) = pi*(sum(sum(V_para .* (V_para.^2 + V_perp.^2) .* f_hat .* (V_perp.*dv_para.*dv_perp))));
    end
end

% returns maxwellian functions for each of the Nx spatial nodes
function [f] = maxwellian(n, v_perp, v_para, u_para, T, R)
    Nx = numel(n);
    Nv = size(v_perp);
    f = zeros(Nv(1), Nv(2), Nx);
    for i = 1:Nx
        f(:, :, i) = (n(i)/((2*pi*R*T(i))^(3/2)))*exp(-((v_para-u_para(i)).^2 + v_perp.^2)/(2*R*T(i)));
    end
end

function [n0, u_para0, T_ae0] = get_boundaries(x_min, x_max)
    x_mid = x_min + (x_max - x_min)/2;
    n0 = @(x) 0.36*tanh(0.05*(x-x_mid)) + 0.64;
    u_para0 = @(x) -1.1154*tanh(0.05*(x-x_mid)) + 1.983;
    T_ae0 = @(x) 0.4424*tanh(0.05*(x-x_mid)) + 0.5576; % Ta = Te at boundary!
end


% ------- LoMaC Truncation -------
function [Vr, S, Vz, rank] = LoMaC(Vr, S, Vz, Rmat, Zmat, rvals, zvals, tolerance, rhoM, JzM, kappaM, wr, wz, c, w_norm_1_squared, w_norm_v_squared, w_norm_v2_squared)
    % LoMaC Truncates given maxwellian (assumed Low-Rank) to given tolerance while conserving
    % macroscopic quantities.

    Nr = numel(rvals); Nz = numel(zvals);
    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);

    % Step 1: Integrate to calculate macro quantities
    p = 2*pi*dr*dz*sum(sum(((Vr) * S * (Vz)') .* Rmat));
    J = 2*pi*dr*dz*sum(sum(((Vr) * S * (Vz.*zvals)') .* Rmat));
    k = pi*dr*dz*sum(sum((((Vr.*(rvals.^2)) * S * Vz') + (Vr* S * ((Vz.*(zvals.^2))'))) .* Rmat));

    % Step 2: Scale by maxwellian to ensure inner product is well defined
    % (f -> 0 as v -> infinity)
    
    f1_proj_S_mtx11 = (p / w_norm_1_squared) - ((2*k - c*p)*c / w_norm_v2_squared);
    f1_proj_S_mtx12 = (J / w_norm_v_squared);
    f1_proj_S_mtx13 = ((2*k - c*p) / w_norm_v2_squared);

    proj_basis_r = wr.*[ones(Nr, 1), rvals.^2];
    proj_basis_z = wz.*[ones(Nz, 1), zvals, zvals.^2];
    f1_proj_S_mtx   = [f1_proj_S_mtx11, f1_proj_S_mtx12, f1_proj_S_mtx13;
                       f1_proj_S_mtx13,               0,               0];

    % f2 = f - f1 (do it via SVD)
    f2_U = [Vr, proj_basis_r];
    f2_S = blkdiag(S, -f1_proj_S_mtx);
    f2_V = [Vz, proj_basis_z];

    % QR factorize
    [f2_Vr, f2_S, f2_Vz, ~] = truncate(f2_U, f2_S, f2_V, rvals, tolerance);

    f2 = f2_Vr*f2_S*f2_Vz';

    % compute Pn(Te(f)) to ensure moments are kept
    trun_f2_proj_S_mtx11 = 2*pi*dr*dz*((sum(sum(f2.*Rmat)) / w_norm_1_squared) - (c*sum(sum(f2.*(Rmat.^2 + Zmat.^2 - c) .* Rmat)) / w_norm_v2_squared));
    trun_f2_proj_S_mtx12 = 2*pi*dr*dz*((sum(sum(f2.*Rmat.*Zmat))) / w_norm_v_squared);
    trun_f2_proj_S_mtx13 = 2*pi*dr*dz*((sum(sum(f2.*(Rmat.^2 + Zmat.^2 - c) .* Rmat)) / w_norm_v2_squared));

    trun_f2_proj_S_mtx   = [trun_f2_proj_S_mtx11, trun_f2_proj_S_mtx12, trun_f2_proj_S_mtx13;
                    trun_f2_proj_S_mtx13,            0,            0];

    % compute fM
    fM_proj_S_mtx11 = (rhoM / w_norm_1_squared) - ((2*kappaM - c.*rhoM)*c / w_norm_v2_squared);
    fM_proj_S_mtx12 = (JzM / w_norm_v_squared);
    fM_proj_S_mtx13 = ((2*kappaM - c*rhoM) / w_norm_v2_squared);

    fM_proj_S_mtx   = [fM_proj_S_mtx11, fM_proj_S_mtx12, fM_proj_S_mtx13;
                       fM_proj_S_mtx13,               0,              0];

    f_mass_S = fM_proj_S_mtx - trun_f2_proj_S_mtx;

    [Vr, S, Vz, rank] = truncate([proj_basis_r, f2_Vr], blkdiag(f_mass_S, f2_S), [proj_basis_z, f2_Vz], rvals, 1e-14);


    % % Nakao comparison
    % w1 = exp(-rvals.^2); w2 = exp(-zvals.^2);                    %w=w1*w2
    % c1 = 2*pi*(dr*sum(w1.*rvals))*(dz*sum(w2));                         %||1||^2, <.,.>_w
    % c2 = 2*pi*(dr*sum(w1.*rvals))*(dz*sum(zvals.^2.*w2));               %||vz||^2, <.,.>_w
    % c = (dr*sum(rvals.^2.*w1.*rvals))/(dr*sum(w1.*rvals)) + (dz*sum(zvals.^2.*w2))/(dz*sum(w2));
    % c3 = 2*pi*dr*dz*sum(sum((Rmat.^2+Zmat.^2-c).^2.*exp(-Rmat.^2-Zmat.^2).*Rmat)); %||vr^2+vz^2-c||^2, <.,.>_w
    % 
    % S_f2P = 2*pi*dr*dz*[sum(sum(f2.*Rmat))/c1 - c*sum(sum(f2.*(Rmat.^2+Zmat.^2-c).*Rmat))/c3, sum(sum(f2.*Zmat.*Rmat))/c2, sum(sum(f2.*(Rmat.^2+Zmat.^2-c).*Rmat))/c3;...
    %     sum(sum(f2.*(Rmat.^2+Zmat.^2-c).*Rmat))/c3, 0, 0];
    % 
    % S_fM = [rhoM/c1 - ((2*kappaM-c*rhoM)/c3)*c, JzM/c2, (2*kappaM-c*rhoM)/c3;...
    %     (2*kappaM-c*rhoM)/c3, 0, 0];
    % 
    % % S_fM = S_fM - S_f2P;
    % sum(sum(fM_proj_S_mtx - S_fM))


end

function [Vr, S, Vz, rank] = truncate(Vr_aug, S_aug, Vz_aug, rvals, tolerance)
    [Qr, Rr] = qr2(Vr_aug, rvals); [Qz, Rz] = qr(Vz_aug, 0);
    [U, Sigma, V] = svd(Rr*S_aug*Rz', 0); 
    rank = find(diag(Sigma) > tolerance, 1, 'last');
    Vr = Qr*U(:, 1:rank);
    S = Sigma(1:rank, 1:rank);
    Vz = Qz*V(:, 1:rank);
end

function [Vr, S, Vz, rank] = truncate_svd(Vr, S, Vz, tolerance)
    [U, Sigma, V] = svd(S, 0);
    rank = find(diag(Sigma) > tolerance, 1, 'last');
    if (sum(rank) == 0)
        rank = 1;
    end
    Vr = Vr*U(:, 1:rank);
    S = Sigma(1:rank, 1:rank);
    Vz = Vz*V(:, 1:rank);
end

function [] = PlotF(f_vals_low_rank, Xmat, Rmat, Zmat2, Nz)
    Nx = size(Xmat, 1);
    F = zeros(Nx, Nz);
    for i = 1:Nx
        f = f_vals_low_rank{i,1}*f_vals_low_rank{i,2}*f_vals_low_rank{i,3}';
        F(i,:) = 2*pi*sum(f.*Rmat, 1);
    end

    figure(9); clf;
    surf(Xmat, Zmat2, F); shading interp;
    % contour(Xmat, Zmat2, F, 30,'linewidth',1.2); shading interp;
    xlabel('x');ylabel('v_{||}');view(2);colorbar;
    drawnow;
end















function [Vr_nn,S_nn,Vz_nn,r_nn] = LoMaC_trun(Vr_hat,S_hat,Vz_hat,dr,dz,rvals,zvals,Rmat,Zmat,Nr,Nz,w1,w2,tol,c1,c2,c3,c,rhoM,JzM,kappaM)
fhat = Vr_hat*S_hat*Vz_hat';
rhoH = 2*pi*dr*dz*sum(sum(fhat.*Rmat));
JzH = 2*pi*dr*dz*sum(sum(Zmat.*fhat.*Rmat));
kappaH = 0.5*2*pi*dr*dz*sum(sum((Rmat.^2+Zmat.^2).*fhat.*Rmat));
% Compute f1
Vr_f1 = w1.*[ones(Nr,1),rvals.^2];
Vz_f1 = w2.*[ones(Nz,1),zvals,zvals.^2];
S_f1 = [rhoH/c1 - ((2*kappaH-c*rhoH)/c3)*c, JzH/c2, (2*kappaH-c*rhoH)/c3;...
        (2*kappaH-c*rhoH)/c3, 0, 0];
% Compute and truncate f2=f-f1
[Qr,Rr] = qr2([Vr_hat,Vr_f1],rvals);
[Qz,Rz] = qr([Vz_hat,Vz_f1],0);
[U,S,V] = svd(Rr*blkdiag(S_hat,-S_f1)*Rz',0);
r_f2 = find(diag(S)>tol,1,'last');
Vr_f2 = Qr*U(:,1:r_f2);
Vz_f2 = Qz*V(:,1:r_f2);
S_f2 = S(1:r_f2,1:r_f2);
% Compute PN(trun(f2)), to subtract later to ensure zero moments.
f2 = Vr_f2*S_f2*Vz_f2';
S_f2P = 2*pi*dr*dz*[sum(sum(f2.*Rmat))/c1 - c*sum(sum(f2.*(Rmat.^2+Zmat.^2-c).*Rmat))/c3, sum(sum(f2.*Zmat.*Rmat))/c2, sum(sum(f2.*(Rmat.^2+Zmat.^2-c).*Rmat))/c3;...
        sum(sum(f2.*(Rmat.^2+Zmat.^2-c).*Rmat))/c3, 0, 0];
% Compute fM
% Vr_fM = Vr_f1;
% Vz_fM = Vz_f1;
S_fM = [rhoM/c1 - ((2*kappaM-c*rhoM)/c3)*c, JzM/c2, (2*kappaM-c*rhoM)/c3;...
        (2*kappaM-c*rhoM)/c3, 0, 0];
% Compute fM-PN(trun(f2)).
% Since they share the same O.N. basis, only add(subtract) core matrices.
% Redefine fM = fM-PN(trun(f2))
S_fM = S_fM - S_f2P;
% Final solution
[Qr,Rr] = qr2([Vr_f1,Vr_f2],rvals);
[Qz,Rz] = qr([Vz_f1,Vz_f2],0);
[U,S,V] = svd(Rr*blkdiag(S_fM,S_f2)*Rz',0);
r_nn = find(diag(S)>1.0e-14,1,'last'); %Ensure square
Vr_nn = Qr*U(:,1:r_nn);
Vz_nn = Qz*V(:,1:r_nn);
S_nn = S(1:r_nn,1:r_nn);
end

function [n_out,u_out,Ti_out,Te_out,nuhat,Shat,Qhat,uhat,nTehat,KTxhat] = QN(f_vals_low_rank,nvals,uvals,Tivals,Tevals,q,qe,m,me,Nx,dx,dt,leftBC,rightBC,dvperp,dvpar,vpar,Vperp,Vpar,n_IC,u_IC,Te_IC)
sz = size(Vperp); %size of array in velocity space


% Compute numerical flux with upwinding
nuhat = zeros(Nx+1,1);     %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
Shat = zeros(Nx+1,1);      %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
Qhat = zeros(Nx+1,1);      %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
nTehat = zeros(Nx+1,1);    %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
KTxhat = zeros(Nx+1,1);    %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
pos_neg = find(vpar>0,1);  %first entry for which vpar is positive


% At x_{1/2}=0
fleft = leftBC;                                  %at x_0=0-dx/2, ghost point
fright = f_vals_low_rank{1,1}*f_vals_low_rank{1,2}*f_vals_low_rank{1,3}';   %at x_1
fhat = zeros(sz(1),sz(2));
fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
nuhat(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
Shat(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
Qhat(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));
for i = 1:Nx-1 % At x_{i+1/2}, interior cell boundaries
    fleft = f_vals_low_rank{i,1}*f_vals_low_rank{i,2}*f_vals_low_rank{i,3}';         %at x_i
    fright = f_vals_low_rank{i+1,1}*f_vals_low_rank{i+1,2}*f_vals_low_rank{i+1,3}';  %at x_{i+1}
    fhat = zeros(sz(1),sz(2));
    fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
    fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
    nuhat(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
    Shat(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
    Qhat(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));
end
% At x_{Nx+1/2}=200
fleft = f_vals_low_rank{Nx,1}*f_vals_low_rank{Nx,2}*f_vals_low_rank{Nx,3}'; %at x_Nx
fright = rightBC;                                %at x_{Nx+1}=200+dx/2, ghost point
fhat = zeros(sz(1),sz(2));
fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
nuhat(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
Shat(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
Qhat(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));


uhat = ([u_IC(0-dx/2);uvals] + [uvals;u_IC(200+dx/2)])/2; %u_{i+1/2}=(u_i+u_{i+1})/2. Size (Nx+1)x1 for Nx+1 cell boundaries
% At x_{1/2}=0
if uhat(1)>0
    nTehat(1) = n_IC(0-dx/2)*Te_IC(0-dx/2);
else
    nTehat(1) = nvals(1)*Tevals(1);
end
for i = 1:Nx-1 % At x_{i+1/2}, interior cell boundaries
    if uhat(i+1)>0
        nTehat(i+1) = nvals(i)*Tevals(i);
    else
        nTehat(i+1) = nvals(i+1)*Tevals(i+1);
    end
end
% At x_{Nx+1/2}=200
if uhat(Nx+1)>0
    nTehat(Nx+1) = nvals(Nx)*Tevals(Nx);
else
    nTehat(Nx+1) = n_IC(200+dx/2)*Te_IC(200+dx/2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


err = 1;
k = 0; %number of iterations


n_tnn = nvals - (dt/dx)*(nuhat(2:end)-nuhat(1:end-1));
n_tnn_right = [n_tnn(2:Nx);n_IC(200+dx/2)]; %with BC for ghost cell
n_tnn_left = [n_IC(0-dx/2);n_tnn(1:Nx-1)];  %with BC for ghost cell


u_k = uvals; %k=0
nu_k = n_tnn.*u_k; %k=0
nU_k = (3*n_tnn.*Tivals/m + n_tnn.*uvals.^2)/2; %k=0
Te_k = Tevals; %k=0
x_k = [nu_k;nU_k;Te_k];


Te_k_right = [Te_k(2:Nx);Te_IC(200+dx/2)]; %at x_2, x_3, ..., x_N, x_{N+1} (ghost cell)
Te_k_left = [Te_IC(0-dx/2);Te_k(1:Nx-1)];  %at x_0 (ghost cell), x_1, ..., x_N
Rk_nu = nu_k - nvals.*uvals + (dt/dx)*(Shat(2:end)-Shat(1:end-1)) - (dt/(2*dx))*(q/(m*qe))*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
Rk_nU = nU_k - (3*nvals.*Tivals/m + nvals.*uvals.^2)/2 + (dt/dx)*(Qhat(2:end)-Qhat(1:end-1)) - ((3*dt*sqrt(2)*sqrt(me))/m^2)*((n_tnn.^2.*Te_k - (m/3)*(2*n_tnn.*nU_k-nu_k.^2))./(Te_k.^(1.5)))...
    - (dt/(2*dx))*(q/qe)*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
Rk_Te = n_tnn.*Te_k - nvals.*Tevals + ((5*dt)/(3*dx))*(uhat(2:end).*nTehat(2:end)-uhat(1:end-1).*nTehat(1:end-1)) - (dt/(3*dx))*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)...
    - (dt/(3*dx^2))*(3.2/sqrt(2*me))*((Te_k.^(2.5)+Te_k_right.^(2.5)).*(Te_k_right-Te_k) - (Te_k_left.^(2.5)+Te_k.^(2.5)).*(Te_k-Te_k_left))...
    - ((2*dt*sqrt(2*me))/m)*((m/3)*(2*n_tnn.*nU_k - nu_k.^2) - n_tnn.^2.*Te_k)./(Te_k.^(1.5));
Rk = -[Rk_nu;Rk_nU;Rk_Te];


tol = min(5e-12,max(abs(Rk))*5e-10);
%tol = 5.0e-4;


while err > tol
    Pk_nunu = eye(Nx,Nx);
    Pk_nunU = zeros(Nx,Nx);
    Pk_nuTe = zeros(Nx,Nx);
        Pk_nuTe(1:end-1,2:end) = Pk_nuTe(1:end-1,2:end) -(dt/(2*dx))*(q/(m*qe))*diag(n_tnn(2:end)); %superdiagonal
        Pk_nuTe(2:end,1:end-1) = Pk_nuTe(2:end,1:end-1) +(dt/(2*dx))*(q/(m*qe))*diag(n_tnn(1:end-1)); %subdiagonal


    Pk_nUnu = -((dt*q)/(2*dx*qe))*diag((n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)./n_tnn);
    Pk_nUnU = eye(Nx,Nx);
    Pk_nUTe = zeros(Nx,Nx);
        Pk_nUTe(1:end-1,2:end) = Pk_nUTe(1:end-1,2:end) -(dt/(2*dx))*(q/qe)*diag(nu_k(1:end-1).*n_tnn(2:end)./n_tnn(1:end-1)); %superdiagonal
        Pk_nUTe(2:end,1:end-1) = Pk_nUTe(2:end,1:end-1) +(dt/(2*dx))*(q/qe)*diag(nu_k(2:end).*n_tnn(1:end-1)./n_tnn(2:end)); %subdiagonal
        Pk_nUTe = Pk_nUTe - (3*dt*sqrt(2*me)/m^2)*diag( -0.5*(n_tnn.^2)./(Te_k.^(1.5)) - (m/3)*(-3/2)*(2*n_tnn.*nU_k-nu_k.^2)./(Te_k.^(2.5)) ); %diagonal


    Pk_Tenu = -(dt/(3*dx))*diag((n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)./n_tnn);
    Pk_TenU = zeros(Nx,Nx);
    Pk_TeTe = zeros(Nx,Nx);
        Pk_TeTe(1:end-1,2:end) = Pk_TeTe(1:end-1,2:end) -(dt/(3*dx))*diag(nu_k(1:end-1).*n_tnn(2:end)./n_tnn(1:end-1))...
            - (dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k(2:end).^(1.5).*(Te_k(2:end)-Te_k(1:end-1)) + (Te_k(1:end-1).^(2.5) + Te_k(2:end).^(2.5))); %superdiagonal
        Pk_TeTe(2:end,1:end-1) = Pk_TeTe(2:end,1:end-1) +(dt/(3*dx))*diag(nu_k(2:end).*n_tnn(1:end-1)./n_tnn(2:end))...
            + (dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k(1:end-1).^(1.5).*(Te_k(2:end)-Te_k(1:end-1)) - (Te_k(1:end-1).^(2.5) + Te_k(2:end).^(2.5))); %subdiagonal
        Pk_TeTe = Pk_TeTe + diag(n_tnn) - (dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k.^(1.5).*(Te_k_right - 2*Te_k + Te_k_left) - (Te_k_right.^(2.5) + 2*Te_k.^(2.5) + Te_k_left.^(2.5)))...
            - (2*dt*sqrt(2*me)/m)*diag(-1.5*(m/3)*((2*n_tnn.*nU_k-nu_k.^2)./(Te_k.^(2.5))) + 0.5*(n_tnn.^2)./(Te_k.^(1.5))); %diagonal


    Pk = [Pk_nunu, Pk_nunU, Pk_nuTe; Pk_nUnu, Pk_nUnU, Pk_nUTe; Pk_Tenu, Pk_TenU, Pk_TeTe];


    dx_kk = Pk\Rk;
    x_kk = x_k + dx_kk;


    nu_kk = x_kk(1:Nx);
    nU_kk = x_kk(Nx+1:2*Nx);
    Te_kk = x_kk(2*Nx+1:3*Nx);


    nu_k = nu_kk;
    nU_k = nU_kk;
    Te_k = Te_kk;
    x_k = x_kk;




    Te_k_right = [Te_k(2:Nx);Te_IC(200+dx/2)]; %at x_2, x_3, ..., x_N, x_{N+1} (ghost cell)
    Te_k_left = [Te_IC(0-dx/2);Te_k(1:Nx-1)];  %at x_0 (ghost cell), x_1, ..., x_N
    Rk_nu = nu_k - nvals.*uvals + (dt/dx)*(Shat(2:end)-Shat(1:end-1)) - (dt/(2*dx))*(q/(m*qe))*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
    Rk_nU = nU_k - (3*nvals.*Tivals/m + nvals.*uvals.^2)/2 + (dt/dx)*(Qhat(2:end)-Qhat(1:end-1)) - ((3*dt*sqrt(2)*sqrt(me))/m^2)*((n_tnn.^2.*Te_k - (m/3)*(2*n_tnn.*nU_k-nu_k.^2))./(Te_k.^(1.5)))...
        - (dt/(2*dx))*(q/qe)*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
    Rk_Te = n_tnn.*Te_k - nvals.*Tevals + ((5*dt)/(3*dx))*(uhat(2:end).*nTehat(2:end)-uhat(1:end-1).*nTehat(1:end-1)) - (dt/(3*dx))*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)...
        - (dt/(3*dx^2))*(3.2/sqrt(2*me))*((Te_k.^(2.5)+Te_k_right.^(2.5)).*(Te_k_right-Te_k) - (Te_k_left.^(2.5)+Te_k.^(2.5)).*(Te_k-Te_k_left))...
        - ((2*dt*sqrt(2*me))/m)*((m/3)*(2*n_tnn.*nU_k - nu_k.^2) - n_tnn.^2.*Te_k)./(Te_k.^(1.5));
    Rk = -[Rk_nu;Rk_nU;Rk_Te];


    err = max(abs(Rk));%norm(Rk,2)/Rn_err;


    k = k+1;
    if k==500
        disp('Completed 500 iterations');
        err
        return
    end    
end
%disp(k)
n_out = n_tnn;
u_out = nu_k./n_tnn;
Ti_out = (m/3)*(2*nU_k./n_tnn - (nu_k.^2)./(n_tnn.^2));
Te_out = Te_k;


% Compute thermal diffusion (at time t^{n+1}) for energy balance law.
Te_all = [Te_IC(0-dx/2);Te_k;Te_IC(200+dx/2)];
KTxhat = (3.2/(2*sqrt(2*me)*dx))*(Te_all(1:end-1).^(2.5)+Te_all(2:end).^(2.5)).*(Te_all(2:end)-Te_all(1:end-1));
end


