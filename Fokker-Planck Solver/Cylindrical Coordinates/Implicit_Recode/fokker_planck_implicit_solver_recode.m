clc; clear variables; close all;

% initial rank
r0 = 10;
% time-stepping method: 1=B.Euler, 2=DIRK2, 3=DIRK3
method = '2';
tolerance = 1e-6;

% mesh parameters
vmin = 0; vmax = 14;
zmin = -16; zmax = 16;
Nr = 100; Nz = 100;
tf = 15;

[Rmat, Zmat, dr, dz] = GetRZ(vmin, vmax, zmin, zmax, Nr, Nz);
rvals = Rmat(:, 1);
zvals = Zmat(1, :)';

% initial conditions
f_M = @(vr,vz,n,ur,uz,T,R) (n/(2*pi*R*T)^(3/2))*exp(-(vr.^2+(vz-uz).^2)/(2*R*T)); %assumes ur=0!!!

n1 = 2;
u1r = 0;
u1z = -0.5;
T1 = 2;
n2 = 1;
u2r = 0;
u2z = 0.9;
T2 = 1;
R = 1/6;

n = n1+n2;
ur = n1*u1r/n + n2*u2r/n;
uz = n1*u1z/n + n2*u2z/n;
T = n1*T1/n + n2*T2/n + (n1*(u1r^2+u1z^2) + n2*(u2r^2+u2z^2))/(3*R*n) - ((n1*u1r+n2*u2r)^2 + (n1*u1z+n2*u2z)^2)/(3*R*(n)^2);
D = R*T;

f0 = @(vr,vz) f_M(vr,vz,n1,u1r,u1z,T1,R) + f_M(vr,vz,n2,u2r,u2z,T2,R); % IC
f = f0(Rmat,Zmat);

% discrete moments of f0
rho0 = sum(sum(f.*Rmat))*2*pi*dr*dz;
Jz0  = sum(sum(f.*Rmat.*Zmat))*2*pi*dr*dz;
kappa0 = sum(sum(f.*Rmat.*((Rmat.^2 + Zmat.^2)/2)))*2*pi*dr*dz;

% f_inf = f_M(Rmat,Zmat,n,ur,uz,T,R); % Equilibrium solution
% f_inf = ((sum(sum(f.*Rmat))*2*pi*dr*dz)/(2*pi*dr*dz*sum(sum(f_inf.*Rmat))))*f_inf;
[f_inf] = QCM(rho0, Jz0, kappa0, R, Rmat, Zmat);

% moments at equilibrium
rhoM = 2*pi*dr*dz*sum(sum(f_inf.*Rmat));
JzM = 2*pi*dr*dz*sum(sum(Zmat.*f_inf.*Rmat));
kappaM = pi*dr*dz*sum(sum((Rmat.^2+Zmat.^2).*f_inf.*Rmat));
% rhoM = rho0;
% JzM = Jz0;
% kappaM = kappa0;

lambdavals = 1;%(0.2:0.1:6)';
errors = zeros(numel(lambdavals), 3);
methods = ['1', '2'];
for x = 1:1
method = methods(x);

for k = 1:numel(lambdavals)
dt = lambdavals(k)/((1/dr) + (1/dz));
tvals = (0:dt:tf)';
if tvals(end) ~= tf
    tvals = [tvals; tf];
end
Nt = numel(tvals);

Ar = @(u, w, t) w; %cell centers
Br = @(u, w, t) w.*(w - u); %evaluated on cell boundaries
Cr = @(u, w, t) D*w; %evaluated cell boundaries
Az = @(u, w, t) w.^0;
Bz = @(u, w, t) w - u;
Cz = @(u, w, t) D*w.^0;

% init bases
[Vr, S, Vz] = svd2(f, rvals);
r0 = min(r0, size(Vr, 2));
Vr = Vr(:, 1:r0); S = S(1:r0, 1:r0); Vz = Vz(:, 1:r0);
f = Vr*S*Vz';

% store rank, mass, momentum, energy, l1 decay, etc...
l1 = zeros(Nt, 1);
mass = zeros(Nt, 1);
uzvals = zeros(Nt, 1);
E = zeros(Nt, 1);
relative_entropy = zeros(Nt, 1);
min_vals = zeros(Nt, 1);
ranks = zeros(Nt, 1);

l1(1) = 2*pi*dr*dz*sum(sum(abs(Rmat .* (f - f_inf))));
mass(1) = 2*pi*dr*dz*sum(sum(Rmat .* f));
uzvals(1) = 2*pi*dr*dz*sum(sum(f .* Rmat .* Zmat));
E(1) = pi*dr*dz*sum(sum(f .* (Rmat.^2 + Zmat.^2) .* Rmat));
relative_entropy(1) = 2*pi*dr*dz*sum(sum(Rmat .* f.*(log(f./f_inf))));
min_vals(1) = min(min(f));
ranks(1) = r0;

% time-stepping loop
for n = 2:Nt
    tval = tvals(n);
    dt = tval - tvals(n-1);
    switch(method)
        case '1'
            [Vr, S, Vz, rank] = BackwardEulerTimestep(Vr, S, Vz, ur, uz, dt, tval, rvals, zvals, Rmat, Zmat, Ar, Az, Br, Bz, Cr, Cz, tolerance, rhoM, JzM, kappaM);
        case '2'
            [Vr, S, Vz, rank] = DIRK2Timestep(Vr, S, Vz, ur, uz, dt, tval, rvals, zvals, Rmat, Zmat, Ar, Az, Br, Bz, Cr, Cz, tolerance, rhoM, JzM, kappaM);
        case '3'
            [Vr, S, Vz, rank] = DIRK3Timestep(Vr, S, Vz, ur, uz, dt, tval, rvals, zvals, Rmat, Zmat, Ar, Az, Br, Bz, Cr, Cz, tolerance, rhoM, JzM, kappaM);
    end

    f = Vr*S*Vz';
    
    l1(n) = 2*pi*dr*dz*sum(sum(abs(Rmat .* (f - f_inf))));
    mass(n) = 2*pi*dr*dz*sum(sum(Rmat .* f));    
    uzvals(n) = 2*pi*dr*dz*sum(sum(f .* Rmat .* Zmat));
    E(n) = pi*dr*dz*sum(sum(f .* (Rmat.^2 + Zmat.^2) .* Rmat));
    relative_entropy(n) = 2*pi*dr*dz*sum(sum(Rmat .* f.*(log(f./f_inf))));
    min_vals(n) = min(min(f));
    ranks(n) = rank;
end

errors(k, x) = 2*pi*dr*dz*sum((sum(Rmat .* abs(f - f_inf)))); % L1 error

end
end

% figure(8); clf;
% loglog(lambdavals, errors(:, 1), 'b-', 'LineWidth', 1.5); hold on;
% loglog(lambdavals, 5e-10*(lambdavals .^ 1), 'b--', 'LineWidth', 1.5);
% loglog(lambdavals, errors(:, 2), 'r-', 'LineWidth', 1.5); hold on;
% loglog(lambdavals, 5e-10*(lambdavals .^ 2), 'r--', 'LineWidth', 1.5);
% xlabel('\lambda'); ylabel('Accuracy'); title('Accuracy plot');
% legend('Backward Euler', '$\mathcal{O}(1)$', 'DIRK2','$\mathcal{O}(2)$', 'interpreter', 'latex');
% return

figure(1); clf; surf(Rmat, Zmat, f);
colorbar; shading interp;
legend(sprintf('N_r = %s', num2str(Nr, 3)), 'Location','northwest');
xlabel('V_r'); ylabel('V_z'); zlabel('U'); title([sprintf('Backward Euler approximation of 0D2V Fokker-Planck system at time %s', num2str(tf, 4))]);

figure(2); clf; surf(Rmat, Zmat, f_inf);
colorbar; shading interp;
xlabel('V_r'); ylabel('V_z'); zlabel('f(V_r, V_z, t)'); title([sprintf('f_{exact} at time t=%s', num2str(tf, 4))]);
 
% l1 decay
figure(3); clf; semilogy(tvals, l1, 'black-', 'LineWidth', 1.5);
xlabel('t'); ylabel('L_1(f(V_r, V_z))'); title('L_1 decay of numerical solution over time');

figure(8); clf; plot(tvals(2:end), abs(uzvals(2:end)-uzvals(1)), 'red-', 'LineWidth', 1.5);
xlabel('t'); ylabel('Absolute error (Uz)'); title('Absolute error of bulk velocity over time');
figure(9); clf; plot(tvals(2:end), abs(E(2:end)-E(1))/E(1), 'red-', 'LineWidth', 1.5);
xlabel('t'); ylabel('Relative error (Energy)'); title('Relative energy of numerical solution over time');

% Positivity
figure(4); clf; plot(tvals, min_vals, 'green-', 'LineWidth', 1.5);
xlabel('t'); ylabel('min(f(V_r, V_z))'); title('Minimum values of numerical solution over time');

% Relative entropy
figure(5); clf; semilogy(tvals, relative_entropy, 'magenta-', 'LineWidth', 1.5);
xlabel('t'); ylabel('Relative entropy'); title('Relative entropy of numerical solution over time');

% Mass
figure(6); clf; plot(tvals(2:end), abs(mass(2:end)-mass(1))/mass(1), 'red-', 'LineWidth', 1.5);
xlabel('t'); ylabel('relative mass'); title('Relative mass of numerical solution over time');

% Rank plot
figure(7); clf;
plot(tvals, ranks, 'black-', 'LineWidth', 1.5); hold on;
% plot(ranks2(:, 1), ranks2(:, 2), 'blue-', 'LineWidth', 1.5);
% plot(ranks3(:, 1), ranks3(:, 2), 'green-', 'LineWidth', 1.5);
xlabel('time'); ylabel('rank'); title('Rank plot over time');
legend('Backward Euler', 'RK2', 'RK3');





%%%%%% FUNCTIONS %%%%%%
function [Flux] = GetFlux(A, B, C, u, xvals, t, dx)
    N = numel(xvals);
    A = A(u, xvals, t);
    B = B(u, xvals(1:N-1) + dx/2, t);
    C = C(u, xvals(1:N-1) + dx/2, t);
    
    w = dx*B./(C + 1e-14);
    delta = (1 ./ w) - (1 ./ (exp(w) - 1));

    F1 = -((1/dx)*C - delta.*B);
    F2 = (1 - delta).*B + (1/dx).*C;

    F_pos = spdiags([F1;0], 0, N, N) + spdiags([0;F2], 1, N, N);
    F_neg = spdiags([0; F2], 0, N, N) + spdiags(F1, -1, N, N);

    Flux = diag(1./A)*(1/dx)*(F_pos - F_neg);
end

function [Vr, S, Vz, rank] = BackwardEulerTimestep(Vr0, S0, Vz0, ur, uz, dt, tval, rvals, zvals, Rmat, Zmat, Ar, Az, Br, Bz, Cr, Cz, tolerance, rhoM, JzM, kappaM)

    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);
    Nr = numel(rvals);
    Nz = numel(zvals);

    % ASK ABOUT T_NN VS T_N!!!
    Fr1 = GetFlux(Ar, Br, Cr, ur, rvals, tval, dr);
    Fz1 = GetFlux(Az, Bz, Cz, uz, zvals, tval, dz);

    Vr0_star = Vr0;
    Vz0_star = Vz0;

    K0 = Vr0*S0;
    L0 = Vz0*S0';

    K1 = sylvester(eye(Nr) - (dt*Fr1), -dt*(Fz1*Vz0_star)'*Vz0_star, K0);
    L1 = sylvester(eye(Nz) - (dt*Fz1), -dt*(Fr1*Vr0_star)'*(rvals.*Vr0_star), L0);

    [Vr1_ddagger, ~] = qr2(K1, rvals);
    [Vz1_ddagger, ~] = qr(L1, 0);

    [Vr1_hat, Vz1_hat] = reduced_augmentation([Vr1_ddagger, Vr0], [Vz1_ddagger, Vz0], rvals);

    S1_hat = sylvester((speye(size(Vr1_hat, 2)) - (dt*((rvals .* Vr1_hat)')*(Fr1*Vr1_hat))), -dt*(Fz1*Vz1_hat)'*Vz1_hat, ((rvals .* Vr1_hat)'*Vr0)*S0*((Vz0')*Vz1_hat));
    % [Vr, S, Vz, rank] = truncate_svd(Vr1_hat, S1_hat, Vz1_hat, tolerance);
    [Vr, S, Vz, rank] = LoMaC(Vr1_hat, S1_hat, Vz1_hat, Rmat, Zmat, rvals, zvals, tolerance, rhoM, JzM, kappaM);
end

function [Vr, S, Vz, rank] = DIRK2Timestep(Vr0, S0, Vz0, ur, uz, dt, tval, rvals, zvals, Rmat, Zmat, Ar, Az, Br, Bz, Cr, Cz, tolerance, rhoM, JzM, kappaM)
    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);
    Nr = numel(rvals);
    Nz = numel(zvals);

    Fr1 = GetFlux(Ar, Br, Cr, ur, rvals, tval, dr);
    Fz1 = GetFlux(Az, Bz, Cz, uz, zvals, tval, dz);

    nu = 1-(sqrt(2)/2);

    % Stage 1: Backward Euler
    [Vr1, S1, Vz1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, ur, uz, nu*dt, tval, rvals, zvals, Rmat, Zmat, Ar, Az, Br, Bz, Cr, Cz, tolerance, rhoM, JzM, kappaM);

    W0 = (Vr0*S0*(Vz0')) + ((1-nu)*dt*(((Fr1*(Vr1)*S1*(Vz1')) + ((Vr1)*S1*((Fz1*Vz1)')))));

    % Reduced Augmentation
    % Predict V_dagger using B. Euler for second stage
    [Vr1_dagger, ~, Vz1_dagger, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, ur, uz, dt, tval, rvals, zvals, Rmat, Zmat, Ar, Az, Br, Bz, Cr, Cz, tolerance, rhoM, JzM, kappaM);
    [Vr_star, Vz_star] = reduced_augmentation([Vr1_dagger, Vr1, Vr0], [Vz1_dagger, Vz1, Vz0], rvals);

    % Stage 2: KLS Steps
    Fr2 = GetFlux(Ar, Br, Cr, ur, rvals, tval+dt, dr);
    Fz2 = GetFlux(Az, Bz, Cz, uz, zvals, tval+dt, dz);
   
    % K/L-Step
    K1 = sylvester(eye(Nr) - (nu*dt*Fr2), -nu*dt*(Fz2*Vz_star)'*Vz_star, W0*Vz_star);
    L1 = sylvester(eye(Nz) - (nu*dt*Fz2), -nu*dt*(Fr2*Vr_star)'*(rvals .* Vr_star), W0'*(rvals .* Vr_star));

    % Get bases
    [Vr_ddagger, ~] = qr2(K1, rvals); [Vz_ddagger, ~] = qr(L1, 0);

    % Reduced Augmentation
    [Vr1_hat, Vz1_hat] = reduced_augmentation([Vr_ddagger, Vr1, Vr0], [Vz_ddagger, Vz1, Vz0], rvals);

    % S-Step
    S1_hat = sylvester(eye(size(Vr1_hat, 2)) - (nu*dt*((rvals .* Vr1_hat)')*Fr2*Vr1_hat), -nu*dt*(Fz2*Vz1_hat)'*Vz1_hat, ((rvals .* Vr1_hat)')*W0*Vz1_hat);
    % [Vr, S, Vz, rank] = truncate_svd(Vr, S, Vz, tolerance);
    [Vr, S, Vz, rank] = LoMaC(Vr1_hat, S1_hat, Vz1_hat, Rmat, Zmat, rvals, zvals, tolerance, rhoM, JzM, kappaM);
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

function [Vr, Vz] = reduced_augmentation(Vr_aug, Vz_aug, rvals)
    tolerance = 1e-12;
    [Qr, Rr] = qr2(Vr_aug, rvals);
    [Qz, Rz] = qr(Vz_aug, 0);
    [Ur, Sr, ~] = svd(Rr, 0);
    [Uz, Sz, ~] = svd(Rz, 0);
    rr = find(diag(Sr) > tolerance, 1, 'last');
    rz = find(diag(Sz) > tolerance, 1, 'last');
    rank = max(rr, rz);
    rank = min(rank, min(size(Qr, 2), size(Qz, 2)));
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

function [rmat, zmat, dr, dz] = GetRZ(vmin, vmax, zmin, zmax, Nv, Nz)
    rvals = linspace(vmin, vmax, Nv+1)';
    zvals = linspace(zmin, zmax, Nz+1)';
    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);
    rmid = rvals(1:end-1) + (dr/2);
    zmid = zvals(1:end-1) + (dz/2);
    [rmat, zmat] = meshgrid(rmid, zmid);
    rmat = rmat';
    zmat = zmat';
end


% ------- LoMaC Truncation -------
function [Vr, S, Vz, rank] = LoMaC(Vr, S, Vz, Rmat, Zmat, rvals, zvals, tolerance, rhoM, JzM, kappaM)
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
    wr = exp(-(rvals.^2));
    wz = exp(-(zvals.^2));

    % Step 3: Orthogonal projection
    % bases: 1, v, v.^2 - c
    c = (dr*sum(rvals.^2.*wr.*rvals))/(dr*sum(wr.*rvals)) + (dz*sum(zvals.^2.*wz))/(dz*sum(wz));

    w_norm_1_squared = 2*pi*dr*dz*sum(rvals .* wr)*sum(wz);
    w_norm_v_squared = 2*pi*dr*dz*sum(rvals .* wr)*sum(zvals.^2 .* wz);
    w_norm_v2_squared = 2*pi*dr*dz*sum(sum((Rmat.^2 + Zmat.^2 - c).^2 .* exp(-Rmat.^2 - Zmat.^2) .* Rmat));
    
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