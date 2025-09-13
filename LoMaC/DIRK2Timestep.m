% Completes 1 DIRK2 timestep, updating the bases Vx0, S0, Vy0 --> Vx, S, Vy
% Inputs: bases: Vx0, S0, Vy0
%         timestep size: dt
%         differentiation matrices: Dxx, Dyy
%         final truncation tolerance: tolerance
% Outputs: updated vases Vx, S, Vy
%          rank: r1


function [Vr, S, Vz, r1] = DIRK2Timestep(Vr0, S0, Vz0, dt, tval, R, Z, rvals, zvals, Br, Bz, Drr, Dzz, tolerance, f_inf)

    Fr1 = GetRadialFluxSPCC(tval, rvals, Br, Drr);
    Fz1 = GetZFluxSPCC(tval, zvals, Bz, Dzz);

    nu = 1-(sqrt(2)/2);

    % Stage 1: Backward Euler
    [Vr1, S1, Vz1, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, nu*dt, tval, R, Z, rvals, zvals, Br, Bz, Drr, Dzz, tolerance, f_inf);

    W0 = (Vr0*S0*(Vz0')) + ((1-nu)*dt*(((Fr1*(Vr1)*S1*(Vz1')) + ((Vr1)*S1*((Fz1*Vz1)')))));

    % Reduced Augmentation
    % Predict V_dagger using B. Euler for second stage
    [Vr1_dagger, ~, Vz1_dagger, ~] = BackwardEulerTimestep(Vr0, S0, Vz0, dt, tval, R, Z, rvals, zvals, Br, Bz, Drr, Dzz, tolerance, f_inf);
    [Vr_star, Vz_star] = reduced_augmentation([Vr1_dagger, Vr1, Vr0], [Vz1_dagger, Vz1, Vz0], rvals);

    % Stage 2: KLS Steps
    Fr2 = GetRadialFluxSPCC(tval + dt, rvals, Br, Drr);
    Fz2 = GetZFluxSPCC(tval + dt, zvals, Bz, Dzz);
   
    % K/L-Step
    K1 = sylvester(eye(size(Fr2)) - (nu*dt*Fr2), -nu*dt*(Fz2*Vz_star)'*Vz_star, W0*Vz_star);
    L1 = sylvester(eye(size(Fz2)) - (nu*dt*Fz2), -nu*dt*(Fr2*Vr_star)'*(rvals .* Vr_star), (W0')*(rvals .* Vr_star));

    % Get bases
    [Vr_ddagger, ~] = qr2(K1, rvals); [Vz_ddagger, ~] = qr(L1, 0);

    % Reduced Augmentation
    [Vr, Vz] = reduced_augmentation([Vr_ddagger, Vr1, Vr0], [Vz_ddagger, Vz1, Vz0], rvals);

    % S-Step
    S = sylvester(eye(size(Vr, 2)) - (nu*dt*((rvals .* Vr)')*Fr2*Vr), -nu*dt*(Fz2*Vz)'*Vz, ((rvals .* Vr)')*W0*Vz);
    % [Vr, S, Vz, r1] = truncate_svd(Vr, S, Vz, tolerance);
    [Vr, S, Vz, r1] = LoMaC(Vr, S, Vz, R, Z, rvals, zvals, tolerance, f_inf);
end