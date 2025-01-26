% Completes 1 DIRK3 timestep, updating the bases Vx0, S0, Vy0 --> Vx3, S3, Vy3
% Inputs: bases: Vx0, S0, Vy0
%         timestep size: dt
%         differentiation matrices: Dxx, Dyy
%         final truncation tolerance: tolerance
% Outputs: updated vases Vx3, S3, Vy3
%          rank: r3

function [Vr3, S3, Vz3, r3] = DIRK3(Vr0, S0, Vz0, rvals, dt, Drr, Dzz, tolerance)
    
    % RK butcher table values
    nu = 0.435866521508459;
    beta1 = -(3/2)*(nu^2) + (4*nu) - (1/4);
    beta2 = (3/2)*(nu^2) - (5*nu) + (5/4);
    Nrr = size(Drr);
    Nzz = size(Dzz);
    
    % Stage 1: Backward Euler
    [Vr1, S1, Vz1, ~] = BackwardEuler(Vr0, S0, Vz0, rvals, nu*dt, Drr, Dzz, tolerance);
    [Vr_dagger1, ~, Vz_dagger1, ~] = BackwardEuler(Vr0, S0, Vz0, rvals, ((1+nu)/2)*dt, Drr, Dzz, tolerance);

    Y1 = (((Drr*Vr1*S1*(Vz1')) + (Vr1*S1*((Dzz*Vz1)'))));
    W1 = (Vr0*S0*(Vz0')) + (((1-nu)/2)*dt*Y1);

    % Reduced Augmentation
    [Vr_star1, Vz_star1] = reduced_augmentation([Vr_dagger1, Vr1, Vr0], [Vz_dagger1, Vz1, Vz0], rvals);

    % Stage 2:
    % K/L-Step
    K2 = sylvester(eye(Nrr) - (nu*dt*Drr), -nu*dt*(Dzz*Vz_star1)'*Vz_star1, W1*Vz_star1);
    L2 = sylvester(eye(Nzz) - (nu*dt*Dzz), -nu*dt*(Drr*(rvals .* Vr_star1))'*Vr_star1, (W1')*(rvals .* Vr_star1));

    % Get bases
    [Vr_ddagger2, ~] = qr2(K2, rvals); [Vz_ddagger2, ~] = qr(L2, 0);

    % Reduced Augmentation
    [Vr2, Vz2] = reduced_augmentation([Vr_ddagger2, Vr1, Vr0], [Vz_ddagger2, Vz1, Vz0], rvals);

    % S-Step
    S2 = sylvester(eye(size(Vr2, 2)) - (nu*dt*(rvals .* Vr2)'*Drr*Vr2), -nu*dt*(Dzz*Vz2)'*Vz2, ((rvals .* Vr2)')*W1*Vz2);
    [Vr2, S2, Vz2, ~] = truncate_svd(Vr2, S2, Vz2, tolerance);

    % Stage 3:
    % Predict V_dagger using B. Euler
    [Vr_dagger3, ~, Vz_dagger3, ~] = BackwardEuler(Vr0, S0, Vz0, rvals, dt, Drr, Dzz, tolerance);
    Y2 = (((Drr*Vr2*S2*(Vz2')) + (Vr2*S2*((Dzz*Vz2)'))));
    W2 = (Vr0*S0*(Vz0')) + (beta1*dt*Y1) + (beta2*dt*Y2);
      
    % Reduced augmentation
    [Vr_star3, Vz_star3] = reduced_augmentation([Vr_dagger3, Vr2, Vr1, Vr0], [Vz_dagger3, Vz2, Vz1, Vz0], rvals);

    % K/L-Step
    K3 = sylvester(eye(Nrr) - (nu*dt*Drr), -nu*dt*(Dzz*Vz_star3)'*Vz_star3, W2*Vz_star3);
    L3 = sylvester(eye(Nzz) - (nu*dt*Dzz), -nu*dt*(Drr*Vr_star3)'*(rvals .* Vr_star3), (W2')*(rvals .* Vr_star3));

    % Get bases
    [Vr_ddagger3, ~] = qr2(K3, rvals); [Vz_ddagger3, ~] = qr(L3, 0);

    % Reduced Augmentation
    [Vr3, Vz3] = reduced_augmentation([Vr_ddagger3, Vr2, Vr1, Vr0], [Vz_ddagger3, Vz2, Vz1, Vz0], rvals);

    % S-Step
    S3 = sylvester(eye(size(Vr3, 2)) - (nu*dt*((rvals .* Vr3)')*Drr*Vr3), -nu*dt*(Dzz*Vz3)'*Vz3, ((rvals .* Vr3)')*W2*Vz3);
    [Vr3, S3, Vz3, r3] = truncate_svd(Vr3, S3, Vz3, tolerance);

end