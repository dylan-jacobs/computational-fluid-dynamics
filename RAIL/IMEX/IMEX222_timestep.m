% IMEX(2, 2, 2) timestep

function [Vx, S, Vy, r1] = IMEX222_timestep(Vx0, S0, Vy0, tn, dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance)

    nu = 1-(sqrt(2)/2);
    delta = 1-(1/(2*nu));
    N = size(Dxx);
    
    % Stage 1: Backward Euler
    [Vx1, S1, Vy1, ~] = IMEX111_timestep(Vx0, S0, Vy0, tn, nu*dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance);

    F1 = -(((Dx*(A{1, 1}(tn).*Vx0))*(A{1, 2}(tn)*S0)*((A{1, 3}(tn).*Vy0)')) + ((A{2, 1}(tn).*Vx0)*(A{2, 2}(tn)*S0)*((Dy*(A{2, 3}(tn).*Vy0))')));
    W1 = (Vx0*S0*(Vy0')) + (dt*(nu*F1 + nu*(phi{1, 1}(tn)*phi{1, 2}(tn+(nu*dt))*(phi{1, 3}(tn)')))) + ((1-nu)*dt*(((Dxx*Vx1*S1*(Vy1')) + (Vx1*S1*((Dyy*Vy1)')))));

    % Reduced Augmentation
    % Predict V_dagger using IMEX111 for second stage
    [Vx_dagger, ~, Vy1_dagger, ~] = IMEX111_timestep(Vx0, S0, Vy0, tn, dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance);
    [Vx_star, Vy_star] = reduced_augmentation([Vx_dagger, Vx1, Vx0], [Vy1_dagger, Vy1, Vy0]);

    % Stage 2: KLS Steps
    % K/L-Step
    K1 = sylvester(eye(N) - (nu*dt*Dxx), -nu*dt*(Dyy*Vy_star)'*Vy_star, W1*Vy_star);
    L1 = sylvester(eye(N) - (nu*dt*Dyy), -nu*dt*(Dxx*Vx_star)'*Vx_star, (W1')*Vx_star);

    % Get bases
    [Vx_ddagger, ~] = qr(K1, 0); [Vy_ddagger, ~] = qr(L1, 0);

    % Reduced Augmentation
    [Vx, Vy] = reduced_augmentation([Vx_ddagger, Vx1, Vx0], [Vy_ddagger, Vx1, Vy0]);

    % S-Step
    S = sylvester(eye(size(Vx, 2)) - (nu*dt*(Vx')*Dxx*Vx), -nu*dt*(Dyy*Vy)'*Vy, (Vx')*W1*Vy);
    [Vx, S, Vy, r1] = truncate_svd(Vx, S, Vy, tolerance);
end