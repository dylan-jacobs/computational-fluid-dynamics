% Completes 1 DIRK2 timestep, updating the bases Vx0, S0, Vy0 --> Vx, S, Vy
% Inputs: bases: Vx0, S0, Vy0
%         timestep size: dt
%         differentiation matrices: Dxx, Dyy
%         final truncation tolerance: tolerance
% Outputs: updated vases Vx, S, Vy
%          rank: r1


function [Vx, S, Vy, r1] = RAIL_DIRK2_timestep(Vx0, S0, Vy0, dt, Dxx, Dyy, tolerance)
    
    nu = 1-(sqrt(2)/2);
    N = size(Dxx);
    
    % Stage 1: Backward Euler
    [Vx1, S1, Vy1, ~] = RAIL_B_Euler_timestep(Vx0, S0, Vy0, nu*dt, Dxx, Dyy, tolerance);

    W0 = (Vx0*S0*(Vy0')) + ((1-nu)*dt*(((Dxx*Vx1*S1*(Vy1')) + (Vx1*S1*((Dyy*Vy1)')))));

    % Reduced Augmentation
    % Predict V_dagger using B. Euler for second stage
    [Vx_dagger, ~, Vy1_dagger, ~] = RAIL_B_Euler_timestep(Vx0, S0, Vy0, dt, Dxx, Dyy, tolerance);
    [Vx_star, Vy_star] = reduced_augmentation([Vx_dagger, Vx1, Vx0], [Vy1_dagger, Vy1, Vy0]);

    % Stage 2: KLS Steps
    % K/L-Step
    K1 = sylvester(eye(N) - (nu*dt*Dxx), -nu*dt*(Dyy*Vy_star)'*Vy_star, W0*Vy_star);
    L1 = sylvester(eye(N) - (nu*dt*Dyy), -nu*dt*(Dxx*Vx_star)'*Vx_star, (W0')*Vx_star);

    % Get bases
    [Vx_ddagger, ~] = qr(K1, 0); [Vy_ddagger, ~] = qr(L1, 0);

    % Reduced Augmentation
    [Vx, Vy] = reduced_augmentation([Vx_ddagger, Vx1, Vx0], [Vy_ddagger, Vx1, Vy0]);

    % S-Step
    S = sylvester(eye(size(Vx, 2)) - (nu*dt*(Vx')*Dxx*Vx), -nu*dt*(Dyy*Vy)'*Vy, (Vx')*W0*Vy);
    [Vx, S, Vy, r1] = truncate_svd(Vx, S, Vy, tolerance);
end