% Completes 1 DIRK3 timestep, updating the bases Vx0, S0, Vy0 --> Vx3, S3, Vy3
% Inputs: bases: Vx0, S0, Vy0
%         timestep size: dt
%         differentiation matrices: Dxx, Dyy
%         final truncation tolerance: tolerance
% Outputs: updated vases Vx3, S3, Vy3
%          rank: r3

function [Vx3, S3, Vy3, r3] = RAIL_DIRK3_timestep(Vx0, S0, Vy0, dt, Dxx, Dyy, tolerance)
    
    % RK butcher table values
    nu = 0.435866521508459;
    beta1 = -(3/2)*(nu^2) + (4*nu) - (1/4);
    beta2 = (3/2)*(nu^2) - (5*nu) + (5/4);
    N = size(Dxx);
    
    % Stage 1: Backward Euler
    [Vx1, S1, Vy1, ~] = RAIL_B_Euler_timestep(Vx0, S0, Vy0, nu*dt, Dxx, Dyy, tolerance);
    [Vx_dagger1, ~, Vy_dagger1, ~] = RAIL_B_Euler_timestep(Vx0, S0, Vy0, ((1+nu)/2)*dt, Dxx, Dyy, tolerance);

    Y1 = (((Dxx*Vx1*S1*(Vy1')) + (Vx1*S1*((Dyy*Vy1)'))));
    W1 = (Vx0*S0*(Vy0')) + (((1-nu)/2)*dt*Y1);

    % Reduced Augmentation
    [Vx_star1, Vy_star1] = reduced_augmentation([Vx_dagger1, Vx1, Vx0], [Vy_dagger1, Vy1, Vy0]);

    % Stage 2:
    % K/L-Step
    K2 = sylvester(eye(N) - (nu*dt*Dxx), -nu*dt*(Dyy*Vy_star1)'*Vy_star1, W1*Vy_star1);
    L2 = sylvester(eye(N) - (nu*dt*Dyy), -nu*dt*(Dxx*Vx_star1)'*Vx_star1, (W1')*Vx_star1);

    % Get bases
    [Vx_ddagger2, ~] = qr(K2, 0); [Vy_ddagger2, ~] = qr(L2, 0);

    % Reduced Augmentation
    [Vx2, Vy2] = reduced_augmentation([Vx_ddagger2, Vx1, Vx0], [Vy_ddagger2, Vx1, Vy0]);

    % S-Step
    S2 = sylvester(eye(size(Vx2, 2)) - (nu*dt*(Vx2')*Dxx*Vx2), -nu*dt*(Dyy*Vy2)'*Vy2, (Vx2')*W1*Vy2);
    [Vx2, S2, Vy2, ~] = truncate_svd(Vx2, S2, Vy2, tolerance);

    % Stage 3:
    % Predict V_dagger using B. Euler
    [Vx_dagger3, ~, Vy_dagger3, ~] = RAIL_B_Euler_timestep(Vx2, S2, Vy2, dt, Dxx, Dyy, tolerance);
    Y2 = (((Dxx*Vx2*S2*(Vy2')) + (Vx2*S2*((Dyy*Vy2)'))));
    W2 = (Vx0*S0*(Vy0')) + (beta1*dt*Y1) + (beta2*dt*Y2);
      
    % Reduced augmentation
    [Vx_star3, Vy_star3] = reduced_augmentation([Vx_dagger3, Vx2, Vx1, Vx0], [Vy_dagger3, Vy2, Vy1, Vy0]);

    % K/L-Step
    K3 = sylvester(eye(N) - (nu*dt*Dxx), -nu*dt*(Dyy*Vy_star3)'*Vy_star3, W2*Vy_star3);
    L3 = sylvester(eye(N) - (nu*dt*Dyy), -nu*dt*(Dxx*Vx_star3)'*Vx_star3, (W2')*Vx_star3);

    % Get bases
    [Vx_ddagger3, ~] = qr(K3, 0); [Vy_ddagger3, ~] = qr(L3, 0);

    % Reduced Augmentation
    [Vx3, Vy3] = reduced_augmentation([Vx_ddagger3, Vx2, Vx1, Vx0], [Vy_ddagger3, Vx2, Vx1, Vy0]);

    % S-Step
    S3 = sylvester(eye(size(Vx3, 2)) - (nu*dt*(Vx3')*Dxx*Vx3), -nu*dt*(Dyy*Vy3)'*Vy3, (Vx3')*W2*Vy3);
    [Vx3, S3, Vy3, r3] = truncate_svd(Vx3, S3, Vy3, tolerance);

end