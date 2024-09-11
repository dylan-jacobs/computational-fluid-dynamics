% Completes 1 B. Euler timestep, updating the bases Vx0, S0, Vy0 --> Vx, S, Vy
% Inputs: bases: Vx0, S0, Vy0
%         timestep size: dt
%         differentiation matrices: Dxx, Dyy
%         final truncation tolerance: tolerance
% Outputs: updated vases Vx, S, Vy
%          rank: r1

function [Vx, S, Vy, r1] = RAIL_B_Euler_timestep(Vx0, S0, Vy0, dt, Dxx, Dyy, tolerance)
    N = size(Dxx);

    Vx_star = Vx0; Vy_star = Vy0;
    K = Vx0*S0;
    L = Vy0*(S0');

    % (I - dtDxx)Kn+1 + Kn+1(-dt(DyyVy_star)'Vy_star = Kn
    K1 = sylvester(eye(N) - (dt*Dxx), -dt*(Dyy*Vy_star)'*Vy_star, K);
    L1 = sylvester(eye(N) - (dt*Dyy), -dt*(Dxx*Vx_star)'*Vx_star, L);

    [Vx_ddagger, ~] = qr(K1, 0); [Vy_ddagger, ~] = qr(L1, 0);

    [Vx, Vy] = reduced_augmentation([Vx_ddagger, Vx0], [Vy_ddagger, Vy0]);
    S = sylvester(eye(size(Vx, 2)) - (dt*(Vx')*Dxx*Vx), -dt*(Dyy*Vy)'*Vy, (Vx')*Vx0*S0*(Vy0')*Vy);
    [Vx, S, Vy, r1] = truncate_svd(Vx, S, Vy, tolerance);
end

















