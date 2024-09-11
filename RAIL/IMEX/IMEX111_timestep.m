% IMEX(1, 1, 1) timestep


function [Vx, S, Vy, r1] = IMEX111_timestep(Vx0, S0, Vy0, tn, dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance)
    
    N = size(Dxx);

    Vx_star = Vx0; Vy_star = Vy0;

    F0 = -(((Dx*(A{1, 1}(tn).*Vx_star))*(A{1, 2}(tn)*S0)*((A{1, 3}(tn).*Vy_star)')) + ((A{2, 1}(tn).*Vx_star)*(A{2, 2}(tn)*S0)*((Dy*(A{2, 3}(tn).*Vy_star))')));
    W0 = (Vx0*S0*(Vy0')) + (dt*(F0 + (phi{1, 1}(tn)*phi{1, 2}(tn+dt)*(phi{1, 3}(tn)'))));

    K1 = sylvester(eye(N) - (dt*Dxx), -dt*(Dyy*Vy_star)'*Vy_star, W0*Vy_star);
    L1 = sylvester(eye(N) - (dt*Dyy), -dt*(Dxx*Vx_star)'*Vx_star, (W0')*Vx_star);

    [Vx_ddagger, ~] = qr(K1, 0); [Vy_ddagger, ~] = qr(L1, 0);

    [Vx, Vy] = reduced_augmentation([Vx_ddagger, Vx0], [Vy_ddagger, Vy0]);
    S = sylvester(eye(size(Vx, 2)) - (dt*(Vx')*Dxx*Vx), -dt*(Dyy*Vy)'*Vy, (Vx')*W0*Vy);
    [Vx, S, Vy, r1] = truncate_svd(Vx, S, Vy, tolerance);
end