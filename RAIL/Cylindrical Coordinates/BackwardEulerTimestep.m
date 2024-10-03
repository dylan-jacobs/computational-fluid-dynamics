function [Vr, S, Vz, rank] = BackwardEulerTimestep(Vr0, S0, Vz0, rvals, dt, Drr, Dzz, tolerance)

    N = size(Drr);

    Vr0_star = Vr0; Vz0_star = Vz0;
   
    K0 = Vr0*S0;
    L0 = Vz0*(S0')*(Vr0')*(rvals .* Vr0_star);

    K1 = sylvester(eye(N) - (dt*Drr), ((-dt)*((Dzz*Vz0_star)')*Vz0_star), K0);
    L1 = sylvester(eye(N) - (dt*Dzz), ((-dt)*((Drr*Vr0_star)')*(rvals.*Vr0_star)), L0);

    Vr1_ddagger = qr2(K1, rvals);
    [Vz1_ddagger, ~] = qr(L1, 0);

    [Vr, Vz] = reduced_augmentation([Vr1_ddagger, Vr0], [Vz1_ddagger, Vz0]);
    S = sylvester(eye(size(Vr, 2)) - (dt*(Vr')*Drr*(Vr)), -dt*(Dzz*Vz)'*Vz, (Vr)'*Vr0*S0*(Vz0')*Vz);
    [Vr, S, Vz, rank] = truncate_svd(Vr, S, Vz, tolerance);
end



function [V] = qr2(V0, rvals)
    [V, ~] = qr((rvals.^0.5) .* V0, 0);
end