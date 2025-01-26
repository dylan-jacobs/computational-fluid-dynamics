function [Vr, S, Vz, rank] = BackwardEuler(Vr0, S0, Vz0, rvals, dt, Drr, Dzz, tolerance)

    Nr = size(Drr);
    Nz = size(Dzz);

    Vr0_star = Vr0; Vz0_star = Vz0;
   
    K0 = Vr0*S0;
    L0 = Vz0*(S0');

    K1 = sylvester(eye(Nr) - (dt*Drr), ((-dt)*((Dzz*Vz0_star)')*Vz0_star), K0);
    L1 = sylvester(eye(Nz) - (dt*Dzz), ((-dt)*((Drr*Vr0_star)')*(rvals.*Vr0_star)), L0);

    [Vr1_ddagger, ~] = qr2(K1, rvals);
    [Vz1_ddagger, ~] = qr(L1, 0);

    [Vr, Vz] = reduced_augmentation([Vr1_ddagger, Vr0], [Vz1_ddagger, Vz0], rvals);
    S = sylvester(eye(size(Vr, 2)) - (dt*((rvals .* Vr)')*(Drr*Vr)), -dt*(Dzz*Vz)'*Vz, ((rvals .* Vr)'*Vr0)*S0*((Vz0')*Vz));
    [Vr, S, Vz, rank] = truncate_svd(Vr, S, Vz, tolerance);
end