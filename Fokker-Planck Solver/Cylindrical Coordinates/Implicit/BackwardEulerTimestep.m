% Facilitates single forward Euler timestep for the 2D Fokker-Planck system

function [Vr, S, Vz, rank] = BackwardEulerTimestep(Vr0, S0, Vz0, dt, tval, rvals, zvals, Br, Bz, Drr, Dzz, tolerance)
    Fr = GetRadialFluxSPCC(tval, rvals, Br, Drr);
    Fz = GetZFluxSPCC(tval, zvals, Bz, Dzz);
    rank = 1;
    Nr = size(Fr);
    Nz = size(Fz);

    Vr0_star = Vr0; Vz0_star = Vz0;
   
    K0 = Vr0*S0;
    L0 = Vz0*(S0');

    K1 = sylvester(eye(Nr) - (dt*Fr), ((-dt)*((Fz*Vz0_star)')*Vz0_star), K0);
    L1 = sylvester(eye(Nz) - (dt*Fz), ((-dt)*((Fr*Vr0_star)')*(rvals.*Vr0_star)), L0);

    [Vr1_ddagger, ~] = qr2(K1, rvals);
    [Vz1_ddagger, ~] = qr(L1, 0);

    % Vr = Vr1_ddagger;
    % Vz = Vz1_ddagger;
    [Vr, Vz] = reduced_augmentation([Vr1_ddagger, Vr0], [Vz1_ddagger, Vz0], rvals);
    S = sylvester(full(speye(size(Vr, 2)) - (dt*((rvals .* Vr)')*(Fr*Vr))), -dt*(Fz*Vz)'*Vz, ((rvals .* Vr)'*Vr0)*S0*((Vz0')*Vz));
    [Vr, S, Vz, rank] = truncate_svd(Vr, S, Vz, tolerance);
end












