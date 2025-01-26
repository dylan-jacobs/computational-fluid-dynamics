% Facilitates single forward Euler timestep for the 1D Fokker-Planck system

function [Vr2, S2, Vz2, rank] = RK2Timestep(Vr0, S0, Vz0, dt, tval, rvals, zvals, Br, Bz, Drr, Dzz, tolerance)
    
    [Vr1, S1, Vz1, ~] = ForwardEulerTimestep(Vr0, S0, Vz0, dt, tval, rvals, zvals, Br, Bz, Drr, Dzz, tolerance);
    % ----- First Stage -----

    % Fr0 = GetRadialFluxSPCC(tval, rvals, Br, Drr);
    % Fz0 = GetZFluxSPCC(tval, zvals, Bz, Dzz);
    % 
    % Vr0_star = Vr0; Vz0_star = Vz0;
    % 
    % K0 = Vr0*S0;
    % L0 = Vz0*(S0');
    % 
    % K1 = K0 + (dt*((Fr0*K0) + (K0*(Fz0*Vz0_star)'*Vz0_star)));
    % L1 = L0 + (dt*((Fz0*L0) + (L0*(Fr0*Vr0_star)'*(rvals.*Vr0_star))));
    % 
    % [Vr1_ddagger, ~] = qr2(K1, rvals);
    % [Vz1_ddagger, ~] = qr(L1, 0);
    % 
    % [Vr1, Vz1] = reduced_augmentation([Vr1_ddagger, Vr0], [Vz1_ddagger, Vz0], rvals);
    % S01 = ((rvals .* Vr1)'*Vr0)*S0*((Vz0')*Vz1);
    % S1 = S01 + dt*( (rvals.*Vr1)'*(Fr0*Vr1)*S01 + (S01*(Fz0*Vz1)'*Vz1) );

    % ------- Second Stage -------
    Fr0 = GetRadialFluxSPCC(tval, rvals, Br, Drr);
    Fz0 = GetZFluxSPCC(tval, zvals, Bz, Dzz);
    Fr1 = GetRadialFluxSPCC(tval + dt, rvals, Br, Drr);
    Fz1 = GetZFluxSPCC(tval + dt, zvals, Bz, Dzz);

    [Vr1_star, Vz1_star] = reduced_augmentation([Vr1, Vr0], [Vz1, Vz0], rvals);
    K01 = Vr0*S0*(Vz0')*Vz1_star;
    L01 = (Vz0*S0')*(Vr0'*(rvals .* Vr1_star));
    K1 = Vr1*S1*(Vz1')*Vz1_star;
    L1 = Vz1*(S1')*(Vr1')*(rvals .* Vr1_star);

    K2 = K01 + (0.5*dt*((Fr0*K01) + (K01*(Fz0*Vz1_star)'*Vz1_star) + (Fr1*K1) + (K1*(Fz1*Vz1_star)'*Vz1_star)));
    L2 = L01 + (0.5*dt*((Fz0*L01) + (L01*(Fr0*Vr1_star)'*(rvals.*Vr1_star) + (Fz1*L1) + (L1*(Fr1*Vr1_star)'*(rvals.*Vr1_star)))));

    [Vr2_ddagger, ~] = qr2(K2, rvals);
    [Vz2_ddagger, ~] = qr(L2, 0);

    [Vr2, Vz2] = reduced_augmentation([Vr2_ddagger, Vr1, Vr0], [Vz2_ddagger, Vz1, Vz0], rvals);
    S02 = ((rvals .* Vr2)'*Vr0)*S0*((Vz0')*Vz2);
    S1 = ((rvals .* Vr2)'*Vr1)*S1*((Vz1')*Vz2);
    S2 = S02 + 0.5*dt*((rvals.*Vr2)'*(Fr0*Vr2)*S02 + (S02*(Fz0*Vz2)'*Vz2) + (rvals.*Vr2)'*(Fr1*Vr2)*S1 + (S1*(Fz1*Vz2)'*Vz2) );
    [Vr2, S2, Vz2, rank] = truncate_svd(Vr2, S2, Vz2, tolerance);



end












