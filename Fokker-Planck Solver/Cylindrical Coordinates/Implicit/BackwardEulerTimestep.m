% Facilitates single forward Euler timestep for the 2D Fokker-Planck system

function [Vr, S, Vz, rank] = BackwardEulerTimestep(Vr0, S0, Vz0, dt, tval, R, Z, rvals, zvals, Br, Bz, Drr, Dzz, tolerance, f_inf)
    rank = 1;
    Fr = GetRadialFluxSPCC(tval, rvals, Br, Drr);
    Fz = GetZFluxSPCC(tval, zvals, Bz, Dzz);
    Ar = @(w) w;
    Az = @(w) w.^0;
    Br = @(w, u, t) w.*(w - u);
    Bz = @(w, u, t) (w - u);
    Cr = @(w, t) 3.*w;
    Cz = @(w, t) 3;
    Fr = getF(tval, 0, rvals, rvals(2)-rvals(1), Ar, Br, Cr, size(rvals, 1));
    Fz = getF(tval, 0, zvals, zvals(2) - zvals(1), Az, Bz, Cz, size(rvals, 2));

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
    % [Vr, S, Vz, rank] = LoMaC(Vr, S, Vz, R, Z, rvals, zvals, tolerance, f_inf);
end



function F = getF(t,u,wvals,dw,A,B,C,N)
    % u = u(t) bulk velocity at time t.
    % wvals = center-centered grid
    
    
    Avals = A(wvals);
    Bvals = B(wvals(1:N-1)+dw/2,u,t); %B_{j+1/2} for j=1,2,...,N-1
    Cvals = C(wvals(1:N-1)+dw/2,t); %C_{j+1/2} for j=1,2,...,N-1
    pvals = dw*Bvals./Cvals + 1.0e-14;
    deltavals = 1./pvals - 1./(exp(pvals)-1);
    
    
    V1 = -((1/dw)*Cvals - deltavals.*Bvals);
    V2 = (1-deltavals).*Bvals + (1/dw)*Cvals;
    
    
    F_plus = spdiags([0;V2],1,N,N) + spdiags([V1;0],0,N,N);
    F_neg = spdiags(V1,-1,N,N) + spdiags([0;V2],0,N,N);
    
    
    F = diag(1./Avals)*(1/dw)*(F_plus-F_neg);
end









