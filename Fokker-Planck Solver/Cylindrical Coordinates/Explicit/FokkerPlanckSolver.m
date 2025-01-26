% Solves the 1D Fokker-Planck system using specified time discretization
% method


% type = time-discretization to use
% U    = initial condition
% dt   = timestep
% Nr   = number of points in x
% tf   = final time 
% interval = spatial interval
% C    = advection function
% D    = diffusive function
function [f, data] = FokkerPlanckSolver(type, f, dt, Nr, Nz, tf, interval, B, D, plotTimeSteps, f_inf, tolerance)
    [R, Z, dr, dz] = GetRZ(Nr, Nz, interval);

    rvals = R(:, 1);
    zvals = Z(1, :)';
    tvals = (0:dt:tf)';
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end
    l1 = zeros(numel(tvals), 1); l1(1) = dr*dz*sum(sum(abs(rvals .* (f-f_inf))));
    positivity = zeros(numel(tvals), 1); positivity(1) = min(min(f));
    relative_entropy = zeros(numel(tvals), 1); relative_entropy(1) = dr*dz*sum(sum(rvals .* (f.*(log(f./f_inf)))));
    mass = zeros(numel(tvals), 1); mass(1) = dr*dz*sum(sum(rvals .* f));
    ranks = zeros(numel(tvals), 1);

    [Vr, S, Vz] = svd2(f, rvals);
    r0 = min(10, size(Vr, 2));

    % truncate solution to pre-set initial rank r0
    Vr = Vr(:, 1:r0);
    S = S(1:r0, 1:r0);
    Vz = Vz(:, 1:r0);
    ranks(1) = r0;

    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);
        switch (type)
            case '1'
                [Vr, S, Vz, ranks(n)] = ForwardEulerTimestep(Vr, S, Vz, dt, tvals(n), rvals, zvals, B, B, D, D, tolerance);
            case '2'
                [Vr, S, Vz, ranks(n)] = RK2Timestep(Vr, S, Vz, dt, tvals(n), rvals, zvals, B, B, D, D, tolerance);

        end
        f = Vr*S*(Vz');
        l1(n) = dr*dz*sum(sum(abs(rvals .* (f-f_inf))));
        positivity(n) = min(min(f));
        relative_entropy(n) = dr*dz*sum(sum(rvals .* (f.*(log(f./f_inf)))));
        mass(n) = dr*dz*sum(sum(rvals .* f));

        if plotTimeSteps
            figure(1); clf; surf(zvals, rvals, f);
            legend(sprintf('N_r = %s', num2str(Nr, 3)), 'Location','northwest');
            colorbar;
            shading interp;
            xlabel('V_z'); ylabel('V_r'); zlabel('f(V_x, V_z)'); title([sprintf('Forward Euler approximation of 1D Fokker-Planck system at time %s', num2str(tf, 4))]);
        end
        
    end
    data = [l1, positivity, relative_entropy, mass, tvals, ranks];
end