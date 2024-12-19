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
function [f, data] = FokkerPlanckSolver(type, f, dt, Nr, Nz, tf, interval, B, D, plotTimeSteps, f_exact, tolerance)
    [R, Z, dr, dz] = GetRZ(Nr, Nz, interval);

    rvals = R(:, 1);
    zvals = Z(1, :)';
    tvals = (0:dt:tf)';
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end
    l1 = zeros(numel(tvals), 1); l1(1) = dr*dz*sum(sum(abs(f)));
    positivity = zeros(numel(tvals), 1); positivity(1) = min(min(f));
    relative_entropy = zeros(numel(tvals), 1); relative_entropy(1) = dr*dz*sum(sum(f.*(log(f./f_exact))));
    mass = zeros(numel(tvals), 1); mass(1) = dr*dz*sum(sum(f));

    [Vr0, S0, Vz0] = svd2(f, rvals);

    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);
        switch (type)
            case '1'
                [f] = ForwardEulerTimestep(Vr0, S0, Vz0, dt, tvals(n), rvals, zvals, B, B, D, D, tolerance);
        end
        size(f)
        size(f_exact)
        dr*dz*sum(sum(abs(f-f_exact)))
        l1(n) = dr*dz*sum(sum(abs(f-f_exact)));
        positivity(n) = min(min(f));
        relative_entropy(n) = dr*dz*sum(sum(f.*(log(f./f_exact))));
        mass(n) = dr*dz*sum(sum(f));

        if plotTimeSteps
            figure(1); clf; surf(rvals, zvals, f);
            legend(sprintf('N_r = %s', num2str(Nr, 3)), 'Location','northwest');
            xlabel('V_r'); ylabel('f(V_r, V_z)'); title([sprintf('Forward Euler approximation of 1D Fokker-Planck system at time %s', num2str(tf, 4))]);
        end
        return
    end
    data = [l1, positivity, relative_entropy, mass, tvals];
end