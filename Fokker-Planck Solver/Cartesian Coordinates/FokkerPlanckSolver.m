% Solves the 1D Fokker-Planck system using specified time discretization
% method


% type = time-discretization to use
% U    = initial condition
% dt   = timestep
% Nx   = number of points in x
% tf   = final time 
% interval = spatial interval
% C    = advection function
% D    = diffusive function
function [f, data] = FokkerPlanckSolver(type, f, dt, Nx, tf, interval, B, D, plotTimeSteps, f_exact)
    [xvals, dx] = GetXVals(Nx, interval);

    tvals = (0:dt:tf)';
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end

    l1 = zeros(numel(tvals), 1); l1(1) = dx*sum(abs(f));
    positivity = zeros(numel(tvals), 1); positivity(1) = min(f);
    relative_entropy = zeros(numel(tvals), 1); relative_entropy(1) = dx*sum(f.*(log(f./f_exact)));
    mass = zeros(numel(tvals), 1); mass(1) = dx*sum(f);

    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);
        switch (type)
            case '1'
                [f] = ForwardEulerTimestep(f, dt, dx, tvals(n), xvals, B, D);
        end

        l1(n) = dx*sum(abs(f-f_exact));
        positivity(n) = min(f);
        relative_entropy(n) = dx*sum(f.*(log(f./f_exact)));
        mass(n) = dx*sum(f);

        if plotTimeSteps
            figure(1); clf; plot(xvals, f);
            legend(sprintf('N_x = %s', num2str(Nx, 3)), 'Location','northwest');
            xlabel('V_x'); ylabel('f(V_x)'); title([sprintf('Forward Euler approximation of 1D Fokker-Planck system at time %s', num2str(tf, 4))]);
        end
    end
    data = [l1, positivity, relative_entropy, mass, tvals];
end