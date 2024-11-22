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
function [f] = FokkerPlanckSolver(type, f, dt, Nx, tf, interval, C, D)
    xvals = linspace(interval(1), interval(2), Nx + 1);
    dx = xvals(2) - xvals(1);

    tvals = (0:dt:tf)';
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end
    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);
        switch (type)
            case '1'
                [f] = ForwardEulerTimestep(f, dt, dx, xvals, C(xvals, tvals(n)), D(xvals));
        end

        figure(1); clf; plot(xvals, f);
        legend(sprintf('N_x = %s', num2str(Nx, 3)), 'Location','northwest');
        xlabel('V_x'); ylabel('f(V_x)'); title([sprintf('Forward Euler approximation of 1D Fokker-Planck system at time %s', num2str(tf, 4))]);

    end
end