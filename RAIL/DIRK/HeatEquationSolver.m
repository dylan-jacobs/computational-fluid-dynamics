function [U, ranks] = HeatEquationSolver(type, U, dt, Nx, Ny, tf, interval, d1, d2, tolerance)
    dx = (interval(2) - interval(1)) / (Nx);
    dy = (interval(4) - interval(3)) / (Ny);

    tvals = (0:dt:tf)';
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end

    [Vx, S, Vy] = svd(U);
    ranks = zeros(numel(tvals), 2);
    [~, ~, Dxx, Dyy] = spectral_matrix(interval(2) - interval(1), Nx, Ny, dx, dy, d1, d2);

    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);
        switch (type)
            case '1'
                [Vx, S, Vy, ranks(n, 2)] = RAIL_B_Euler_timestep(Vx, S, Vy, dt, Dxx, Dyy, tolerance);
            case '2'
                [Vx, S, Vy, ranks(n, 2)] = RAIL_DIRK2_timestep(Vx, S, Vy, dt, Dxx, Dyy, tolerance);
            case '3'
                [Vx, S, Vy, ranks(n, 2)] = RAIL_DIRK3_timestep(Vx, S, Vy, dt, Dxx, Dyy, tolerance);
        end
    end
    U = Vx*S*(Vy');
    ranks(:, 1) = tvals;
end