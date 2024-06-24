% Use WENO-5, Lax-Friedrichs Flux and various temporal discretizations (SSP
% RK1/RK2/RK3) to solve 2D conservation law equations
% Spatial mesh size = N+1, N cells in both X and Y directions
% tf = final time
% Flux function f in X
% Flux function g in Y
% u_t + f(u)_x + g(u)_y = 0


function [u] = Finite_Differences_2D(discretizationType, Nx, Ny, lambda, interval, tf, f, g, u0, alpha, beta)
    xvals = linspace(interval(1), interval(2), Nx+1)';
    yvals = linspace(interval(3), interval(4), Ny+1)';
    dx = xvals(2) - xvals(1);
    dy = yvals(2) - yvals(1);
    xmid = xvals(1:end-1) + dx/2;
    ymid = yvals(1:end-1) + dy/2;
    dt = lambda/((alpha/dx) + (beta/dy));

    [X, Y] = meshgrid(xmid, ymid); X = X'; Y = Y';

    tvals = (0:dt:tf)'; 
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end

    u = u0(X, Y);
    for n = 2:numel(tvals)
        t0 = tvals(n-1);
        for i = 1:Ny
            % update x
            u(:, i) = Time_Discretization(discretizationType, u(:, i), t0, dt/2, X(:, i), Y(:, i), alpha, f);
        end
        for j = 1:Nx
            % update y
            u(j, :) = Time_Discretization(discretizationType, u(j, :)', t0, dt, X(j, :)', Y(j, :)', beta, g);
        end
        for i = 1:Ny
            % update x
            u(:, i) = Time_Discretization(discretizationType, u(:, i), t0 + (dt/2), dt/2, X(:, i), Y(:, i), alpha, f);
        end
    end

    % figure; clf; surf(X, Y, u);
    % colorbar;
    % shading interp; % removes gridlines
    % legend(sprintf('Nx = Ny = %s', num2str(Nx, 3)));
    % xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('2D WENO+%s', discretizationType), sprintf(' approximation at time %s', num2str(tf, 2))]);

end





