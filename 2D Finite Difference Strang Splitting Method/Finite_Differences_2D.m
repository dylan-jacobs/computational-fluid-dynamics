% Use WENO-5, Lax-Friedrichs Flux and various temporal discretizations (SSP
% RK1/RK2/RK3) to solve 2D conservation law equations
% Spatial mesh size = N+1, N cells in both X and Y directions
% tf = final time
% Flux function f in X
% Flux function g in Y
% u_t + f(u)_x + g(u)_y = 0


function [u] = Finite_Differences_2D(discretizationType, Nx, Ny, lambda, interval, tf, f, g, u0, alpha, beta)
    [X, Y, dx, dy] = GetXY(Nx, Ny, interval);
    dt = lambda/((alpha/dx) + (beta/dy)); % CFL condition

    tvals = (0:dt:tf)'; 
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end

    u = u0(X, Y);
    for n = 2:numel(tvals)
        t0 = tvals(n-1);
        dt = tvals(n) - tvals(n-1);
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
    % shading flat; % removes gridlines
    % legend(sprintf('Nx = Ny = %s', num2str(Nx, 3)), 'Location','northwest');
    % xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('2D WENO+%s', discretizationType), sprintf(' approximation at time %s', num2str(tf, 4))]);

end





