% Use WENO-5, Lax-Friedrichs Flux and various temporal discretizations (SSP
% RK1/RK2/RK3) to solve 2D conservation law equations
% Spatial mesh size = N+1, N cells in both X and Y directions
% tf = final time
% Flux function f in X
% Flux function g in Y
% u_t + f(u)_x + g(u)_y = 0


function [u] = Finite_Differences_2D_Nonsplitting(discretizationType, Nx, Ny, lambda, interval, tf, f, g, u0, alpha, beta, plot)
    [X, Y, dx, dy] = GetXY(Nx, Ny, interval);
    dt = lambda/((alpha/dx) + (beta/dy)); % CFL condition

    tvals = (0:dt:tf)'; 
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end

    u = u0(X, Y);
    figure; surf(X, Y, u);
    colorbar; shading flat;
    for n = 2:numel(tvals)
        t0 = tvals(n-1);
        dt = tvals(n) - tvals(n-1); % DO NOT FORGET TO INCLUDE THIS!!!!!
        u = Time_Discretization_Nonsplitting(discretizationType, u, t0, dt, X, Y, alpha, beta, f, g);
    end

    if plot
        figure; clf; surf(X, Y, u);
        colorbar;
        shading flat; % removes gridlines
        legend(sprintf('Nx = Ny = %s', num2str(Nx, 3)), 'Location','northwest');
        xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('2D WENO+%s', discretizationType), sprintf(' approximation at time %s', num2str(tf, 4))]);
    end
end





