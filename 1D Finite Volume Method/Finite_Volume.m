% Use WENO-5, Lax-Friedrichs Flux and various temporal discretizations (SSP
% RK1/RK2/RK3) to solve conservation law equations
% Spatial mesh size = NX+1, Nx cells
% tf = final time
% Flux function f
% u_t + (f(u))_x = 0


function [u] = Finite_Volume(discretizationType, Nx, lambda, interval, tf, f, u0, alpha)
    xvals = linspace(interval(1), interval(2), Nx + 1)';
    dx = xvals(2) - xvals(1);
    dt = lambda*dx;
    tvals = (0:dt:tf)'; 
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end

    % initialize u to cell averages
    u = zeros(numel(Nx), 1);
    for i = 1:Nx
        u(i, 1) = (1/dx)*GLQuad(u0, xvals(i), xvals(i+1), 7);
    end

    u = Time_Discretization(discretizationType, u, tvals, xvals, alpha, f);

end





