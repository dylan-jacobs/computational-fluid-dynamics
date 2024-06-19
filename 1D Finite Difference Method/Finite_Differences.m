% Use WENO-5, Lax-Friedrichs Flux and various temporal discretizations (SSP
% RK1/RK2/RK3) to solve conservation law equations
% Spatial mesh size = NX+1, Nx cells
% tf = final time
% Flux function f
% u_t + (f(u))_x = 0


function [u] = Finite_Differences(discretizationType, Nx, lambda, interval, tf, f, u0, alpha)
    xvals = linspace(interval(1), interval(2), Nx + 1)';
    dx = xvals(2) - xvals(1);
    xmid = xvals(1:end-1) + dx/2;
    dt = lambda*dx;
    tvals = (0:dt:tf)'; 
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end

    u = u0(xmid);
    u = Time_Discretization(discretizationType, u, tvals, xvals, alpha, f);

end





