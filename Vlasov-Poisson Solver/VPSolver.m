% Use WENO-5, Lax-Friedrichs Flux and various temporal discretizations (SSP
% RK1/RK2/RK3) to solve Vlasov-Poisson Equation
% Spatial mesh size = N+1, N cells in both X and Y directions
% Inputs:
%   discretizationType - temporal discretization method (RK1, RK2 ... RK4)
%   Nx, Nv - x, v mesh-size, respectively
%   lambda - dt = lambda*dx
%   interval - [xmin, xmax, vmin, vmax]
%   tf - final time
%   f0 - initial condition
% Outputs:
%   f_matrix - solution f at all times
%   EF - electric field at all times
%   mass
%   L1, L2 - norms of numerical solution f at all times
%   energy, entropy
%   tvals


function [f_matrix, EF, mass, L1, L2, energy, entropy, tvals] = VPSolver(discretizationType, Nx, Nv, lambda, interval, tf, f0)

    [X, V, dx, dv] = GetXY(Nx, Nv, interval);
    dt = lambda.*dx; % CFL condition

    tvals = (0:dt:tf)'; 
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end

    % initialize vectors
    Nt = numel(tvals);
    f_matrix = zeros(Nx, Nv, Nt);
    EF = zeros(Nt, 1);
    mass = zeros(Nt, 1);
    L1 = zeros(Nt, 1);
    L2 = zeros(Nt, 1);
    energy = zeros(Nt, 1);
    entropy = zeros(Nt, 1);

    f = f0(X, V);
    f_matrix(:, :, 1) = f;
    [EF(1), mass(1), L1(1), L2(1), energy(1), entropy(1)] = quantities(f, V, interval(2) - interval(1), Nx, dx, dv);
    for n = 2:numel(tvals)
        t0 = tvals(n-1)
        dt = tvals(n) - tvals(n-1); % DO NOT FORGET TO INCLUDE THIS!!!!!
        f = Time_Discretization_Nonsplitting(discretizationType, f, t0, dt, X, V);

        % truncate???
        f = LoMaC(f, V(1, :), dv, 1e-4);

        f_matrix(:, :, n) = f;

        [EF(n), mass(n), L1(n), L2(n), energy(n), entropy(n)] = quantities(f, V, interval(2) - interval(1), Nx, dx, dv);
        % if mod(n, 10) == 0
        %     figure(8); clf; surf(X, V, f_matrix(:, :, n));
        %     colorbar;
        %     shading flat; % removes gridlines
        %     legend(sprintf('N_x = %s, N_v = %s', num2str(Nx, 3), num2str(Nv, 3)), 'Location','northwest');
        %     xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('2D WENO+%s', discretizationType), sprintf(' approximation at time %s', num2str(tvals(n), 4))]);
        %     view(2); % bird's eye view
        %     xlim([interval(1), interval(2)]);
        % end
    end
end





