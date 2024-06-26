% Choose a type of time discretization: Forward Euler (SSP-RK1), SSP-RK2,
% SSP-RK3
% f = flux function

function [u] = Time_Discretization_Nonsplitting(type, u, t0, dt, X, Y, alpha, beta, f, g)
    dx = X(2, 1) - X(1, 1);
    dy = Y(1, 2) - Y(1, 1);

    [F0, G0] = GetFlux(X, Y, u, t0, alpha, beta, f, g); % F is horizontal flux, G is vertical

    switch type
        case 'RK1'
            u = u + (dt)*((F0/dx) + (G0/dy));
        case 'RK2'
            % Stage 1
            u1 = u + (dt)*((F0/dx) + (G0/dy)); % at time t1 = tn+dt
    
            % Stage 2
            t1 = t0 + dt;
            [F1, G1] = GetFlux(X, Y, u1, t1, alpha, beta, f, g);

            u = u + ((dt/(2*dx))*(F0 + F1)) + ((dt/(2*dy))*(G0 + G1));
        case 'RK3'
            % Stage 1
            u1 = u + (dt)*((F0/dx) + (G0/dy));
    
            % Stage 2
            t1 = t0 + dt;
            [F1, G1] = GetFlux(X, Y, u1, t1, alpha, beta, f, g);

            u2 = u + ((dt/(4*dx))*(F0 + F1)) + ((dt/(4*dy))*(G0 + G1)); % at time t2 = tn + dt

            % Stage 3
            t2 = t0 + dt/2;
            [F2, G2] = GetFlux(X, Y, u2, t2, alpha, beta, f, g);

            u = u + ((dt/(dx))*((F0/6) + (F1/6) + (F2*2/3))) + ((dt/(dy))*((G0/6) + (G1/6) + (G2*2/3))); 
            
        otherwise % RK4
            % Stage 1
            u1 = u + (dt/2)*((F0/dx) + (G0/dy));
    
            % Stage 2
            t1 = t0 + dt/2;
            [F1, G1] = GetFlux(X, Y, u1, t1, alpha, beta, f, g);

            u2 = u + ((dt/(2*dx))*(F1)) + ((dt/(2*dy))*(G1)); % at time t2 = tn + dt

            % Stage 3
            t2 = t0 + dt/2;
            [F2, G2] = GetFlux(X, Y, u2, t2, alpha, beta, f, g);

            u3 = u + ((dt/(dx))*(F2)) + ((dt/(dy))*(G2)); 

            % Stage 4
            t3 = t0 + dt;
            [F3, G3] = GetFlux(X, Y, u3, t3, alpha, beta, f, g);

            u = u + ((dt/(6*dx))*(F0 + 2*F1 + 2*F2 + F3)) + ((dt/(6*dy))*(G0 + 2*G1 + 2*G2 + G3)); 
    end
end

function [F, G] = GetFlux(X, Y, u, t, alpha, beta, f, g)
    F = zeros(size(X));
    G = zeros(size(Y));
    
    % first iterate to find F
    for i = 1:size(X, 1)
        xmid = X(:, i); ymid = Y(:, i); % hold y constant, vary x
        f_pos = (f(u(:, i), xmid, ymid, t) + (alpha.*u(:, i)))/2;
        f_neg = (f(u(:, i), xmid, ymid, t) - (alpha.*u(:, i)))/2;
        
        F(:, i) = ComputeFluxDifference(f_pos, f_neg);
    end

    % iterate to find G
    for j = 1:size(Y, 1)
        xmid = X(j, :)'; ymid = Y(j, :)'; % hold x constant, vary y
        g_pos = (g(u(j, :)', xmid, ymid, t) + (beta.*u(j, :)'))/2; 
        g_neg = (g(u(j, :)', xmid, ymid, t) - (beta.*u(j, :)'))/2;
        
        G(j, :) = ComputeFluxDifference(g_pos, g_neg);
    end
end

function [F] = ComputeFluxDifference(f_pos, f_neg)
    flux_left_lim = WENO(f_pos, 1); % left limit, right boundary
    flux_right_lim = WENO(f_neg, 0); % right limit, left boundary
    flux_right_lim = [flux_right_lim(2:end); flux_right_lim(1)];
    f_hat_pos = flux_left_lim + flux_right_lim;
    f_hat_neg = [f_hat_pos(end); f_hat_pos(1:end-1)];
    F = -(f_hat_pos - f_hat_neg);
end





























