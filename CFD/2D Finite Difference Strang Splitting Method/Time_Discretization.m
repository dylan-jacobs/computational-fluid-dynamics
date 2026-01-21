% Choose a type of time discretization: Forward Euler (SSP-RK1), SSP-RK2,
% SSP-RK3
% f = flux function

function [u] = Time_Discretization(type, u0, t0, dt, xmid, ymid, alpha, f)
    u = u0;
    dx = xmid(2) - xmid(1);
    if dx == 0
        dx = ymid(2) - ymid(1);
    end

    F0 = GetFlux(xmid, ymid, u, t0, alpha, f);

    switch type
        case 'RK1'
            u = u + (dt/dx)*F0;
        case 'RK2'
            % Stage 1
            u1 = u + (dt/(dx))*F0; % at time t1 = tn+dt
    
            % Stage 2
            t1 = t0 + dt;
            F1 = GetFlux(xmid, ymid, u1, t1, alpha, f);

            u = u + (dt/(2*dx))*(F0 + F1);
        case 'RK3'
            % Stage 1
            u1 = u + (dt/(dx))*F0; % at time t1 = tn+dt
    
            % Stage 2
            t1 = t0 + dt;
            F1 = GetFlux(xmid, ymid, u1, t1, alpha, f);

            u2 = u + (dt/(4*dx))*(F0 + F1); % at time t2 = tn + dt

            % Stage 3
            t2 = t0 + dt/2;
            F2 = GetFlux(xmid, ymid, u2, t2, alpha, f);

            u = u + (dt/(dx))*((F0/6) + (F1/6) + (F2*2/3)); 
            
        otherwise % RK4
            % Stage 1
            u1 = u + (dt/(2*dx))*F0; % at time t1 = tn+dt
    
            % Stage 2
            t1 = t0 + dt/2;
            F1 = GetFlux(xmid, ymid, u1, t1, alpha, f);

            u2 = u + (dt/(2*dx))*(F1); % at time t2 = tn + dt

            % Stage 3
            t2 = t0 + dt/2;
            F2 = GetFlux(xmid, ymid, u2, t2, alpha, f);

            u3 = u + (dt/(dx))*(F2); 

            % Stage 4
            t3 = t0 + dt;
            F3 = GetFlux(xmid, ymid, u3, t3, alpha, f);

            u = u + (dt/(6*dx))*(F0 + 2*F1 + 2*F2 + F3); 

    end
    % figure(2); clf;
    % xmid = xvals(1:end-1) + dx/2;
    % plot(xmid, u, 'b-', 'LineWidth', 1.5);
end


function [F] = GetFlux(X, Y, u, t, alpha, f) % can reuse for both f(u)|dx|alpha and g(u)|dy|beta
    f_pos = (f(u, X, Y, t) + (alpha.*u))/2;
    f_neg = (f(u, X, Y, t) - (alpha.*u))/2;
    
    flux_left_lim = WENO(f_pos, 1); % left limit, right boundary
    flux_right_lim = WENO(f_neg, 0); % right limit, left boundary
    flux_right_lim = [flux_right_lim(2:end); flux_right_lim(1)];
    f_hat_pos = flux_left_lim + flux_right_lim;
    f_hat_neg = [f_hat_pos(end); f_hat_pos(1:end-1)];
    F = -(f_hat_pos - f_hat_neg);
end
