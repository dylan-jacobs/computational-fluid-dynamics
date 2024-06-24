% Choose a type of time discretization: Forward Euler (SSP-RK1), SSP-RK2,
% SSP-RK3
% f = flux function

function [u] = Time_Discretization(type, u0, t0, dt, xmid, ymid, alpha, f)
    u = u0;
    dx = xmid(2) - xmid(1);
    if dx == 0
        dx = ymid(2) - ymid(1);
    end
    f_pos = (f(u, xmid, ymid, t0) + (alpha.*u))/2;
    f_neg = (f(u, xmid, ymid, t0) - (alpha.*u))/2;
    
    flux_left_lim = WENO(f_pos, 1); % left limit, right boundary
    flux_right_lim = WENO(f_neg, 0); % right limit, left boundary
    flux_right_lim = [flux_right_lim(2:end); flux_right_lim(1)];
    f_hat_pos = flux_left_lim + flux_right_lim;
    f_hat_neg = [f_hat_pos(end); f_hat_pos(1:end-1)];
    F0 = -(f_hat_pos - f_hat_neg);

    switch type
        case 'RK1'
            u = u + (dt/dx)*F0;
        case 'RK2'
            % Stage 1
            u1 = u + (dt/(dx))*F0; % at time t1 = tn+dt
    
            % Stage 2
            t1 = t0 + dt;

            f_pos = (f(u1, xmid, ymid, t1) + (alpha*u1))/2;
            f_neg = (f(u1, xmid, ymid, t1) - (alpha*u1))/2;
            
            flux_left_lim = WENO(f_pos, 1); % left limit, right boundary
            flux_right_lim = WENO(f_neg, 0); % right limit, left boundary
            flux_right_lim = [flux_right_lim(2:end); flux_right_lim(1)];
            f_hat_pos = flux_left_lim + flux_right_lim;
            f_hat_neg = [f_hat_pos(end); f_hat_pos(1:end-1)];
            F1 = -(f_hat_pos - f_hat_neg);

            u = u + (dt/(2*dx))*(F0 + F1);
        otherwise % RK3
            % Stage 1
            u1 = u + (dt/(dx))*F0; % at time t1 = tn+dt
    
            % Stage 2
            t1 = t0 + dt;

            f_pos = (f(u1, xmid, ymid, t1) + (alpha*u1))/2;
            f_neg = (f(u1, xmid, ymid, t1) - (alpha*u1))/2;
            
            flux_left_lim = WENO(f_pos, 1); % left limit, right boundary
            flux_right_lim = WENO(f_neg, 0); % right limit, left boundary
            flux_right_lim = [flux_right_lim(2:end); flux_right_lim(1)];
            f_hat_pos = flux_left_lim + flux_right_lim;
            f_hat_neg = [f_hat_pos(end); f_hat_pos(1:end-1)];
            F1 = -(f_hat_pos - f_hat_neg);
            u2 = u + (dt/(4*dx))*(F0 + F1);

            % Stage 3
            t2 = t0 + dt/2;

            f_pos = (f(u2, xmid, ymid, t2) + (alpha*u2))/2;
            f_neg = (f(u2, xmid, ymid, t2) - (alpha*u2))/2;

            flux_left_lim = WENO(f_pos, 1); % left limit, right boundary
            flux_right_lim = WENO(f_neg, 0); % right limit, left boundary
            flux_right_lim = [flux_right_lim(2:end); flux_right_lim(1)];
            f_hat_pos = flux_left_lim + flux_right_lim;
            f_hat_neg = [f_hat_pos(end); f_hat_pos(1:end-1)];
            F2 = -(f_hat_pos - f_hat_neg);

            u = u + (dt/(dx))*((F0/6) + (F1/6) + (F2*2/3)); 

    end
    % figure(2); clf;
    % xmid = xvals(1:end-1) + dx/2;
    % plot(xmid, u, 'b-', 'LineWidth', 1.5);
end
