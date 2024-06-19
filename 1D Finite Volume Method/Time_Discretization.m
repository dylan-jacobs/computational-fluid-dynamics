% Choose a type of time discretization: Forward Euler (SSP-RK1), SSP-RK2,
% SSP-RK3
% f = flux function

function [u] = Time_Discretization(type, u0, tvals, xvals, alpha, f)
    dx = xvals(2) - xvals(1);
    xmid = xvals(1:end-1) + dx/2;
    u = u0;

    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);
        t0 = tvals(n-1);
        leftlim = WENO(u, 1); % left limit, right boundary (1 = left limit at right boundary, 0 = right limit at left boundary)
        rightlim = WENO(u, 0);
        rightlim = [rightlim(2:end); rightlim(1)]; % right limit, right boundary
        [flux_pos, flux_neg] = LF_Flux(@(u) f(xvals(2:end), t0, u), rightlim, leftlim, alpha);
        F0 = -(flux_pos - flux_neg);
         
        switch type
            case 'RK1'
                u = u + (dt/dx)*F0;
            case 'RK2'
                % Stage 1
                u1 = u + (dt/(dx))*F0; % at time t1 = tn+dt
        
                % Stage 2
                t1 = t0 + dt;
                leftlim = WENO(u1, 1); % left limit, right boundary (1 = left limit at right boundary, 0 = right limit at left boundary)
                rightlim = WENO(u1, 0);
                rightlim = [rightlim(2:end); rightlim(1)]; % right limit, right boundary
                [flux_pos, flux_neg] = LF_Flux(@(v) f(xvals(2:end), t1, v), rightlim, leftlim, alpha);
                F1 = -(flux_pos - flux_neg);
                u = u + (dt/(2*dx))*(F0 + F1); 
            otherwise
                % Stage 1
                u1 = u + (dt/(dx))*F0; % at time t1 = tn+dt
        
                % Stage 2
                t1 = t0 + dt;
                leftlim = WENO(u1, 1); % left limit, right boundary (1 = left limit at right boundary, 0 = right limit at left boundary)
                rightlim = WENO(u1, 0);
                rightlim = [rightlim(2:end); rightlim(1)]; % right limit, right boundary
                [flux_pos, flux_neg] = LF_Flux(@(v) f(xvals(2:end), t1, v), rightlim, leftlim, alpha);
                F1 = -(flux_pos - flux_neg);
                u2 = u + (dt/(4*dx))*(F0 + F1); 

                % Stage 3
                t2 = t0 + dt/2;
                leftlim = WENO(u2, 1); % left limit, right boundary (1 = left limit at right boundary, 0 = right limit at left boundary)
                rightlim = WENO(u2, 0);
                rightlim = [rightlim(2:end); rightlim(1)]; % right limit, right boundary
                [flux_pos, flux_neg] = LF_Flux(@(v) f(xvals(2:end), t2, v), rightlim, leftlim, alpha);
                F2 = -(flux_pos - flux_neg);
                u = u + (dt/(dx))*((F0/6) + (F1/6) + (F2*2/3)); 

        end
        % figure(2); clf;
        % xmid = xvals(1:end-1) + dx/2;
        % plot(xmid, u, 'b-', 'LineWidth', 1.5);
    end
end
