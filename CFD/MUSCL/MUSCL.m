% Approximates nonlinear PDE given spatial mesh size Nx, interval I,
% initial condition u0, final time tf, exact solution u_exact, and flux
% function f
function [] = MUSCL(xvals, u0, tf, u_exact, f)
    dx = xvals(2) - xvals(1);
    xmid = xvals(1:end-1) + dx/2;
    dt = dx.^2; % CFL condition
    tvals = (0:dt:tf)';
    if tvals(end)~=tf
        tvals = [tvals; tf];
    end
    
    u_bar = (0.5)*((u0(xvals(2:end)) + u0(xvals(1:end-1)))); % initial condition

    fig_idx = numel(findobj('type', 'figure')) + 1;
    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);
        u_mid = u_bar;
        u_pos = [u_bar(2:end) ; u_bar(1)];
        u_neg = [u_bar(end); u_bar(1:end-1)];

        u1_tilde_pos = 0.5*(u_mid - u_neg);
        u2_tilde_pos = 0.5*(u_pos - u_mid);

        f_pos = f(u_mid + minmod(u1_tilde_pos, u2_tilde_pos));
        f_neg = [f_pos(end); f_pos(1:end-1)];

        u_bar = u_bar - (dt/(dx))*(f_pos - f_neg);
        
        figure(fig_idx); clf;
        plot(xmid, u_exact, 'black--'); hold on;
        plot(xmid, u_bar, 'b-', 'LineWidth', 1.5);
        ylim([-0.2, 2.2]);
    end
    title(sprintf('2-Point MUSCL Scheme at t=%d', tf));
    xlabel('x'); ylabel('u'); legend('Exact', sprintf('Nx=%d', numel(xvals)-1));
end



















