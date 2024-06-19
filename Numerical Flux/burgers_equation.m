
% Approximates Burgers' Equation given spatial mesh size Nx, interval I,
% initial condition u0, final time tf, alpha, and exact solution u_exact
function [] = burgers_equation(Nx, I, u0, tf, alpha, u_exact)
    xvals = linspace(I(1), I(2), Nx + 1)';
    dx = xvals(2) - xvals(1);
    xmid = xvals(1:end-1) + dx/2;
    dt = 0.5*dx; % CFL condition: dt = 0.5dx
    tvals = (0:dt:tf)';
    if tvals(end)~=tf
        tvals = [tvals; tf];
    end
    
    u_bar = (0.5)*((u0(xvals(2:end)) + u0(xvals(1:end-1)))); % initial condition

    fig_idx = numel(findobj('type', 'figure')) + 1;
    for n = 2:numel(tvals)
        u_mid = u_bar;
        u_pos = [u_bar(2:end) ; u_bar(1)];
        f_pos = 0.25*(u_mid.^2 + u_pos.^2) + (alpha/2)*(u_mid - u_pos);
        f_neg = [f_pos(end); f_pos(1:end-1)];
    
        u_bar = u_bar - (dt/dx)*sign(alpha)*(f_pos - f_neg);
        figure(fig_idx); clf;
        plot(xmid, u_exact(xmid), 'black--'); hold on;
        plot(xvals(2:end), u_bar, 'b-', 'LineWidth', 1.5);
        ylim([-1.2, 1.2]);
    end
    title(sprintf('Lax-Friedrichs Method at t=%d', tf));
    xlabel('x'); ylabel('u'); legend('Exact', sprintf('Lax-Friedrichs Approximation Nx=%d', Nx));
end



















