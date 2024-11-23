% Facilitates single forward Euler timestep for the 1D Fokker-Planck system

function [f] = ForwardEulerTimestep(f, dt, dx, tval, xvals, C, D)
    D = D(xvals);
    N = numel(xvals);
    lambdavals = 0.5*(C([xvals(2:end); xvals(end)+dx], tval).^2) - 0.5*(C([xvals(1)-dx; xvals(1:end-1)], tval).^2);
    deltavals = (1./lambdavals) + (1./(1-exp(lambdavals)));
    C_tilde = D.*lambdavals./dx;
    
    C_tilde_matrix = diag(C_tilde);
    f_tilde_matrix = gallery('tridiag', N, zeros(1, N-1), deltavals, 1-deltavals(1:end-1));
    D_matrix       = gallery('tridiag', N, zeros(1, N-1), -D, D(1:end-1));

    F_matrix = (C_tilde_matrix*f_tilde_matrix) + ((1/dx)*D_matrix);

    F_pos = F_matrix*f;
    F_neg = [0; F_pos(1:end-1)];

    F_pos(end) = 0;

    f = f + (dt/dx)*(F_pos - F_neg);
end