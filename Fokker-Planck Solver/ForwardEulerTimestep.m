% Facilitates single forward Euler timestep for the 1D Fokker-Planck system

function [f] = ForwardEulerTimestep(f, dt, dx, xvals, C, D)

    N = numel(xvals);
    lambdavals = 0.5*(C.^2);
    deltavals = (1./lambdavals) + (1./(1-exp(lambdavals)));
    C_tilde = D.*lambdavals./dx;
    
    C_tilde_matrix = diag(C_tilde);
    f_tilde_matrix = gallery('tridiag', N, zeros(1, N-1), deltavals, 1-deltavals(1:end-1));
    D_matrix       = gallery('tridiag', N, zeros(1, N-1), -D, D(1:end-1));

    F_matrix = (C_tilde_matrix*f_tilde_matrix) + ((1/dx).*D_matrix);

    F_pos = F_matrix*f;
    F_neg = [0; F_pos(1:end-1)];

    F_pos(1) = 0; F_pos(end) = 0;
    F_neg(1) = 0; F_neg(end) = 0;

    f = f + (dt/dx).*(F_pos - F_neg);
end