% Facilitates single forward Euler timestep for the 1D Fokker-Planck system

function [f] = ForwardEulerTimestep(f, dt, dx, tval, xvals, B, D)
    N = numel(xvals);
    [weights, nodes] = GetWeightsAndNodes(7);
       
    D = D(xvals + (dx/2));
    lambdavals = zeros(N, 1);
    for i = 1:N-1
        lambdavals(i) = GaussLegendre(xvals(i), xvals(i+1), weights, nodes, @(x) B(x, tval)./D(i));
    end
    lambdavals(N) = GaussLegendre(xvals(N), xvals(N) + dx, weights, nodes, @(x) B(x, tval)./D(N));
    C_tilde = (D./dx).*lambdavals;
    lambdavals = lambdavals + 1E-14; % prevent zero division
    deltavals = (1./lambdavals) + (1./(1-exp(lambdavals)));

    C_tilde_matrix = diag(C_tilde);
    f_tilde_matrix = gallery('tridiag', N, zeros(1, N-1), deltavals, 1-deltavals(1:end-1));
    D_matrix       = gallery('tridiag', N, zeros(1, N-1), -D, D(1:end-1));
    F_matrix = (C_tilde_matrix*f_tilde_matrix) + ((1/dx)*D_matrix);
    
    F_pos = F_matrix*f;
    
    F_pos(end) = 0; % zero-flux boundary conditions!
    F_neg = [0; F_pos(1:end-1)]; 
    
    f = f + ((dt/dx)*(F_pos - F_neg));
end












