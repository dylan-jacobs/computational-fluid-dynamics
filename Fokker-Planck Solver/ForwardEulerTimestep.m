% Facilitates single forward Euler timestep for the 1D Fokker-Planck system

function [f] = ForwardEulerTimestep(f, dt, dx, tval, xvals, B, D)
    N = numel(xvals);
    [weights, nodes] = GetWeightsAndNodes(7);
    
    % lambdavals = zeros(N, 1);
    % % lambdavals(1) = GaussLegendre(xvals(1) - dx, xvals(1), weights, nodes, @(x) C(x, tval));
    % for i = 1:N-1
    %     lambdavals(i) = GaussLegendre(xvals(i), xvals(i+1), weights, nodes, @(x) B(x, tval));
    % end
    % lambdavals(N) = GaussLegendre(xvals(N), xvals(end) + dx, weights, nodes, @(x) B(x, tval));
    % a = lambdavals(1:5);
    % lambdavals = 0.5*(B([xvals(2:end); xvals(end)+dx], tval).^2) - 0.5*(B(xvals, tval).^2);
    % b = lambdavals(1:5);
    % lambdavals = lambdavals + 1E-14;
    % deltavals = (1./lambdavals) + (1./(1-exp(lambdavals)));
    % 
    % D = D(xvals + (dx/2));
    % 
    % C_tilde = D.*lambdavals./dx;
    % 
    % C_tilde_matrix = diag(C_tilde);
    % f_tilde_matrix = gallery('tridiag', N, zeros(1, N-1), deltavals, 1-deltavals(1:end-1));
    % D_matrix       = gallery('tridiag', N, zeros(1, N-1), -D, D(1:end-1));
    % F_matrix = (C_tilde_matrix*f_tilde_matrix) + ((1/dx)*D_matrix);
    % 
    % F_pos = F_matrix*f;
    % F_neg = [0; F_pos(1:end-1)];
    % 
    % F_pos(end) = 0;
    
    D = D(xvals + (dx/2));
    C_tilde = zeros(N, 1);
    for i = 1:N-1
        C_tilde(i) = GaussLegendre(xvals(i), xvals(i+1), weights, nodes, @(x) B(x, tval));
    end
    C_tilde(N) = GaussLegendre(xvals(N), xvals(N) + dx, weights, nodes, @(x) B(x, tval));
    C_tilde = C_tilde + 1E-14;
    C_tilde = (D./dx).*C_tilde;
    lambdavals = (dx./D).*C_tilde;
    deltavals = (1./lambdavals) + (1./(1-exp(lambdavals)));

    C_tilde_matrix = diag(C_tilde);
    f_tilde_matrix = gallery('tridiag', N, zeros(1, N-1), deltavals, 1-deltavals(1:end-1));
    D_matrix       = gallery('tridiag', N, zeros(1, N-1), -D, D(1:end-1));
    F_matrix = (C_tilde_matrix*f_tilde_matrix) + ((1/dx)*D_matrix);
    
    F_pos = F_matrix*f;


    F_neg = [0; F_pos(1:end-1)];
    F_pos(end) = 0;
    
    f = f + ((dt/dx)*(F_pos - F_neg));
end