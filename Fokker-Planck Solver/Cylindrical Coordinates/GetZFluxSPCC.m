function [Flux] = GetZFluxSPCC(tval, xvals, B, D)
    dx = xvals(2) - xvals(1);
    N = numel(xvals);

    [weights, nodes] = GetWeightsAndNodes(7);

    D = D(xvals + (dx/2));
    lambdavals = zeros(N, 1);
    for i = 1:N-1
        lambdavals(i) = GaussLegendre(xvals(i), xvals(i+1), weights, nodes, @(x) B(x, tval)./D(i));
    end
    lambdavals(N) = GaussLegendre(xvals(N), xvals(N) + dx, weights, nodes, @(x) B(x, tval)./D(N));
    lambdavals = lambdavals + 1E-12; % prevent zero division
    C_tilde = (D./dx).*lambdavals;
    deltavals = (1./lambdavals) + (1./(1-exp(lambdavals)));
    C_tilde_matrix = diag(C_tilde);
    f_tilde_matrix = gallery('tridiag', N, zeros(1, N-1), deltavals, 1-deltavals(1:end-1));
    D_matrix       = gallery('tridiag', N, zeros(1, N-1), -D, D(1:end-1));
    F_pos = (C_tilde_matrix*f_tilde_matrix) + ((1/dx).*D_matrix);
    F_pos(end, :) = zeros(1, N);
    F_neg = [zeros(1, N); F_pos(1:end-1, :)];
    Flux = (1/dx)*(F_pos - F_neg);
end