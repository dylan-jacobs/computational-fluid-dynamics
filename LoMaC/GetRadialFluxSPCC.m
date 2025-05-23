function [Flux] = GetRadialFluxSPCC(tval, rvals, B, D)
    dr = rvals(2) - rvals(1);
    N = numel(rvals);

    [weights, nodes] = GetWeightsAndNodes(7);
    D = D(rvals + (dr/2));
    lambdavals = zeros(N, 1);
    for i = 1:N-1
        lambdavals(i) = GaussLegendre(rvals(i), rvals(i+1), weights, nodes, @(x) B(x, tval)./D(i));
    end
    lambdavals(N) = GaussLegendre(rvals(N), rvals(N) + dr, weights, nodes, @(x) B(x, tval)./D(N));
    lambdavals = lambdavals + 1E-14; % prevent zero division
    C_tilde = (D./dr).*lambdavals;
    deltavals = (1./lambdavals) + (1./(1-exp(lambdavals)));
    C_tilde_matrix = spdiags(C_tilde, 0, N, N);
    f_tilde_matrix = gallery('tridiag', N, zeros(1, N-1), deltavals, 1-deltavals(1:end-1));
    D_matrix       = gallery('tridiag', N, zeros(1, N-1), -D, D(1:end-1));

    F_pos = (C_tilde_matrix*f_tilde_matrix) + ((1/dr).*D_matrix);
    F_pos(end, :) = zeros(1, N);
    F_neg = [zeros(1, N); F_pos(1:end-1, :)];

    r_inv = spdiags(1./rvals, 0, N, N);
    r_pos = spdiags(rvals + (dr/2), 0, N, N);
    r_neg = spdiags(rvals - (dr/2), 0, N, N);
    Flux = (1/dr)*r_inv*((r_pos*F_pos) - (r_neg*F_neg));
end