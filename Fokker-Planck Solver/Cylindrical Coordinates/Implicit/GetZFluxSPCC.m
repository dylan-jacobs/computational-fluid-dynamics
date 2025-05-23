function [Flux] = GetZFluxSPCC(tval, zvals, B, D)
    dz = zvals(2) - zvals(1);
    N = numel(zvals);

    [weights, nodes] = GetWeightsAndNodes(7);

    D = D(zvals + (dz/2));
    lambdavals = zeros(N, 1);
    for i = 1:N-1
        lambdavals(i) = GaussLegendre(zvals(i), zvals(i+1), weights, nodes, @(x) B(x, tval)./D(i));
    end
    lambdavals(N) = GaussLegendre(zvals(N), zvals(N) + dz, weights, nodes, @(x) B(x, tval)./D(N));
    lambdavals = lambdavals + 1E-12; % prevent zero division
    C_tilde = (D./dz).*lambdavals;
    deltavals = (1./lambdavals) + (1./(1-exp(lambdavals)));
    C_tilde_matrix = spdiags(C_tilde, 0, N, N);
    f_tilde_matrix = gallery('tridiag', N, zeros(1, N-1), deltavals, 1-deltavals(1:end-1));
    D_matrix       = gallery('tridiag', N, zeros(1, N-1), -D, D(1:end-1));
    F_pos = (C_tilde_matrix*f_tilde_matrix) + ((1/dz).*D_matrix);
    F_pos(end, :) = zeros(1, N);
    F_neg = [zeros(1, N); F_pos(1:end-1, :)];
    Flux = (1/dz)*(F_pos - F_neg);
end