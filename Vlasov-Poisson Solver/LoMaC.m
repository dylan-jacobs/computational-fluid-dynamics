function [f_truncated] = LoMaC(f, vvals, dv, tolerance)
%LoMaC Truncates given maxwellian (assumed Low-Rank) while conserving
%macroscopic quantities.
%   Detailed explanation goes here

    % Step 1: Integrate to calculate macro quantities
    p = dv*sum(f);
    J = dv*sum(f .* vvals);
    k = dv*0.5*sum(f .* (vvals.^2));

    % Step 2: Scale by maxwellian to ensure inner product is well defined
    % (f -> 0 and v -> infinity)
    w = exp(-(vvals.^2)/2);
    f_tilde = f ./ w;

    % Step 3: Orthogonal projection
    % bases:
    % 1, v, v.^2 - c
    w_norm = @(x) sqrt(dv*sum(x.^2 .* w));
    c = dv*sum((vvals.^2) .* w) / w_norm(1);
    f1 = ((p ./ w_norm(1)).*w) + ((J ./ w_norm(vvals)) .* w .* vvals) + (((2*k - c.*p) ./ w_norm((vvals.^2) - c)) .* w .* ((vvals.^2) - c));
    f2 = f - f1;

    % Step 4: Perform SVD truncation on non-mass part of f (f2)
    [U, S, V] = svd(f2./sqrt(w), 0);

    rank = find(diag(S) > tolerance, 1, 'last');
    if (sum(rank) == 0)
        rank = 1;
    end
    rank = min(rank, min(size(U, 2), size(V, 2)));

    U = U(:, 1:rank);
    S = S(1:rank, 1:rank);
    V = V(:, 1:rank);

    f2_truncated = sqrt(w).*(U*S*V');

    f_truncated = f1 + f2_truncated;

end