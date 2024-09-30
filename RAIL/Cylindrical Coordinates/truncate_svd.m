% truncates Vx, S, Vy to given tolerance epsilon
% returns updated Vx, S, Vy and rank r1

function [Vr, S, Vz, rank] = truncate_svd(Vr, S, Vz, tolerance)  
    [U, sigma, V] = svd(S);
    rank = find(diag(sigma) > tolerance, 1, 'last');
    Vr = Vr*U(:, 1:rank);
    S = sigma(1:rank, 1:rank);
    Vz = Vz*V(:, 1:rank);
end