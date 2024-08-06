% truncates Vx, S, Vy to given tolerance epsilon
% returns updated Vx, S, Vy and rank r1

function [Vx, S, Vy, r1] = truncate_svd(Vx, S, Vy, tolerance)  
    [U, sigma, V] = svd(S);
    r1 = find(diag(sigma) > tolerance, 1, 'last');
    Vx = Vx*U(:, 1:r1);
    S = sigma(1:r1, 1:r1);
    Vy = Vy*V(:, 1:r1);
end