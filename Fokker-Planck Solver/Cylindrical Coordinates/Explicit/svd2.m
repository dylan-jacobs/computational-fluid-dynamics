function [U, S, V] = svd2(A, rvals)
    [U, S, V] = svd(sqrt(rvals) .* A, 0);
    U = (1./sqrt(rvals)) .* U;
end