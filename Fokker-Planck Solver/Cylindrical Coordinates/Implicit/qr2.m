function [Q, R] = qr2(V0, rvals)
    [Q, R] = qr(sqrt(rvals) .* V0, 0);
    Q = (1./sqrt(rvals)).*Q;
end