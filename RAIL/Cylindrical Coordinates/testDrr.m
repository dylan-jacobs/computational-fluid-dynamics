



Nrvals = [10,20,40,80,160];
L1vals = zeros(numel(Nrvals),1);

for k = 1:numel(Nrvals)

    Nr = Nrvals(k);
    rvals = linspace(0,3,Nr+1)';
    dr = rvals(2)-rvals(1);
    rvals = rvals(1:Nr)+dr/2;

    Drr = gallery('tridiag', Nr, rvals(1:end-1)+(dr/2), -2*rvals, rvals(1:end-1)+(dr/2));
    Drr(1, 1) = -(rvals(1) + (dr/2));
    Drr = (1./(dr^2)) .* (diag(1./rvals) * Drr);

    Drru = Drr*exp(-4*rvals.^2);

    Drru_exact = -8*exp(-4*rvals.^2).*(2 - 8*rvals.^2);

    L1vals(k) = dr*sum(abs(Drru - Drru_exact));
end

disp('Order = ')
disp(log2(L1vals(1:end-1) ./ (L1vals(2:end))))











