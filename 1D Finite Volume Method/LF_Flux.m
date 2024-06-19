% computes the Lax-Friedrichs Flux
% f_hat = 0.5(f(u-) + f(u+) - alpha(u+ - u-))
% alpha = max(f'(u))

function [flux_pos, flux_neg] = LF_Flux(f, rightlim, leftlim, alpha)
    flux_pos = 0.5*(f(leftlim) + f(rightlim) - (alpha*(rightlim - leftlim)));
    flux_neg = [flux_pos(end); flux_pos(1:end-1)];
end


























