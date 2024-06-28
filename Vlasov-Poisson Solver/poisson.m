% AUTHOR: J. Nakao
% PURPOSE: Poisson solver for Vlasov-Poisson system (1d1v).
% 
% 
function E = poisson(rho, Nx, xLength)
% Inputs:
    % rho - charge density
    % Nx - mesh size
    % xlength - period wrt x
% Outputs:
    % E - electric field
    % Emax - Linf norm of electric field

% fft mode
k_rescale = 2*pi/xLength*[0:Nx/2-1 -Nx/2:-1]';
% coefficient for E with periodic b.c.
E_mode = -fft(rho)*sqrt(-1)./k_rescale;
E_mode(1) = 0;
% E via ifft from E_mode
E = real(ifft(E_mode));
end