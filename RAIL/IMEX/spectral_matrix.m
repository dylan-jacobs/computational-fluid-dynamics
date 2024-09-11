% Spectrally accurate differentiation matrix [Trefethen, 2000]
function [Dx, Dy, Dxx, Dyy] = spectral_matrix(L, Nx, Ny, dx, dy, d1, d2)
    % First derivative matrices Dx, Dy if we ever need them
    columnx = [0 0.5*(-1).^(1:Nx-1).*cot((1:Nx-1)*(2*pi*dx/(2*L)))];
    Dx = (2*pi/L)*toeplitz(columnx,columnx([1 Nx:-1:2]));
    columny = [0 0.5*(-1).^(1:Ny-1).*cot((1:Ny-1)*(2*pi*dy/(2*L)))];
    Dy = (2*pi/L)*toeplitz(columny,columny([1 Ny:-1:2]));

    % Second derivative matrices Dxx, Dyy
    Dxx = (2*pi/L)^2*toeplitz([-1/(3*(2*dx/L)^2)-1/6 ...
    .5*(-1).^(2:Nx)./sin((2*pi*dx/L)*(1:Nx-1)/2).^2]);
    
    Dyy = (2*pi/L)^2*toeplitz([-1/(3*(2*dy/L)^2)-1/6 ...
    .5*(-1).^(2:Ny)./sin((2*pi*dy/L)*(1:Ny-1)/2).^2]);
    Dxx = d1*Dxx; Dyy = d2*Dyy;
end