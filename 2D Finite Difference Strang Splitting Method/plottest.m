Nx = 100;
Ny = 100;
xvals = linspace(0, 2*pi, Nx+1)';
yvals = linspace(0, 2*pi, Ny+1)';
[X, Y] = meshgrid(xvals, yvals);
X = X'; Y = Y';
u = sin(X+Y);
figure(1); clf; surf(X, Y, u); axis([0, 2*pi, 0, 2*pi, -1.1, 1.1] );
colorbar;
shading interp; % removes gridlines
% figure(1); view(2); % bird's eye view