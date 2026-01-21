function [X, Y, dx, dy] = GetXY(Nx, Ny, interval)
    xvals = linspace(interval(1), interval(2), Nx+1)';
    yvals = linspace(interval(3), interval(4), Ny+1)';
    dx = xvals(2) - xvals(1);
    dy = yvals(2) - yvals(1);
    xmid = xvals(1:end-1) + dx/2;
    ymid = yvals(1:end-1) + dy/2;
    [X, Y] = meshgrid(xmid, ymid); X = X'; Y = Y';
end