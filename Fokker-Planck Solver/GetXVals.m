function [xmid, dx] = GetXVals(Nx, interval)
    xvals = linspace(interval(1), interval(2), Nx+1)';
    dx = xvals(2) - xvals(1);
    xmid = xvals(1:end-1) + dx/2;
end