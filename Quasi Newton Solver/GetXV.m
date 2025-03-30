function [X, R, Z, dx, dr, dz] = GetXV(Nx, N_vx, N_vz, interval)
    xvals = linspace(interval(1), interval(2), Nx+1)';
    rvals = linspace(interval(3), interval(4), N_vx+1)';
    zvals = linspace(interval(5), interval(6), N_vz+1)';
    dx = xvals(2) - xvals(1);
    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);
    xmid = xvals(1:end-1) + dx/2; % centered mesh
    rmid = rvals(1:end-1) + dr/2; % centered mesh
    zmid = zvals(1:end-1) + dz/2; % centered mesh
    [R, Z] = meshgrid(rmid, zmid); 
    X = xmid; R = R'; Z = Z';
end