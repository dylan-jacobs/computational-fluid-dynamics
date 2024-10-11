function [R, Z, dr, dz] = GetRZ(Nr, Nz, interval)
    rvals = linspace(interval(1), interval(2), Nr+1)';
    zvals = linspace(interval(3), interval(4), Nz+1)';
    dr = rvals(2) - rvals(1);
    dz = zvals(2) - zvals(1);
    rmid = rvals(1:end) + dr/2; % centered mesh
    zmid = zvals(2:end) + dz/2; % centered mesh
    [R, Z] = meshgrid(rmid, zmid); R = R'; Z = Z';
end