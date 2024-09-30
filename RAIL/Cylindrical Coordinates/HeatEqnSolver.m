function [U, ranks] = HeatEqnSolver(type, U, dt, Nr, Nz, tf, interval, d1, d2, tolerance, r0)
    [R, Z, dr, dz] = GetRZ(Nr, Nz, interval);

    tvals = (0:dt:tf)';
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end

    Drr = (1./(dr^2)) .* (diag(1./R) .* gallery('tridiag', Nr, R(2:end, 1)-(dr/2), -2*R(1:end, 1), R(1:end-1, 1)+(dr/2)));
    Drr(1, 1) = -(R(1) + (dr/2));
    Dzz = (1./(dz^2)) .* gallery('tridiag', Nz, 1, -2, 1);

    [Vr, S, Vz] = svd(U);
    ranks = zeros(numel(tvals), 2);
    ranks(:, 1) = tvals;
    r0 = min(r0, size(Vr, 2));

    % truncate solution to pre-set initial rank r0
    Vr = Vr(:, 1:r0);
    S = S(1:r0, 1:r0);
    Vz = Vz(:, 1:r0);
    ranks(1) = r0;
  
    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);

        switch type
            case '1'
                [Vr, S, Vz, ~] = BackwardEulerTimestep(Vr, S, Vz, R, dt, Drr, Dzz, tolerance);
        end

        % if mod(n, 1) == 0
        %     figure(1); clf; surf(R, Z, U);
        %     colorbar;
        %     shading flat; % removes gridlines
        %     xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('RAIL approximation at time %s', num2str(tf, 4))]);
        %     view(3);
        % end
    end
    U = Vr*S*(Vz');
end