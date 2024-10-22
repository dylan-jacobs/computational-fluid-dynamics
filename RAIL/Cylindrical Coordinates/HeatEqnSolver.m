function [U, ranks] = HeatEqnSolver(type, U, dt, Nr, Nz, tf, interval, tolerance, r0)
    [R, Z, dr, dz] = GetRZ(Nr, Nz, interval);

    rvals = R(:, 1); % get one column of r values
    tvals = (0:dt:tf)';
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end

    Drr = gallery('tridiag', Nr, rvals(1:end-1)+(dr/2), -2*rvals, rvals(1:end-1)+(dr/2));
    Drr(1, 1) = -(rvals(1) + (dr/2));
    Drr(end, end-1) = ((1/3) * (rvals(end) + (dr/2))) + (rvals(end) - (dr/2));
    Drr(end, end) =  ((-3) * (rvals(end) + (dr/2))) - (rvals(end) - (dr/2));
    Drr = (1/(dr^2)) * (diag(1./rvals) * Drr);


    % Dzz = gallery('tridiag', Nz, 1, -2, 1); % centered nodes
    Dzz = (2*pi/interval(4))^2*toeplitz([-1/(3*(2*dz/interval(4))^2)-1/6 ...
  .5*(-1).^(2:Nz)./sin((2*pi*dz/interval(4))*(1:Nz-1)/2).^2]);

    [Vr, S, Vz] = svd2(U, rvals);
    ranks = zeros(numel(tvals), 2);
    ranks(:, 1) = tvals;
    r0 = min(r0, size(Vr, 2));

    % truncate solution to pre-set initial rank r0
    Vr = Vr(:, 1:r0);
    S = S(1:r0, 1:r0);
    Vz = Vz(:, 1:r0);
    ranks(1, 2) = r0;
  
    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);

        switch type
            case '1'
                [Vr, S, Vz, ranks(n, 2)] = BackwardEuler(Vr, S, Vz, rvals, dt, Drr, Dzz, tolerance);
            case '2'
                [Vr, S, Vz, ranks(n, 2)] = DIRK2(Vr, S, Vz, rvals, dt, Drr, Dzz, tolerance);
            case '3'
                [Vr, S, Vz, ranks(n, 2)] = DIRK3(Vr, S, Vz, rvals, dt, Drr, Dzz, tolerance);
        end

        % if mod(n, 1) == 0
        %     figure(1); clf; surf(R, Z, Vr*S*(Vz'));
        %     colorbar;
        %     shading flat; % removes gridlines
        %     xlabel('R'); ylabel('Z'); zlabel('U(R, Z)'); title([sprintf('RAIL approximation at time %s', num2str(tf, 4))]);
        %     view(3);
        % end
    end
    U = Vr*S*(Vz');
end

function [U, S, V] = svd2(A, rvals)
    [U, S, V] = svd(sqrt(rvals) .* A, 0);
    U = (1./sqrt(rvals)) .* U;
end