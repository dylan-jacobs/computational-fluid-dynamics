% Completes approximation using IMEX(1, 1, 1)

function [U, ranks] = IMEX(type, U, dt, Nx, Ny, tf, interval, A, d1, d2, phi, tolerance)
    [X, Y, dx, dy] = GetXY(Nx, Ny, interval);

    % Discretize in space
    xvals = X(:, 1); yvals = Y(1, :)';
    A{1, 1} = @(t) A{1, 1}(xvals);
    A{1, 3} = @(t) A{1, 3}(yvals);
   
    A{2, 1} = @(t) A{2, 1}(xvals);
    A{2, 3} = @(t) A{2, 3}(yvals);
    phi{1, 1} = @(t) phi{1, 1}(xvals, t);
    phi{1, 3} = @(t) phi{1, 3}(yvals, t);  

    tvals = (0:dt:tf)';
    if tvals(end) ~= tf
        tvals = [tvals; tf];
    end

    [Vx, S, Vy] = svd(U);
    ranks = zeros(numel(tvals), 2);
    
    [Dx, Dy, Dxx, Dyy] = spectral_matrix(interval(2) - interval(1), Nx, Ny, dx, dy, d1, d2);

    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);

        switch type
            case '111'
                [Vx, S, Vy, ranks(n, 2)] = IMEX111_timestep(Vx, S, Vy, tvals(n-1), dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance);
            case '222'
                [Vx, S, Vy, ranks(n, 2)] = IMEX222_timestep(Vx, S, Vy, tvals(n-1), dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance);
            case '443'
                [Vx, S, Vy, ranks(n, 2)] = IMEX443_timestep(Vx, S, Vy, tvals(n-1), dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance);
  
        end

        % if mod(n, 1) == 0
        %     figure(1); clf; surf(X, Y, Vx*S*(Vy'));
        %     colorbar;
        %     shading flat; % removes gridlines
        %     legend(sprintf('N_x = %s, N_y = %s', num2str(Nx, 3), num2str(Ny, 3)), 'Location','northwest');
        %     xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('RAIL approximation at time %s', num2str(tf, 4))]);
        %     view(3);
        % end
    end
    U = Vx*S*(Vy');
    ranks(:, 1) = tvals;
end



