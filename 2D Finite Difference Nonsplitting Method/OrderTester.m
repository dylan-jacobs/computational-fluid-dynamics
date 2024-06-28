% Compute order for 2D Finite Difference Problems for Conservation Laws

function [output_table] = OrderTester(discretizationType, lambda, tf, u0, interval, f, g, u_exact_eqn, alpha, beta)
    Nxvals = [20, 40, 80, 160, 320]';   
    errors = zeros(numel(Nxvals), 2); % L1, L2 errors
    
    for i = 1:size(errors, 1)
        Nx = Nxvals(i);
        Ny = Nx; % uniform mesh
        
        % compute exact solution
        [X, Y, dx, dy] = GetXY(Nx, Ny, interval);
        u_exact = u_exact_eqn(X, Y);
        
        u = Finite_Differences_2D_Nonsplitting(discretizationType, Nx, Ny, lambda, interval, tf, f, g, u0, alpha, beta, 0);
        errors(i, 1) = dx*dy*(sum(sum(abs(u - u_exact)))); % L1 error
        errors(i, 2) = sqrt(dx*dy*sum(sum((u - u_exact).^2))); % L2 error
    end

    figure; clf; surf(X, Y, u); hold on; surf(X, Y, u_exact);
    colorbar;
    shading flat; % removes gridlines
    legend(sprintf('Nx = Ny = %s', num2str(Nx, 3)), 'Location', 'northwest');
    xlabel('X'); ylabel('Y'); zlabel('U(X, Y)'); title([sprintf('2D WENO+%s', discretizationType), sprintf(' approximation at time %s', num2str(tf, 2))]);

    L1_error = errors(:, 1);
    L2_error = errors(:, 2);
    L1_order = [0; log2(errors(1:end-1, 1) ./ errors(2:end, 1))];
    L2_order = [0; log2(errors(1:end-1, 2) ./ errors(2:end, 2))];
    
    output_table = table(L1_error, L1_order, L2_error, L2_order);
end



