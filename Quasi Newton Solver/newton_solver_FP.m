% Given current solution f, solves for the macroscopic parameters of the Fokker-Planck equation 
% mass (n), momentum (nu), kinetic energy (nU), temperature (Te)
% f = [v_para, v_perp, space]

% TEST ON 4.4.2!!!

function [n, nu, nU, T] = newton_solver_FP(f, n0, u_para0, U0, Te0, dt, dv_para, dv_perp, v_para, v_perp, qa, qe, ma, me, x_end)

    Nx = numel(n0); % max val of i, j
    vx = size(f, 2); % velocity dimens

    % get boundary conditions
    [n0, ~, T_ae] = get_boundaries(); % Ta = Te at boundary --> T_ae

    Te0_pos = [Te0(2:end), T_ae(x_end)]; Te0_neg = [T_ae(0), Te0(1:end-1)];
    
    % first, compute fluxes via summation
    f_hat = get_f_hat(f, v_para, x_end);
    nu_hat = zeros(Nx+1, vx, vx);
    S_hat = zeros(Nx+1, vx, vx);
    Q_hat = zeros(Nx+1, vx, vx);
    nTe_hat = zeros(Nx+1, vx, vx);
    for i = 1:Nx+1
        nu_hat(i) = 2*pi*(sum(sum(v_para(i) .* f_hat(:, :, i) .* (v_perp(i).*dv_para.*dv_perp))));
    
        S_hat(i) = 2*pi*(sum(sum(v_para(i).^2 .* f_hat(:, :, i) .* (v_perp(i).*dv_para.*dv_perp))));
    
        Q_hat(i) = pi*(sum(sum(v_para(i).^3 .* f_hat(:, :, i) .* (v_perp(i).*dv_para.*dv_perp))));
    
        nTe_i = n0(i)*Te0(i);
        nTe_i1 = n0(i+1)*Te0(i+1);
        nTe_hat(i) = (1./u_para0) .* ((nTe_i)*(u_para0(i) > 0) + nTe_i1.*(u_para0(i) <= 0)); % more upwinding
    end

    % shift bounds to get flux pos/neg
    nu_hat_pos = nu_hat(2:end); nu_hat_neg = [0, nu_hat(1:end-1)];
    S_hat_pos = S_hat(2:end); S_hat_neg = S_hat(1:end-1);
    Q_hat_pos = Q_hat(2:end); Q_hat_neg = Q_hat(1:end-1);
    nTe_hat_pos = nTe_hat(2:end); nTe_hat_neg = nTe_hat(1:end-1);

    % explicitly find n via Forward Euler
    n = n0 - (dt/dx)*(nu_hat_pos - nu_hat_neg); % now find n_i+1, n_i-1
    n_pos = [n(2:end), n0(x_end)]; n_neg = [n0(0), n(1:end-1)];

    % init y_vec, R_norm    
    y = [n0.*u_para0; n0.*U0; Te0];

    kappa_pos = (3.2/(2*sqrt(2*me)))*((Te0_pos.^(5/2) + Te0.^(5/2)));
    kappa_neg = (3.2/(2*sqrt(2*me)))*((Te0.^(5/2) + Te0_neg.^(5/2)));
    R_norm = 1;
    tol = min(5e-12);

    while (R_norm > tol)
        R1 = y(1) - (n0.*u_para0) + (dt/dx).*(S_hat_pos - S_hat_neg) - ((dt*qa)/(2*dx*qe*ma)).*((n_pos.*Te0_pos) - (n_neg.*Te0_neg));
        R2 = y(2) - (n0.*U0) + (dt/dx).*(Q_hat_pos - Q_hat_neg) - ((dt*qa*u_para0)./(2*dx*qe*ma)).*((n_pos.*Te0_pos) - (n_neg.*Te0_neg)) - (((dt.*3.*sqrt(2*me).*(n.^2))./((ma.^2).*(Te0.^(3/2)))) .* (Te0 - (ma/3).*(2*U0 - (u_para0.^2))));
        R3 = (y(3)./n) - (Te0) + ((5*dt)/(3*dx)).*(nTe_hat_pos - nTe_hat_neg) - (((dt*u_para0)./(2*dx)).*((n_pos.*Te0_pos) - (n_neg.*Te0_neg))) - (((2*dt)/(3*dx.^2)) .* ((kappa_pos.*(Te0_pos - Te0)) - (kappa_neg.*(Te0 - Te0_neg)))) - (((dt.*2.*sqrt(2*me).*(n.^2))./(ma.*(Te0.^(3/2)))) .* (-Te0 + (ma/3).*(2*U0 - (u_para0.^2))));
        R = [R1; R2; R3];

        % define partial derivatives
        nu_nu = diag(ones(Nx));
        nu_nU = 0;
        nu_Te = gallery('tridiag', (dt*qa*n_neg)/(2*dx*qe*ma), 0, -(dt*qa*n_pos)/(2*dx*qe*ma));
    
        nU_nu = diag( ((-dt*qa)./(2*dx*qe*n*ma))*((n_pos.*Te0_pos) - (n_neg.*Te0_neg)) );
        nU_nU = diag(ones(Nx));
        nU_Te = gallery('tridiag', -(dt.*qa.*u_para0.*n_pos)./(2*dx*qe*ma), (3*dt*sqrt(2*me).*(n.^2))./(2*ma.^2.*(Te0.^(3/2))), -(dt.*qa.*u_para0.*n_neg)/(2*dx*qe*ma) );
    
        Te_nu = diag(-(dt.*((n_pos.*Te0_pos) - (n_neg.*Te0_neg)))./(3*dx*n));
        Te_nU = 0;
        Te_Te_left = ((dt.*(u_para0.*n0).*n_pos)./(3*dx.*n)) + (((dt*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((-5.*(Te0_neg.^(3/2))/2).*(Te0 - Te0_neg)) + (Te0.^(5/2) + Te0_neg.^(5/2))));
        Te_Te_mid = n - (((dt*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((-5.*(Te0.^(3/2))/2).*(Te0_pos - Te0)) - (Te0_pos.^(5/2) + Te0.^(5/2)) -((5.*(Te0.^(3/2))/2).*(Te0 - Te0_neg)) - (Te0.^(5/2) + Te0_neg.^(5/2)))) - (dt*sqrt(2*me).*(n.^2))/(ma*(Te0.^(3/2)));
        Te_Te_right = ((-dt.*(u_para0.*n0).*n_pos)./(3*dx.*n)) - (((dt*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((-5.*(Te0_pos.^(3/2))/2).*(Te0_pos - Te0)) + (Te0_pos.^(5/2) + Te0.^(5/2))));
        Te_Te = gallery('tridiag', Te_Te_left, Te_Te_mid, Te_Te_right);
    
        P = [nu_nu, nu_nU, nu_Te;
             nU_nu, nU_nU, nU_Te;
             Te_nu, Te_nU, Te_Te]; % holy block matrix

        dy = -P\R; % solve for delta y at time l+1
        y = y + dy;
        
        tol = min(5e-12, 5e-10 * max(R)); % update tolerance to ensure we always are less than starting norm
    end

    nu = y(1);
    nU = y(2);
    T = y(3);
end

function [f_hat] = get_f_hat(f, v_para, x_end)

    f_hat = zeros(size(f));
    [n0, u_para0, T0] = get_boundaries();

    f0 = maxwell(n0(0), u_para0(0), T0(0));
    f_end = maxwell(n0(x_end), u_para0(x_end), T0(x_end));
    f = [f0, f, f_end];

    for i = 1:numel(v_para)
        f_hat(:, :, i) = (v_para(i) > 0)*f(:, :, i) + (v_para(i) <= 0)*f(:, :, i+1);
    end

end

function [n0, u_para0, T_ae0] = get_boundaries()
    n0 = @(x) 0.36*tanh(0.05*(x-100)) + 0.64;
    u_para0 = @(x) -1.1154*tanh(0.05*(x-100)) + 1.983;
    T_ae0 = 0.4424*tanh(0.05*(x-100)) + 0.5576; % Ta = Te at boundary!
end