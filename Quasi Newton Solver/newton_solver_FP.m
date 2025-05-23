% Given current solution f, solves for the macroscopic parameters of the Fokker-Planck equation 
% mass (n), momentum (nu), kinetic energy (nU), temperature (Te)
% f in [v_para, v_perp, space]

% TESTED ON 4.4.2

function [n, nu, nU, T_e] = newton_solver_FP(f, n0, u_para0, U0, T_e0, dt, dx, dv_para, dv_perp, v_para, v_perp, qa, qe, ma, me, R_const, x_min, x_max)

    Nx = numel(n0); % max val of i, j
    x_ghost_left = x_min - dx/2; % left ghost cell
    x_ghost_right = x_max + dx/2; % right ghost cell

    % get boundary conditions
    [n0_boundary, u_para0_boundary, T_ae_boundary] = get_boundaries(); % Ta = Te at boundary --> T_ae
    n0_min = n0_boundary(x_ghost_left); n0_max = n0_boundary(x_ghost_right);
    u_para0_min = u_para0_boundary(x_ghost_left); u_para0_max = u_para0_boundary(x_ghost_right);
    T_ae_min = T_ae_boundary(x_ghost_left); T_ae_max = T_ae_boundary(x_ghost_right);

    n0_pos = [n0; n0_max]; n0_neg = [n0_min; n0];
    u_para0_half_nodes = (u_para0(1:end-1) + u_para0(2:end))/2; u_para0_half_nodes = [u_para0_min; u_para0_half_nodes; u_para0_max];
    Te0_pos = [T_e0(2:end); T_ae_max]; Te0_neg = [T_ae_min; T_e0(1:end-1)];

    % first, compute fluxes via summation
    f_hat = get_f_hat(f, v_para, v_perp, R_const, x_ghost_left, x_ghost_right);

    % ---- Compute f_hat using f -----
    % f0 = maxwellian(n0_min, v_para, v_perp, u_para0_min, T_ae_min, R_const);
    % f_end = maxwellian(n0_max, v_para, v_perp, u_para0_max, T_ae_max, R_const);
    % f_hat_construction = cat(3, f0, f, f_end);
    % ------------------------

    nu_hat = zeros(Nx+1, 1);
    S_hat = zeros(Nx+1, 1);
    Q_hat = zeros(Nx+1, 1);
    for i = 1:Nx+1
        % f_hat = zeros(size(f_hat_construction, 1), size(f_hat_construction, 2));
        % v_para_split_idx = find(v_para > 0, 1); % upwinding!
        % f_hat(1:v_para_split_idx-1, :) = f_hat_construction(1:v_para_split_idx-1, :, i+1);
        % f_hat(v_para_split_idx:end, :) = f_hat_construction(v_para_split_idx:end, :, i);
  
        nu_hat(i) = 2*pi*(sum(sum(v_para .* f_hat(:, :, i) .* (v_perp.*dv_para.*dv_perp))));
    
        S_hat(i) = 2*pi*(sum(sum(v_para.^2 .* f_hat(:, :, i) .* (v_perp.*dv_para.*dv_perp))));
    
        Q_hat(i) = pi*(sum(sum(v_para .* (v_para.^2 + v_perp.^2) .* f_hat(:, :, i) .* (v_perp.*dv_para.*dv_perp))));
    end
    nTe_hat = ((n0_neg.*[T_ae_min; T_e0]).*(u_para0_half_nodes > 0) + (n0_pos.*[T_e0; T_ae_max]).*(u_para0_half_nodes <= 0)); % upwinding
    % shift bounds to get flux pos/neg
    nu_hat_pos = nu_hat(2:end); nu_hat_neg = nu_hat(1:end-1);

    S_hat_pos = S_hat(2:end); S_hat_neg = S_hat(1:end-1);
    Q_hat_pos = Q_hat(2:end); Q_hat_neg = Q_hat(1:end-1);
    nTe_hat_pos = nTe_hat(2:end); nTe_hat_neg = nTe_hat(1:end-1);

    % explicitly find n via Forward Euler
    n = n0 - (dt/dx)*(nu_hat_pos - nu_hat_neg); % now find n_i+1, n_i-1
    n_pos = [n(2:end); n0_max]; n_neg = [n0_min; n(1:end-1)];

    % ---- init y_vec, R_norm ----  
    y = [n0.*u_para0; n0.*U0; T_e0];

    nu = y(1:Nx); u = nu ./ n;
    nU = y(Nx+1:2*Nx);  U = nU./n;
    T_e = y(2*Nx+1:end);
    Te_pos = [T_e(2:end); T_ae_max]; Te_neg = [T_ae_min; T_e(1:end-1)];
    nTe_pos = n_pos.*Te_pos; nTe_neg = n_neg.*Te_neg;


    kappa_pos = (3.2/(2*sqrt(2*me)))*((Te_pos.^(5/2) + T_e.^(5/2)));
    kappa_neg = (3.2/(2*sqrt(2*me)))*((T_e.^(5/2) + Te_neg.^(5/2))); 
    kappa_half_nodes = (kappa_pos(1:end-1) + kappa_neg(2:end))/2; kappa_half_nodes = [(3.2/(2*sqrt(2*me)))*((T_ae_min.^(5/2) + T_e(1).^(5/2))); kappa_half_nodes; (3.2/(2*sqrt(2*me)))*((T_ae_max.^(5/2) + T_e(end).^(5/2)));];
    kappa_pos = kappa_half_nodes(2:end);
    kappa_neg = kappa_half_nodes(1:end-1);

    R1 = nu - (n0.*u_para0) + (dt/dx).*(S_hat_pos - S_hat_neg) - ((dt*qa)/(2*dx*qe*ma)).*((n_pos.*Te_pos) - (n_neg.*Te_neg));
    R2 = nU - (n0.*U0) + (dt/dx).*(Q_hat_pos - Q_hat_neg) - (((dt*qa*nu)./(2*dx*qe*n)).*((n_pos.*Te_pos) - (n_neg.*Te_neg))) - (((dt.*3.*sqrt(2*me))./((ma.^2).*(T_e.^(3/2)))) .* (((n.^2).*T_e) - ((ma/3).*((2.*n.*nU) - (nu.^2)))));
    R3 = (n0.*T_e) - (n0.*T_e0) + ((5*dt)/(3*dx)).*((u_para0_half_nodes(2:end).*nTe_hat_pos) - (u_para0_half_nodes(1:end-1).*nTe_hat_neg)) - (((dt*u)./(3*dx)).*(nTe_pos - nTe_neg)) - (((2*dt)/(3*dx.^2)) .* ((kappa_pos.*(Te_pos - T_e)) - (kappa_neg.*(T_e - Te_neg)))) - (((dt.*2.*sqrt(2*me))./(ma.*(T_e.^(3/2)))) .* (((ma/3).*((2*n.*nU) - (nu.^2))) - ((n.^2).*T_e)));
    R = [R1; R2; R3];

    
    tol = min(5e-12, max(abs(R))*5e-10); % ensure we don't get worse!
   
    while (max(abs(R)) > tol)
        % define partial derivatives of residual
        nu_nu = spdiags(ones(Nx), 0, Nx, Nx);
        nu_nU = spdiags(zeros(Nx),0, Nx, Nx);
        nu_Te = gallery('tridiag', (dt*qa*n(1:end-1))/(2*dx*qe*ma), zeros(Nx,1), -(dt*qa*n(2:end))/(2*dx*qe*ma));
    
        nU_nu = diag( ((-dt*qa)./(2*dx.*qe.*n)).*((n_pos.*Te_pos) - (n_neg.*Te_neg)) );
        nU_nU = spdiags(ones(Nx), 0, Nx, Nx);
        % nU_Te_mid = (-(3*dt*sqrt(2*me))./(ma.^2)).*( (-(n.^2)./(2.*T_e.^(3/2))) + (ma/2)*((2*n.*y(Nx+1:2*Nx)) - (y(1:Nx).^2))./(T_e.^(5/2)) );
        nU_Te_mid = ((3*dt*sqrt(2*me).*(n.^2))./(2*(ma.^2))).*((1./(T_e.^(3/2))  + (ma./(T_e.^(5/2))).*((2*U) - (u.^2)) ));
        nU_Te = gallery('tridiag', (dt.*qa.*nu(2:end).*(n(1:end-1)./n(2:end)))./(2*dx*qe), nU_Te_mid, -(dt.*qa.*nu(1:end-1).*(n(2:end)./n(1:end-1)))/(2*dx*qe) );
    
        Te_nu = diag(-(dt.*((n_pos.*Te_pos) - (n_neg.*Te_neg)))./(3*dx*n));
        Te_nU = spdiags(zeros(Nx), 0, Nx, Nx);
        Te_Te_left = ((dt.*nu.*n_neg)./(3*dx.*n)) + (((dt*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(Te_neg.^(3/2))/2).*(T_e - Te_neg)) - (T_e.^(5/2) + Te_neg.^(5/2)))); Te_Te_left = Te_Te_left(2:end);
        % Te_Te_mid = n - (((dt*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(T_e.^(3/2))./2).*(Te_pos - T_e)) - (Te_pos.^(5/2) + T_e.^(5/2)) - ((5.*(T_e.^(3/2))/2).*(T_e - Te_neg)) - (T_e.^(5/2) + Te_neg.^(5/2)))) - ((2*dt*sqrt(2*me)./(ma)).*(((-ma./(3*T_e.^(5/2))).*(2*n.*(nU) - (nu.^2)) + ((n.^2)./(2*T_e.^(3/2))))));
        Te_Te_mid = n - (((dt*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(T_e.^(3/2))./2).*(Te_pos - T_e)) - (Te_pos.^(5/2) + T_e.^(5/2)) - ((5.*(T_e.^(3/2))/2).*(T_e - Te_neg)) - (T_e.^(5/2) + Te_neg.^(5/2)))) - ((2*dt*sqrt(2*me).*(n.^2)./(2*ma)).*(((-ma./(T_e.^(5/2))).*(2*U - (u.^2)) + (1./(T_e.^(3/2))))));
        Te_Te_right = ((-dt.*nu.*n_pos)./(3*dx.*n)) - (((dt*3.2)./(3*(dx.^2).*sqrt(2*me))).*(((5.*(Te_pos.^(3/2))/2).*(Te_pos - T_e)) + (Te_pos.^(5/2) + T_e.^(5/2)))); Te_Te_right = Te_Te_right(1:end-1);
        Te_Te = gallery('tridiag', Te_Te_left, Te_Te_mid, Te_Te_right);

        P = [nu_nu, nu_nU, nu_Te;
             nU_nu, nU_nU, nU_Te;
             Te_nu, Te_nU, Te_Te]; % holy block matrix
     
        dy = -P\R; % solve for delta y at time l+1
        y = y + dy;

        % update y_vec stuff
        nu = y(1:Nx); u = nu ./ n;
        nU = y(Nx+1:2*Nx); U = nU ./ n;
        T_e = y(2*Nx+1:end);
        Te_pos = [T_e(2:end); T_ae_max]; Te_neg = [T_ae_min; T_e(1:end-1)];
        nTe_pos = n_pos.*Te_pos; nTe_neg = n_neg.*Te_neg;

        kappa_pos = (3.2/(2*sqrt(2*me)))*((Te_pos.^(5/2) + T_e.^(5/2)));
        kappa_neg = (3.2/(2*sqrt(2*me)))*((T_e.^(5/2) + Te_neg.^(5/2)));
        kappa_half_nodes = (kappa_pos(1:end-1) + kappa_neg(2:end))/2; kappa_half_nodes = [(3.2/(2*sqrt(2*me)))*((T_ae_min.^(5/2) + T_e(1).^(5/2))); kappa_half_nodes; (3.2/(2*sqrt(2*me)))*((T_ae_max.^(5/2) + T_e(end).^(5/2)))];
        kappa_pos = kappa_half_nodes(2:end);
        kappa_neg = kappa_half_nodes(1:end-1);
        
        R1 = nu - (n0.*u_para0) + (dt/dx).*(S_hat_pos - S_hat_neg) - ((dt*qa)/(2*dx*qe*ma)).*((n_pos.*Te_pos) - (n_neg.*Te_neg));
        R2 = nU - (n0.*U0) + (dt/dx).*(Q_hat_pos - Q_hat_neg) - (((dt*qa*nu)./(2*dx*qe*n)).*((n_pos.*Te_pos) - (n_neg.*Te_neg))) - (((dt.*3.*sqrt(2*me))./((ma.^2).*(T_e.^(3/2)))) .* (((n.^2).*T_e) - ((ma/3).*((2.*n.*nU) - (nu.^2)))));
        R3 = (n.*T_e) - (n0.*T_e0) + ((5*dt)/(3*dx)).*((u_para0_half_nodes(2:end).*nTe_hat_pos) - (u_para0_half_nodes(1:end-1).*nTe_hat_neg)) - (((dt*u)./(3*dx)).*(nTe_pos - nTe_neg)) - (((2*dt)/(3*dx.^2)) .* ((kappa_pos.*(Te_pos - T_e)) - (kappa_neg.*(T_e - Te_neg)))) - (((dt.*2.*sqrt(2*me))./(ma.*(T_e.^(3/2)))) .* (((ma/3).*((2*n.*nU) - (nu.^2))) - ((n.^2).*T_e)));
        R = [R1; R2; R3];
    end
  
    % parameters to return
    nu = (y(1:Nx));
    nU = y(Nx+1:2*Nx);
    T_e = y(2*Nx+1:end);

    % should delete: attempt to force all macros to BCs artificially
    % nu(1) = n0_min*u_para0_min;
    % nu(end) = n0_max*u_para0_max;
    % U_min = 0.5*((3*T_ae_min/ma) + (u_para0_min.^2));
    % U_max = 0.5*((3*T_ae_max/ma) + (u_para0_max.^2));
    % nU(1) = n0_min*U_min;
    % nU(end) = n0_max*U_max;
    % T_e(1) = T_ae_min;
    % T_e(end) = T_ae_max;
end

function [f_hat] = get_f_hat(f, v_para, v_perp, R_const, x_min, x_max)

    [n0, u_para0, T0] = get_boundaries();

    f0 = maxwellian(n0(x_min), v_para, v_perp, u_para0(x_min), T0(x_min), R_const);
    f_end = maxwellian(n0(x_max), v_para, v_perp, u_para0(x_max), T0(x_max), R_const);
    
    f = cat(3, f0, f, f_end);
    f_size = size(f);
    f_hat = zeros(f_size(1), f_size(2), f_size(3)-1);
    v_para_split_idx = find(v_para(1, :) > 0, 1); % upwinding!
    for i = 1:size(f, 3)-1 % loop through spatial nodes
        f_hat(:, 1:v_para_split_idx-1, i) = f(:, 1:v_para_split_idx-1, i+1);
        f_hat(:, v_para_split_idx:end, i) = f(:, v_para_split_idx:end, i);
    end
end
