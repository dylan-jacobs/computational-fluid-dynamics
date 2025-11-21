function [f_inf, n_inf, uz_inf, T_inf] = QCM(rho0, Jz0, kappa0, R, Rmat, Zmat)

    dr = Rmat(2, 1) - Rmat(1, 1);
    dz = Zmat(1, 2) - Zmat(2, 1);

    tol = 1e-15;

    Mk = [pi^(3/2); 0; 3]; % init guess for moments at equilibrium
    n_k = Mk(1);
    u_para_k = Mk(2);
    T_k = Mk(3);

    Exp = exp(-((Zmat-u_para_k).^2 + Rmat.^2)/(2*R*T_k));
    Rk_n = rho0 - sum(sum((n_k/((2*pi*R*T_k)^(3/2)))*Exp.*Rmat))*2*pi*dr*dz;
    Rk_u = Jz0  - sum(sum((n_k/((2*pi*R*T_k)^(3/2)))*Exp.*Rmat.*Zmat))*2*pi*dr*dz;
    Rk_T = kappa0   - sum(sum((n_k/((2*pi*R*T_k)^(3/2)))*Exp.*Rmat.*((Rmat.^2 + Zmat.^2)/2)))*2*pi*dr*dz;

    Rk = [Rk_n; Rk_u; Rk_T];  

    iters = 0;

    while norm(Rk, 2) > tol
        Exp = exp(-((Zmat-u_para_k).^2 + Rmat.^2)/(2*R*T_k));
        Jk_nn = -sum(sum((1/((2*pi*R*T_k)^(3/2)))*Exp.*Rmat))*2*pi*dr*dz;
        Jk_nu = -sum(sum((n_k/((2*pi*R*T_k)^(3/2)))*Exp.*((Zmat - u_para_k)/(R*T_k)).*Rmat))*2*pi*dr*dz;
        Jk_nT = -sum(sum((n_k/((2*pi*R)^(3/2)))*Exp.*((-1.5*T_k^(-5/2)) + (((Rmat.^2) + (Zmat - u_para_k).^2)/(2*R*T_k^(7/2)))).*Rmat))*2*pi*dr*dz;
    
        Jk_un = -sum(sum((1/((2*pi*R*T_k)^(3/2)))*Exp.*Rmat.*Zmat))*2*pi*dr*dz;
        Jk_uu = -sum(sum((n_k/((2*pi*R*T_k)^(3/2)))*Exp.*((Zmat - u_para_k)/(R*T_k)).*Rmat.*Zmat))*2*pi*dr*dz;
        Jk_uT = -sum(sum((n_k/((2*pi*R)^(3/2)))*Exp.*((-1.5*T_k^(-5/2)) + (((Rmat.^2) + (Zmat - u_para_k).^2)/(2*R*T_k^(7/2)))).*Rmat.*Zmat))*2*pi*dr*dz;
        
        Jk_Tn = -sum(sum((1/((2*pi*R.*T_k)^(3/2)))*Exp.*Rmat.*(Rmat.^2 + Zmat.^2)/2))*2*pi*dr*dz;
        Jk_Tu = -sum(sum((n_k/((2*pi*R.*T_k)^(3/2)))*Exp.*((Zmat - u_para_k)/(R*T_k)).*Rmat.*(Rmat.^2 + Zmat.^2)/2))*2*pi*dr*dz;
        Jk_TT = -sum(sum((n_k/((2*pi*R)^(3/2)))*Exp.*((-1.5*T_k^(-5/2)) + (((Rmat.^2) + (Zmat - u_para_k).^2)/(2*R*T_k^(7/2)))).*Rmat.*(Rmat.^2 + Zmat.^2)/2))*2*pi*dr*dz;

        J = [Jk_nn, Jk_nu, Jk_nT;
             Jk_un, Jk_uu, Jk_uT;
             Jk_Tn, Jk_Tu, Jk_TT];

        dMk = -J\Rk;
        Mk = Mk + dMk;

        n_k = Mk(1);
        u_para_k = Mk(2); 
        T_k = Mk(3);

        Exp = exp(-((Zmat-u_para_k).^2 + Rmat.^2)/(2*R*T_k));
        Rk_n = rho0 - sum(sum((n_k/((2*pi*R*T_k)^(3/2)))*Exp.*Rmat))*2*pi*dr*dz;
        Rk_u = Jz0  - sum(sum((n_k/((2*pi*R*T_k)^(3/2)))*Exp.*Rmat.*Zmat))*2*pi*dr*dz;
        Rk_T = kappa0   - sum(sum((n_k/((2*pi*R*T_k)^(3/2)))*Exp.*Rmat.*((Rmat.^2 + Zmat.^2)/2)))*2*pi*dr*dz;

        Rk = [Rk_n; Rk_u; Rk_T];

        % disp('Error: ')
        % disp(norm(Rk, 2));
        iters = iters + 1;
        if iters > 100
            break
        end
    end

    f_inf = (n_k/((2*pi*R*T_k)^(3/2)))*exp(-((Zmat-u_para_k).^2 + Rmat.^2)/(2*R*T_k));
    n_inf = n_k;
    uz_inf = u_para_k;
    T_inf = T_k;

end