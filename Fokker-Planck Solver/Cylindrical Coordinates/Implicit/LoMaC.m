%% ------- LoMaC Truncation -------

function [Vr, S, Vz, rank] = LoMaC(Vr, S, Vz, R, Z, rvals, zvals, tolerance, f_inf)
% LoMaC Truncates given maxwellian (assumed Low-Rank) to given tolerance while conserving
% macroscopic quantities.

    Nr = numel(rvals); Nz = numel(zvals);
    dvr = rvals(2) - rvals(1);
    dvz = zvals(2) - zvals(1);

    rhoM = 2*pi*dvr*dvz*sum(sum(f_inf.*R));
    JzM = 2*pi*dvr*dvz*sum(sum(Z.*f_inf.*R));
    kappaM = pi*dvr*dvz*sum(sum((R.^2+Z.^2).*f_inf.*R));

    % Step 1: Integrate to calculate macro quantities
    p = 2*pi*sum(sum(((Vr) * S * (Vz)') .* (R .* dvr .* dvz)));
    J = 2*pi*sum(sum(((Vr) * S * (Vz.*zvals)') .* (R .* dvr .* dvz)));
    k = pi*sum(sum((((Vr.*(rvals.^2)) * S * Vz') + (Vr* S * ((Vz.*(zvals.^2))'))) .* (R .* dvr .* dvz)));

    % Step 2: Scale by maxwellian to ensure inner product is well defined
    % (f -> 0 as v -> infinity)
    wr = exp(-(rvals.^2));
    wz = exp(-(zvals.^2));

    % Step 3: Orthogonal projection
    % bases: 1, v, v.^2 - c
    c = (dvr*sum(rvals.^2.*wr.*rvals))/(dvr*sum(wr.*rvals)) + (dvz*sum(zvals.^2.*wz))/(dvz*sum(wz));

    w_norm_1_squared = (2*pi*dvr*dvz*sum(rvals .* wr)*sum(wz));
    w_norm_v_squared = (2*pi*dvr*dvz*sum(rvals .* wr)*sum(zvals.^2 .* wz));
    w_norm_v2_squared = (2*pi*dvr*dvz*sum(sum((R.^2 + Z.^2 - c).^2 .* exp(-R.^2 - Z.^2) .* R)));
    
    f1_proj_S_mtx11 = (p ./ w_norm_1_squared) - ((2*k - c.*p).*c ./ w_norm_v2_squared);
    f1_proj_S_mtx12 = (J ./ w_norm_v_squared);
    f1_proj_S_mtx13 = ((2*k - c.*p) ./ w_norm_v2_squared);

    proj_basis_r = wr.*[ones(Nr, 1), rvals.^2];
    proj_basis_z = wz.*[ones(Nz, 1), zvals, zvals.^2];
    f1_proj_S_mtx   = [f1_proj_S_mtx11, f1_proj_S_mtx12, f1_proj_S_mtx13;
                       f1_proj_S_mtx13,               0,               0];

    % f2 = f - f1 (do it via SVD)
    f2_U = [Vr, proj_basis_r];
    f2_S = blkdiag(S, -f1_proj_S_mtx);
    f2_V = [Vz, proj_basis_z];

    % QR factorize
    [f2_Vr, f2_S, f2_Vz, ~] = truncate(f2_U, f2_S, f2_V, rvals, tolerance);

    f2 = f2_Vr*f2_S*f2_Vz';

    % compute Pn(Te(f)) to ensure moments are kept
    trun_f2_proj_S_mtx11 = 2*pi*dvr*dvz*((sum(sum(f2.*R)) ./ w_norm_1_squared) - (c*sum(sum(f2.*(R.^2 + Z.^2 - c) .* R)) ./ w_norm_v2_squared));
    trun_f2_proj_S_mtx12 = 2*pi*dvr*dvz*((sum(sum(f2.*R.*Z))) ./ w_norm_v_squared);
    trun_f2_proj_S_mtx13 = 2*pi*dvr*dvz*((sum(sum(f2.*(R.^2 + Z.^2 - c) .* R)) ./ w_norm_v2_squared));

    trun_f2_proj_S_mtx   = [trun_f2_proj_S_mtx11, trun_f2_proj_S_mtx12, trun_f2_proj_S_mtx13;
                    trun_f2_proj_S_mtx13,            0,            0];

    % compute fM
    fM_proj_S_mtx11 = (rhoM ./ w_norm_1_squared) - ((2*kappaM - c.*rhoM).*c ./ w_norm_v2_squared);
    fM_proj_S_mtx12 = (JzM ./ w_norm_v_squared);
    fM_proj_S_mtx13 = ((2*kappaM - c.*rhoM) ./ w_norm_v2_squared);

    fM_proj_S_mtx   = [fM_proj_S_mtx11, fM_proj_S_mtx12, fM_proj_S_mtx13;
                       fM_proj_S_mtx13,               0,              0];

    f_mass_S = fM_proj_S_mtx - trun_f2_proj_S_mtx;

    [Vr, S, Vz, rank] = truncate([proj_basis_r, f2_Vr], blkdiag(f_mass_S, f2_S), [proj_basis_z, f2_Vz], rvals, 1e-14);
end

function [Vr, S, Vz, rank] = truncate(Vr_aug, S_aug, Vz_aug, rvals, tolerance)
    [Qr, Rr] = qr2(Vr_aug, rvals); [Qz, Rz] = qr(Vz_aug, 0);
    [U, Sigma, V] = svd(Rr*S_aug*(Rz'), 0); 
    rank = find(diag(Sigma) > tolerance, 1, 'last');
    Vr = Qr*U(:, 1:rank);
    S = Sigma(1:rank, 1:rank);
    Vz = Qz*V(:, 1:rank);
end