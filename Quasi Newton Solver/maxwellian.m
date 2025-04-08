function [f] = maxwellian(n, v_para, v_perp, u_para, T, R)
    Nx = numel(n);
    Nv = size(v_para);
    f = zeros(Nv(1), Nv(2), Nx);
    
    for i = 1:Nx
        f(:, :, i) = (n(i)./((2*pi*R.*T(i)).^(3/2))).*exp(-((v_para-u_para(i)).^2 + v_perp.^2)./(2*R*T(i)));
    end
end