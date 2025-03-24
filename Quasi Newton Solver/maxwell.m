function [f] = maxwell(n, v_para, v_perp, u_para, T, R)
    f = (n./(2*pi*R.*T).^(3/2)).*exp(-((v_para-u_para).^2 + v_perp.^2)./(2*R*T));
end