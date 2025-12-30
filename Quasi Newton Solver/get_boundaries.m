function [n0, u_para0, T_ae0] = get_boundaries()
    n0 = @(x) 0.36*tanh(0.05*(x-100)) + 0.64;
    u_para0 = @(x) -1.1154*tanh(0.05*(x-100)) + 1.983;
    T_ae0 = @(x) 0.4424*tanh(0.05*(x-100)) + 0.5576; % Ta = Te at boundary
end