% reconstructs function f using WENO-5, 

function [u_boundaries] = WENO(u, leftright) 
    epsilon = 1e-6;
    c = [11/6, -7/6, 1/3;
         1/3, 5/6, -1/6;
        -1/6, 5/6, 1/3;
         1/3, -7/6, 11/6];
    
    c = c(1+leftright:end-(1-leftright), :); % shift c nodes left by 1 if approximating left boundaries
    d = [0.3*leftright + 0.1*(1-leftright), 0.6, 0.1*leftright + 0.3*(1-leftright)]; % do same shift here

    upos = [u(2:end); u(1)];
    upospos = [upos(2:end); upos(1)];
    uneg = [u(end); u(1:end-1)];
    unegneg = [uneg(end); uneg(1:end-1)];

    % smoothness indicator (lower -> smoother)
    beta = [(13/12)*(u - 2*upos + upospos).^2 + (1/4)*(3*u - 4*upos + upospos).^2, (13/12)*(uneg - 2*u + upos).^2 + (1/4)*(uneg - upos).^2, (13/12)*(unegneg - 2*uneg + u).^2 + (1/4)*(unegneg - 4*uneg + 3*u).^2];

    % weights depend on smoothness
    alpha = d ./ ((epsilon + beta).^2);
    omega = alpha ./ (sum(alpha, 2)); % should sum to 1

    u_mtx1 = [u, upos, upospos];
    u_mtx2 = [uneg, u, upos];
    u_mtx3 = [unegneg, uneg, u];

    u1 = sum(c(1, :).*u_mtx1, 2);
    u2 = sum(c(2, :).*u_mtx2, 2);
    u3 = sum(c(3, :).*u_mtx3, 2);

    u_mtx = [u1, u2, u3];

    % omega = repmat(d, [Nx, 1]); % turn on to see oscillations
    u_boundaries = sum(omega.*u_mtx, 2);
end
