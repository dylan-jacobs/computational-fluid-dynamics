% WENO
% 3 3-point stencils (k=3)
% 5th Order accuracy

clear variables; close all; clc;

leftright = 0; % 1 = right, 0 = left;

f = @(x) sin(x);
% f = @(x) (x>=-1) - (x>=1);
Nxvals = [40, 80, 160, 320, 640]'; % Nx + 1 points, Nx cells

epsilon = 1e-6;
c = [11/6, -7/6, 1/3;
     1/3, 5/6, -1/6;
    -1/6, 5/6, 1/3;
     1/3, -7/6, 11/6];

c = c(1+leftright:end-(1-leftright), :); % shift nodes left by 1 if approximating left boundaries
d = [0.3*leftright + 0.1*(1-leftright), 0.6, 0.1*leftright + 0.3*(1-leftright)]; % do same shift here

errors = zeros(numel(Nxvals), 2);

for k = 1:size(errors, 1)
    Nx = Nxvals(k);
    xvals = linspace(-pi, pi, Nx + 1)'; 
    dx = xvals(2) - xvals(1);

    u_exact = f(xvals(1+leftright:end-(1-leftright))); % exact point values at left/right boundaries
    u = zeros(Nx, 1); % initialize u as "exact" cell averages
    for j = 1:Nx
        u(j, 1) = (1/dx)*GLquad(f,xvals(j),xvals(j+1),7);
    end

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
    u_approx = sum(omega.*u_mtx, 2);

    errors(k, 1) = dx*sum(abs(u_approx - u_exact)); % L1 error
    errors(k, 2) = sqrt(dx*sum((u_approx - u_exact).^2)); % L2 error
end

figure(1); clf;
plot(xvals(2:end), u_approx, '-b', 'LineWidth', 1.5); hold on;
plot(xvals(2:end), u_exact, '--black');
legend(sprintf('WENO without weighting Nx = %d', Nx), 'Exact');
xlabel('x'); ylabel('u'); title('WENO No Weights Solution Curve');

L1_error = errors(:, 1);
L2_error = errors(:, 2);
L1_order = [0; log2(errors(1:end-1, 1) ./ errors(2:end, 1))];
L2_order = [0; log2(errors(1:end-1, 2) ./ errors(2:end, 2))];

table(L1_error, L1_order, L2_error, L2_order)


%% FUNCTIONS
function output = GLquad(f,a,b,M)
    output = 0;
    [nodes,weights] = GLnodesweights(M);
    for i = 1:M
        output = output + weights(i)*f(((b-a)/2)*nodes(i)+(b+a)/2);
    end
    output = output*(b-a)/2;
end

function [nodes,weights] = GLnodesweights(M)
if M == 2
    nodes = [-1/sqrt(3), 1/sqrt(3)];
    weights = [1, 1];
elseif M == 3
    nodes = [-sqrt(3/5), 0, sqrt(3/5)];
    weights = [5/9, 8/9, 5/9];
elseif M == 4
    nodes = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053];
    weights = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454];
elseif M == 5
    nodes = [-0.906179845938664, -0.538469310105683, 0, 0.538469310105683, 0.906179845938664];
    weights = [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189];
elseif M == 6
    nodes = [-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152];
    weights = [0.171324492379170, 0.360761573048139, 0.4679139345726910, 0.4679139345726910, 0.360761573048139, 0.171324492379170];
elseif M == 7
    nodes = [-0.949107912342759, -0.741531185599394, -0.405845151377397, 0, 0.405845151377397, 0.741531185599394, 0.949107912342759];
    weights = [0.129484966168870, 0.279705391489277, 0.3818300505051189, 0.4179591836734694, 0.3818300505051189, 0.279705391489277, 0.129484966168870];
end
end




