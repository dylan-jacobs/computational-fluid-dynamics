clear variables; close all; clc;

% TEST SPATIAL ACCURACY OF MUSCL
Nxvals = [50, 100, 200, 400, 800]'; % Nx + 1 points, Nx intervals
errors = zeros(numel(Nxvals), 2);

tf = 0.5; % final time --> choose BEFORE breaking time
a = -1;
f = @(u) 0.5*(u.^2);
%
f = @(u) a*u;
u0 = @(x) sin(x); % initial condition

for k = 1:size(errors, 1)
    Nx = Nxvals(k);
    xvals = linspace(-pi, pi, Nx+1)';
    dx = xvals(2) - xvals(1);
    xmid = xvals(1:end-1) + dx/2;
    dt = (dx^2);

    % u_exact = U_Exact(xvals, tf); %EXACT SOLN COMPUTED BELOW

    tvals = (0:dt:tf)';
    if tvals(end)~=tf
        tvals = [tvals; tf];
    end
    
    % initialize u as array of cell averages
    u_bar = zeros(numel(Nx), 1);
    for j = 1:Nx
        u_bar(j, 1) = (1/dx)*GLquad(u0,xvals(j),xvals(j+1),7);
    end

    for n = 2:numel(tvals)
        dt = tvals(n) - tvals(n-1);
        u_mid = u_bar;
        u_pos = [u_bar(2:end) ; u_bar(1)];
        u_pos_pos = [u_bar(3:end); u_bar(1); u_bar(2)];
        u_neg = [u_bar(end); u_bar(1:end-1)];

        u1_tilde_upwind = 0.5*(u_mid - u_neg);
        u2_tilde_upwind = 0.5*(u_pos - u_mid);

        u1_tilde_downwind = 0.5*(u_mid - u_pos);
        u2_tilde_downwind = 0.5*(u_neg - u_mid);

        f_pos_upwind = f(u_mid + minmod(u1_tilde_upwind, u2_tilde_upwind));
        f_pos_downwind = f(u_mid + minmod(u1_tilde_downwind, u2_tilde_downwind));

        if (a >= 0)
            f_pos = f_pos_upwind; 
            f_neg = [f_pos(end); f_pos(1:end-1)];
        else
            f_neg = f_pos_downwind;
            f_pos = [f_neg(2:end); f_neg(1)];
        end

        u_bar = u_bar - (dt/(dx))*(f_pos - f_neg);
    end

    u_linadv = @(x1) u0(x1+(-a*tf));
    u_exact = zeros(numel(Nx), 1);
    for j = 1:Nx
        u_exact(j, 1) = (1/dx)*GLquad(u_linadv, xvals(j), xvals(j+1), 7);
    end

    figure(1); clf;
    plot(xvals(1:end-1), u_exact, 'black--'); hold on;
    plot(xvals(1:end-1), u_bar, 'b-', 'LineWidth', 1.5);
    ylim([-1.2, 1.2]);
    title(sprintf('2-Point MUSCL Scheme at t=%s', num2str(tf, 2)));
    xlabel('x'); ylabel('u'); legend('Exact', sprintf('Nx=%d', numel(xvals)-1));


    errors(k, 1) = dx*sum(abs(u_bar - u_exact)); % L1 error
    errors(k, 2) = sqrt(dx*sum((u_bar - u_exact).^2)); % L2 error
end

L1_error = errors(:, 1);
L2_error = errors(:, 2);
L1_order = [0; log2(errors(1:end-1, 1) ./ errors(2:end, 1))];
L2_order = [0; log2(errors(1:end-1, 2) ./ errors(2:end, 2))];

table(L1_error, L1_order, L2_error, L2_order)




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










