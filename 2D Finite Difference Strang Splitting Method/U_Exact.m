% compute u_exact cell averages

function [u_exact] = U_Exact(xvals, tf)
    u_exact = zeros(numel(xvals)-1, 1);
    dx = xvals(2) - xvals(1);    
    for i = 1:numel(xvals)-1
        sub_a = xvals(1) + (dx * (i-1));
        sub_b = xvals(1) + (dx * (i));
        [weights, nodes] = GetWeightsAndNodes(7);
        sub_approx = GaussLegendreCustom(sub_a, sub_b, weights, nodes, tf);
        u_exact(i) = sub_approx/dx;
    end
end

function [approximation] = GaussLegendreCustom(a, b, weights, nodes, tf)
    % interval [a, b]
    % N intervals, N+1 nodes
    % f = function

    approximation = 0;
    dXi = (b-a) / 2;
    for i = 1:numel(nodes)
        node = nodes(i);
        Xi = node * (b-a) / 2;

        % For Burgers' Equation
        f = @(y) y-(sin((Xi + (b+a)/2)-(tf*y)));
        u_exact = fsolve(f, Xi + (b+a)/2, optimset('Display', 'off'));

        % For Linear Advection
        % u_exact = sin((Xi + (b+a)/2)-tf);
        approximation = approximation + (weights(i) * u_exact);
    end
    approximation = dXi * approximation;
end






