% AUTHOR: J. Nakao
% PURPOSE: Compute the physical quantities.
% "A high order multi-dimensional characteristic tracing strategy for the
% Vlasov–Poisson system" by Qiu and Russo, J. Sci. Comput. (2017).
% 
% 
function [EF,mass,L1, L2,energy,entropy] = quantities(f,V,xlength,Nx,dx,dv)
% Inputs:
    % f - dependent variable at time t^n
    % V - velocity
    % Nx - mesh wrt x
    % dx = mesh size wrt x
    % dv - mesh size wrt v
% Outputs:
    % EF - norm of electric field
    % mass - integrate f over x and v
    % energy - integrate f*v^2 over x and v + integrate E^2 over x
    % entropy - integrate f*log(f) over x and v
rho = (dv*sum(f, 2))-1;
E = poisson(rho,Nx,xlength);

%EF = max(abs(E)); %Linf norm
EF = sqrt(dx*sum(real(E).^2)); %L2 norm
mass = dx*dv*sum(sum(f)); %Riemann sum
L1 = dx*dv*sum(sum(abs(f))); % L1 norm
L2 = sqrt(dx*dv*sum(sum((f.^2)))); %L2 norm
energy = dx*dv*sum(sum(f.*V.^2)) + dx*sum(E.^2); %Riemann sum
entropy = dx*dv*sum(sum(-f.*log(f)));

end