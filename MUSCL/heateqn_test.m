% 1D heat equation
% u_t = u_xx, 0<x<1
% Homo. Dirichlet BCs, u(0,t)=u(1,t)=0.
% IC: u(x,0) = sin(2*pi*x)
clear all;close all;

Nx = 100; %number of cells
xvals = linspace(0,1,Nx+1)'; %Nx+1 grid points
dx = xvals(2)-xvals(1);

lambdavals = 0.01:0.01:10; %dt = lambda*dx^2
L1err = zeros(numel(lambdavals),1);

for k = 1:numel(lambdavals)
Tf = 0.1;
lambda = lambdavals(k);
%dt = lambda*dx^2; %f.Euler: dt<0.5dx^2
dt = lambda*dx; %b.Euler
tvals = 0:dt:Tf;
if tvals(end)~=Tf
    tvals = [tvals, Tf];
end
Nt = numel(tvals);

% Build differentiation matrix
Dxx = (1/dx^2)*gallery('tridiag',Nx-1,1,-2,1); %only interior nodes

% Initial condition
u = sin(2*pi*xvals); %IC
u = u(2:end-1); %only interior nodes

for n = 2:Nt
    dtn = tvals(n)-tvals(n-1);
    %u = (speye(Nx-1,Nx-1) + dtn*Dxx)*u; %f.Euler
    u = (speye(Nx-1,Nx-1) - dtn*Dxx)\u; %b.Euler
end

uexact = exp(-4*pi^2*Tf)*sin(2*pi*xvals);
uexact = uexact(2:end-1); %only interior nodes
L1err(k) = dx*sum(abs(u-uexact));
end

figure(1);clf;loglog(lambdavals,L1err,'blue','linewidth',1.5);
hold on;loglog(lambdavals,0.01*lambdavals.^1,'black-.','linewidth',1.5);













