% Compute 2D Finite Differences Problems for Conservation Laws
clear variables; close all; clc;

tf = 0.5; % final time
a = 1;
b = 2;
lambda = 0.95; % dt = lambda*dx
interval = [0, 2*pi, 0, 2*pi]; % xmin, xmax, ymin, ymax
discretizationType = 'RK1'; % time discretization type (RK1, RK2, RK3)

alpha = 1;
beta = 2;

%% Test 1 - Linear Advection
f = @(u, x, y, t) (a*u); 
g = @(u, x, y, t) (b*u);
u0 = @(x, y) sin(x + y);
u_exact = @(x, y) sin(x + y - (tf*(a + b)));

OrderTester(discretizationType, lambda, tf, u0, interval, f, g, u_exact, alpha, beta)



%% Test 2 - Variable Coefficient 1
interval = [-pi, pi, -pi, pi];
f = @(u, x, y, t) (-y.*u); 
g = @(u, x, y, t) (x.*u);
u0 = @(x, y) exp(-3*(x.^2 + y.^2));
u_exact = @(x, y) u0(x, y);

OrderTester(discretizationType, lambda, tf, u0, interval, f, g, u_exact, alpha, beta)
return
%% Test 3 - Variable Coefficient 2 (Testing WENO on Discontinuities)
interval = [-2, 2, -2, 2];
tf = 2*pi;
f = @(u, x, y, t) (-y.*u); 
g = @(u, x, y, t) (x.*u);
u0 = @(x, y) ((x>=-1) - (x>=1)).*((y>=-1) - (y>=1)); % 2D step function centered at origin

u = Finite_Differences_2D(discretizationType, 100, 100, lambda, interval, tf, f, g, u0, alpha, alpha);
