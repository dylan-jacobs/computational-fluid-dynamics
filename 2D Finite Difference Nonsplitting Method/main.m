% Compute 2D Finite Differences Problems for Conservation Laws
clear variables; close all; clc;

tf = 0.5; % final time
lambda = 0.95; % dt = lambda*CFL condition
interval = [0, 2*pi, 0, 2*pi]; % xmin, xmax, ymin, ymax
discretizationType = 'RK4'; % time discretization type (RK1, RK2, RK3)


%% Test 1 - Linear Advection
a = 1;
b = 2;
f = @(u, x, y, t) (a*u); 
g = @(u, x, y, t) (b*u);
u0 = @(x, y) sin(x + y);
u_exact = @(x, y) sin(x + y - (tf*(a + b)));
alpha = 1; beta = 2;

OrderTester(discretizationType, lambda, tf, u0, interval, f, g, u_exact, alpha, beta)

%% Test 2 - Rigidbody Rotation
interval = [-pi, pi, -pi, pi];
f = @(u, x, y, t) (-y.*u); 
g = @(u, x, y, t) (x.*u);
u0 = @(x, y) exp(-3*(x.^2 + y.^2));
u_exact = @(x, y) u0(x, y);
alpha = pi; beta = pi;

OrderTester(discretizationType, lambda, tf, u0, interval, f, g, u_exact, alpha, beta)


% %% Test 3 - WENO on Discontinuous 2D Step Function)
% interval = [-pi, pi, -pi, pi];
% tf = 2*pi;
% f = @(u, x, y, t) (-y.*u); 
% g = @(u, x, y, t) (x.*u);
% u0 = @(x, y) ((x>=-1) - (x>=1)).*((y>=-1) - (y>=1)); % 2D step function centered at origin
% alpha = 1; beta = 1;
% 
% u = Finite_Differences_2D(discretizationType, 100, 100, lambda, interval, tf, f, g, u0, alpha, beta);

%% Test 4 - Swirling Deformation 1
interval = [-pi, pi, -pi, pi];
tf = 1.5;
gt = @(t) cos(pi.*t./tf).*pi;
f = @(u, x, y, t) -(cos(x./2).^2).*(sin(y).*gt(t).*u);
g = @(u, x, y, t) (sin(x).*(cos(y./2).^2).*gt(t).*u);

u0 = @(x, y) cos_bell(x, y);
u_exact = @(x, y) u0(x, y);
alpha = pi; beta = pi;

OrderTester(discretizationType, lambda, tf, u0, interval, f, g, u_exact, alpha, beta)

%% Test 5 - Swirling Deformation 2: WENO Discontinuity Test
interval = [-pi, pi, -pi, pi];
tf = 5*pi;
gt = @(t) 1;
f = @(u, x, y, t) -(cos(x./2).^2).*(sin(y).*gt(t).*u);
g = @(u, x, y, t) (sin(x).*(cos(y./2).^2).*gt(t).*u);

u0 = @(x, y) swirl(x, y);
u_exact = @(x, y) u0(x, y);
alpha = 1; beta = 1;

u = Finite_Differences_2D_Nonsplitting(discretizationType, 100, 100, lambda, interval, tf, f, g, u0, alpha, beta);



function [output] = cos_bell(x, y)
    rb0 = 0.3*pi;
    rb = sqrt((x-(0.3*pi)).^2 + ((y).^2));
    output = (rb < rb0).*(rb.*(cos(rb.*pi./(2.*rb0))).^6);
end

function [output] = swirl(x, y)    
    rb0 = 8*pi/5;
    rb = sqrt((x-(pi)).^2 + ((y-pi).^2));
    output = (rb < rb0).*1;
end






























