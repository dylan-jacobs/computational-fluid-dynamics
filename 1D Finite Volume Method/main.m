% Compute 1D Finite Volume Problems for Conservation Laws
clear variables; close all; clc;

tf = 0.3;
a = 1;
lambda = 0.3; % dt = lambda*dx
discretizationType = 'RK2';

% Test 1 - Burger's Equation
f = @(x, t, u) (a*(u.^2)/2);
u0 = @(x) sin(x);

u_exact_eqn = @(x) (@(y) y-(sin(x-(a*tf*y))));
u_exact = @(x) fsolve(u_exact_eqn(x), x, optimset('Display', 'off'));

OrderTester(discretizationType, lambda, tf, u0, f, u_exact)

% Test 2 - Linear Advection
f = @(x, t, u) (a*u); 
u0 = @(x) sin(x);
u_exact = @(xvals) u0(xvals-(a*tf));

OrderTester(discretizationType, lambda, tf, u0, f, u_exact)

%% Test 3 - Variable Coefficient 1
f = @(x, t, u) a.*sin(x).*u;    
u0 = @(x) 1;
u_exact = @(x) sin(2*atan(exp(a*-tf).*tan(x./2))) ./ sin(x);

OrderTester(discretizationType, lambda, tf, u0, f, u_exact)

%% Test 4 - Variable Coefficient 2
f = @(x, t, u) a*u./(1+t);
u0 = @(x) exp(-5*(x-pi).^2);
u_exact = @(x) exp(-5*(x-a*log(tf+1)-pi).^2);

OrderTester(discretizationType, lambda, tf, u0, f, u_exact)

