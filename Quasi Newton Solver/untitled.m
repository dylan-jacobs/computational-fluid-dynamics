clc;
syms y(x)

Dy = diff(y);

ode = 66.5*diff(y, x, 2) == -557.7*diff(y, x, 1) - (40000*y)+1000;
% ode = diff(y, x, 2) == -4*diff(y, x, 1) - (8*y)+1;

c1 = y(0) == 0;
c2 = Dy(0) == 0;

ySol(x) = dsolve(ode, [c1, c2]);
ySol = simplify(ySol);

xvals = linspace(0, 1.5, 1000);

plot(xvals, ySol(xvals))




