clc; clear all; close all;

x = [0, 1, 2, 3]';
y = [0.2, 2, 1.2, 3]';

A = [ones(4, 1), exp(x), x.^2, sin(x)]
coeffs = A\y

xvals = linspace(0, 5, 100);
yvals = coeffs(1) + (coeffs(2)*exp(xvals)) + (coeffs(3)*(xvals.^2)) + (coeffs(4)*sin(xvals));

scatter(x, y); hold on;

plot(xvals, yvals)