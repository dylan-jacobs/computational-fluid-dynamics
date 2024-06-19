% Linear reconstruction inverse solver
clear variables; close all; clc;
syms A B C D E uneg2 uneg1 u u1 u2

Amtx3 = [49/12, -2, 1;
     13/12, -1, 1;
     1/12,   0, 1];
Amtx3 = [49/12, 2, 1;
     13/12, 1, 1;
     1/12,   0, 1];

uvec3 = [u2; u1; u]; 

soln3 = linsolve(Amtx3, uvec3)

Amtx5 = [1441/80, -17/2, 49/12, -2, 1;
          121/80,  -5/4, 13/12, -1, 1;
            1/80,     0,  1/12,  0, 1;
          121/80,   5/4, 13/12,  1, 1;
         1441/80,  17/2, 49/12,  2, 1];

uvec5 = [uneg2; uneg1; u; u1; u2]; 

soln5 = linsolve(Amtx5, uvec5)

Amtx2 = [1, 1;
         2, 1];
Amtx2 = [0, 1;
         1, 1];
uvec2 = [u1; u2];
soln2 = linsolve(Amtx2, uvec2)


% Ax = u2 - u1
% b = 2u1 - u2

% (u2 - u1)(x - xj) + (2u1 - u2) = 0.5*(u_pos_pos - u_pos) + 2u_pos - u_pos_pos = 0.5*(u_pos - u_pos_pos)

% (u2 - u1)(x - xj) + u1 = 0.5*(u_mid - u_pos)














