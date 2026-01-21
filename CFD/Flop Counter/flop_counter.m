% Verify flop count for SVD and QR factorization

% %% How to use SVD, QR
% A = rand(10, 4);
% 
% % Full QR
% [Q, R] = qr(A);
% size(Q)
% 
% % Reduced QR
% [Q, R] = qr(A, 0); % economy QR
% size(Q)
% 
% % Full SVD
% [U, S, V] = svd(A);
% size(S)
% 
% % Reduced SVD
% [U, S, V] = svd(A, 0);
% size(S)
% 
% 
% norm(A - U*S*V', 2)


%% Test flop count using tic-toc
clear variables; close all; clc; 

Nvals = [40, 80, 160, 320, 640, 1280, 2560]';
Nt = numel(Nvals);
tvals = zeros(Nt, 5); % 4 tests: SVD, QR, QR (skinny), matrix multiplication associate property speed reduction (2 cols)

for k = 1:Nt
    N = Nvals(k);
    A = rand(N, N);
    B = rand(N, N/10); % skinny

    % SVD
    tic
    [U, S, V] = svd(A, 0);
    tvals(k, 1) = toc;

    % QR
    tic
    [Q, R] = qr(A, 0);
    tvals(k, 2) = toc;

    % Skinny QR 
    tic
    [Q, R] = qr(B, 0);
    tvals(k, 3) = toc;
  
    % Matrix Multiplication Associate Property Test
    D = A; E = rand(N, N); F = B;
    tic
    (D*E)*(F);
    tvals(k, 4) = toc;

    tic
    (D)*(E*F);
    tvals(k, 5) = toc;
end

svd = tvals(:, 1);
qr = tvals(:, 2);
qrSkinny = tvals(:, 3);
matrixMult1 = tvals(:, 4);
matrixMult2 = tvals(:, 5);
svdOrder = [0; log2(svd(2:end) ./ svd(1:end-1))];
qrOrder = [0; log2(qr(2:end) ./ qr(1:end-1))];
qrSkinnyOrder = [0; log2(qrSkinny(2:end) ./ qrSkinny(1:end-1))];
m1Order = [0; log2(matrixMult1(2:end) ./ matrixMult1(1:end-1))];
m2Order = [0; log2(matrixMult2(2:end) ./ matrixMult2(1:end-1))];

output_table = table(svd, svdOrder, qr, qrOrder, qrSkinny, qrSkinnyOrder, matrixMult1, m1Order, matrixMult2, m2Order)




















