
function [Vx, Vy] = reduced_augmentation(Vx_aug, Vy_aug)
    tolerance = 1e-12;
    [Qx, Rx] = qr(Vx_aug, 0); [Qy, Ry] = qr(Vy_aug, 0);
    [Ux, SigmaX, ~] = svd(Rx, 0); 
    [Uy, SigmaY, ~] = svd(Ry, 0); 
    rx = find(diag(SigmaX) > tolerance, 1, 'last');
    ry = find(diag(SigmaY) > tolerance, 1, 'last');
    R = max(rx, ry);
    Vx = Qx*Ux(:, 1:R); Vy = Qy*Uy(:, 1:R);
end


