
function [Vr, Vz] = reduced_augmentation(Vr_aug, Vz_aug, rvals)
    tolerance = 1e-14;
    [Qr, Rr] = qr2(Vr_aug, rvals); [Qz, Rz] = qr(Vz_aug, 0);
    [Ur, SigmaR, ~] = svd(Rr, 0); 
    [Uz, SigmaZ, ~] = svd(Rz, 0); 
    rr = find(diag(SigmaR) > tolerance, 1, 'last');
    rz = find(diag(SigmaZ) > tolerance, 1, 'last');
    r = max(rr, rz);
    r = min(r, min(size(Qz, 2), size(Qr, 2)));
    Vr = Qr*Ur(:, 1:r);
    Vz = Qz*Uz(:, 1:r);
end


