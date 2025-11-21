function [Vr_nn,S_nn,Vz_nn,r_nn] = LoMaC_trun(Vr_hat,S_hat,Vz_hat,dr,dz,rvals,zvals,Rmat,Zmat,Nr,Nz,w1,w2,tol,c1,c2,c3,c,rhoM,JzM,kappaM)
fhat = Vr_hat*S_hat*Vz_hat';
rhoH = 2*pi*dr*dz*sum(sum(fhat.*Rmat));
JzH = 2*pi*dr*dz*sum(sum(Zmat.*fhat.*Rmat));
kappaH = 0.5*2*pi*dr*dz*sum(sum((Rmat.^2+Zmat.^2).*fhat.*Rmat));
% Compute f1
Vr_f1 = w1.*[ones(Nr,1),rvals.^2];
Vz_f1 = w2.*[ones(Nz,1),zvals,zvals.^2];
S_f1 = [rhoH/c1 - ((2*kappaH-c*rhoH)/c3)*c, JzH/c2, (2*kappaH-c*rhoH)/c3;...
        (2*kappaH-c*rhoH)/c3, 0, 0];
% Compute and truncate f2=f-f1
[Qr,Rr] = qr2([Vr_hat,Vr_f1],rvals);
[Qz,Rz] = qr([Vz_hat,Vz_f1],0);
[U,S,V] = svd(Rr*blkdiag(S_hat,-S_f1)*Rz',0);
r_f2 = find(diag(S)>tol,1,'last');
Vr_f2 = Qr*U(:,1:r_f2);
Vz_f2 = Qz*V(:,1:r_f2);
S_f2 = S(1:r_f2,1:r_f2);
% Compute PN(trun(f2)), to subtract later to ensure zero moments.
f2 = Vr_f2*S_f2*Vz_f2';
S_f2P = 2*pi*dr*dz*[sum(sum(f2.*Rmat))/c1 - c*sum(sum(f2.*(Rmat.^2+Zmat.^2-c).*Rmat))/c3, sum(sum(f2.*Zmat.*Rmat))/c2, sum(sum(f2.*(Rmat.^2+Zmat.^2-c).*Rmat))/c3;...
        sum(sum(f2.*(Rmat.^2+Zmat.^2-c).*Rmat))/c3, 0, 0];
% Compute fM
% Vr_fM = Vr_f1;
% Vz_fM = Vz_f1;
S_fM = [rhoM/c1 - ((2*kappaM-c*rhoM)/c3)*c, JzM/c2, (2*kappaM-c*rhoM)/c3;...
        (2*kappaM-c*rhoM)/c3, 0, 0];
% Compute fM-PN(trun(f2)).
% Since they share the same O.N. basis, only add(subtract) core matrices.
% Redefine fM = fM-PN(trun(f2))
S_fM = S_fM - S_f2P;
% Final solution
[Qr,Rr] = qr2([Vr_f1,Vr_f2],rvals);
[Qz,Rz] = qr([Vz_f1,Vz_f2],0);
[U,S,V] = svd(Rr*blkdiag(S_fM,S_f2)*Rz',0);
r_nn = find(diag(S)>1.0e-14,1,'last'); %Ensure square
Vr_nn = Qr*U(:,1:r_nn);
Vz_nn = Qz*V(:,1:r_nn);
S_nn = S(1:r_nn,1:r_nn);
end
function [Q, R] = qr2(X, rvals)
    [Q, R] = qr(sqrt(rvals) .* X, 0);
    Q = Q ./ sqrt(rvals);
end