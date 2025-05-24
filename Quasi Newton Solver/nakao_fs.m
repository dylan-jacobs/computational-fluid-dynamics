% clear all;close all;
tic
Nx = 80;
xvals = linspace(0,200,Nx+1)';
dx = xvals(2)-xvals(1);
xvals = xvals(1:Nx)+dx/2; %cell-centered grid

q = 1;         %ion charge
qe = -1;       %electron charge
m = 1;         %ion mass
me = 1/1836;   %electron mass

n_IC = @(x) 0.36*tanh(0.05*(x-100)) + 0.64;
u_IC = @(x) -1.1154*tanh(0.05*(x-100))+1.983;
Ti_IC = @(x) 0.4424*tanh(0.05*(x-100))+0.5576;
Te_IC = Ti_IC; %only for initial condition

Nvperp = 50; Nvpar = 100;
vperp = linspace(0,8,Nvperp+1)'; vpar = linspace(-8,10,Nvpar+1)';
dvperp = vperp(2)-vperp(1); dvpar = vpar(2)-vpar(1);
vperp = vperp(1:Nvperp)+dvperp/2; vpar = vpar(1:Nvpar)+dvpar/2;
[Vperp,Vpar] = meshgrid(vperp,vpar); Vperp = Vperp'; Vpar = Vpar';
R = 1; %gas constant(!!!)
maxwell = @(n,u,T) (n/(2*pi*R*T)^(3/2))*exp(-((Vpar-u).^2+Vperp.^2)/(2*R*T));
leftBC = maxwell(n_IC(0-dx/2),u_IC(0-dx/2),Ti_IC(0-dx/2));        %maxwellian at x=0-dx/2 (ghost point; Dirichlet and Neumann BC)
rightBC = maxwell(n_IC(200+dx/2),u_IC(200+dx/2),Ti_IC(200+dx/2)); %maxwellian at x=200+dx/2 (ghost point; Dirichlet and Neumann BC)
% !!! Might need to normalize the distribution functions by their numerical number densities

nvals = n_IC(xvals);          %number density
uvals = u_IC(xvals);          %bulk velocity (in parallel component)
Tivals = Ti_IC(xvals);        %ion temperature
Tevals = Te_IC(xvals);        %electron temperature

Tf = 100;
tvals = [0,0.005:0.3:Tf];
if tvals(end) ~= Tf
    tvals = [tvals,Tf];
end
Nt = numel(tvals);


for n = 2:Nt
dtn = tvals(n)-tvals(n-1);
tn = tvals(n-1); tnn = tvals(n);
disp(num2str(tn))
[nvals,uvals,Tivals,Tevals] = QN(nvals,uvals,Tivals,Tevals,q,qe,m,me,Nx,dx,dtn,maxwell,leftBC,rightBC,dvperp,dvpar,vperp,vpar,Vperp,Vpar,n_IC,u_IC,Ti_IC,Te_IC);
end

figure(1);clf;plot(xvals,nvals,'r-','linewidth',1.5);
hold on;plot(xvals,uvals/3.0984,'b-','linewidth',1.5);
hold on;plot(xvals,Tivals,'k-','linewidth',1.5);
hold on;plot(xvals,Tevals,'g-','linewidth',1.5);
hold on;xlabel('x');
legend('n','u/u0','Ti','Te');

toc

%% Functions

function [n_out,u_out,Ti_out,Te_out] = QN(nvals,uvals,Tivals,Tevals,q,qe,m,me,Nx,dx,dt,maxwell,leftBC,rightBC,dvperp,dvpar,vperp,vpar,Vperp,Vpar,n_IC,u_IC,Ti_IC,Te_IC)
sz = size(Vperp); %size of array in velocity space

% Compute numerical flux with upwinding
nuhat = zeros(Nx+1,1);     %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
Shat = zeros(Nx+1,1);      %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
Qhat = zeros(Nx+1,1);      %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
nTehat = zeros(Nx+1,1);    %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
pos_neg = find(vpar>0,1);  %first entry for which vpar is positive

% At x_{1/2}=0
fleft = leftBC;                                  %at x_0=0-dx/2, ghost point
fright = maxwell(nvals(1),uvals(1),Tivals(1));   %at x_1
fhat = zeros(sz(1),sz(2));
fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
nuhat(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
Shat(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
Qhat(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));
for i = 1:Nx-1 % At x_{i+1/2}, interior cell boundaries
    fleft = maxwell(nvals(i),uvals(i),Tivals(i));         %at x_i
    fright = maxwell(nvals(i+1),uvals(i+1),Tivals(i+1));  %at x_{i+1}
    fhat = zeros(sz(1),sz(2));
    fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
    fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
    nuhat(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
    Shat(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
    Qhat(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));
end
% At x_{Nx+1/2}=200
fleft = maxwell(nvals(Nx),uvals(Nx),Tivals(Nx)); %at x_Nx
fright = rightBC;                                %at x_{Nx+1}=200+dx/2, ghost point
fhat = zeros(sz(1),sz(2));
fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
nuhat(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
Shat(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
Qhat(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));
uhat = ([u_IC(0-dx/2);uvals] + [uvals;u_IC(200+dx/2)])/2; %u_{i+1/2}=(u_i+u_{i+1})/2. Size (Nx+1)x1 for Nx+1 cell boundaries
% At x_{1/2}=0
if uhat(1)>0
    nTehat(1) = n_IC(0-dx/2)*Te_IC(0-dx/2);
else
    nTehat(1) = nvals(1)*Tevals(1);
end
for i = 1:Nx-1 % At x_{i+1/2}, interior cell boundaries
    if uhat(i+1)>0
        nTehat(i+1) = nvals(i)*Tevals(i);
    else
        nTehat(i+1) = nvals(i+1)*Tevals(i+1);
    end
end
% At x_{Nx+1/2}=200
if uhat(Nx+1)>0
    nTehat(Nx+1) = nvals(Nx)*Tevals(Nx);
else
    nTehat(Nx+1) = n_IC(200+dx/2)*Te_IC(200+dx/2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol = 5.0e-4;
err = 1;
k = 0; %number of iterations

n_tnn = nvals - (dt/dx)*(nuhat(2:end)-nuhat(1:end-1));
n_tnn_right = [n_tnn(2:Nx);n_IC(200+dx/2)]; %with BC for ghost cell
n_tnn_left = [n_IC(0-dx/2);n_tnn(1:Nx-1)];  %with BC for ghost cell

u_k = uvals; %k=0
nu_k = n_tnn.*u_k; %k=0
nU_k = (3*n_tnn.*Tivals/m + n_tnn.*uvals.^2)/2; %k=0
Te_k = Tevals; %k=0
x_k = [nu_k;nU_k;Te_k];

Te_k_right = [Te_k(2:Nx);Te_IC(200+dx/2)]; %at x_2, x_3, ..., x_N, x_{N+1} (ghost cell)
Te_k_left = [Te_IC(0-dx/2);Te_k(1:Nx-1)];  %at x_0 (ghost cell), x_1, ..., x_N
Rk_nu = nu_k - nvals.*uvals + (dt/dx)*(Shat(2:end)-Shat(1:end-1)) - (dt/(2*dx))*(q/(m*qe))*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
Rk_nU = nU_k - (3*nvals.*Tivals/m + nvals.*uvals.^2)/2 + (dt/dx)*(Qhat(2:end)-Qhat(1:end-1)) - ((3*dt*sqrt(2)*sqrt(me))/m^2)*((n_tnn.^2.*Te_k - (m/3)*(2*n_tnn.*nU_k-nu_k.^2))./(Te_k.^(1.5)))...
    - (dt/(2*dx))*(q/qe)*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
Rk_Te = n_tnn.*Te_k - nvals.*Tevals + ((5*dt)/(3*dx))*(uhat(2:end).*nTehat(2:end)-uhat(1:end-1).*nTehat(1:end-1)) - (dt/(3*dx))*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)...
    - (dt/(3*dx^2))*(3.2/sqrt(2*me))*((Te_k.^(2.5)+Te_k_right.^(2.5)).*(Te_k_right-Te_k) - (Te_k_left.^(2.5)+Te_k.^(2.5)).*(Te_k-Te_k_left))...
    - ((2*dt*sqrt(2*me))/m)*((m/3)*(2*n_tnn.*nU_k - nu_k.^2) - n_tnn.^2.*Te_k)./(Te_k.^(1.5));
Rk = -[Rk_nu;Rk_nU;Rk_Te];

Rn_err = norm(Rk,2);

while err > tol
    Pk_nunu = eye(Nx,Nx);
    Pk_nunU = zeros(Nx,Nx);
    Pk_nuTe = zeros(Nx,Nx);
        Pk_nuTe(1:end-1,2:end) = Pk_nuTe(1:end-1,2:end) -(dt/(2*dx))*(q/(m*qe))*diag(n_tnn(2:end)); %superdiagonal
        Pk_nuTe(2:end,1:end-1) = Pk_nuTe(2:end,1:end-1) +(dt/(2*dx))*(q/(m*qe))*diag(n_tnn(1:end-1)); %subdiagonal

    Pk_nUnu = -((dt*q)/(2*dx*qe))*diag((n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)./n_tnn);
    Pk_nUnU = eye(Nx,Nx);
    Pk_nUTe = zeros(Nx,Nx);
        Pk_nUTe(1:end-1,2:end) = Pk_nUTe(1:end-1,2:end) -(dt/(2*dx))*(q/qe)*diag(nu_k(1:end-1).*n_tnn(2:end)./n_tnn(1:end-1)); %superdiagonal
        Pk_nUTe(2:end,1:end-1) = Pk_nUTe(2:end,1:end-1) +(dt/(2*dx))*(q/qe)*diag(nu_k(2:end).*n_tnn(1:end-1)./n_tnn(2:end)); %subdiagonal
        Pk_nUTe = Pk_nUTe - (3*dt*sqrt(2*me)/m^2)*diag( -0.5*(n_tnn.^2)./(Te_k.^(1.5)) - (m/3)*(-3/2)*(2*n_tnn.*nU_k-nu_k.^2)./(Te_k.^(2.5)) ); %diagonal

    Pk_Tenu = -(dt/(3*dx))*diag((n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)./n_tnn);
    Pk_TenU = zeros(Nx,Nx);
    Pk_TeTe = zeros(Nx,Nx);
        Pk_TeTe(1:end-1,2:end) = Pk_TeTe(1:end-1,2:end) -(dt/(3*dx))*diag(nu_k(1:end-1).*n_tnn(2:end)./n_tnn(1:end-1))...
            - (dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k(2:end).^(1.5).*(Te_k(2:end)-Te_k(1:end-1)) + (Te_k(1:end-1).^(2.5) + Te_k(2:end).^(2.5))); %superdiagonal
        Pk_TeTe(2:end,1:end-1) = Pk_TeTe(2:end,1:end-1) +(dt/(3*dx))*diag(nu_k(2:end).*n_tnn(1:end-1)./n_tnn(2:end))...
            + (dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k(1:end-1).^(1.5).*(Te_k(2:end)-Te_k(1:end-1)) - (Te_k(1:end-1).^(2.5) + Te_k(2:end).^(2.5))); %subdiagonal
        Pk_TeTe = Pk_TeTe + diag(n_tnn) - (dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k.^(1.5).*(Te_k_right - 2*Te_k + Te_k_left) - (Te_k_right.^(2.5) + 2*Te_k.^(2.5) + Te_k_left.^(2.5)))...
            - (2*dt*sqrt(2*me)/m)*diag(-1.5*(m/3)*((2*n_tnn.*nU_k-nu_k.^2)./(Te_k.^(2.5))) + 0.5*(n_tnn.^2)./(Te_k.^(1.5))); %diagonal

    Pk = [Pk_nunu, Pk_nunU, Pk_nuTe; Pk_nUnu, Pk_nUnU, Pk_nUTe; Pk_Tenu, Pk_TenU, Pk_TeTe];

    dx_kk = Pk\Rk;
    x_kk = x_k + dx_kk;

    nu_kk = x_kk(1:Nx);
    nU_kk = x_kk(Nx+1:2*Nx);
    Te_kk = x_kk(2*Nx+1:3*Nx);

    nu_k = nu_kk;
    nU_k = nU_kk;
    Te_k = Te_kk;
    x_k = x_kk;


    Te_k_right = [Te_k(2:Nx);Te_IC(200+dx/2)]; %at x_2, x_3, ..., x_N, x_{N+1} (ghost cell)
    Te_k_left = [Te_IC(0-dx/2);Te_k(1:Nx-1)];  %at x_0 (ghost cell), x_1, ..., x_N
    Rk_nu = nu_k - nvals.*uvals + (dt/dx)*(Shat(2:end)-Shat(1:end-1)) - (dt/(2*dx))*(q/(m*qe))*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
    Rk_nU = nU_k - (3*nvals.*Tivals/m + nvals.*uvals.^2)/2 + (dt/dx)*(Qhat(2:end)-Qhat(1:end-1)) - ((3*dt*sqrt(2)*sqrt(me))/m^2)*((n_tnn.^2.*Te_k - (m/3)*(2*n_tnn.*nU_k-nu_k.^2))./(Te_k.^(1.5)))...
        - (dt/(2*dx))*(q/qe)*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
    Rk_Te = n_tnn.*Te_k - nvals.*Tevals + ((5*dt)/(3*dx))*(uhat(2:end).*nTehat(2:end)-uhat(1:end-1).*nTehat(1:end-1)) - (dt/(3*dx))*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)...
        - (dt/(3*dx^2))*(3.2/sqrt(2*me))*((Te_k.^(2.5)+Te_k_right.^(2.5)).*(Te_k_right-Te_k) - (Te_k_left.^(2.5)+Te_k.^(2.5)).*(Te_k-Te_k_left))...
        - ((2*dt*sqrt(2*me))/m)*((m/3)*(2*n_tnn.*nU_k - nu_k.^2) - n_tnn.^2.*Te_k)./(Te_k.^(1.5));
    Rk = -[Rk_nu;Rk_nU;Rk_Te];

    err = norm(Rk,2)/Rn_err;

    k = k+1;
    if k==500
        disp('Completed 500 iterations');
        err
        return
    end    
end

n_out = n_tnn;
u_out = nu_k./n_tnn;
Ti_out = (m/3)*(2*nU_k./n_tnn - (nu_k.^2)./(n_tnn.^2));
Te_out = Te_k;
end



% function [n_out,u_out,Ti_out,Te_out] = QN(nvals,uvals,Tivals,Tevals,q,qe,m,me,Nx,dx,dt,maxwell,leftBC,rightBC,dvperp,dvpar,vperp,vpar,Vperp,Vpar,n_IC,u_IC,Ti_IC,Te_IC)
% sz = size(Vperp); %size of array in velocity space
% 
% % Compute numerical flux with upwinding
% nuhat = zeros(Nx+1,1);     %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
% Shat = zeros(Nx+1,1);      %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
% Qhat = zeros(Nx+1,1);      %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
% nTehat = zeros(Nx+1,1);    %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
% pos_neg = find(vpar>0,1);  %first entry for which vpar is positive
% 
% % At x_{1/2}=0
% fleft = leftBC;                                  %at x_0=0-dx/2, ghost point
% fright = maxwell(nvals(1),uvals(1),Tivals(1));   %at x_1
% fhat = zeros(sz(1),sz(2));
% fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
% fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
% nuhat(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
% Shat(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
% Qhat(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^3.*fhat.*Vperp/2));
% for i = 1:Nx-1 % At x_{i+1/2}, interior cell boundaries
%     fleft = maxwell(nvals(i),uvals(i),Tivals(i));         %at x_i
%     fright = maxwell(nvals(i+1),uvals(i+1),Tivals(i+1));  %at x_{i+1}
%     fhat = zeros(sz(1),sz(2));
%     fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
%     fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
%     nuhat(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
%     Shat(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
%     Qhat(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^3.*fhat.*Vperp/2));
% end
% % At x_{Nx+1/2}=200
% fleft = maxwell(nvals(Nx),uvals(Nx),Tivals(Nx)); %at x_Nx
% fright = rightBC;                                %at x_{Nx+1}=200+dx/2, ghost point
% fhat = zeros(sz(1),sz(2));
% fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
% fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
% nuhat(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
% Shat(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
% Qhat(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^3.*fhat.*Vperp/2));
% 
% uhat = ([u_IC(0-dx/2);uvals] + [uvals;u_IC(100+dx/2)])/2; %u_{i+1/2}=(u_i+u_{i+1})/2. Size (Nx+1)x1 for Nx+1 cell boundaries
% % At x_{1/2}=0
% if uhat(1)>0
%     nTehat(1) = n_IC(0-dx/2)*Te_IC(0-dx/2);
% else
%     nTehat(1) = nvals(1)*Tevals(1);
% end
% for i = 1:Nx-1 % At x_{i+1/2}, interior cell boundaries
%     if uhat(i+1)>0
%         nTehat(i+1) = nvals(i)*Tevals(i);
%     else
%         nTehat(i+1) = nvals(i+1)*Tevals(i+1);
%     end
% end
% % At x_{Nx+1/2}=200
% if uhat(Nx+1)>0
%     nTehat(Nx+1) = nvals(Nx)*Tevals(Nx);
% else
%     nTehat(Nx+1) = n_IC(200+dx/2)*Te_IC(200+dx/2);
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% tol = 5.0e-4;
% err = 1;
% k = 0; %number of iterations
% 
% n_tnn = nvals - (dt/dx)*(nuhat(2:end)-nuhat(1:end-1));
% n_tnn_right = [n_tnn(2:Nx);n_IC(200+dx/2)]; %with BC for ghost cell
% n_tnn_left = [n_IC(0-dx/2);n_tnn(1:Nx-1)];  %with BC for ghost cell
% 
% u_k = uvals; %k=0
% nu_k = n_tnn.*u_k; %k=0
% nU_k = (3*n_tnn.*Tivals/m + n_tnn.*uvals.^2)/2; %k=0
% Te_k = Tevals; %k=0
% x_k = [nu_k;nU_k;Te_k];
% 
% Te_k_right = [Te_k(2:Nx);Te_IC(200+dx/2)]; %at x_2, x_3, ..., x_N, x_{N+1} (ghost cell)
% Te_k_left = [Te_IC(0-dx/2);Te_k(1:Nx-1)];  %at x_0 (ghost cell), x_1, ..., x_N
% Rk_nu = nu_k - nvals.*uvals + (dt/dx)*(Shat(2:end)-Shat(1:end-1)) - (dt/(2*dx))*(q/(m*qe))*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
% Rk_nU = nU_k - (3*nvals.*Tivals/m + nvals.*uvals.^2)/2 + (dt/dx)*(Qhat(2:end)-Qhat(1:end-1)) - ((3*dt*sqrt(me))/m^2)*((n_tnn.^2.*Te_k - (m/3)*(2*n_tnn.*nU_k-nu_k.^2))./(Te_k.^(1.5)))...
%     - (dt/(2*dx))*(q/qe)*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
% Rk_Te = n_tnn.*Te_k - nvals.*Tevals + ((5*dt)/(3*dx))*(uhat(2:end).*nTehat(2:end)-uhat(1:end-1).*nTehat(1:end-1)) - (dt/(3*dx))*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)...
%     - (dt/(3*dx^2))*(3.2/sqrt(me))*((Te_k.^(2.5)+Te_k_right.^(2.5)).*(Te_k_right-Te_k) - (Te_k_left.^(2.5)+Te_k.^(2.5)).*(Te_k-Te_k_left))...
%     - ((2*dt*sqrt(me))/m)*((m/3)*(2*n_tnn.*nU_k - nu_k.^2) - n_tnn.^2.*Te_k)./(Te_k.^(1.5));
% Rk = -[Rk_nu;Rk_nU;Rk_Te];
% 
% Rn_err = norm(Rk,2);
% 
% while err > tol
%     Pk_nunu = eye(Nx,Nx);
%     Pk_nunU = zeros(Nx,Nx);
%     Pk_nuTe = zeros(Nx,Nx);
%         Pk_nuTe(1:end-1,2:end) = Pk_nuTe(1:end-1,2:end) -(dt/(2*dx))*(q/(m*qe))*diag(n_tnn(2:end)); %superdiagonal
%         Pk_nuTe(2:end,1:end-1) = Pk_nuTe(2:end,1:end-1) +(dt/(2*dx))*(q/(m*qe))*diag(n_tnn(1:end-1)); %subdiagonal
% 
%     Pk_nUnu = -((dt*q)/(2*dx*qe))*diag((n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)./n_tnn);
%     Pk_nUnU = eye(Nx,Nx);
%     Pk_nUTe = zeros(Nx,Nx);
%         Pk_nUTe(1:end-1,2:end) = Pk_nUTe(1:end-1,2:end) -(dt/(2*dx))*(q/qe)*diag(nu_k(1:end-1).*n_tnn(2:end)./n_tnn(1:end-1)); %superdiagonal
%         Pk_nUTe(2:end,1:end-1) = Pk_nUTe(2:end,1:end-1) +(dt/(2*dx))*(q/qe)*diag(nu_k(2:end).*n_tnn(1:end-1)./n_tnn(2:end)); %subdiagonal
%         Pk_nUTe = Pk_nUTe - (3*dt*sqrt(2*me)/m^2)*diag( -0.5*(n_tnn.^2)./(Te_k.^(1.5)) - (m/3)*(-3/2)*(2*n_tnn.*nU_k-nu_k.^2)./(Te_k.^(2.5)) ); %diagonal
% 
%     Pk_Tenu = -(dt/(3*dx))*diag((n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)./n_tnn);
%     Pk_TenU = zeros(Nx,Nx);
%     Pk_TeTe = zeros(Nx,Nx);
%         Pk_TeTe(1:end-1,2:end) = Pk_TeTe(1:end-1,2:end) -(dt/(3*dx))*diag(nu_k(1:end-1).*n_tnn(2:end)./n_tnn(1:end-1))...
%             - (dt/(3*dx^2))*(3.2/sqrt(me))*diag(2.5*Te_k(2:end).^(1.5).*(Te_k(2:end)-Te_k(1:end-1)) + (Te_k(1:end-1).^(2.5) + Te_k(2:end).^(2.5))); %superdiagonal
%         Pk_TeTe(2:end,1:end-1) = Pk_TeTe(2:end,1:end-1) +(dt/(3*dx))*diag(nu_k(2:end).*n_tnn(1:end-1)./n_tnn(2:end))...
%             + (dt/(3*dx^2))*(3.2/sqrt(me))*diag(2.5*Te_k(1:end-1).^(1.5).*(Te_k(2:end)-Te_k(1:end-1)) - (Te_k(1:end-1).^(2.5) + Te_k(2:end).^(2.5))); %subdiagonal
%         Pk_TeTe = Pk_TeTe + diag(n_tnn) - (dt/(3*dx^2))*(3.2/sqrt(me))*diag(2.5*Te_k.^(1.5).*(Te_k_right - 2*Te_k + Te_k_left) - (Te_k_right.^(2.5) + 2*Te_k.^(2.5) + Te_k_left.^(2.5)))...
%             - (2*dt*sqrt(me)/m)*diag(-1.5*(m/3)*((2*n_tnn.*nU_k-nu_k.^2)./(Te_k.^(2.5))) + 0.5*(n_tnn.^2)./(Te_k.^(1.5))); %diagonal
% 
%     Pk = [Pk_nunu, Pk_nunU, Pk_nuTe; Pk_nUnu, Pk_nUnU, Pk_nUTe; Pk_Tenu, Pk_TenU, Pk_TeTe];
% 
%     dx_kk = Pk\Rk;
%     x_kk = x_k + dx_kk;
% 
%     nu_kk = x_kk(1:Nx);
%     nU_kk = x_kk(Nx+1:2*Nx);
%     Te_kk = x_kk(2*Nx+1:3*Nx);
% 
%     nu_k = nu_kk;
%     nU_k = nU_kk;
%     Te_k = Te_kk;
%     x_k = x_kk;
% 
% 
%     Te_k_right = [Te_k(2:Nx);Te_IC(200+dx/2)]; %at x_2, x_3, ..., x_N, x_{N+1} (ghost cell)
%     Te_k_left = [Te_IC(0-dx/2);Te_k(1:Nx-1)];  %at x_0 (ghost cell), x_1, ..., x_N
%     Rk_nu = nu_k - nvals.*uvals + (dt/dx)*(Shat(2:end)-Shat(1:end-1)) - (dt/(2*dx))*(q/(m*qe))*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
%     Rk_nU = nU_k - (3*nvals.*Tivals/m + nvals.*uvals.^2)/2 + (dt/dx)*(Qhat(2:end)-Qhat(1:end-1)) - ((3*dt*sqrt(me))/m^2)*((n_tnn.^2.*Te_k - (m/3)*(2*n_tnn.*nU_k-nu_k.^2))./(Te_k.^(1.5)))...
%         - (dt/(2*dx))*(q/qe)*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
%     Rk_Te = n_tnn.*Te_k - nvals.*Tevals + ((5*dt)/(3*dx))*(uhat(2:end).*nTehat(2:end)-uhat(1:end-1).*nTehat(1:end-1)) - (dt/(3*dx))*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)...
%         - (dt/(3*dx^2))*(3.2/sqrt(me))*((Te_k.^(2.5)+Te_k_right.^(2.5)).*(Te_k_right-Te_k) - (Te_k_left.^(2.5)+Te_k.^(2.5)).*(Te_k-Te_k_left))...
%         - ((2*dt*sqrt(me))/m)*((m/3)*(2*n_tnn.*nU_k - nu_k.^2) - n_tnn.^2.*Te_k)./(Te_k.^(1.5));
%     Rk = -[Rk_nu;Rk_nU;Rk_Te];
% 
%     err = norm(Rk,2)/Rn_err;
% 
%     k = k+1;
%     if k==500
%         disp('Completed 500 iterations');
%         err
%         return
%     end    
% end
% 
% n_out = n_tnn;
% u_out = nu_k./n_tnn;
% Ti_out = (m/3)*(2*nU_k./n_tnn - (nu_k.^2)./(n_tnn.^2));
% Te_out = Te_k;
% end
















