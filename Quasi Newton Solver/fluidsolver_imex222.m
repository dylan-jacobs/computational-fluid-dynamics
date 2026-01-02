clear all;close all;

% Nx=100, Nvperp=50, Nvpar=100, Lv=8, Tf=15, dt=0.001
% vals = load('fluid_refsoln.mat');
% vals = vals.soln;

DTvals = [0.3];
% DTvals = [0.4,0.2,0.1,0.05,0.025,0.0125];
ERR = zeros(numel(DTvals),4);
for pp = 1:numel(DTvals)

Nx = 100;
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


Nvperp = 50; Nvpar = 100; Lv = 8;
vperp = linspace(0,Lv,Nvperp+1)'; vpar = linspace(-Lv,Lv,Nvpar+1)';
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

Tf = 1;
tvals = [0,0.005:DTvals(pp):Tf];% tvals = sort([tvals,25,50,75,100,125,150,175,200,225,250,275,300,325,350,375]);
if tvals(end) ~= Tf
    tvals = [tvals,Tf];
end
Nt = numel(tvals);


% Global quantities
totalmass = zeros(Nt,1);
totalmomentum = zeros(Nt,1);
totalenergy = zeros(Nt,1);

totalmass(1) = dx*m*sum(nvals);
totalmomentum(1) = dx*m*sum(nvals.*uvals);
totalenergy(1) = dx*m*sum(0.5*((3/m)*nvals.*Tivals + nvals.*uvals.^2) + (3/2)*nvals.*Tevals);


tic
for n = 2:Nt
dtn = tvals(n)-tvals(n-1);
tn = tvals(n-1); tnn = tvals(n);
disp(num2str(tnn))
[nvals,uvals,Tivals,Tevals] = QN(nvals,uvals,Tivals,Tevals,q,qe,m,me,Nx,dx,dtn,maxwell,leftBC,rightBC,dvperp,dvpar,vperp,vpar,Vperp,Vpar,n_IC,u_IC,Ti_IC,Te_IC);
figure(2);clf;plot(xvals,nvals,'r-','linewidth',1.5);
hold on;plot(xvals,uvals/3.0984,'b-','linewidth',1.5);
hold on;plot(xvals,Tivals,'k-','linewidth',1.5);
hold on;plot(xvals,Tevals,'g-','linewidth',1.5);
hold on;xlabel('x');title(['Numerical solution at time t=',num2str(tnn)]);
legend('n','u/u0','Ti','Te');
axis([0,200,0,1.5]);

end
toc
figure(1);clf;plot(xvals,nvals,'r-','linewidth',1.5);
hold on;plot(xvals,uvals/3.0984,'b-','linewidth',1.5);
hold on;plot(xvals,Tivals,'k-','linewidth',1.5);
hold on;plot(xvals,Tevals,'g-','linewidth',1.5);
hold on;xlabel('x');
legend('n','u/u0','Ti','Te');
axis([0,200,0,1.5]);

ERR(pp,1) = dx*sum(abs(vals(:,1)-nvals));
ERR(pp,2) = dx*sum(abs(vals(:,2)-uvals));
ERR(pp,3) = dx*sum(abs(vals(:,3)-Tivals));
ERR(pp,4) = dx*sum(abs(vals(:,4)-Tevals));
end


disp(ERR);
disp(log2(ERR(1:end-1,:)./ERR(2:end,:)));


%% Functions

function [n_out,u_out,Ti_out,Te_out] = QN(nvals0,uvals0,Tivals0,Tevals0,q,qe,m,me,Nx,dx,dt,maxwell,leftBC,rightBC,dvperp,dvpar,vperp,vpar,Vperp,Vpar,n_IC,u_IC,Ti_IC,Te_IC)
gamma = 1 - sqrt(2)/2;
delta = 1 - 1/(2*gamma);
sz = size(Vperp); %size of array in velocity space

%%%%%%%%%%%%%%%%%%%
% % % STAGE 1 % % %
%%%%%%%%%%%%%%%%%%%

% Compute numerical flux with upwinding
nuhat0 = zeros(Nx+1,1);     %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
Shat0 = zeros(Nx+1,1);      %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
Qhat0 = zeros(Nx+1,1);      %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
nTehat0 = zeros(Nx+1,1);    %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
pos_neg = find(vpar>0,1);  %first entry for which vpar is positive

% At x_{1/2}=0
fleft = leftBC;                                  %at x_0=0-dx/2, ghost point
fright = maxwell(nvals0(1),uvals0(1),Tivals0(1));   %at x_1
fhat = zeros(sz(1),sz(2));
fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
nuhat0(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
Shat0(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
Qhat0(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));
for i = 1:Nx-1 % At x_{i+1/2}, interior cell boundaries
    fleft = maxwell(nvals0(i),uvals0(i),Tivals0(i));         %at x_i
    fright = maxwell(nvals0(i+1),uvals0(i+1),Tivals0(i+1));  %at x_{i+1}
    fhat = zeros(sz(1),sz(2));
    fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
    fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
    nuhat0(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
    Shat0(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
    Qhat0(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));
end
% At x_{Nx+1/2}=200
fleft = maxwell(nvals0(Nx),uvals0(Nx),Tivals0(Nx)); %at x_Nx
fright = rightBC;                                %at x_{Nx+1}=200+dx/2, ghost point
fhat = zeros(sz(1),sz(2));
fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
nuhat0(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
Shat0(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
Qhat0(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));

uhat0 = ([u_IC(0-dx/2);uvals0] + [uvals0;u_IC(200+dx/2)])/2; %u_{i+1/2}=(u_i+u_{i+1})/2. Size (Nx+1)x1 for Nx+1 cell boundaries
% At x_{1/2}=0
if uhat0(1)>0
    nTehat0(1) = n_IC(0-dx/2)*Te_IC(0-dx/2);
else
    nTehat0(1) = nvals0(1)*Tevals0(1);
end
for i = 1:Nx-1 % At x_{i+1/2}, interior cell boundaries
    if uhat0(i+1)>0
        nTehat0(i+1) = nvals0(i)*Tevals0(i);
    else
        nTehat0(i+1) = nvals0(i+1)*Tevals0(i+1);
    end
end
% At x_{Nx+1/2}=200
if uhat0(Nx+1)>0
    nTehat0(Nx+1) = nvals0(Nx)*Tevals0(Nx);
else
    nTehat0(Nx+1) = n_IC(200+dx/2)*Te_IC(200+dx/2);
end

err = 1;
k = 0; %number of iterations

n_t1 = nvals0 - (gamma*dt/dx)*(nuhat0(2:end)-nuhat0(1:end-1));
n_t1_right = [n_t1(2:Nx);n_IC(200+dx/2)]; %with BC for ghost cell
n_t1_left = [n_IC(0-dx/2);n_t1(1:Nx-1)];  %with BC for ghost cell

u_k = uvals0; %k=0
nu_k = n_t1.*u_k; %k=0
nU_k = (3*n_t1.*Tivals0/m + n_t1.*uvals0.^2)/2; %k=0
Te_k = Tevals0; %k=0
x_k = [nu_k;nU_k;Te_k];

Te_k_right = [Te_k(2:Nx);Te_IC(200+dx/2)]; %at x_2, x_3, ..., x_N, x_{N+1} (ghost cell)
Te_k_left = [Te_IC(0-dx/2);Te_k(1:Nx-1)];  %at x_0 (ghost cell), x_1, ..., x_N
Rk_nu = nu_k - nvals0.*uvals0 + (gamma*dt/dx)*(Shat0(2:end)-Shat0(1:end-1)) - (gamma*dt/(2*dx))*(q/(m*qe))*(n_t1_right.*Te_k_right - n_t1_left.*Te_k_left);
Rk_nU = nU_k - (3*nvals0.*Tivals0/m + nvals0.*uvals0.^2)/2 + (gamma*dt/dx)*(Qhat0(2:end)-Qhat0(1:end-1)) - ((3*gamma*dt*sqrt(2)*sqrt(me))/m^2)*((n_t1.^2.*Te_k - (m/3)*(2*n_t1.*nU_k-nu_k.^2))./(Te_k.^(1.5)))...
    - (gamma*dt/(2*dx))*(q/qe)*(nu_k./n_t1).*(n_t1_right.*Te_k_right - n_t1_left.*Te_k_left);
Rk_Te = n_t1.*Te_k - nvals0.*Tevals0 + ((5*gamma*dt)/(3*dx))*(uhat0(2:end).*nTehat0(2:end)-uhat0(1:end-1).*nTehat0(1:end-1)) - (gamma*dt/(3*dx))*(nu_k./n_t1).*(n_t1_right.*Te_k_right - n_t1_left.*Te_k_left)...
    - (gamma*dt/(3*dx^2))*(3.2/sqrt(2*me))*((Te_k.^(2.5)+Te_k_right.^(2.5)).*(Te_k_right-Te_k) - (Te_k_left.^(2.5)+Te_k.^(2.5)).*(Te_k-Te_k_left))...
    - ((2*gamma*dt*sqrt(2*me))/m)*((m/3)*(2*n_t1.*nU_k - nu_k.^2) - n_t1.^2.*Te_k)./(Te_k.^(1.5));
Rk = -[Rk_nu;Rk_nU;Rk_Te];

tol = min(5e-12,max(abs(Rk))*5e-10);
%tol = 5.0e-4;

while err > tol
    Pk_nunu = eye(Nx,Nx);
    Pk_nunU = zeros(Nx,Nx);
    Pk_nuTe = zeros(Nx,Nx);
        Pk_nuTe(1:end-1,2:end) = Pk_nuTe(1:end-1,2:end) -(gamma*dt/(2*dx))*(q/(m*qe))*diag(n_t1(2:end)); %superdiagonal
        Pk_nuTe(2:end,1:end-1) = Pk_nuTe(2:end,1:end-1) +(gamma*dt/(2*dx))*(q/(m*qe))*diag(n_t1(1:end-1)); %subdiagonal

    Pk_nUnu = -((gamma*dt*q)/(2*dx*qe))*diag((n_t1_right.*Te_k_right - n_t1_left.*Te_k_left)./n_t1);
    Pk_nUnU = eye(Nx,Nx);
    Pk_nUTe = zeros(Nx,Nx);
        Pk_nUTe(1:end-1,2:end) = Pk_nUTe(1:end-1,2:end) -(gamma*dt/(2*dx))*(q/qe)*diag(nu_k(1:end-1).*n_t1(2:end)./n_t1(1:end-1)); %superdiagonal
        Pk_nUTe(2:end,1:end-1) = Pk_nUTe(2:end,1:end-1) +(gamma*dt/(2*dx))*(q/qe)*diag(nu_k(2:end).*n_t1(1:end-1)./n_t1(2:end)); %subdiagonal
        Pk_nUTe = Pk_nUTe - (3*gamma*dt*sqrt(2*me)/m^2)*diag( -0.5*(n_t1.^2)./(Te_k.^(1.5)) - (m/3)*(-3/2)*(2*n_t1.*nU_k-nu_k.^2)./(Te_k.^(2.5)) ); %diagonal

    Pk_Tenu = -(gamma*dt/(3*dx))*diag((n_t1_right.*Te_k_right - n_t1_left.*Te_k_left)./n_t1);
    Pk_TenU = zeros(Nx,Nx);
    Pk_TeTe = zeros(Nx,Nx);
        Pk_TeTe(1:end-1,2:end) = Pk_TeTe(1:end-1,2:end) -(gamma*dt/(3*dx))*diag(nu_k(1:end-1).*n_t1(2:end)./n_t1(1:end-1))...
            - (gamma*dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k(2:end).^(1.5).*(Te_k(2:end)-Te_k(1:end-1)) + (Te_k(1:end-1).^(2.5) + Te_k(2:end).^(2.5))); %superdiagonal
        Pk_TeTe(2:end,1:end-1) = Pk_TeTe(2:end,1:end-1) +(gamma*dt/(3*dx))*diag(nu_k(2:end).*n_t1(1:end-1)./n_t1(2:end))...
            + (gamma*dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k(1:end-1).^(1.5).*(Te_k(2:end)-Te_k(1:end-1)) - (Te_k(1:end-1).^(2.5) + Te_k(2:end).^(2.5))); %subdiagonal
        Pk_TeTe = Pk_TeTe + diag(n_t1) - (gamma*dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k.^(1.5).*(Te_k_right - 2*Te_k + Te_k_left) - (Te_k_right.^(2.5) + 2*Te_k.^(2.5) + Te_k_left.^(2.5)))...
            - (2*gamma*dt*sqrt(2*me)/m)*diag(-1.5*(m/3)*((2*n_t1.*nU_k-nu_k.^2)./(Te_k.^(2.5))) + 0.5*(n_t1.^2)./(Te_k.^(1.5))); %diagonal

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
    Rk_nu = nu_k - nvals0.*uvals0 + (gamma*dt/dx)*(Shat0(2:end)-Shat0(1:end-1)) - (gamma*dt/(2*dx))*(q/(m*qe))*(n_t1_right.*Te_k_right - n_t1_left.*Te_k_left);
    Rk_nU = nU_k - (3*nvals0.*Tivals0/m + nvals0.*uvals0.^2)/2 + (gamma*dt/dx)*(Qhat0(2:end)-Qhat0(1:end-1)) - ((3*gamma*dt*sqrt(2)*sqrt(me))/m^2)*((n_t1.^2.*Te_k - (m/3)*(2*n_t1.*nU_k-nu_k.^2))./(Te_k.^(1.5)))...
        - (gamma*dt/(2*dx))*(q/qe)*(nu_k./n_t1).*(n_t1_right.*Te_k_right - n_t1_left.*Te_k_left);
    Rk_Te = n_t1.*Te_k - nvals0.*Tevals0 + ((5*gamma*dt)/(3*dx))*(uhat0(2:end).*nTehat0(2:end)-uhat0(1:end-1).*nTehat0(1:end-1)) - (gamma*dt/(3*dx))*(nu_k./n_t1).*(n_t1_right.*Te_k_right - n_t1_left.*Te_k_left)...
        - (gamma*dt/(3*dx^2))*(3.2/sqrt(2*me))*((Te_k.^(2.5)+Te_k_right.^(2.5)).*(Te_k_right-Te_k) - (Te_k_left.^(2.5)+Te_k.^(2.5)).*(Te_k-Te_k_left))...
        - ((2*gamma*dt*sqrt(2*me))/m)*((m/3)*(2*n_t1.*nU_k - nu_k.^2) - n_t1.^2.*Te_k)./(Te_k.^(1.5));
    Rk = -[Rk_nu;Rk_nU;Rk_Te];

    err = max(abs(Rk));%norm(Rk,2)/Rn_err;

    k = k+1;
    if k==500
        disp('Completed 500 iterations - stage1');
        err
        return
    end    
end
%disp(k)
nvals1 = n_t1;
uvals1 = nu_k./n_t1;
Tivals1 = (m/3)*(2*nU_k./n_t1 - (nu_k.^2)./(n_t1.^2));
Tevals1 = Te_k;

%%%%%%%%%%%%%%%%%%%
% % % STAGE 2 % % %
%%%%%%%%%%%%%%%%%%%

u_1 = uvals1; 
nu_1 = n_t1.*u_1; 
nU_1 = (3*n_t1.*Tivals1/m + n_t1.*uvals1.^2)/2; 
Te_1 = Tevals1; 
Te_1_right = [Te_1(2:Nx);Te_IC(200+dx/2)]; %at x_2, x_3, ..., x_N, x_{N+1} (ghost cell)
Te_1_left = [Te_IC(0-dx/2);Te_1(1:Nx-1)];  %at x_0 (ghost cell), x_1, ..., x_N

% Compute numerical flux with upwinding
nuhat1 = zeros(Nx+1,1);     %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
Shat1 = zeros(Nx+1,1);      %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
Qhat1 = zeros(Nx+1,1);      %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
nTehat1 = zeros(Nx+1,1);    %at cell boundaries: 0=x1-dx/2 < x1+dx/2 < ... < xN-dx/2 < xN+dx/2=200
pos_neg = find(vpar>0,1);  %first entry for which vpar is positive

% At x_{1/2}=0
fleft = leftBC;                                  %at x_0=0-dx/2, ghost point
fright = maxwell(nvals1(1),uvals1(1),Tivals1(1));   %at x_1
fhat = zeros(sz(1),sz(2));
fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
nuhat1(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
Shat1(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
Qhat1(1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));
for i = 1:Nx-1 % At x_{i+1/2}, interior cell boundaries
    fleft = maxwell(nvals1(i),uvals1(i),Tivals1(i));         %at x_i
    fright = maxwell(nvals1(i+1),uvals1(i+1),Tivals1(i+1));  %at x_{i+1}
    fhat = zeros(sz(1),sz(2));
    fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
    fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
    nuhat1(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
    Shat1(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
    Qhat1(i+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));
end
% At x_{Nx+1/2}=200
fleft = maxwell(nvals1(Nx),uvals1(Nx),Tivals1(Nx)); %at x_Nx
fright = rightBC;                                %at x_{Nx+1}=200+dx/2, ghost point
fhat = zeros(sz(1),sz(2));
fhat(:,1:pos_neg-1) = fright(:,1:pos_neg-1);
fhat(:,pos_neg:end) = fleft(:,pos_neg:end);
nuhat1(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*fhat.*Vperp));
Shat1(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.^2.*fhat.*Vperp));
Qhat1(Nx+1) = 2*pi*dvperp*dvpar*sum(sum(Vpar.*(Vpar.^2+Vperp.^2).*fhat.*Vperp/2));

uhat1 = ([u_IC(0-dx/2);uvals1] + [uvals1;u_IC(200+dx/2)])/2; %u_{i+1/2}=(u_i+u_{i+1})/2. Size (Nx+1)x1 for Nx+1 cell boundaries
% At x_{1/2}=0
if uhat1(1)>0
    nTehat1(1) = n_IC(0-dx/2)*Te_IC(0-dx/2);
else
    nTehat1(1) = nvals1(1)*Tevals1(1);
end
for i = 1:Nx-1 % At x_{i+1/2}, interior cell boundaries
    if uhat1(i+1)>0
        nTehat1(i+1) = nvals1(i)*Tevals1(i);
    else
        nTehat1(i+1) = nvals1(i+1)*Tevals1(i+1);
    end
end
% At x_{Nx+1/2}=200
if uhat1(Nx+1)>0
    nTehat1(Nx+1) = nvals1(Nx)*Tevals1(Nx);
else
    nTehat1(Nx+1) = n_IC(200+dx/2)*Te_IC(200+dx/2);
end

err = 1;
k = 0; %number of iterations

n_tnn = nvals0 - (delta*dt/dx)*(nuhat0(2:end)-nuhat0(1:end-1)) - ((1-delta)*dt/dx)*(nuhat1(2:end)-nuhat1(1:end-1));
n_tnn_right = [n_tnn(2:Nx);n_IC(200+dx/2)]; %with BC for ghost cell
n_tnn_left = [n_IC(0-dx/2);n_tnn(1:Nx-1)];  %with BC for ghost cell

u_k = uvals0; %k=0
nu_k = n_tnn.*u_k; %k=0
nU_k = (3*n_tnn.*Tivals0/m + n_tnn.*uvals0.^2)/2; %k=0
Te_k = Tevals0; %k=0
x_k = [nu_k;nU_k;Te_k];

Te_k_right = [Te_k(2:Nx);Te_IC(200+dx/2)]; %at x_2, x_3, ..., x_N, x_{N+1} (ghost cell)
Te_k_left = [Te_IC(0-dx/2);Te_k(1:Nx-1)];  %at x_0 (ghost cell), x_1, ..., x_N
Rk_nu = nu_k - nvals0.*uvals0...
    + (delta*dt/dx)*(Shat0(2:end)-Shat0(1:end-1)) ...
    + ((1-delta)*dt/dx)*(Shat1(2:end)-Shat1(1:end-1))...
    - ((1-gamma)*dt/(2*dx))*(q/(m*qe))*(n_t1_right.*Te_1_right - n_t1_left.*Te_1_left)...
    - (gamma*dt/(2*dx))*(q/(m*qe))*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
Rk_nU = nU_k - (3*nvals0.*Tivals0/m + nvals0.*uvals0.^2)/2 ...
    + (delta*dt/dx)*(Qhat0(2:end)-Qhat0(1:end-1)) ...
    + ((1-delta)*dt/dx)*(Qhat1(2:end)-Qhat1(1:end-1))...
    - ((3*(1-gamma)*dt*sqrt(2)*sqrt(me))/m^2)*((n_t1.^2.*Te_1 - (m/3)*(2*n_t1.*nU_1-nu_1.^2))./(Te_1.^(1.5)))...
    - ((1-gamma)*dt/(2*dx))*(q/qe)*(nu_1./n_t1).*(n_t1_right.*Te_1_right - n_t1_left.*Te_1_left)...
    - ((3*gamma*dt*sqrt(2)*sqrt(me))/m^2)*((n_tnn.^2.*Te_k - (m/3)*(2*n_tnn.*nU_k-nu_k.^2))./(Te_k.^(1.5)))...
    - (gamma*dt/(2*dx))*(q/qe)*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
Rk_Te = n_tnn.*Te_k - nvals0.*Tevals0 ...
    + ((5*delta*dt)/(3*dx))*(uhat0(2:end).*nTehat0(2:end)-uhat0(1:end-1).*nTehat0(1:end-1)) ...
    + ((5*(1-delta)*dt)/(3*dx))*(uhat1(2:end).*nTehat1(2:end)-uhat1(1:end-1).*nTehat1(1:end-1)) ...
    - ((1-gamma)*dt/(3*dx))*(nu_1./n_t1).*(n_t1_right.*Te_1_right - n_t1_left.*Te_1_left)...
    - ((1-gamma)*dt/(3*dx^2))*(3.2/sqrt(2*me))*((Te_1.^(2.5)+Te_1_right.^(2.5)).*(Te_1_right-Te_1) - (Te_1_left.^(2.5)+Te_1.^(2.5)).*(Te_1-Te_1_left))...
    - ((2*(1-gamma)*dt*sqrt(2*me))/m)*((m/3)*(2*n_t1.*nU_1 - nu_1.^2) - n_t1.^2.*Te_1)./(Te_1.^(1.5))...
    - (gamma*dt/(3*dx))*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)...
    - (gamma*dt/(3*dx^2))*(3.2/sqrt(2*me))*((Te_k.^(2.5)+Te_k_right.^(2.5)).*(Te_k_right-Te_k) - (Te_k_left.^(2.5)+Te_k.^(2.5)).*(Te_k-Te_k_left))...
    - ((2*gamma*dt*sqrt(2*me))/m)*((m/3)*(2*n_tnn.*nU_k - nu_k.^2) - n_tnn.^2.*Te_k)./(Te_k.^(1.5));
Rk = -[Rk_nu;Rk_nU;Rk_Te];

tol = min(5e-12,max(abs(Rk))*5e-10);
%tol = 5.0e-4;

while err > tol
    Pk_nunu = eye(Nx,Nx);
    Pk_nunU = zeros(Nx,Nx);
    Pk_nuTe = zeros(Nx,Nx);
        Pk_nuTe(1:end-1,2:end) = Pk_nuTe(1:end-1,2:end) -(gamma*dt/(2*dx))*(q/(m*qe))*diag(n_tnn(2:end)); %superdiagonal
        Pk_nuTe(2:end,1:end-1) = Pk_nuTe(2:end,1:end-1) +(gamma*dt/(2*dx))*(q/(m*qe))*diag(n_tnn(1:end-1)); %subdiagonal

    Pk_nUnu = -((gamma*dt*q)/(2*dx*qe))*diag((n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)./n_tnn);
    Pk_nUnU = eye(Nx,Nx);
    Pk_nUTe = zeros(Nx,Nx);
        Pk_nUTe(1:end-1,2:end) = Pk_nUTe(1:end-1,2:end) -(gamma*dt/(2*dx))*(q/qe)*diag(nu_k(1:end-1).*n_tnn(2:end)./n_tnn(1:end-1)); %superdiagonal
        Pk_nUTe(2:end,1:end-1) = Pk_nUTe(2:end,1:end-1) +(gamma*dt/(2*dx))*(q/qe)*diag(nu_k(2:end).*n_tnn(1:end-1)./n_tnn(2:end)); %subdiagonal
        Pk_nUTe = Pk_nUTe - (3*gamma*dt*sqrt(2*me)/m^2)*diag( -0.5*(n_tnn.^2)./(Te_k.^(1.5)) - (m/3)*(-3/2)*(2*n_tnn.*nU_k-nu_k.^2)./(Te_k.^(2.5)) ); %diagonal

    Pk_Tenu = -(gamma*dt/(3*dx))*diag((n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)./n_tnn);
    Pk_TenU = zeros(Nx,Nx);
    Pk_TeTe = zeros(Nx,Nx);
        Pk_TeTe(1:end-1,2:end) = Pk_TeTe(1:end-1,2:end) -(gamma*dt/(3*dx))*diag(nu_k(1:end-1).*n_tnn(2:end)./n_tnn(1:end-1))...
            - (gamma*dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k(2:end).^(1.5).*(Te_k(2:end)-Te_k(1:end-1)) + (Te_k(1:end-1).^(2.5) + Te_k(2:end).^(2.5))); %superdiagonal
        Pk_TeTe(2:end,1:end-1) = Pk_TeTe(2:end,1:end-1) +(gamma*dt/(3*dx))*diag(nu_k(2:end).*n_tnn(1:end-1)./n_tnn(2:end))...
            + (gamma*dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k(1:end-1).^(1.5).*(Te_k(2:end)-Te_k(1:end-1)) - (Te_k(1:end-1).^(2.5) + Te_k(2:end).^(2.5))); %subdiagonal
        Pk_TeTe = Pk_TeTe + diag(n_tnn) - (gamma*dt/(3*dx^2))*(3.2/sqrt(2*me))*diag(2.5*Te_k.^(1.5).*(Te_k_right - 2*Te_k + Te_k_left) - (Te_k_right.^(2.5) + 2*Te_k.^(2.5) + Te_k_left.^(2.5)))...
            - (2*gamma*dt*sqrt(2*me)/m)*diag(-1.5*(m/3)*((2*n_tnn.*nU_k-nu_k.^2)./(Te_k.^(2.5))) + 0.5*(n_tnn.^2)./(Te_k.^(1.5))); %diagonal

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
    Rk_nu = nu_k - nvals0.*uvals0...
        + (delta*dt/dx)*(Shat0(2:end)-Shat0(1:end-1)) ...
        + ((1-delta)*dt/dx)*(Shat1(2:end)-Shat1(1:end-1))...
        - ((1-gamma)*dt/(2*dx))*(q/(m*qe))*(n_t1_right.*Te_1_right - n_t1_left.*Te_1_left)...
        - (gamma*dt/(2*dx))*(q/(m*qe))*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
    Rk_nU = nU_k - (3*nvals0.*Tivals0/m + nvals0.*uvals0.^2)/2 ...
        + (delta*dt/dx)*(Qhat0(2:end)-Qhat0(1:end-1)) ...
        + ((1-delta)*dt/dx)*(Qhat1(2:end)-Qhat1(1:end-1))...
        - ((3*(1-gamma)*dt*sqrt(2)*sqrt(me))/m^2)*((n_t1.^2.*Te_1 - (m/3)*(2*n_t1.*nU_1-nu_1.^2))./(Te_1.^(1.5)))...
        - ((1-gamma)*dt/(2*dx))*(q/qe)*(nu_1./n_t1).*(n_t1_right.*Te_1_right - n_t1_left.*Te_1_left)...
        - ((3*gamma*dt*sqrt(2)*sqrt(me))/m^2)*((n_tnn.^2.*Te_k - (m/3)*(2*n_tnn.*nU_k-nu_k.^2))./(Te_k.^(1.5)))...
        - (gamma*dt/(2*dx))*(q/qe)*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left);
    Rk_Te = n_tnn.*Te_k - nvals0.*Tevals0 ...
        + ((5*delta*dt)/(3*dx))*(uhat0(2:end).*nTehat0(2:end)-uhat0(1:end-1).*nTehat0(1:end-1)) ...
        + ((5*(1-delta)*dt)/(3*dx))*(uhat1(2:end).*nTehat1(2:end)-uhat1(1:end-1).*nTehat1(1:end-1)) ...
        - ((1-gamma)*dt/(3*dx))*(nu_1./n_t1).*(n_t1_right.*Te_1_right - n_t1_left.*Te_1_left)...
        - ((1-gamma)*dt/(3*dx^2))*(3.2/sqrt(2*me))*((Te_1.^(2.5)+Te_1_right.^(2.5)).*(Te_1_right-Te_1) - (Te_1_left.^(2.5)+Te_1.^(2.5)).*(Te_1-Te_1_left))...
        - ((2*(1-gamma)*dt*sqrt(2*me))/m)*((m/3)*(2*n_t1.*nU_1 - nu_1.^2) - n_t1.^2.*Te_1)./(Te_1.^(1.5))...
        - (gamma*dt/(3*dx))*(nu_k./n_tnn).*(n_tnn_right.*Te_k_right - n_tnn_left.*Te_k_left)...
        - (gamma*dt/(3*dx^2))*(3.2/sqrt(2*me))*((Te_k.^(2.5)+Te_k_right.^(2.5)).*(Te_k_right-Te_k) - (Te_k_left.^(2.5)+Te_k.^(2.5)).*(Te_k-Te_k_left))...
        - ((2*gamma*dt*sqrt(2*me))/m)*((m/3)*(2*n_tnn.*nU_k - nu_k.^2) - n_tnn.^2.*Te_k)./(Te_k.^(1.5));
    Rk = -[Rk_nu;Rk_nU;Rk_Te];

    err = max(abs(Rk));%norm(Rk,2)/Rn_err;

    k = k+1;
    if k==500
        disp('Completed 500 iterations - stage2');
        err
        return
    end    
end
%disp(k)
n_out = n_tnn;
u_out = nu_k./n_tnn;
Ti_out = (m/3)*(2*nU_k./n_tnn - (nu_k.^2)./(n_tnn.^2));
Te_out = Te_k;
end















