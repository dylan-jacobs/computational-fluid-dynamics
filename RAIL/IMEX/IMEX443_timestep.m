% IMEX(4, 4, 3) timestep

function [Vx4, S4, Vy4, r1] = IMEX443_timestep(Vx0, S0, Vy0, tn, dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance)

    N = size(Dxx);

    Y0_hat = -(((Dx*(A{1, 1}(tn).*Vx0))*(A{1, 2}(tn).*S0)*((A{1, 3}(tn).*Vy0)')) + ((A{2, 1}(tn).*Vx0)*(A{2, 2}(tn).*S0)*((Dy*(A{2, 3}(tn).*Vy0))')));
    
    % ------------- Stage 1: (IMEX111) -------------
    [Vx1, S1, Vy1, ~] = IMEX111_timestep(Vx0, S0, Vy0, tn, (1/2)*dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance);

    % Compute/store Y1, Y1_hat for next step
    Y1_hat = -(((Dx*(A{1, 1}(tn+((1/2)*dt)).*Vx1))*(A{1, 2}(tn+((1/2)*dt)).*S1)*((A{1, 3}(tn+((1/2)*dt)).*Vy1)')) + ((A{2, 1}(tn+((1/2)*dt)).*Vx1)*(A{2, 2}(tn+((1/2)*dt)).*S1)*((Dy*(A{2, 3}(tn+((1/2)*dt)).*Vy1))')));
    Y1 = ((((Dxx*Vx1*S1*(Vy1')) + (Vx1*S1*((Dyy*Vy1)'))))) + (phi{1, 1}(tn+((1/2)*dt))*phi{1, 2}(tn+((1/2)*dt))*(phi{1, 3}(tn+((1/2)*dt))'));


    % ------------- Stage 2: -------------

    % Reduced Augmentation
    % Predict V_dagger using IMEX111 for second stage
    [Vx2_dagger, ~, Vy2_dagger, ~] = IMEX111_timestep(Vx0, S0, Vy0, tn, (2/3)*dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance);
    [Vx_star2, Vy_star2] = reduced_augmentation([Vx2_dagger, Vx1, Vx0], [Vy2_dagger, Vy1, Vy0]);

    W1 = (Vx0*S0*(Vy0')) + (dt*(1/2)*(phi{1, 1}(tn+((2/3)*dt))*phi{1, 2}(tn+((2/3)*dt))*(phi{1, 3}(tn+((2/3)*dt))'))) + (dt*( ((11/18)*Y0_hat) + ((1/18)*Y1_hat) + ((1/6)*Y1) ));

    % K/L-Step
    K2 = sylvester(eye(N) - ((1/2)*dt*Dxx), -(1/2)*dt*(Dyy*Vy_star2)'*Vy_star2, W1*Vy_star2);
    L2 = sylvester(eye(N) - ((1/2)*dt*Dyy), -(1/2)*dt*(Dxx*Vx_star2)'*Vx_star2, (W1')*Vx_star2);

    % Get bases
    [Vx_ddagger2, ~] = qr(K2, 0); [Vy_ddagger2, ~] = qr(L2, 0);

    % Reduced Augmentation
    [Vx2, Vy2] = reduced_augmentation([Vx_ddagger2, Vx1, Vx0], [Vy_ddagger2, Vy1, Vy0]);

    % S-Step
    S2 = sylvester(eye(size(Vx2, 2)) - ((1/2)*dt*(Vx2')*Dxx*Vx2), -(1/2)*dt*(Dyy*Vy2)'*Vy2, (Vx2')*W1*Vy2);
    [Vx2, S2, Vy2, ~] = truncate_svd(Vx2, S2, Vy2, tolerance);

    % Compute/store Y2, Y2_hat for next step
    Y2_hat = -(((Dx*(A{1, 1}(tn+((2/3)*dt)).*Vx2))*(A{1, 2}(tn+((2/3)*dt)).*S2)*((A{1, 3}(tn+((2/3)*dt)).*Vy2)')) + ((A{2, 1}(tn+((2/3)*dt)).*Vx2)*(A{2, 2}(tn+((2/3)*dt)).*S2)*((Dy*(A{2, 3}(tn+((2/3)*dt)).*Vy2))')));
    Y2 = ((((Dxx*Vx2*S2*(Vy2')) + (Vx2*S2*((Dyy*Vy2)'))))) + (phi{1, 1}(tn+((2/3)*dt))*phi{1, 2}(tn+((2/3)*dt))*(phi{1, 3}(tn+((2/3)*dt))'));


    % ------------- Stage 3: -------------

    % Reduced Augmentation
    % Predict V_dagger using IMEX111 for second stage
    [Vx3_dagger, ~, Vy3_dagger, ~] = IMEX111_timestep(Vx0, S0, Vy0, tn, (1/2)*dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance);
    [Vx_star3, Vy_star3] = reduced_augmentation([Vx3_dagger, Vx2, Vx1, Vx0], [Vy3_dagger, Vy2, Vy1, Vy0]);

    W2 = (Vx0*S0*(Vy0')) + (dt*(1/2)*(phi{1, 1}(tn+((1/2)*dt))*phi{1, 2}(tn+((1/2)*dt))*(phi{1, 3}(tn+((1/2)*dt))'))) + (dt*( ((5/6)*Y0_hat) + ((-5/6)*Y1_hat) + ((1/2)*Y2_hat) + ((-1/2)*Y1) + ((1/2)*Y2)));

    % K/L-Step
    K3 = sylvester(eye(N) - ((1/2)*dt*Dxx), -(1/2)*dt*(Dyy*Vy_star3)'*Vy_star3, W2*Vy_star3);
    L3 = sylvester(eye(N) - ((1/2)*dt*Dyy), -(1/2)*dt*(Dxx*Vx_star3)'*Vx_star3, (W2')*Vx_star3);

    % Get bases
    [Vx_ddagger3, ~] = qr(K3, 0); [Vy_ddagger3, ~] = qr(L3, 0);

    % Reduced Augmentation
    [Vx3, Vy3] = reduced_augmentation([Vx_ddagger3, Vx2, Vx1, Vx0], [Vy_ddagger3, Vy2, Vy1, Vy0]);

    % S-Step
    S3 = sylvester(eye(size(Vx3, 2)) - ((1/2)*dt*(Vx3')*Dxx*Vx3), -(1/2)*dt*(Dyy*Vy3)'*Vy3, (Vx3')*W2*Vy3);
    [Vx3, S3, Vy3, ~] = truncate_svd(Vx3, S3, Vy3, tolerance);

    % Compute/store Y2, Y2_hat for next step
    Y3_hat = -(((Dx*(A{1, 1}(tn+((1/2)*dt)).*Vx3))*(A{1, 2}(tn+((1/2)*dt)).*S3)*((A{1, 3}(tn+((1/2)*dt)).*Vy3)')) + ((A{2, 1}(tn+((1/2)*dt)).*Vx3)*(A{2, 2}(tn+((1/2)*dt)).*S3)*((Dy*(A{2, 3}(tn+((1/2)*dt)).*Vy3))')));
    Y3 = ((((Dxx*Vx3*S3*(Vy3')) + (Vx3*S3*((Dyy*Vy3)'))))) + (phi{1, 1}(tn+((1/2)*dt))*phi{1, 2}(tn+((1/2)*dt))*(phi{1, 3}(tn+((1/2)*dt))'));

    
    
    
    % ------------- Stage 4: -------------

    % Reduced Augmentation
    % Predict V_dagger using IMEX111 for second stage
    [Vx4_dagger, ~, Vy4_dagger, ~] = IMEX111_timestep(Vx0, S0, Vy0, tn, dt, A, Dx, Dy, Dxx, Dyy, phi, tolerance);
    [Vx_star4, Vy_star4] = reduced_augmentation([Vx4_dagger, Vx3, Vx2, Vx1, Vx0], [Vy4_dagger, Vy3, Vy2, Vy1, Vy0]);

    W3 = (Vx0*S0*(Vy0')) + ((1/2)*dt*(phi{1, 1}(tn+(dt))*phi{1, 2}(tn+(dt))*(phi{1, 3}(tn+(dt))'))) + (dt*( ((1/4)*Y0_hat) + ((7/4)*Y1_hat) + ((3/4)*Y2_hat) + ((-7/4)*Y3_hat) + ((3/2)*Y1) + ((-3/2)*Y2) + ((1/2)*Y3) ) );

    % K/L-Step
    K4 = sylvester(eye(N) - ((1/2)*dt*Dxx), -(1/2)*dt*(Dyy*Vy_star4)'*Vy_star4, W3*Vy_star4);
    L4 = sylvester(eye(N) - ((1/2)*dt*Dyy), -(1/2)*dt*(Dxx*Vx_star4)'*Vx_star4, (W3')*Vx_star4);

    % Get bases
    [Vx_ddagger4, ~] = qr(K4, 0); [Vy_ddagger4, ~] = qr(L4, 0);

    % Reduced Augmentation
    [Vx4, Vy4] = reduced_augmentation([Vx_ddagger4, Vx3, Vx2, Vx1, Vx0], [Vy_ddagger4, Vy3,  Vy2, Vy1, Vy0]);

    % S-Step
    S4 = sylvester(eye(size(Vx4, 2)) - ((1/2)*dt*(Vx4')*Dxx*Vx4), -(1/2)*dt*(Dyy*Vy4)'*Vy4, (Vx4')*W3*Vy4);
    [Vx4, S4, Vy4, r1] = truncate_svd(Vx4, S4, Vy4, tolerance);
end