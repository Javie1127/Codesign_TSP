function [fdcomm_op, radar_op, cov_op] = tsp_UL_subgradient(fdcomm, radar_comm, radar, cov, ii,k)
% UL subgradient method 
% update P_ui within the sub-gradient method
tu_max = fdcomm.tu_max;
t = 1;
lambda_ui_k_0 = 10;
mu_ui_k_0 = 10;
lambda_ui_k_t = lambda_ui_k_0;
mu_ui_k_t = mu_ui_k_0;
% mu_ui_k_t = fdcomm.mu_UL(ii,k);
% lambda_ui_k_t = fdcomm.lambda_UL(ii,k);
%% Initialization 
fdcomm_temp = fdcomm; % to track 
radar_temp = radar;
cov_temp = cov;
fdcomm_op = fdcomm;
radar_op = radar;
cov_op = cov;
R_UL = fdcomm.R_UL; % UL minimum rate
%% Compute the initial Xi_mse

Xi_min = fdcomm_op.xi_UL(k)+fdcomm_op.xi_DL(k)+sum(radar_op.Xi_radar_Nr);
% lambda_iu_k_t = fdcomm.lambda_UL(ii,k);
P_U_i_max = fdcomm.UL_power(ii);
d_UL_i = fdcomm.ULstream_num(ii);
HiB = fdcomm.ULchannels{ii};
Xi = zeros(tu_max,1);
Xi_op = zeros(tu_max,1);
mu_i_k = zeros(tu_max,1);
lambda_ui_k = zeros(tu_max,1);
while t <= tu_max
    % fdcomm tracks the optimal results
    % fdcomm_temp tracks the instantaneous lambda, Piu,
    Xi_t = fdcomm_temp.xi_UL(k)+fdcomm_temp.xi_DL(k)+sum(radar_temp.Xi_radar_Nr);
    P_ui_k_t = fdcomm_temp.ULprecoders{ii,k}; 
    R_in_ui = cov_temp.in_UL{ii,k};
    R_ui_k_t = abs(log2(det((eye(d_UL_i)+P_ui_k_t'*HiB'/R_in_ui*HiB*P_ui_k_t))));
    g_t = abs(trace(P_ui_k_t*P_ui_k_t'))-P_U_i_max;
    q_t = R_UL-R_ui_k_t;
    switch fdcomm.step_size_rules
        case 'Square_summable'
            epsilon_i_k_t = 1/t;
            beta_i_k_t = 1/t;
        case 'Nonsummable_diminishing'
            epsilon_i_k_t = (1/sqrt(t))/(norm(q_t));
            beta_i_k_t = (1/t)/(norm(g_t)+1e-3*P_U_i_max);
        case 'Polyak'
            gamma = 0.5;
            epsilon_i_k_t = (Xi_t-Xi_min+gamma^t)/(norm(q_t))^2;
            beta_i_k_t = (Xi_t-Xi_min+gamma^t)/(norm(g_t))^2;
    end
    tilde_P_ui_k = fdcomm_op.ULprecoders{ii,k};
    
    %% update lambda and mu
    %epsilon_i_k_t = (1/t)/norm(P_U_i_max-abs(trace(P_ui_k_t*P_ui_k_t')));
    lambda_i_k_new = lambda_ui_k_t + beta_i_k_t*(g_t);
    %lambda_ui_k_t = max(min(lambda_i_k_new,lambda_ui_k_t),0);
    lambda_ui_k_t = max((lambda_i_k_new),0);
    mu_i_k_new = mu_ui_k_t+epsilon_i_k_t*(q_t); 
    mu_ui_k_t = max(mu_i_k_new,0);
    %mu_ui_k_t = max(min(mu_i_k_new,mu_ui_k_t),0);
    fdcomm_temp.mu_UL(ii,k) = mu_ui_k_t;
    fdcomm_temp.lambda_UL(ii,k) = lambda_ui_k_t;
    %% Update PiB with new mu_i_k_t and lambda_i_k_t
    fdcomm_temp = tsp_UL_precoders(k, ii, fdcomm_temp, cov_temp,tilde_P_ui_k);
    %%  update the MMSE matrix E_ui_k and Xi_MSE
    cov_temp = update_cov_ui_k(fdcomm_temp,radar_temp,cov_temp,radar_comm,P_ui_k_t,ii,k);
    [fdcomm_temp,radar_temp] = update_Xi_WMMSE(fdcomm_temp, radar_temp,cov_temp);
    Xi_t = fdcomm_temp.Xi_WMMSE_total;
    Xi(t) = Xi_t;
    lambda_ui_k(t) = lambda_ui_k_t;
    mu_i_k(t) = mu_ui_k_t;
    P_ui_k_temp =  fdcomm_temp.ULprecoders{ii,k};
    R_ui_k_temp = abs(log2(det(eye(d_UL_i)+nearestSPD(P_ui_k_temp'*HiB'/(R_in_ui)*HiB*P_ui_k_temp))));
    if Xi_t < Xi_min && (R_UL-R_ui_k_temp)<=0 && abs(trace(P_ui_k_temp*P_ui_k_temp')) <= P_U_i_max    
        Xi_min = Xi_t;
        Xi_op(t) = Xi_t;
        fdcomm_op = fdcomm_temp;
        radar_op = radar_temp;
        cov_op = cov_temp;
    else
        Xi_op(t) = Xi_min; 
    end
    t = t+1;
end

% cov = tsp_covmat_rev(fdcomm, radar, radar_comm,cov_temp);
% %% Update MSE matrix
% [fdcomm, radar] = tsp_mse(fdcomm, radar, cov, k);
% %% Update mean square error 
% fdcomm = Xi_comm_k(fdcomm, k);
% radar = Xi_radar(radar);
end

