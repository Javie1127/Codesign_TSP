function [fdcomm, radar, cov] = tsp_UL_subgradient_power(fdcomm, radar_comm, radar, cov, ii,k)
% UL subgradient method 
% update P_ui within the sub-gradient method
tu_max = fdcomm.tu_max;
t = 1;
lambda_ui_k_t = 1;
% mu_i_k_t = fdcomm.mu_UL(ii,k);
mu_ui_k_t = 1;
%% Initialization 
fdcomm_temp = fdcomm; % to track 
radar_temp = radar;
R_UL = fdcomm.R_UL; % UL minimum rate
%% Compute the initial Xi_mse
fdcomm_temp = Xi_comm_k(fdcomm_temp, k);
radar_temp = Xi_radar(radar_temp);
cov_temp = cov;
Xi_min = fdcomm_temp.Xi_UL(k) + fdcomm_temp.Xi_DL(k) + sum(radar_temp.Xi_r);
% Xi_k = fdcomm.Xi_UL(k)+fdcomm.Xi_DL(k) + sum(radar.Xi_r);
% lambda_iu_k_t = fdcomm.lambda_UL(ii,k);
P_U_i_max = fdcomm.ULpower(ii);
d_UL_i = fdcomm.ULstream_num(ii);
HiB = fdcomm.ULchannels{ii};
Xi = zeros(tu_max,1);
Xi_op = zeros(tu_max,1);
lambda_op = fdcomm.lambda_UL(ii,k);
mu_op = fdcomm.mu_UL(ii,k);
P_ui_k_op = fdcomm.ULprecoders{ii,k};
lambda_ui_k = zeros(tu_max,1);
mu_i_k = zeros(tu_max,1);
U_ui_k = fdcomm.UL_WMMSE_RX{ii,k};
while t <= tu_max
    % fdcomm tracks the optimal results
    % fdcomm_temp tracks the instantaneous lambda, Piu,
    beta_i_k_t = 1/t;
    epsilon_i_k_t = 1/t;
    P_ui_k_t = fdcomm_temp.ULprecoders{ii,k}; 
    R_in_ui = cov_temp.in_UL{ii,k};
    R_ui_k_t = real(log2(det((eye(d_UL_i)+P_ui_k_t'*HiB'/R_in_ui*HiB*P_ui_k_t))));
    %% update lambda and mu
    lambda_i_k_new = lambda_ui_k_t + beta_i_k_t*(real(trace(P_ui_k_t*P_ui_k_t'))-P_U_i_max);
    lambda_ui_k_t = max(lambda_i_k_new,0);
    mu_i_k_new = mu_ui_k_t+epsilon_i_k_t*(R_UL-R_ui_k_t); 
    mu_ui_k_t = max(mu_i_k_new,0);
    fdcomm_temp.mu_UL(ii,k) = mu_ui_k_t;
    fdcomm_temp.lambda_UL(ii,k) = lambda_ui_k_t;
    %% Update PiB with new mu_i_k_t and lambda_i_k_t
    fdcomm_temp = tsp_UL_precoders(k, ii, fdcomm_temp, cov_temp);
    %%  update the MMSE matrix E_ui_k
    cov_temp = tsp_covmat(fdcomm_temp, radar_temp, radar_comm);
    P_ui_k_temp = fdcomm_temp.ULprecoders{ii,k};
    R_UL_i_k = cov_temp.total_UL{ii,k};
    fdcomm_temp.UL_MMSE_nop{ii,k} = ...
        eye(fdcomm_temp.ULstream_num(ii))-2*U_ui_k*HiB*P_ui_k_temp + nearestSPD(U_ui_k*R_UL_i_k*U_ui_k');
    %% update Xi_MSE
    fdcomm_temp = Xi_comm_k(fdcomm_temp, k);
%     radar_temp = Xi_radar(radar_temp);
%     Xi_temp = fdcomm_temp.Xi_UL(k)+fdcomm_temp.Xi_DL(k) + sum(radar_temp.Xi_r);
    Xi_t = fdcomm_temp.Xi_UL(k)+fdcomm_temp.Xi_DL(k) + sum(radar_temp.Xi_r);
    Xi(t) = Xi_t;
    lambda_ui_k(t) = lambda_ui_k_t;
    mu_i_k(t) = mu_ui_k_t;
    R_in_ui_temp = cov_temp.in_UL{ii,k};
    R_ui_k_temp = real(log2(det(eye(d_UL_i)+(P_ui_k_temp'*HiB'/R_in_ui_temp*HiB*P_ui_k_temp))));
    if Xi_t < Xi_min && (R_UL-R_ui_k_temp)<=0 && real(trace(P_ui_k_temp*P_ui_k_temp')) <= P_U_i_max    
        lambda_op = lambda_ui_k_t;
        mu_op = mu_ui_k_t;
        P_ui_k_op = P_ui_k_temp;
        Xi_min = Xi_t;
        Xi_op(t) = Xi_t;
    else
        Xi_op(t) = Xi_min; 
    end
    t = t+1;
end

fdcomm.lambda_UL(ii,k) = lambda_op;
fdcomm.mu_UL(ii,k) = mu_op;
fdcomm.ULprecoders{ii,k} = P_ui_k_op;
cov = tsp_covmat(fdcomm, radar, radar_comm);
%% Update MSE matrix
[fdcomm, radar] = tsp_mse(fdcomm, radar, cov, k);
%% Update mean square error 
fdcomm = Xi_comm_k(fdcomm, k);
radar = Xi_radar(radar);
end

