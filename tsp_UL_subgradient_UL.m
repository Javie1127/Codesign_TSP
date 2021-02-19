function [ULcomm, radar_UL, cov_UL] = tsp_UL_subgradient_UL(ULcomm, radar_comm_UL, radar_UL, cov_UL, ii,k)
% UL subgradient method 
% update P_ui within the sub-gradient method
tu_max = ULcomm.tu_max;
t = 1;
lambda_ui_k_t = 1;
% mu_i_k_t = fdcomm.mu_UL(ii,k);
mu_ui_k_t = 1;
%% Initialization 
ULcomm_temp = ULcomm; % to track 
radar_UL_temp = radar_UL;
R_UL = ULcomm.R_UL; % UL minimum rate
%% Compute the initial Xi_mse
ULcomm_temp = Xi_comm_k_UL(ULcomm_temp, k);
radar_UL_temp = Xi_radar(radar_UL_temp);
cov_temp = cov_UL;
Xi_min = ULcomm_temp.Xi_UL(k) + sum(radar_UL_temp.Xi_r);
% Xi_k = fdcomm.Xi_UL(k)+fdcomm.Xi_DL(k) + sum(radar.Xi_r);
% lambda_iu_k_t = fdcomm.lambda_UL(ii,k);
P_U_i_max = ULcomm.ULpower(ii);
d_UL_i = ULcomm.ULstream_num(ii);
HiB = ULcomm.ULchannels{ii};
Xi = zeros(tu_max,1);
Xi_op = zeros(tu_max,1);
lambda_op = ULcomm.lambda_UL(ii,k);
mu_op = ULcomm.mu_UL(ii,k);
P_ui_k_op = ULcomm.ULprecoders{ii,k};
lambda_ui_k = zeros(tu_max,1);
mu_i_k = zeros(tu_max,1);
U_ui_k = ULcomm.UL_WMMSE_RX{ii,k};
while t <= tu_max
    % fdcomm tracks the optimal results
    % fdcomm_temp tracks the instantaneous lambda, Piu,
    beta_i_k_t = 1/t;
    epsilon_i_k_t = 1/t;
    P_ui_k_t = ULcomm_temp.ULprecoders{ii,k}; 
    R_in_ui = cov_temp.in_UL{ii,k};
    R_ui_k_t = real(log2(det((eye(d_UL_i)+nearestSPD(P_ui_k_t'*HiB'/R_in_ui*HiB*P_ui_k_t)))));
    %% update lambda and mu
    lambda_i_k_new = lambda_ui_k_t + beta_i_k_t*(real(trace(P_ui_k_t*P_ui_k_t'))-P_U_i_max);
    lambda_ui_k_t = max(lambda_i_k_new,0);
    mu_i_k_new = mu_ui_k_t+epsilon_i_k_t*(R_UL-R_ui_k_t); 
    mu_ui_k_t = max(mu_i_k_new,0);
    ULcomm_temp.mu_UL(ii,k) = mu_ui_k_t;
    ULcomm_temp.lambda_UL(ii,k) = lambda_ui_k_t;
    %% Update PiB with new mu_i_k_t and lambda_i_k_t
    ULcomm_temp = tsp_UL_precoders_UL(k, ii, ULcomm_temp, cov_temp);
    %%  update the MMSE matrix E_ui_k
    cov_temp = tsp_covmat_UL(ULcomm_temp, radar_UL_temp, radar_comm_UL);
    P_ui_k_temp = ULcomm_temp.ULprecoders{ii,k};
    R_UL_i_k = cov_temp.total_UL{ii,k};
    ULcomm_temp.UL_MMSE_nop{ii,k} = ...
        eye(ULcomm_temp.ULstream_num(ii))-2*U_ui_k*HiB*P_ui_k_temp + (U_ui_k*R_UL_i_k*U_ui_k');
    %% update Xi_MSE
    ULcomm_temp = Xi_comm_k_UL(ULcomm_temp, k);
%     radar_temp = Xi_radar(radar_temp);
%     Xi_temp = fdcomm_temp.Xi_UL(k)+fdcomm_temp.Xi_DL(k) + sum(radar_temp.Xi_r);
    Xi_t = ULcomm_temp.Xi_UL(k) + sum(radar_UL_temp.Xi_r);
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

ULcomm.lambda_UL(ii,k) = lambda_op;
ULcomm.mu_UL(ii,k) = mu_op;
ULcomm.ULprecoders{ii,k} = P_ui_k_op;
cov_UL = tsp_covmat_UL(ULcomm, radar_UL, radar_comm_UL);
%% Update MSE matrix
[ULcomm, radar_UL] = tsp_mse_UL(ULcomm, radar_UL, cov_UL, k);
%% Update mean square error 
ULcomm = Xi_comm_k_UL(ULcomm, k);
radar_UL = Xi_radar(radar_UL);
end

