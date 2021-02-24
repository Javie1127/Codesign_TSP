function [DLcomm, radar_DL, cov_DL] = tsp_DL_subgradient_DL_wc(DLcomm, radar_comm_DL, radar_DL, cov_DL, jj, k)
%------------------Algorithm 1 (WMMMSE minimization)-------------------------------------------        
%% initialization 
J = DLcomm.DL_num;
td_max = DLcomm.td_max;
P_B_max = DLcomm.DLpower;
t = 1;
R_DL = DLcomm.R_DL; % UL minimum rate
DLcomm_temp = DLcomm; 
cov_DL_temp = cov_DL;
% radar_temp = radar;
% Xi_k = fdcomm.Xi_UL(k)+fdcomm.Xi_DL(k) + sum(radar.Xi_r);
% lambda_d_k_t = fdcomm.lambda_DL(k);
% mu_dj_k_t = fdcomm.mu_DL(jj,k);
lambda_d_k_t = 1;
mu_dj_k_t = 1;
d_DL_j = DLcomm.DLstream_num(jj);
HBj = DLcomm.DLchannels{jj};
% R_in_dj = cov.in_DL{jj,k};
% R_DL_j_k = cov.total_DL{jj,k};
lambda_opt = 1;
mu_op = 1;
P_dj_k_op = DLcomm.DLprecoders{jj,k}; 
lambda_j_k = zeros(t,1);
mu_j_k = zeros(t,1);
Xi = zeros(td_max,1);
Xi_op = zeros(td_max,1);
DLcomm = DLcomm_temp;
radar_temp = radar_DL;
Xi_min = DLcomm_temp.Xi_DL(k) + sum(radar_temp.Xi_r);
U_dj_k = DLcomm_temp.DL_WMMSE_RX{jj,k};
while t <= td_max
    beta_d_k_t = 1/t;
    epsilon_j_k_t = 1/sqrt(t);
    % update lambda_k_t
    P_dg_k_sum = 0;
    for g = 1:J
        P_dj_k = DLcomm_temp.DLprecoders{g,k};
        P_dg_k_sum = real(trace(P_dj_k*P_dj_k')) + P_dg_k_sum;
    end
    lambda_k_d_temp = lambda_d_k_t + beta_d_k_t*(P_dg_k_sum-P_B_max);
    lambda_k_d_temp = max(0,lambda_k_d_temp);
    DLcomm_temp.lambda_DL(k) = lambda_k_d_temp;
    %% update mu_j_k_t
%     E_jd_k_t = fdcomm.DL_MMSE{jj,k};
%     R_jd_k_t = log2(det((E_jd_k_t)^(-1)));
    R_in_dj = cov_DL_temp.in_DL{jj,k};
    P_dj_k_t = DLcomm.DLprecoders{jj,k};
    R_dj_k_t = real(log2(det(eye(d_DL_j)+(P_dj_k_t'*HBj'/R_in_dj*HBj*P_dj_k_t))));
    mu_jd_k_temp = mu_dj_k_t + epsilon_j_k_t*(R_DL-R_dj_k_t);
    mu_dj_k_t = max(mu_jd_k_temp,0);
    DLcomm_temp.mu_DL(jj,k) = mu_dj_k_t;
    %% update P_jd_k
    DLcomm_temp = tsp_DL_precoders_DL(k, jj, DLcomm,cov_DL);
    %% Update the MSE matrix E_dj_k_t
    cov_DL_temp = tsp_covmat_DL_wc(DLcomm_temp, radar_DL, radar_comm_DL);
    P_dj_k_t = DLcomm_temp.DLprecoders{jj,k};
    R_DL_j_k = cov_DL_temp.total_DL{jj,k};
    E_dj_k_t = eye(DLcomm.DLstream_num(jj))-2*U_dj_k*HBj*P_dj_k_t+(U_dj_k*R_DL_j_k*U_dj_k');
    DLcomm_temp.DL_MMSE_nop{jj,k} = E_dj_k_t;
    %% update Xi_MSE
    DLcomm_temp = Xi_comm_k_DL(DLcomm_temp, k);
%     radar_temp = Xi_radar(radar_temp);
    Xi_t = DLcomm_temp.Xi_DL(k) + sum(radar_temp.Xi_r);
    Xi(t) = Xi_t;
    lambda_j_k(t) = lambda_k_d_temp;
    mu_j_k(t) = mu_dj_k_t;
    power_temp = real(trace((P_dj_k_t*P_dj_k_t'))) + P_dg_k_sum;
    R_in_dj_temp = cov_DL_temp.in_DL{jj,k};
    R_dj_k_t_temp = real(log2(det(eye(d_DL_j)+(P_dj_k_t'*HBj'/R_in_dj_temp*HBj*P_dj_k_t))));
    if Xi_t < Xi_min && power_temp <= P_B_max && R_dj_k_t_temp >= R_DL
        lambda_opt = lambda_d_k_t;
        mu_op = mu_dj_k_t;
        Xi_min = Xi_t; 
        Xi_op(t) = Xi_t; 
        P_dj_k_op = P_dj_k_t;
        DLcomm.Xi_UL(k) = DLcomm_temp.Xi_UL(k);
        DLcomm.Xi_DL(k) = DLcomm_temp.Xi_DL(k);
    else
        Xi_op(t) = Xi_min;
    end
    t = t+1;
end
DLcomm.lambda_DL(jj,k) = lambda_opt;
DLcomm.mu_DL(jj,k) = mu_op;
DLcomm.DLprecoders{jj,k} = P_dj_k_op;
cov_DL = tsp_covmat_DL_wc(DLcomm, radar_DL, radar_comm_DL);
%% Update MSE matrix
[DLcomm, radar_DL] = tsp_mse_DL(DLcomm, radar_DL, cov_DL, k);
%% Update Xi_mse
DLcomm = Xi_comm_k_DL(DLcomm, k);
radar_DL = Xi_radar(radar_DL);
end

