function [fdcomm_op, radar_op, cov_op] = tsp_DL_subgradient(fdcomm, radar_comm, radar, cov, jj, k)
%------------------Algorithm 1 (WMMMSE minimization)-------------------------------------------        
%% initialization 
J = fdcomm.DL_num;
td_max = fdcomm.td_max;
P_B_max = fdcomm.BS_power;
t = 1;
R_DL = fdcomm.R_DL; % UL minimum rate
fdcomm_temp = fdcomm; 
radar_temp = radar;
cov_temp = cov;
fdcomm_op = fdcomm;
radar_op = radar;
cov_op = cov;
% lambda_d_k_t = fdcomm.lambda_DL(k);
% mu_dj_k_t = fdcomm.mu_DL(jj,k);
lambda_d_k_t = 10;
mu_dj_k_t = 50;
d_DL_j = fdcomm.DLstream_num(jj);
HBj = fdcomm.DLchannels{jj};
% R_in_dj = cov.in_DL{jj,k};
% R_DL_j_k = cov.total_DL{jj,k};
lambda_j_k = zeros(t,1);
mu_j_k = zeros(t,1);
Xi = zeros(td_max,1);
Xi_op = zeros(td_max,1);

% Xi_min = fdcomm_temp.Xi_UL(k) + fdcomm_temp.Xi_DL(k) + sum(radar_temp.Xi_r);
% U_dj_k = fdcomm_temp.DL_WMMSE_RX{jj,k};
Xi_min = fdcomm_temp.xi_UL(k)+fdcomm_temp.xi_DL(k)+sum(radar.Xi_radar_Nr);
P_dg_k_sum = 0;
for g = 1:J
    P_dg_k = fdcomm_temp.DLprecoders{g,k};
    if g ~= jj
        P_dg_k_sum = abs(trace((P_dg_k*P_dg_k'))) + P_dg_k_sum;
    end
end
while t <= td_max
    Xi_t = fdcomm_temp.Xi_WMMSE_total;
    tilde_P_dj_k = fdcomm_op.DLprecoders{jj,k};
    P_dj_k_t = fdcomm_temp.DLprecoders{jj,k};
    P_B_t = P_dg_k_sum + abs(trace((P_dj_k_t*P_dj_k_t')));
    g_t = P_B_t-P_B_max;
    R_in_dj = cov_temp.in_DL{jj,k};
    R_dj_k_t = abs(log2(det(eye(d_DL_j)+nearestSPD(P_dj_k_t'*HBj'/R_in_dj*HBj*P_dj_k_t))));
    q_t = R_DL-R_dj_k_t;
    switch fdcomm.step_size_rules
        case 'Square_summable'
            epsilon_j_k_t = 1/t;
            beta_d_k_t = 1/(t);
        case 'Nonsummable_diminishing'
            epsilon_j_k_t = (1/sqrt(t))/(norm(q_t));
            beta_d_k_t = (1/sqrt(t))/(norm(g_t)+1e-3*P_B_max);
        case 'Polyak'
            gamma = 0.5;
            epsilon_j_k_t = (Xi_t-Xi_min+gamma^t)/(norm(q_t))^2;
            beta_d_k_t = (Xi_t-Xi_min+gamma^t)/(norm(g_t))^2;        
    end
    % update lambda_k_t
    
    lambda_k_d_temp = lambda_d_k_t + beta_d_k_t*(P_B_t-P_B_max);
    %lambda_k_d_temp = max(0,min(lambda_k_d_temp,lambda_d_k_t));
    lambda_d_k_t = max(0,lambda_k_d_temp);
    fdcomm_temp.lambda_DL(k) = lambda_d_k_t;
    %% update mu_j_k_t
    
    mu_jd_k_temp = mu_dj_k_t + epsilon_j_k_t*q_t;
    mu_dj_k_t = max(mu_jd_k_temp,0);
    fdcomm_temp.mu_DL(jj,k) = mu_dj_k_t;
    %% update P_jd_k
    fdcomm_temp = tsp_DL_precoders(k, jj, fdcomm_temp, cov_temp,tilde_P_dj_k);
%     fdcomm_temp = Comm_MMSE(fdcomm_temp,radar,cov_temp); 
%     radar_temp = radar_MMSE(radar_temp,cov_temp);
    %% Update the MSE matrix E_dj_k_t
%     cov_temp = tsp_covmat_rev(fdcomm_temp, radar, radar_comm,cov);
%     P_dj_k_t = fdcomm_temp.DLprecoders{jj,k};
%     R_DL_j_k = cov_temp.total_DL{jj,k};
%     E_dj_k_t = eye(fdcomm.DLstream_num(jj))-2*U_dj_k*HBj*P_dj_k_t+(U_dj_k*R_DL_j_k*U_dj_k');
%     fdcomm_temp.DL_MMSE_nop{jj,k} = E_dj_k_t;
    cov_temp = update_cov_dj_k(fdcomm_temp,radar_temp,cov_temp,radar_comm,P_dj_k_t,jj,k);
    [fdcomm_temp,radar_temp] = update_Xi_WMMSE_k(fdcomm_temp, radar_temp,cov_temp,k);
    Xi_t = fdcomm_temp.Xi_WMMSE_total;
    %% update Xi_MSE
%     fdcomm_temp = Xi_comm_k(fdcomm_temp, k);
% %     radar_temp = Xi_radar(radar_temp);
%     Xi_t = fdcomm_temp.Xi_UL(k)+ fdcomm_temp.Xi_DL(k) + sum(radar_temp.Xi_r);
    Xi(t) = Xi_t;
    lambda_j_k(t) = lambda_d_k_t;
    mu_j_k(t) = mu_dj_k_t;
    P_dj_k_temp = fdcomm_temp.DLprecoders{jj,k};
    power_temp = abs(trace((P_dj_k_temp*P_dj_k_temp'))) + P_dg_k_sum;
    R_dj_k_t_temp = abs(log2(det(eye(d_DL_j)+(P_dj_k_temp'*HBj'/R_in_dj*HBj*P_dj_k_temp))));
    if Xi_t < Xi_min && power_temp <= P_B_max && R_dj_k_t_temp >= R_DL
        Xi_min = Xi_t; 
        Xi_op(t) = Xi_t; 
        fdcomm_op = fdcomm_temp;
        radar_op = radar_temp;
        cov_op = cov_temp;
    else
        Xi_op(t) = Xi_min;
    end
    if Xi_t>10^2*Xi_min
        break
    end
    t = t+1;
end
% cov = tsp_covmat_rev(fdcomm, radar, radar_comm,cov_temp);
% %% Update MSE matrix
% [fdcomm, radar] = tsp_mse(fdcomm, radar, cov, k);
% %% Update Xi_mse
% fdcomm = Xi_comm_k(fdcomm, k);
% radar = Xi_radar(radar);
end

