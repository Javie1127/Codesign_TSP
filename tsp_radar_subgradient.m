function [fdcomm, radar, cov] = tsp_radar_subgradient(fdcomm, radar_comm, radar, cov,k)
% UL subgradient method 
%
t_r_max = radar.t_r_max;
t_r = 1;
Mr = radar.TX;
eta_T = radar.channelgain;
Nr = radar.RX;
K = radar.codelength;
% fdcomm_temp = fdcomm; 
radar_temp = radar;
cov_temp = cov;
fdcomm = Xi_comm_k(fdcomm,k);
radar = Xi_radar(radar);
lambda_opt = radar.lambda(k);
ak_opt = radar.codematrix(k,:).';
% Xi_total = sum(fdcomm.Xi_UL) + sum(fdcomm.Xi_DL)+ sum(radar.Xi_r);
lambda_r_k_t = radar.lambda(k);
%P_rk_max = radar.max_power;
P_rk_max = sum(radar.power/K);
Xi_temp = zeros(t_r_max,1);
while t_r <= t_r_max
    % fdcomm tracks the optimal results
    % fdcomm_temp tracks the instantaneous lambda, Piu,
    beta_i_k_t = 1/t_r; % step number
    a_k_t = radar_temp.codematrix(k,:).';
    % update lambda
    lambda_r_k_new = lambda_r_k_t + beta_i_k_t*(real(a_k_t'*a_k_t)-P_rk_max);
    lambda_r_k_t = max(lambda_r_k_new,0);
    radar_temp.lambda(k) = lambda_r_k_t;
    % Update a_k with lambda_r_k_t
    [radar_temp,cov_temp] = radar_code(fdcomm, radar_temp, radar_comm, cov_temp,k);
    % Update covariance matrices 
    cov_temp = covmat(fdcomm,radar_temp,radar_comm); 
    % update Comm WMMSE receivers and weight matrices
    % fdcomm_temp = Comm_MMSE(fdcomm_temp, radar_temp, cov_temp);
    % radar_temp = radar_MMSE(radar_temp,cov_temp);
    for nr = 1:Nr
         eta_t = eta_T(nr);
         St_nr_t = cov_temp.S_tr(:,:,nr);
         R_in_nr = cov_temp.inr(:,:,nr);
         Urnr = radar_temp.WMMSE_RX{nr,1};
         Ernr = nearestSPD(eta_t^2*(eye(Mr)-2*nearestSPD(Urnr*(St_nr_t*St_nr_t')*Urnr'))+nearestSPD(Urnr*R_in_nr*Urnr'));
         radar_temp.MMSE_nop{nr,1} = Ernr;
     end
    % update Xi_MSE
    % fdcomm_temp = Xi_comm_k(fdcomm_temp, k);
    radar_temp = Xi_radar(radar_temp);
    Xi_temp_t = sum(fdcomm.Xi_UL) + sum(fdcomm.Xi_DL) + sum(radar_temp.Xi_r);
    Xi_temp(t_r) = Xi_temp_t;
    if t_r > 1 && Xi_temp_t < min(Xi_temp(1:t_r-1)) && any(radar_temp.Xi_r)>=0
        ak_temp = radar_temp.codematrix(k,:).';
        if real(trace(ak_temp'*ak_temp)) <= P_rk_max && real(trace(ak_temp'*ak_temp))>=0
            lambda_opt = lambda_r_k_t;
            ak_opt = ak_temp;
            radar.Xi_r = radar_temp.Xi_r;
        end
    end
    % calculate the new Xi_mse
    t_r = t_r+1;
end
radar.lambda(k) = lambda_opt;
radar.codematrix(k,:) = ak_opt.';
cov = covmat(fdcomm, radar, radar_comm);
% fdcomm = Comm_MMSE(fdcomm,radar,cov);
% radar = radar_MMSE(radar,cov);
end

