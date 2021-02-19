function [DLcomm, radar_temp_DL] =...
    tsp_WMMSE_algorithm_DL(DLcomm, radar_DL, radar_comm_DL,cov_DL)
%------------------Algorithm 3 (WMMMSE minimization)-------------------------------------------        
K = radar_DL.codelength;
J = DLcomm.DL_num;
% Xi_min = sum(fdcomm.Xi_UL) + sum(fdcomm.Xi_DL) + sum(radar.Xi_r);
iota_max = radar_DL.iota_max;
iota = 1;
% Xi_mse_op = zeros(iota_max,1);
% Xi_mse = zeros(iota_max,1);
while iota<= iota_max
    disp(iota);
    for k = 1: K
        %% DL precoder
        for jj = 1 : J
            %disp(['jj=',num2str(ii)])
            [DLcomm] = tsp_DL_napla_DL(k, jj, DLcomm, radar_DL, radar_comm_DL, cov_DL);
            [DLcomm,radar_DL,cov_DL] = tsp_DL_subgradient_DL(DLcomm, radar_comm_DL, radar_DL, cov_DL, jj ,k);
        end
        %% radar code 
        radar_temp_DL = tsp_radar_code_DL(DLcomm, radar_DL, radar_comm_DL, cov_DL, k);
    end
%     %% Update the WMMSE receivers and weight matrices
%     cov = tsp_covmat(fdcomm, radar, radar_comm); % update the covariance matrices 
%     fdcomm = tsp_Comm_MMSE(fdcomm, radar, cov);
%     % radar = tsp_radar_MMSE(radar, cov);
%     for k = 1:K
%         fdcomm = Xi_comm_k(fdcomm, k);
%         Xi_mse(iota) = fdcomm.Xi_UL(k)+ fdcomm.Xi_DL(k) + Xi_mse(iota);
%     end
%     % radar = Xi_radar(radar);
%     Xi_mse(iota) = Xi_mse(iota) + sum(radar.Xi_r);
%     if Xi_min < Xi_mse(iota)
% %         fdcomm = fdcomm_op;
% %         radar = radar_op;
%         Xi_mse_op(iota) = Xi_min;
%     else
%         Xi_min = Xi_mse(iota);
%         fdcomm_op = fdcomm;
%         radar_op = radar;
%         Xi_mse_op(iota) = Xi_mse(iota);
%     end
%     % Update optimal covariance matrices, fdcomm, radar 
%     cov_op = tsp_covmat(fdcomm_op,radar_op, radar_comm);
%     fdcomm_op = tsp_Comm_MMSE(fdcomm_op, radar_op, cov_op);
%     radar_op = tsp_radar_MMSE(radar_op, cov_op);
    iota = iota + 1;
end
