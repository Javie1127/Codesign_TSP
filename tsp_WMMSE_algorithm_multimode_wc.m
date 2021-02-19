function [fdcomm, radar_temp] =...
    tsp_WMMSE_algorithm_multimode_wc(fdcomm, radar, radar_comm,cov, dmode)
%------------------Algorithm 3 (WMMMSE minimization)-------------------------------------------        
K = radar.codelength;
I = fdcomm.UL_num;
J = fdcomm.DL_num;
% Xi_min = sum(fdcomm.Xi_UL) + sum(fdcomm.Xi_DL) + sum(radar.Xi_r);
iota_max = radar.iota_max;
iota = 1;
% Xi_mse_op = zeros(iota_max,1);
% Xi_mse = zeros(iota_max,1);
while iota<= iota_max
    disp(iota);
    for k = 1: K
        %% UL precoder
        for ii = 1 : I
            %disp(['ii=',num2str(ii)]);
            [fdcomm] = tsp_UL_napla(k, ii, fdcomm, radar, radar_comm); %% modify the 
            [fdcomm,radar,cov] = tsp_UL_subgradient_wc(fdcomm, radar_comm, radar, cov, ii,k);
        end
        %% DL precoder
        for jj = 1 : J
            %disp(['jj=',num2str(ii)])
            [fdcomm] = tsp_DL_napla_wc(k, jj, fdcomm, radar, radar_comm, cov);
            [fdcomm,radar,cov] = tsp_DL_subgradient_wc(fdcomm, radar_comm, radar, cov, jj ,k);
        end
        %% radar code 
        if strcmp(dmode,'co_design')
            radar_temp = tsp_radar_code(fdcomm, radar, radar_comm, cov, k);
        else 
            radar_temp = radar;
        end
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
