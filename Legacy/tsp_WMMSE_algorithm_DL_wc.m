function [DLcomm, radar_temp] =...
    tsp_WMMSE_algorithm_DL_wc(DLcomm, radar, radar_comm,cov)
%------------------Algorithm 3 (WMMMSE minimization)-------------------------------------------        
K = radar.codelength;
J = DLcomm.DL_num;
% Xi_min = sum(fdcomm.Xi_UL) + sum(fdcomm.Xi_DL) + sum(radar.Xi_r);
iota_max = radar.iota_max;
iota = 1;
% Xi_mse_op = zeros(iota_max,1);
% Xi_mse = zeros(iota_max,1);
while iota<= iota_max
    disp(iota);
    for k = 1: K
        %% DL precoder
        for jj = 1 : J
            %disp(['jj=',num2str(ii)])
            [DLcomm] = tsp_DL_napla_DL_wc(k, jj, DLcomm, radar, radar_comm);
            [DLcomm,radar,cov] = tsp_DL_subgradient_DL_wc(DLcomm, radar_comm, radar, cov, jj ,k);
        end
        %% radar code 
        radar_temp = tsp_radar_code_DL_wc(DLcomm, radar, radar_comm, cov, k);
    end
    iota = iota + 1;
end
