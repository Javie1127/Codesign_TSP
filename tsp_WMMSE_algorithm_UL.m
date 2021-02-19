function [ULcomm, radar_temp_UL] =...
    tsp_WMMSE_algorithm_UL(ULcomm, radar_UL, radar_comm_UL,cov_UL)
%------------------Algorithm 3 (WMMMSE minimization)-------------------------------------------        
K = radar_UL.codelength;
I = ULcomm.UL_num;
% Xi_min = sum(fdcomm.Xi_UL) + sum(fdcomm.Xi_DL) + sum(radar.Xi_r);
iota_max = radar_UL.iota_max;
iota = 1;
% Xi_mse_op = zeros(iota_max,1);
% Xi_mse = zeros(iota_max,1);
while iota<= iota_max
    disp(iota);
    for k = 1: K
        %% UL precoder
        for ii = 1 : I
            %disp(['ii=',num2str(ii)]);
            [ULcomm] = tsp_UL_napla_UL(k, ii, ULcomm, radar_UL, radar_comm_UL); %% modify the 
            [ULcomm,radar_UL,cov_UL] = tsp_UL_subgradient_UL(ULcomm, radar_comm_UL, radar_UL, cov_UL, ii,k);
        end
        %% radar code 
        radar_temp_UL = tsp_radar_code_UL(ULcomm, radar_UL, radar_comm_UL, cov_UL, k);
    end
    iota = iota + 1;
end
