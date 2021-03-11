 function [fdcomm, radar_temp,cov] =...
    tsp_WMMSE_MRMC(fdcomm, radar, radar_comm,cov)
%------------------Algorithm 3 (WMMMSE minimization)-------------------------------------------        
K = radar.codelength;
I = fdcomm.UL_num;
J = fdcomm.DL_num;
% Xi_min = sum(fdcomm.Xi_UL) + sum(fdcomm.Xi_DL) + sum(radar.Xi_r);
iota_max = radar.iota_max;
iota = 1;
% Xi_mse_op = zeros(iota_max,1);
% Xi_mse = zeros(iota_max,1);
[fdcomm,radar] = update_Xi_WMMSE(fdcomm,radar,cov);
while iota<= iota_max
    %disp(iota);
    for k = 1: K
        if ~strcmp(fdcomm.precoder_type,'Uniform')
            %% UL precoder
            for ii = 1 : I
                %disp(['ii=',num2str(ii)]);
                [fdcomm] = tsp_UL_napla_rev(k, ii, fdcomm, radar, radar_comm); 
                [fdcomm,radar,cov] = tsp_UL_subgradient(fdcomm, radar_comm, radar, cov, ii,k);
            end
        end
        if strcmp(fdcomm.precoder_type,'Proposed') || strcmp(fdcomm.precoder_type,'Uniform')
            %% DL precoder
            for jj = 1 : J
                %disp(['jj=',num2str(ii)])
                [fdcomm] = tsp_DL_napla(k, jj, fdcomm, radar, radar_comm);
                [fdcomm,radar,cov] = tsp_DL_subgradient(fdcomm, radar_comm, radar, cov, jj ,k);
            end
        end
        if strcmp(radar.coding_type,'Proposed')
            %% radar code 
            radar_temp = tsp_radar_code(fdcomm, radar, radar_comm, cov, k);
            cov = up_cov_radar(fdcomm,radar_temp,cov,radar_comm);
        else
            radar_temp = radar;
        end
        
    end
    iota = iota + 1;
end
