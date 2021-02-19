function [fdcomm_op, radar_op,cov_op] = tsp_altermating_projection...
(fdcomm, radar,radar_comm, ini)
ell_max = radar.ell_max;
ell = 1;
%% Initialization
if nargin == 3 || strcmp(ini, 'normal_ini')
    [fdcomm_ini, radar_ini, cov_ini] = tsp_ini_normal(radar,fdcomm,radar_comm);
elseif strcmp (ini,'random_ini')
    [fdcomm_ini, radar_ini, cov_ini] = tsp_ini_random(radar,fdcomm, radar_comm);
end
radar_op = radar_ini;
fdcomm_op = fdcomm_ini;
radar = radar_ini;
fdcomm = fdcomm_ini;
cov = cov_ini;
I_total = zeros(ell_max,1);
I_total_op = zeros(ell_max,1);
I_max = sum(fdcomm.MI_UL(:))+sum(fdcomm.MI_DL(:))+ sum(radar.MI_radar);
while ell <= ell_max
    [fdcomm, radar] =...
        tsp_WMMSE_algorithm(fdcomm, radar, radar_comm, cov);
    % Calculate A^star
    radar = tsp_Nearest_PAR(radar);
    % Update covariance matrices
    cov = tsp_covmat(fdcomm,radar,radar_comm);
    % Update linear receivers
    fdcomm = tsp_Comm_MMSE(fdcomm, radar, cov);
    radar = tsp_radar_MMSE(radar, cov);
    I_total(ell) = sum(fdcomm.MI_UL(:))+sum(fdcomm.MI_DL(:))+ sum(radar.MI_radar);
    if I_total(ell) > I_max
        I_total_op(ell) = I_total(ell);
        I_max = I_total(ell);
        radar_op = radar;
        fdcomm_op = fdcomm;
    else
        I_total_op(ell) = I_max;
    end
    ell = ell+1;
end
fdcomm_op.I_total = I_total;
fdcomm_op.I_total_op = I_total_op;
cov_op = tsp_covmat(fdcomm_op,radar_op,radar_comm);
end

