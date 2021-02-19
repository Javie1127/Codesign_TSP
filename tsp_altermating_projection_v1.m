function [fdcomm_op, radar_op,cov_op, radar_ini] = tsp_altermating_projection_v1...
(fdcomm_ini, radar_ini,radar_comm,cov_ini)
ell_max = radar_ini.ell_max;
ell = 1;
radar_op = radar_ini;
fdcomm_op = fdcomm_ini;
radar = radar_ini;
fdcomm = fdcomm_ini;
cov = cov_ini;
I_total = zeros(ell_max,1);
I_total_op = zeros(ell_max,1);
I_DL_op = zeros(ell_max,1);
I_UL_op = zeros(ell_max,1);
I_radar_op = zeros(ell_max,1);
I_UL = fdcomm.alpha_UL.*fdcomm.MI_UL;
I_UL_max = sum(I_UL(:));
I_DL = fdcomm.alpha_DL.*fdcomm.MI_DL;
I_DL_max = sum(I_DL(:));
I_radar_max = sum(radar.alpha_r.*radar.MI_radar);
I_max = I_UL_max+I_DL_max+ I_radar_max;
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
    I_UL_ell = fdcomm.alpha_UL.*fdcomm.MI_UL;
    I_DL_ell = fdcomm.alpha_DL.*fdcomm.MI_DL;
    I_radar_ell = radar.alpha_r.*radar.MI_radar;
    I_total(ell) = sum(I_UL_ell(:))+sum(I_DL_ell(:))+ sum(I_radar_ell(:));
    if I_total(ell) > I_max
        I_total_op(ell) = I_total(ell);
        I_UL_op(ell) = sum(I_UL_ell(:));
        I_DL_op(ell) = sum(I_DL_ell(:));
        I_radar_op(ell) = sum(I_radar_ell);
        I_max = I_total(ell);
        radar_op = radar;
        fdcomm_op = fdcomm;
    else
        I_total_op(ell) = I_max;
        I_UL_op(ell) = I_UL_max;
        I_DL_op(ell) = I_DL_max;
        I_radar_op(ell) = I_radar_max;
    end
    ell = ell+1;
end
fdcomm_op.I_total = I_total;
fdcomm_op.I_total_op = I_total_op;
fdcomm_op.I_UL_op = I_UL_op;
fdcomm_op.I_DL_op = I_DL_op;
fdcomm_op.I_radar_op = I_radar_op;
cov_op = tsp_covmat(fdcomm_op,radar_op,radar_comm);
end

