function [ULcomm_op, radar_op,cov_op] = tsp_altermating_projection_UL...
(ULcomm, radar,radar_comm,cov)
ell_max = radar.ell_max;
ell = 1;
radar_op = radar;
ULcomm_op = ULcomm; 
I_total = zeros(ell_max,1);
I_total_op = zeros(ell_max,1);
I_UL_op = zeros(ell_max,1);
I_radar_op = zeros(ell_max,1);
I_UL = ULcomm.alpha_UL.*ULcomm.MI_UL;
I_UL_max = sum(I_UL(:));
I_radar_max = sum(radar.alpha_r.*radar.MI_radar);
I_max = I_UL_max + I_radar_max;
while ell <= ell_max
    [ULcomm, radar] =...
        tsp_WMMSE_algorithm_UL(ULcomm, radar, radar_comm, cov);
    % Calculate A^star
    radar = tsp_Nearest_PAR(radar);
    % Update covariance matrices
    cov = tsp_covmat_UL(ULcomm,radar,radar_comm);
    % Update linear receivers
    ULcomm = tsp_Comm_MMSE_UL(ULcomm, radar, cov);
    radar = tsp_radar_MMSE(radar, cov);
    I_UL_ell = ULcomm.alpha_UL.*ULcomm.MI_UL;
    I_radar_ell = radar.alpha_r.*radar.MI_radar;
    I_total(ell) = sum(I_UL_ell(:))+ sum(I_radar_ell(:));
    if I_total(ell) > I_max
        I_total_op(ell) = I_total(ell);
        I_UL_op(ell) = sum(I_UL_ell(:));
        I_radar_op(ell) = sum(I_radar_ell);
        I_max = I_total(ell);
        radar_op = radar;
        ULcomm_op = ULcomm;
    else
        I_total_op(ell) = I_max;
        I_UL_op(ell) = I_UL_max;
        I_radar_op(ell) = I_radar_max;
    end
    ell = ell+1;
end
ULcomm_op.I_total = I_total;
ULcomm_op.I_total_op = I_total_op;
ULcomm_op.I_UL_op = I_UL_op;
ULcomm_op.I_radar_op = I_radar_op;
cov_op = tsp_covmat_UL(ULcomm_op,radar_op,radar_comm);
end

