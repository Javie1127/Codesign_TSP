function [DLcomm_op, radar_DL_op,cov_DL_op] = tsp_altermating_projection_DL_wc...
(DLcomm, radar_DL,radar_comm_DL,cov_DL)
ell_max = radar_DL.ell_max;
ell = 1;
radar_DL_op = radar_DL;
DLcomm_op = DLcomm; 
I_total = zeros(ell_max,1);
I_total_op = zeros(ell_max,1);
I_DL_op = zeros(ell_max,1);
I_radar_op = zeros(ell_max,1);
I_DL = DLcomm.alpha_DL.*DLcomm.MI_DL;
I_DL_max = sum(I_DL(:));
I_radar_max = sum(radar_DL.alpha_r.*radar_DL.MI_radar);
I_max = I_DL_max+ I_radar_max;
while ell <= ell_max
    [DLcomm, radar_DL] =...
        tsp_WMMSE_algorithm_DL_wc(DLcomm, radar_DL, radar_comm_DL, cov_DL);
    % Calculate A^star
    radar_DL = tsp_Nearest_PAR(radar_DL);
    % Update covariance matrices
    cov_DL = tsp_covmat_DL_wc(DLcomm,radar_DL,radar_comm_DL);
    % Update linear receivers
    DLcomm = tsp_Comm_MMSE_DL(DLcomm, radar_DL, cov_DL);
    radar_DL = tsp_radar_MMSE(radar_DL, cov_DL);
    I_DL_ell = DLcomm.alpha_DL.*DLcomm.MI_DL;
    I_radar_ell = radar_DL.alpha_r.*radar_DL.MI_radar;
    I_total(ell) = sum(I_DL_ell(:))+ sum(I_radar_ell);
    if I_total(ell) > I_max
        I_total_op(ell) = I_total(ell);
        I_max = I_total(ell);
        I_DL_op(ell) = sum(I_DL_ell(:));
        I_radar_op(ell) = sum(I_radar_ell);
        radar_DL_op = radar_DL;
        DLcomm_op = DLcomm;
    else
        I_total_op(ell) = I_max;
        I_DL_op(ell) = I_DL_max;
        I_radar_op(ell) = I_radar_max;
    end
    ell = ell+1;
end
DLcomm_op.I_total = I_total;
DLcomm_op.I_total_op = I_total_op;
DLcomm_op.I_DL_op = I_DL_op;
DLcomm_op.I_radar_op = I_radar_op;
cov_DL_op = tsp_covmat_DL_wc(DLcomm_op,radar_DL_op,radar_comm_DL);
end

