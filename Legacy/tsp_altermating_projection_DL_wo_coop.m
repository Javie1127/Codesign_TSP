function [DLcomm_op, radar_DL_op,cov_DL_op] = tsp_altermating_projection_DL_wo_coop...
(DLcomm, radar_DL,radar_comm_DL,cov_DL)
ell_max = radar_DL.ell_max;
ell = 1;
radar_DL_op = radar_DL;
DLcomm_op = DLcomm; 
I_total = zeros(ell_max,1);
I_total_op = zeros(ell_max,1);
I_max = sum(DLcomm.MI_DL(:))+ sum(radar_DL.MI_radar);
while ell <= ell_max
    [DLcomm, radar_DL] =...
        tsp_WMMSE_algorithm_DL(DLcomm, radar_DL, radar_comm_DL, cov_DL);
    % Calculate A^star
    radar_DL = tsp_Nearest_PAR(radar_DL);
    % Update covariance matrices
    cov_DL = tsp_covmat_DL(DLcomm,radar_DL,radar_comm_DL);
    % Update linear receivers
    DLcomm = tsp_Comm_MMSE_DL(DLcomm, radar_DL, cov_DL);
    radar_DL = tsp_radar_MMSE(radar_DL, cov_DL);
    I_total(ell) = sum(DLcomm.MI_DL(:))+ sum(radar_DL.MI_radar);
    if I_total(ell) > I_max
        I_total_op(ell) = I_total(ell);
        I_max = I_total(ell);
        radar_DL_op = radar_DL;
        DLcomm_op = DLcomm;
    else
        I_total_op(ell) = I_max;
    end
    ell = ell+1;
end
DLcomm_op.I_total = I_total;
DLcomm_op.I_total_op = I_total_op;
cov_DL_op = tsp_covmat_DL(DLcomm_op,radar_DL_op,radar_comm_DL);
end

