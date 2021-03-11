function [fdcomm_op, radar_op] =...
    tsp_BCD_AP(fdcomm, radar, radar_comm, cov)

%------------------Algorithm 5 (BCD-AP MRMC)-------------------------------------------

ell_max = radar.ell_max;
radar_op = radar;
fdcomm_op = fdcomm;
radar_ell = radar;
fdcomm_ell = fdcomm;
cov_ell = cov;
I_total = zeros(ell_max,1);
I_total_op = zeros(ell_max,1);
I_DL_op = zeros(ell_max,1);
I_UL_op = zeros(ell_max,1);
I_radar_op = zeros(ell_max,1);
I_UL_max = 0;
I_DL_max = 0;
I_radar_max = 0;
I_max = 0;
ell = 1;
%% BCD-AP Iteration
while ell <= ell_max
    disp(ell);
    [fdcomm_ell, radar_ell,cov_ell] =...
        tsp_WMMSE_MRMC(fdcomm_ell, radar_ell, radar_comm, cov_ell);
    if strcmp(radar_ell.coding_type,'Proposed')
        % Calculate A^star
        radar_ell = tsp_Nearest_PAR(radar_ell);
        % Update covariance matrices
        cov_ell = up_cov_radar(fdcomm_ell,radar_ell,cov_ell,radar_comm);
    end
%     cov_ell = tsp_covmat_rev(fdcomm_ell,radar_ell,radar_comm,cov_ell);
    % Update linear receivers
    fdcomm_ell = tsp_Comm_MMSE_rev(fdcomm_ell, radar_ell, cov_ell);
    radar_ell = tsp_radar_MMSE_rev(radar_ell, cov_ell);
    [fdcomm_ell, radar_ell] = tsp_MI(fdcomm_ell,radar_ell,cov_ell);
    I_UL_ell = fdcomm_ell.MI_UL;
    I_DL_ell = fdcomm_ell.MI_DL;
    I_radar_ell = radar_ell.MI_radar;
    I_total(ell) = fdcomm_ell.MI_total;
    if I_total(ell) > I_max
        I_total_op(ell) = I_total(ell);
        I_UL_op(ell) = sum(I_UL_ell(:));
        I_UL_max = I_UL_op(ell);
        I_DL_op(ell) = sum(I_DL_ell(:));
        I_DL_max = I_DL_op(ell);
        I_radar_op(ell) = sum(I_radar_ell);
        I_radar_max = I_radar_op(ell);
        I_max = I_total(ell);
        radar_op = radar_ell;
        fdcomm_op = fdcomm_ell;
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
radar_op.I_radar_op = I_radar_op;
end

