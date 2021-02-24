function [DLcomm] = tsp_Comm_MMSE_DL(DLcomm, radar_DL, cov_DL)
%UComm_MMSE_receiver returns:
%----linear UL/DL MMSE receiver
%----Optimal Weight Matrices
K = radar_DL.codelength;
J = DLcomm.DL_num;
U_DL = cell(J,K);
W_DL = cell(J,K);
E_DL_star = cell(J,K);
E_DL = cell(J,K);
I_DL = zeros(J,K);
xi_DL = cell(K,J);
for k = 1:K
    for jj = 1 : J
        alpha_dj = DLcomm.alpha_DL(jj);
        d_DL_j = DLcomm.DLstream_num(jj);
        HBj = DLcomm.DLchannels{jj};
        R_DL_j_k = cov_DL.total_DL{jj,k};
        P_dj_k = DLcomm.DLprecoders{jj,k};
        U_dj_k = P_dj_k'*HBj'/(R_DL_j_k);
        E_dj_k_star = eye(DLcomm.DLstream_num(jj))-P_dj_k'*HBj'/R_DL_j_k*HBj*P_dj_k;
        W_dj_k = inv(E_dj_k_star);
        E_dj_k = eye(DLcomm.DLstream_num(jj))-2*U_dj_k*HBj*P_dj_k+U_dj_k*R_DL_j_k*U_dj_k';
        R_in_dj = cov_DL.in_DL{jj,k};
        %DL_rate(jj,k) = real(log2(det(nearestSPD(eye(d_DL_j)+nearestSPD(P_jd_k'*HBj'/R_in_dj*HBj*P_jd_k)))));
        I_DL(jj,k) = real(log2(det(eye(d_DL_j)+ (U_dj_k*R_DL_j_k*U_dj_k')/(U_dj_k*R_in_dj*U_dj_k'))));
        W_DL{jj,k} = W_dj_k;
        U_DL{jj,k} = U_dj_k;
        E_DL_star{jj,k} = E_dj_k_star;
        E_DL{jj,k} = E_dj_k;
        %xi_DL{k,jj} = alpha_dj*U_dj_k'*W_dj_k*U_dj_k;
        xi_DL{k,jj} = alpha_dj*U_dj_k'/E_dj_k_star*U_dj_k;
    end
end
DLcomm.DL_WMMSE_RX = U_DL;
DLcomm.DL_weights = W_DL;
DLcomm.DL_MMSE = E_DL_star;
DLcomm.DL_MMSE_nop = E_DL;
%fdcomm.DL_rate = DL_rate;
%fdcomm.UL_rate = UL_rate;
DLcomm.MI_DL = I_DL;
DLcomm.xi_DL = xi_DL;
end

