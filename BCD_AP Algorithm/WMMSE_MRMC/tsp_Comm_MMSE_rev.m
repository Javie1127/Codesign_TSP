function [fdcomm] = tsp_Comm_MMSE_rev(fdcomm, radar, cov)
%UComm_MMSE_receiver returns:
%----linear UL/DL MMSE receiver
%----Optimal Weight Matrices
K = radar.codelength;
%% DL WMMSE 
J = fdcomm.DL_num;
U_DL = cell(J,K);
W_DL = cell(J,K);
E_DL_star = cell(J,K);
E_DL = cell(J,K);
I_DL = zeros(J,K);
xi_DL = cell(K,J);
for k = 1:K
    for jj = 1 : J
        alpha_dj = fdcomm.alpha_DL(jj);
        d_DL_j = fdcomm.DLstream_num(jj);
        HBj = fdcomm.DLchannels{jj};
        R_DL_j_k = cov.total_DL{jj,k};
        P_dj_k = fdcomm.DLprecoders{jj,k};
        U_dj_k = P_dj_k'*HBj'/(R_DL_j_k);
        E_dj_k_star = eye(fdcomm.DLstream_num(jj))-P_dj_k'*HBj'/R_DL_j_k*HBj*P_dj_k;
        W_dj_k = inv(E_dj_k_star);
        E_dj_k = eye(fdcomm.DLstream_num(jj))-2*U_dj_k*HBj*P_dj_k+U_dj_k*R_DL_j_k*U_dj_k';
        R_in_dj = cov.in_DL{jj,k};
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
fdcomm.DL_WMMSE_RX = U_DL;
fdcomm.DL_weights = W_DL;
fdcomm.DL_MMSE = E_DL_star;
fdcomm.DL_MMSE_nop = E_DL;
fdcomm.MI_DL = I_DL;
fdcomm.xi_DL = xi_DL;

%% UL WMMSE
I = fdcomm.UL_num;
Nc = fdcomm.BSRx;
U_UL = cell(I,K);
W_UL = cell(I,K);
E_UL_star = cell(I,K);
E_UL = cell(I,K);
I_UL = zeros(I,K);
xi_UL = zeros(Nc,Nc,K);
for k = 1:K
    xi_UL_k = 0;
    for ii = 1 : I
        alpha_ui = fdcomm.alpha_UL(ii);
        d_UL_i = fdcomm.ULstream_num(ii);
        HiB = fdcomm.ULchannels{ii};
        R_UL_i_k = cov.total_UL{ii,k};
        R_in_i_k = cov.in_UL{ii,k};
        P_ui_k = fdcomm.ULprecoders{ii,k};
        U_ui_k = P_ui_k'*HiB'/(R_UL_i_k);
        E_ui_k_star = eye(fdcomm.ULstream_num(ii))- P_ui_k'*HiB'/R_UL_i_k*HiB*P_ui_k;
        W_ui_k = inv(E_ui_k_star);
        E_ui_k = eye(fdcomm.ULstream_num(ii))-2*U_ui_k*HiB*P_ui_k+ U_ui_k*R_in_i_k*U_ui_k';
        W_UL{ii,k} = W_ui_k;
        U_UL{ii,k} = U_ui_k;
        E_UL_star{ii,k} = E_ui_k_star;
        E_UL{ii,k} = E_ui_k;
%        UL_rate(ii,k) = ...
%             real(log2(det(nearestSPD(eye(d_UL_i)+nearestSPD(P_iB_k'*HiB'/R_in_i_k*HiB*P_iB_k)))));
        I_UL(ii,k) = real(log2(det(eye(d_UL_i)+ (U_ui_k*R_UL_i_k*U_ui_k')/(U_ui_k*R_in_i_k*U_ui_k'))));
        xi_UL_k = xi_UL_k + alpha_ui*U_ui_k'/E_ui_k_star*U_ui_k;
    end
    xi_UL(:,:,k) = xi_UL_k;
end
fdcomm.UL_WMMSE_RX = U_UL;
fdcomm.UL_weights = W_UL;
fdcomm.UL_MMSE = E_UL_star;
fdcomm.UL_MMSE_nop = E_UL;
fdcomm.MI_UL = I_UL;
fdcomm.xi_UL = xi_UL;

end

