function [ULcomm] = tsp_Comm_MMSE_UL(ULcomm, radar_UL, cov_UL)
%UComm_MMSE_receiver returns:
%----linear UL/DL MMSE receiver
%----Optimal Weight Matrices
K = radar_UL.codelength;
I = ULcomm.UL_num;
Nc = ULcomm.BSRX;
U_UL = cell(I,K);
W_UL = cell(I,K);
E_UL_star = cell(I,K);
E_UL = cell(I,K);
I_UL = zeros(I,K);
xi_UL = zeros(Nc,Nc,K);
for k = 1:K
    xi_UL_k = 0;
    for ii = 1 : I
        alpha_ui = ULcomm.alpha_UL(ii);
        d_UL_i = ULcomm.ULstream_num(ii);
        HiB = ULcomm.ULchannels{ii};
        R_UL_i_k = cov_UL.total_UL{ii,k};
        R_in_i_k = cov_UL.in_UL{ii,k};
        P_ui_k = ULcomm.ULprecoders{ii,k};
        U_ui_k = P_ui_k'*HiB'/(R_UL_i_k);
        E_ui_k_star = eye(ULcomm.ULstream_num(ii))- P_ui_k'*HiB'/R_UL_i_k*HiB*P_ui_k;
        W_ui_k = inv(E_ui_k_star);
        E_ui_k = eye(ULcomm.ULstream_num(ii))-2*U_ui_k*HiB*P_ui_k+ U_ui_k*R_in_i_k*U_ui_k';
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
ULcomm.UL_WMMSE_RX = U_UL;
ULcomm.UL_weights = W_UL;
ULcomm.UL_MMSE = E_UL_star;
ULcomm.UL_MMSE_nop = E_UL;
ULcomm.MI_UL = I_UL;
ULcomm.xi_UL = xi_UL;
end

