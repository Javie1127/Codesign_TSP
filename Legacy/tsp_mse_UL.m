function [fdcomm,radar] = tsp_mse_UL(fdcomm, radar, cov,k)
%-------update the instant MSE matrix--------------
I = fdcomm.UL_num;
Nr = radar.RX;
Sigma_tNr = radar.Sigma;
M = size(Sigma_tNr,2);
% UL_rate = zeros(I,K);
% DL_rate = zeros(J,K);
for ii = 1 : I
    HiB = fdcomm.ULchannels{ii};
    R_UL_i_k = cov.total_UL{ii,k};
    R_in_i_k = cov.in_UL{ii,k};
    P_ui_k = fdcomm.ULprecoders{ii,k};
    U_ui_k = fdcomm.UL_WMMSE_RX{ii,k};
    fdcomm.UL_MMSE{ii,k} = eye(fdcomm.ULstream_num(ii))- nearestSPD(P_ui_k'*HiB'/R_UL_i_k*HiB*P_ui_k);
    fdcomm.UL_MMSE_nop{ii,k} = eye(fdcomm.ULstream_num(ii))-2*U_ui_k*HiB*P_ui_k + nearestSPD(U_ui_k*R_in_i_k*U_ui_k');
end
for nr = 1:Nr
    Sigma_tnr = Sigma_tNr(:,:,nr);
    St_nr = cov.S_tr(:,:,nr);
    R_in_nr = cov.inr(:,:,nr);
    R_t_nr = cov.target2radar(:,:,nr);
    Urnr = radar.WMMSE_RX{nr,1};
    R_r_nr = nearestSPD(R_in_nr + R_t_nr);
    Ernr_star = nearestSPD(Sigma_tnr*(eye(M)-nearestSPD(St_nr'/R_r_nr*St_nr*Sigma_tnr)));
    % Ernr = nearestSPD(Sigma_tnr-2*Urnr*St_nr*Sigma_tnr + nearestSPD(Urnr*(St_nr*St_nr')*Urnr')+nearestSPD(Urnr*R_in_nr*Urnr'));
    Ernr = Sigma_tnr-Sigma_tnr*St_nr'*Urnr'-Urnr*St_nr*Sigma_tnr+ nearestSPD(Urnr*St_nr*Sigma_tnr*St_nr'*Urnr')...
        +nearestSPD(Urnr*R_in_nr*Urnr');
    radar.MMSE{nr,1} = Ernr_star;
    radar.MMSE_nop{nr,1} = Ernr;
end
end

