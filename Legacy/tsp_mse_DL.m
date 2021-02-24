function [DLcomm,radar_DL] = tsp_mse_DL(DLcomm, radar_DL, cov_DL,k)
%-------update the instant MSE matrix--------------
J = DLcomm.DL_num;
Nr = radar_DL.RX;
Sigma_tNr = radar_DL.Sigma;
M = size(Sigma_tNr,2);
% UL_rate = zeros(I,K);
% DL_rate = zeros(J,K);
for jj = 1 : J
    U_dj_k = DLcomm.DL_WMMSE_RX{jj,k};
    HBj = DLcomm.DLchannels{jj};
    R_DL_j_k = cov_DL.total_DL{jj,k};
    P_dj_k = DLcomm.DLprecoders{jj,k};
    DLcomm.DL_MMSE{jj,k} = eye(DLcomm.DLstream_num(jj))- nearestSPD(P_dj_k'*HBj'/R_DL_j_k*HBj*P_dj_k);
    DLcomm.DL_MMSE_nop{jj,k} = eye(DLcomm.DLstream_num(jj))-2*U_dj_k*HBj*P_dj_k+nearestSPD(U_dj_k*R_DL_j_k*U_dj_k');
end
for nr = 1:Nr
    Sigma_tnr = Sigma_tNr(:,:,nr);
    St_nr = cov_DL.S_tr(:,:,nr);
    R_in_nr = cov_DL.inr(:,:,nr);
    R_t_nr = cov_DL.target2radar(:,:,nr);
    Urnr = radar_DL.WMMSE_RX{nr,1};
    R_r_nr = nearestSPD(R_in_nr + R_t_nr);
    Ernr_star = nearestSPD(Sigma_tnr*(eye(M)-nearestSPD(St_nr'/R_r_nr*St_nr*Sigma_tnr)));
    % Ernr = nearestSPD(Sigma_tnr-2*Urnr*St_nr*Sigma_tnr + nearestSPD(Urnr*(St_nr*St_nr')*Urnr')+nearestSPD(Urnr*R_in_nr*Urnr'));
    Ernr = Sigma_tnr-Sigma_tnr*St_nr'*Urnr'-Urnr*St_nr*Sigma_tnr+ nearestSPD(Urnr*St_nr*Sigma_tnr*St_nr'*Urnr')...
        +nearestSPD(Urnr*R_in_nr*Urnr');
    radar_DL.MMSE{nr,1} = Ernr_star;
    radar_DL.MMSE_nop{nr,1} = Ernr;
end
end

