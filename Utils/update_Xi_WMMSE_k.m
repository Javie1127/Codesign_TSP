function [fdcomm,radar] = update_Xi_WMMSE_k(fdcomm, radar,cov,k)
J = fdcomm.DL_num;
I = fdcomm.UL_num;

% E_DL_star = fdcomm.DL_MMSE;

Xi_UL_k = 0;
for ii = 1:I
    HiB = fdcomm.ULchannels{ii};
    R_UL_i_k = cov.total_UL{ii,k};
    P_ui_k = fdcomm.ULprecoders{ii,k};
    U_ui_k = fdcomm.UL_WMMSE_RX{ii,k};
    alpha_ui = fdcomm.alpha_UL(ii);
    W_ui_k = fdcomm.UL_weights{ii,k};
%     E_ui_k_star = fdcomm.UL_MMSE{ii,k};
    E_ui_k = eye(fdcomm.ULstream_num(ii))-U_ui_k*HiB*P_ui_k-P_ui_k'*HiB'*U_ui_k'+ (U_ui_k*R_UL_i_k*U_ui_k');
    d_new = eye(size(E_ui_k), 'logical');
    E_ui_k(d_new) = real(diag(E_ui_k));
    %E_ui_k = nearestSPD(E_ui_k);
    fdcomm.UL_MMSE_nop{ii,k} = E_ui_k;
%     fdcomm.Xi_WMMSE_UL{ii,k} = alpha_ui*real(trace(E_ui_k_star\E_ui_k));
    fdcomm.Xi_WMMSE_UL{ii,k} = alpha_ui*real(trace(W_ui_k*E_ui_k));
    Xi_UL_k = Xi_UL_k + fdcomm.Xi_WMMSE_UL{ii,k};
end
Xi_DL_k = 0;
for jj = 1:J
    alpha_dj = fdcomm.alpha_DL(jj);
    HBj = fdcomm.DLchannels{jj};
    U_dj_k = fdcomm.DL_WMMSE_RX{jj,k};
    R_DL_j_k = cov.total_DL{jj,k};
    P_dj_k = fdcomm.DLprecoders{jj,k};
    W_dj_k = fdcomm.DL_weights{jj,k};
    E_dj_k = eye(fdcomm.DLstream_num(jj))-U_dj_k*HBj*P_dj_k-P_dj_k'*HBj'*U_dj_k'+ U_dj_k*R_DL_j_k*U_dj_k';
    d_new = eye(size(E_dj_k), 'logical');
    E_dj_k(d_new) = real(diag(E_dj_k));
    %E_dj_k = nearestSPD(E_dj_k);
    fdcomm.DL_MMSE_nop{jj,k} = E_dj_k;
%     fdcomm.Xi_WMMSE_DL{jj,k} =  alpha_dj*abs(trace(E_DL_star{jj,k}\E_dj_k));
    fdcomm.Xi_WMMSE_DL{jj,k} =  alpha_dj*real(trace(W_dj_k*E_dj_k));
    Xi_DL_k = Xi_DL_k+fdcomm.Xi_WMMSE_DL{jj,k};
end
fdcomm.Xi_UL(k) = Xi_UL_k;
fdcomm.Xi_DL(k) = Xi_DL_k;

%-----------radar WMMSE------------------
Nr = radar.Rx;
Xi_radar_Nr = zeros(Nr,1);
Sigma_t_Nr = radar.Sigma_h_tr;
S_tr = cov.S_tr;
% Er = cell(Nr,1);
Er_star = radar.MMSE;
for nr = 1:Nr
    Sigma_tnr = Sigma_t_Nr(:,:,nr);
    Urnr = radar.WMMSE_RX{nr,1};
    Rr_nr = cov.total_r(:,:,nr);
    Ernr = Sigma_tnr- Urnr*S_tr*Sigma_tnr - Sigma_tnr'*S_tr'*Urnr' + Urnr*Rr_nr*Urnr';
    d_new = eye(size(Ernr), 'logical');
    Ernr(d_new) = real(diag(Ernr));
    %Ernr = nearestSPD(Ernr);
%     Wrnr = radar.WMMSE_weights{nr,1};
    Ernr_star = Er_star{nr,1};
    alpha_nr = radar.alpha_r(nr);
%     ee = Ernr_star\Ernr;
%     Xi_radar_Nr(nr) = alpha_nr*abs(sum(diag(ee)));
    Xi_radar_Nr(nr) = alpha_nr*real(trace(Ernr_star\Ernr));
%     Xi_radar_Nr(nr) = alpha_nr*real(sum(sum(Wrnr.'.*Ernr)));
%     Xi_radar_Nr(nr) = alpha_nr*abs(trace(Wrnr*Ernr));
%     Er{nr,1} = Ernr;
end
% radar.MMSE_nop = Er;
radar.Xi_radar_Nr = Xi_radar_Nr;
fdcomm.Xi_WMMSE_total_k = Xi_UL_k + Xi_DL_k + sum(Xi_radar_Nr);
end