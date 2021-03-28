function [fdcomm,radar] = tsp_MI(fdcomm,radar,cov)
%
%   Detailed explanation goes here
K = radar.codelength;
I = fdcomm.UL_num;
J = fdcomm.DL_num;
I_DL = zeros(J,K);
I_UL = zeros(I,K);
U_UL = fdcomm.UL_WMMSE_RX;
U_DL = fdcomm.DL_WMMSE_RX;
R_DL = cov.total_DL;
R_UL = cov.total_UL;
R_in_UL = cov.in_UL;
R_in_DL = cov.in_DL;
for k = 1:K
    for jj = 1 : J
        d_DL_j = fdcomm.DLstream_num(jj);
        R_in_dj = R_in_DL{jj,k};
        alpha_dj = fdcomm.alpha_DL(jj);
        R_DL_j_k = R_DL{jj,k};
        U_dj_k = U_DL{jj,k};
        I_DL(jj,k) = alpha_dj*abs(log2(det(eye(d_DL_j)+ (U_dj_k*R_DL_j_k*U_dj_k')/(U_dj_k*R_in_dj*U_dj_k'))));
        %I_DL(jj,k) = abs(log2(det(eye(d_DL_j)+ (U_dj_k*R_DL_j_k*U_dj_k')/(U_dj_k*R_in_dj*U_dj_k'))));
    end
    for ii = 1:I
        R_in_i_k = R_in_UL{ii,k};
        alpha_ui = fdcomm.alpha_UL(ii);
        d_UL_i = fdcomm.ULstream_num(ii);
        R_UL_i_k = R_UL{ii,k};
        U_ui_k = U_UL{ii,k};
        I_UL(ii,k) = alpha_ui*abs(log2(det(eye(d_UL_i)+ ...
            (U_ui_k*R_UL_i_k*U_ui_k')/(U_ui_k*R_in_i_k*U_ui_k'))));
    end
end
Nr = radar.Rx;
I_r = zeros(Nr,1); 
%R_r = cov.total_r;
R_in = cov.inr;
Ur = radar.WMMSE_RX;
M = radar.total_Tx;
R_tr = cov.target2radar;
for nr = 1:Nr
    alpha_nr = radar.alpha_r(nr);
    R_in_nr = R_in(:,:,nr);
    %R_r_nr = R_r(:,:,nr);
    R_tr_nr = R_tr(:,:,nr);
    Urnr = Ur{nr,1};
    R_in_tilde = nearestSPD(Urnr*R_in_nr*Urnr');
    %R_r_tilde = nearestSPD(Urnr*R_r_nr*Urnr');
    R_tr_tilde = nearestSPD(Urnr*R_tr_nr*Urnr');
    %term_1 = (R_t_tilde)/(R_in_tilde);
    %ee = eye(size(Sigma_t_Nr,2)) + (R_r_tilde)/(R_in_tilde);
    a1 = chol(R_in_tilde);
    a1 = triu(a1);
    R_in_tilde_inv = a1\(a1'\eye(size(R_in_tilde,1)));
%     ee = abs(log2(det(eye(K*M)+R_tr_tilde/R_in_tilde)));
    ee = abs(log2(det(eye(K*M)+R_tr_tilde*R_in_tilde_inv)));
    I_r_nr= alpha_nr*ee;
    %I_r_nr = alpha_nr*log2(det(R_r_tilde)/det(R_in_tilde));
    I_r(nr)=I_r_nr;
end

fdcomm.MI_DL = I_DL;
fdcomm.MI_UL = I_UL;
radar.MI_radar = I_r;
fdcomm.MI_total = sum(I_DL(:))+sum(I_UL(:))+sum(I_r);
end

