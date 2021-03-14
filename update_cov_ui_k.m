function [cov] = update_cov_ui_k(fdcomm,radar,cov,radar_comm,P_iB_k_old,ii,k)

%% Update UL related covs
R_IB = cov.UL;
H_UL = fdcomm.ULchannels;
H_iB = H_UL{ii};
P_uI = fdcomm.ULprecoders;
P_iB_k = P_uI{ii,k};
R_iB_old = R_IB{ii,k};
R_iB_k = (H_iB*(P_iB_k*P_iB_k')*H_iB');
d_new = eye(size(R_iB_k), 'logical');
R_iB_k(d_new) = abs(diag(R_iB_k));
R_IB{ii,k} = R_iB_k;
I = fdcomm.UL_num;
R_MUI_UL = cov.MUI_UL ;
R_in_UL = cov.in_UL;
R_total_UL = cov.total_UL;
for q = 1:I
    if q~= ii
        R_in_UL{q,k} = R_in_UL{q,k}-R_iB_old+R_IB{ii,k};
        R_MUI_UL{q,k} = R_MUI_UL{q,k}-R_iB_old + R_IB{ii,k};
        R_total_UL{q,k} = R_total_UL{q,k} - R_iB_old+R_IB{ii,k};
    else
        R_total_UL{q,k}= R_total_UL{q,k} - R_iB_old+R_IB{ii,k};
    end
end
cov.MUI_UL = R_MUI_UL;
cov.UL = R_IB;
cov.in_UL = R_in_UL;
cov.total_UL = R_total_UL;
%% Update radar related Covs
S_Ur = cov.S_Ur;
D_UL = fdcomm.ULsymbols;
nu = radar_comm.nu;
Nr = radar.Rx;
K = radar.codelength;
for nr = 1:Nr
    P_ui_k = P_uI{ii,k};
    d_ui_k = D_UL{ii}(:,nu,k);
    S_Ur{k,ii,nr} = P_ui_k*d_ui_k;
end
Sigma_U_Nr = radar_comm.Sigma_U_Nr;
R_U_Nr = cov.U_r;
R_U_Nr_total = cov.UL2radar;
R_in_r = cov.inr;
R_r = cov.total_r;
Nr = radar.Rx;
for nr = 1:Nr
    R_U_nr_temp = 0;
    for ii = 1:I
        R_i_nr = zeros(K,K);
        for m = 1:K
            s_inr_m = S_Ur{m,ii,nr};
           for l = 1:K
               s_inr_l = S_Ur{l,ii,nr};
               R_i_nr(m,l) = abs(trace(s_inr_m*s_inr_l'*Sigma_U_Nr{m,l,ii,nr}));
           end
        end
        %R_i_nr = nearestSPD(R_i_nr);
        R_U_Nr{ii,nr} = R_i_nr;
        R_U_nr_temp = R_U_nr_temp + R_i_nr;
    end
    R_in_r(:,:,nr) = R_in_r(:,:,nr) - R_U_Nr_total(:,:,nr)+R_U_nr_temp;
    R_r(:,:,nr) = R_r(:,:,nr) - R_U_Nr_total(:,:,nr)+R_U_nr_temp;
    R_U_Nr_total(:,:,nr) = R_U_nr_temp;
end
cov.S_Ur = S_Ur;
cov.U_r = R_U_Nr;
cov.UL2radar = R_U_Nr_total;
cov.inr = R_in_r;
cov.total_r = R_r;
%% Update DL related covariance matrices
if fdcomm.DL_num > 0
    H_UL_DL = fdcomm.ULDLchannels;
    R_ULDL = cov.UL2DL;
    R_in_DL = cov.in_DL;
    R_total_DL = cov.total_DL;
    J = fdcomm.DL_num;
    for jj = 1:J
        Hij = H_UL_DL{ii,jj};
        R_ij_old = (Hij*(P_iB_k_old*P_iB_k_old')*Hij');
        d_old =  eye(size(R_ij_old), 'logical');
        R_ij_old(d_old) = abs(diag(R_ij_old));
        R_ij_new = (Hij*(P_iB_k*P_iB_k')*Hij');
        d_new = eye(size(R_ij_new), 'logical');
        R_ij_new(d_new) = abs(diag(R_ij_new));
        R_in_DL{jj,k} = R_in_DL{jj,k}-R_ij_old+ R_ij_new;
        R_ULDL{jj,k} = R_ULDL{jj,k}-R_ij_old+ R_ij_new;
        R_total_DL{jj,k} = R_total_DL{jj,k}-R_ij_old+R_ij_new;
    end
    cov.UL2DL = R_ULDL;
    cov.in_DL = R_in_DL;
    cov.total_DL = R_total_DL;
end
end