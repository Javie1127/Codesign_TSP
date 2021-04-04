function [fdcomm] = tsp_DL_precoders(k, jj, fdcomm,cov,tilde_P_dj_k)
%function [P_iB_k] = UL_precoders(ii, k,...
%   H_UL, H_i_MUI, W_DL_k, W_UL_k, U_DL_k, U_UL_k,Ur,...
%   lambda_UL, mu_DL, mu_UL, napla_PiB_R_DL_k,...,
%    napla_PiB_R_UL_k, napla_PiBk_Xi_r)
% H_DL: downlink channel matrices
% H_BB: Self-interference channel matrix
% W_DL_k: DL UE Weight Matices
% W_UL_k: UL UE Weight Matrices
% Wr : radar Weight matrices; 
% A: K*Mr
% F_rtr Mr*Nr contains the doppler frequencies of each TX-RX pair
% lambda_UL: Lagrange multiplier vector for the UL power constraints
% mu_DL: K*J Lagrange multiplier vector for the DL rate constraints
% mu_UL: K*I Lagrange multiplier vector for the UL rate constraints
% napla_PBj_R_DL_k: Mc*dj*J
% napla_PBj_R_UL_k: Mc*dj*I
I = fdcomm.UL_num;
J = fdcomm.DL_num;
Mc = fdcomm.BSTx;
d_DL = fdcomm.DLstream_num; % number of data streams of the DL UEs
H_BB = fdcomm.BBchannel;
if nargin == 4
    tilde_P_dj_k = fdcomm.DLprecoders{jj,k};
end
%% Update Ad;
Ad = fdcomm.Ad_fixed + fdcomm.lambda_DL(k)*eye(Mc);
%% Derivatives of the UL achivable rate w.r.t P_dj_k
napla_P_dj_k_R_UL_sum = 0;
for q = 1:I
    mu_qu_k = fdcomm.mu_UL(q,k);
%     mu_qu_k =1;
    H_qB = fdcomm.ULchannels{q,1};
    R_in_qu_k = cov.in_UL{q,k};
    P_uq_k = fdcomm.ULprecoders{q,k};
    E_qu_k_star = fdcomm.UL_MMSE{q,k};
    cc = chol(R_in_qu_k);
    qq = qr(cc);
    R_in_qu_k_inv = qq\(qq'\eye(size(R_in_qu_k,1)));
%     term_2 =nearestSPD(H_BB'*(R_in_qu_k\H_qB)*P_uq_k*E_qu_k_star*...
%         P_uq_k'*H_qB'*(R_in_qu_k\H_BB));
    term_2 =nearestSPD(H_BB'*(R_in_qu_k_inv)*H_qB*P_uq_k*E_qu_k_star*...
        P_uq_k'*H_qB'*(R_in_qu_k_inv)*H_BB);
    napla_P_dj_k_R_UL_sum = napla_P_dj_k_R_UL_sum - mu_qu_k*term_2*tilde_P_dj_k;
end
%% Derivatives of the DL achivable rate w.r.t. P_dj_k
napla_P_dj_k_R_DL_sum = 0;
H_Bj = fdcomm.DLchannels{jj,1};
for g = 1:J
    R_in_gd_k = cov.in_DL{g,k};
    mu_gd_k = fdcomm.mu_DL(g,k);
%     mu_gd_k = 0;
    if g == jj
        E_dj_k_star = fdcomm.DL_MMSE{g,k};
       napla_P_dj_k_R_DL_sum = napla_P_dj_k_R_DL_sum + mu_gd_k*(H_Bj'*(R_in_gd_k\H_Bj))*tilde_P_dj_k*E_dj_k_star;
         %(eye(d_DL(g))+tilde_P_dj_k'*H_Bj'/(R_in_gd_k)*H_Bj*tilde_P_dj_k);
%         napla_P_dj_k_R_DL_sum = napla_P_dj_k_R_DL_sum + mu_gd_k*(H_Bj'/(R_in_gd_k)*H_Bj)*tilde_P_dj_k*...
%           pinv(eye(d_DL(g))+tilde_P_dj_k'*H_Bj'/(R_in_gd_k)*H_Bj*tilde_P_dj_k);
    else
        H_Bg = fdcomm.DLchannels{g,1};
        P_dg_k = fdcomm.DLprecoders{g,k};
        E_dg_k_star = fdcomm.DL_MMSE{g,k};
        napla_P_dj_k_R_DL_sum = napla_P_dj_k_R_DL_sum - mu_gd_k*(H_Bg'*(R_in_gd_k\H_Bg))*P_dg_k*E_dg_k_star*...
           P_dg_k'*H_Bg'*(R_in_gd_k\H_Bg)*tilde_P_dj_k;
    end
end
Cd = napla_P_dj_k_R_DL_sum + napla_P_dj_k_R_UL_sum + fdcomm.Cd_fixed;
B_Bm_j_k = fdcomm.B_Bm_j_k;
B_Bt_j_k = fdcomm.B_Bt_j_k;
F_Bt_k = fdcomm.F_Bt_k;
F_Bm_k = fdcomm.F_Bm_k;
D_d_j = d_DL(jj);
P_dj_k = (kron(eye(D_d_j),Ad)+kron(B_Bt_j_k.',F_Bt_k)+kron(B_Bm_j_k.',F_Bm_k))\reshape(Cd,[],1);
P_dj_k = reshape(P_dj_k,Mc,D_d_j);
fdcomm.DLprecoders{jj,k} = P_dj_k;
end
