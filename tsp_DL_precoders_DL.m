function [DLcomm] = tsp_DL_precoders_DL(k, jj, DLcomm,cov_DL)
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
J = DLcomm.DL_num;
Mc = DLcomm.BSTX;
d_DL = DLcomm.DLstream_num; % number of data streams of the DL UEs
tilde_P_dj_k = DLcomm.DLprecoders{jj,k};
%% Update Ad;
Ad = DLcomm.Ad_fixed + DLcomm.lambda_DL(k)*eye(Mc);
%% Derivatives of the DL achivable rate w.r.t. P_dj_k
napla_P_dj_k_R_DL_sum = 0;
H_Bj = DLcomm.DLchannels{jj,1};
for g = 1:J
    R_in_gd_k = cov_DL.in_DL{g,k};
    mu_gd_k = DLcomm.mu_DL(g,k);
    if g == jj
       napla_P_dj_k_R_DL_sum = napla_P_dj_k_R_DL_sum + mu_gd_k*nearestSPD(H_Bj'/R_in_gd_k*H_Bj)*tilde_P_dj_k/...
           nearestSPD(eye(d_DL(g))+tilde_P_dj_k'*H_Bj'/R_in_gd_k*H_Bj*tilde_P_dj_k);
    else
        H_Bg = DLcomm.DLchannels{g,1};
        P_dg_k = DLcomm.DLprecoders{g,k};
        napla_P_dj_k_R_DL_sum = napla_P_dj_k_R_DL_sum - mu_gd_k*nearestSPD(H_Bg'/R_in_gd_k*H_Bg)*tilde_P_dj_k/...
           (eye(d_DL(g))+nearestSPD(P_dg_k'*H_Bg'/R_in_gd_k*H_Bg*P_dg_k))*...
           P_dg_k'*H_Bg'/R_in_gd_k*H_Bg*tilde_P_dj_k;
    end
end
Cd = napla_P_dj_k_R_DL_sum + DLcomm.Cd_fixed;
Bd = DLcomm.Bd_fixed;
P_dj_k = sylvester(Ad,Bd,Cd);
DLcomm.DLprecoders{jj,k} = P_dj_k;
end
