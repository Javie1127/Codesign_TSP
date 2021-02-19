function [fdcomm] = tsp_UL_precoders(k, ii, fdcomm, cov)

I = fdcomm.UL_num;
J = fdcomm.DL_num;
% napla_PiB_R_UL_k = cell(I,1);
% napla_PiB_R_DL_k = cell(J,1);
napla_Pui_R_UL_sum = 0;
napla_Pui_R_DL_sum = 0;
d_UL = fdcomm.ULstream_num; % number of data streams of the UL UE
d_DL = fdcomm.DLstream_num;
H_iB = fdcomm.ULchannels{ii,1};
tilde_P_iu_k = fdcomm.ULprecoders{ii,k}; % being updated
lambda_iu_k = fdcomm.lambda_UL(ii,k);
N_i = fdcomm.UL_UE_Ant(ii);
%% update Au
Au = fdcomm.Au_fixed + lambda_iu_k*eye(N_i);
%% Bu
Bu = fdcomm.Bu_fixed;
%% Derivatives of the UL achivable rate
napla_Pui_R_UL = cell(I,1);
for q = 1:I
    R_in_uq_k = cov.in_UL{q,k};
    H_qB = fdcomm.ULchannels{q,1};
    P_qu_k = fdcomm.ULprecoders{q,k};
    % mu_qu_k = fdcomm.mu_UL(q,k);
    if q == ii
        napla_Pui_R_UL{q,1} = H_qB'/R_in_uq_k*H_qB*P_qu_k/...
          (eye(d_UL(q))+(P_qu_k'*H_qB'/R_in_uq_k*H_qB*tilde_P_iu_k));
    else
        napla_Pui_R_UL{q,1} = -H_qB'/R_in_uq_k*H_qB*P_qu_k/...
           (eye(d_UL(q)) + (P_qu_k'*H_qB'/R_in_uq_k*H_qB*P_qu_k))*...
           P_qu_k'*H_qB'/R_in_uq_k*H_iB*tilde_P_iu_k;
    end
end
%% Derivatives of the DL achivable rate
napla_Pui_R_DL = cell(J,1);
for jj = 1:J
    H_ij = fdcomm.ULDLchannels{ii,jj};
    H_Bj = fdcomm.DLchannels{jj};
    R_in_dj_k = cov.in_DL{jj,k};
    P_dj_k = fdcomm.DLprecoders{jj,k};
    napla_Pui_R_DL{jj,1} = -H_ij'/R_in_dj_k*H_Bj*P_dj_k/...
        (eye(d_DL(jj))+(P_dj_k'*H_Bj'/R_in_dj_k*H_Bj*P_dj_k))*...
        P_dj_k'*H_Bj'/R_in_dj_k*H_ij*tilde_P_iu_k;
end
%% Update Cu
for q = 1:I
    mu_qu_k = fdcomm.mu_UL(q,k);
    if q == ii
        napla_Pui_R_UL_sum = mu_qu_k*napla_Pui_R_UL{q,1}+ napla_Pui_R_UL_sum;
    else
        napla_Pui_R_UL_sum = -mu_qu_k*napla_Pui_R_UL{q,1}+napla_Pui_R_UL_sum;
    end
end
for jj = 1:J
    mu_gd_k = fdcomm.mu_DL(jj,k);
    napla_Pui_R_DL_sum = napla_Pui_R_DL_sum-mu_gd_k*napla_Pui_R_DL{jj,1};
end
Cu = fdcomm.Cu_fixed + napla_Pui_R_UL_sum+napla_Pui_R_DL_sum;
%% Calculate PiB
P_iu_k = sylvester(Au,Bu,Cu);
fdcomm.ULprecoders{ii,k} = P_iu_k;
end
