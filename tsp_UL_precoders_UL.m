function [ULcomm] = tsp_UL_precoders_UL(k, ii, ULcomm, cov_UL)

I = ULcomm.UL_num;
napla_Pui_R_UL_sum = 0;
d_UL = ULcomm.ULstream_num; % number of data streams of the UL UE
H_iB = ULcomm.ULchannels{ii,1};
tilde_P_iu_k = ULcomm.ULprecoders{ii,k}; % being updated
lambda_iu_k = ULcomm.lambda_UL(ii,k);
N_i = ULcomm.UL_UE_Ant(ii);
%% update Au
Au = ULcomm.Au_fixed + lambda_iu_k*eye(N_i);
%% Bu
Bu = ULcomm.Bu_fixed;
%% Derivatives of the UL achivable rate
napla_Pui_R_UL = cell(I,1);
for q = 1:I
    R_in_uq_k = cov_UL.in_UL{q,k};
    H_qB = ULcomm.ULchannels{q,1};
    P_qu_k = ULcomm.ULprecoders{q,k};
    % mu_qu_k = fdcomm.mu_UL(q,k);
    if q == ii
        napla_Pui_R_UL{q,1} = H_qB'/R_in_uq_k*H_qB*P_qu_k/...
          (eye(d_UL(q))+nearestSPD(P_qu_k'*H_qB'/R_in_uq_k*H_qB*tilde_P_iu_k));
    else
        napla_Pui_R_UL{q,1} = -H_qB'/R_in_uq_k*H_qB*P_qu_k/...
           (eye(d_UL(q)) + nearestSPD(P_qu_k'*H_qB'/R_in_uq_k*H_qB*P_qu_k))*...
           P_qu_k'*H_qB'/R_in_uq_k*H_iB*tilde_P_iu_k;
    end
end
%% Update Cu
for q = 1:I
    mu_qu_k = ULcomm.mu_UL(q,k);
    if q == ii
        napla_Pui_R_UL_sum = mu_qu_k*napla_Pui_R_UL{q,1}+ napla_Pui_R_UL_sum;
    else
        napla_Pui_R_UL_sum = -mu_qu_k*napla_Pui_R_UL{q,1}+napla_Pui_R_UL_sum;
    end
end
Cu = ULcomm.Cu_fixed + napla_Pui_R_UL_sum;
%% Calculate PiB
P_iu_k = sylvester(Au,Bu,Cu);
ULcomm.ULprecoders{ii,k} = P_iu_k;
end
