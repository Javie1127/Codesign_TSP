function [fdcomm] = tsp_UL_precoders(k, ii, fdcomm, cov, tilde_P_iu_k)

I = fdcomm.UL_num;
J = fdcomm.DL_num;
% napla_PiB_R_UL_k = cell(I,1);
% napla_PiB_R_DL_k = cell(J,1);
napla_Pui_R_UL_sum = 0;
napla_Pui_R_DL_sum = 0;
d_UL = fdcomm.ULstream_num; % number of data streams of the UL UE
H_iB = fdcomm.ULchannels{ii,1};
if nargin == 4
    tilde_P_iu_k = fdcomm.ULprecoders{ii,k}; % being updated
end
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
    P_qu_k = fdcomm.ULprecoders{q,k};
    %  mu_qu_k = fdcomm.mu_UL(q,k);
    if q == ii
        E_ui_star = fdcomm.UL_MMSE{q,1};
        napla_Pui_R_UL{q,1} = H_iB'*(R_in_uq_k\H_iB)*tilde_P_iu_k*E_ui_star;
%         napla_Pui_R_UL{q,1} = H_iB'/R_in_uq_k*H_iB*tilde_P_iu_k/...
%           (eye(d_UL(q))+(P_qu_k'*H_iB'/R_in_uq_k*H_iB*tilde_P_iu_k));
    else
        H_qB = fdcomm.ULchannels{q,1};
        E_uq_star = fdcomm.UL_MMSE{q,1};
%         napla_Pui_R_UL{q,1} = -H_iB'*(R_in_uq_k\H_qB)*P_qu_k*E_uq_star*P_qu_k'*H_qB'*(R_in_uq_k\H_iB)*tilde_P_iu_k;
%         cc = chol(R_in_uq_k);
%         R_in_uq_k_inv = cc\(cc'\eye(size(R_in_uq_k,1)));
        R_in_uq_k_inv = inv_SVD(R_in_uq_k);
        napla_Pui_R_UL{q,1} = -H_iB'*(R_in_uq_k_inv)*H_qB*P_qu_k*E_uq_star*P_qu_k'*H_qB'*(R_in_uq_k_inv)*H_iB*tilde_P_iu_k;
        %   /(eye(d_UL(q)) + nearestSPD(P_qu_k'*H_qB'/R_in_uq_k*H_qB*P_qu_k))*...
           
    end
end
%% Derivatives of the DL achivable rate
napla_Pui_R_DL = cell(J,1);
for jj = 1:J
    H_ij = fdcomm.ULDLchannels{ii,jj};
    H_Bj = fdcomm.DLchannels{jj};
    R_in_dj_k = cov.in_DL{jj,k};
    P_dj_k = fdcomm.DLprecoders{jj,k};
    E_dj_star_k = fdcomm.DL_MMSE{jj,k};
    term_2 = H_ij'*(R_in_dj_k\H_Bj)*P_dj_k*E_dj_star_k...
        *P_dj_k'*H_Bj'*(R_in_dj_k\H_ij);
    %term_2 = nearestSPD(term_2);
%     napla_Pui_R_DL{jj,1} = -H_ij'/R_in_dj_k*H_Bj*P_dj_k/...
%         (eye(d_DL(jj))+term_1)*...
%         P_dj_k'*H_Bj'/R_in_dj_k*H_ij*tilde_P_iu_k;
    napla_Pui_R_DL{jj,1} = -term_2*tilde_P_iu_k;
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
Fu = fdcomm.Fu;
%% Calculate PiB
% Au_new = Fu\Au;
% Cu_new = Fu\Cu;
p_ui_k = (kron(eye(d_UL(ii)),Au)+kron(Bu.',Fu))\reshape(Cu,[],1);
P_ui_k = reshape(p_ui_k,[],d_UL(ii));
%P_ui_k = sylvester(Au_new,Bu,Cu_new);
fdcomm.ULprecoders{ii,k} = P_ui_k;
end
