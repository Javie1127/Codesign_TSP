function [fdcomm] = tsp_UL_napla_rev(k, ii, fdcomm, radar, radar_comm)

K = radar.codelength;
J = fdcomm.DL_num;
% napla_PiB_R_UL_k = cell(I,1);
% napla_PiB_R_DL_k = cell(J,1);
H_iB = fdcomm.ULchannels{ii,1};
%% Au
xi_UL_k = fdcomm.xi_UL(:,:,k);
xi_DL = fdcomm.xi_DL; 
term_1 = 0;
for g = 1:J
    H_ig = fdcomm.ULDLchannels{ii,g};
    term_1 = term_1 + nearestSPD(H_ig'*xi_DL{k,g}*H_ig);
end
Au_i_k = nearestSPD(H_iB'*xi_UL_k*H_iB)+term_1;
fdcomm.Au_fixed = Au_i_k;
%% Bu
Nr = radar.Rx;
nu = radar_comm.nu;
d_iB = fdcomm.ULsymbols{ii};
Bu_i_k = d_iB(:,nu,k)*d_iB(:,nu,k)';
%% Cu & Fu
Cu_i_k = 0;
Fu_i_k = 0;
Sigma_U_Nr = radar_comm.Sigma_U_Nr;
d_ui_k = d_iB(:,nu,k);
for nr = 1:Nr
    for m = 1:K
        if m == k
            Fu_i_k = radar.xi_r(m,k,nr)*Sigma_U_Nr{k,k,ii,nr}+Fu_i_k;
        else
            d_ui_m = d_iB(:,nu,m);
            Cu_i_k = Cu_i_k - Sigma_U_Nr{m,k,ii,nr}*fdcomm.ULprecoders{ii,k}*d_ui_m*radar.xi_r(m,k,nr)*d_ui_k';
        end
    end
end
fdcomm.Bu_fixed = Bu_i_k;
fdcomm.Cu_fixed = Cu_i_k;
fdcomm.Fu = Fu_i_k;
end
