function [fdcomm] = tsp_UL_napla(k, ii, fdcomm, radar, radar_comm)

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
    term_1 = term_1 + H_ig'*xi_DL{k,g}*H_ig;
end
Au = H_iB'*xi_UL_k*H_iB+term_1;
fdcomm.Au_fixed = Au;
%% Bu
Bu = 0;
Cu = 0;
Nr = radar.RX;
n = radar.CUT_Idx;
d_iB = fdcomm.ULsymbols{ii};
for nr = 1:Nr
    eta_i_nr = radar_comm.UL2rchannelgains(ii,nr);
    for m = 1:K
        if m == k
            d_ui_k_n = d_iB(:,1,m);
            Bu = Bu + eta_i_nr*d_ui_k_n*radar.xi_r(k,k,nr)*d_ui_k_n';
        else
            d_ui_m_n = d_iB(:,n+1,m);
            Cu = Cu - eta_i_nr*d_ui_m_n*radar.xi_r(m,k,nr)*d_ui_m_n';
        end
    end
end
fdcomm.Bu_fixed = Bu;
fdcomm.Cu_fixed = Cu;
end
