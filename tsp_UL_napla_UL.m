function [ULcomm] = tsp_UL_napla_UL(k, ii, ULcomm, radar_UL, radar_comm_UL)

K = radar_UL.codelength;
% napla_PiB_R_UL_k = cell(I,1);
% napla_PiB_R_DL_k = cell(J,1);
H_iB = ULcomm.ULchannels{ii,1};
%% Au
xi_UL_k = ULcomm.xi_UL(:,:,k);
Au = H_iB'*xi_UL_k*H_iB;
ULcomm.Au_fixed = Au;
%% Bu
Bu = 0;
Cu = 0;
Nr = radar_UL.RX;
n = radar_UL.CUT_Idx;
d_iB = ULcomm.ULsymbols{ii};
for nr = 1:Nr
    eta_i_nr = radar_comm_UL.UL2rchannelgains(ii,nr);
    for m = 1:K
        if m == k
            d_ui_k_n = d_iB(:,1,m);
            Bu = Bu + eta_i_nr*d_ui_k_n*radar_UL.xi_r(k,k,nr)*d_ui_k_n';
        else
            d_ui_m_n = d_iB(:,n+1,m);
            Cu = Cu - eta_i_nr*d_ui_m_n*radar_UL.xi_r(m,k,nr)*d_ui_m_n';
        end
    end
end
ULcomm.Bu_fixed = Bu;
ULcomm.Cu_fixed = Cu;
end
