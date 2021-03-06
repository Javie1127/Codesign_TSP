function [DLcomm] = tsp_DL_napla_DL(k, jj, DLcomm, radar_DL, radar_comm_DL, cov_DL)
% I = fdcomm.UL_num;
K = radar_DL.codelength;
J = DLcomm.DL_num;
% napla_PiB_R_UL_k = cell(I,1);
% napla_PiB_R_DL_k = cell(J,1);
%% Ad
xi_DL = DLcomm.xi_DL; 
term_1 = 0;
for g = 1:J
    H_Bg = DLcomm.DLchannels{g,1};
    term_1 = term_1 + H_Bg'*xi_DL{k,g}*H_Bg;
end
Ad = term_1;
DLcomm.Ad_fixed = Ad;
%% Bd
Bd = 0;
Nr = radar_DL.RX;
n = radar_DL.CUT_Idx;
d_DL_j = DLcomm.DLsymbols{jj,1};
for nr = 1:Nr
    eta_Bt_nr = radar_comm_DL.Btrchannelgains(nr);
    eta_Bm_nr = radar_comm_DL.Bmrchannelgains(nr);
    d_dj_k_one = d_DL_j(:,1,k);
    d_dj_k_n = d_DL_j(:,n+1,k);
    Bd = Bd + eta_Bt_nr*d_dj_k_one*radar_DL.xi_r(k,k,nr)*d_dj_k_one'...
                +eta_Bm_nr*d_dj_k_n*radar_DL.xi_r(k,k,nr)*d_dj_k_n';
end
%% Cd_fixed
JB = radar_comm_DL.JB;
term_1 = 0;
term_2 = 0;
term_3 = 0;
d_dj_k_one = d_DL_j(:,1,k);
d_dj_k_n = d_DL_j(:,n+1,k);
for nr = 1:Nr
    Sigma_nr = radar_DL.Sigma(:,:,nr);
    eta_Bm_nr = radar_comm_DL.Bmrchannelgains(nr);
    term_1_temp = 0;
    term_2_temp = 0;
    term_3_temp = 0;
    for m = 1:K
        term_2_tempp = 0;
        st_nr_m = cov_DL.S_tr(m,:,nr).';
        term_1_temp = term_1_temp + st_nr_m*radar_DL.xi_r(m,k,nr)*d_dj_k_one';
        d_dg = DLcomm.DLsymbols{g,1};
        d_dg_m_n = d_dg(:,n+1,m);
        P_dg_m = DLcomm.DLprecoders{g,m};
        if m ~= K
            for g = 1:J
                term_2_tempp = term_2_tempp + P_dg_m*d_dg_m_n*radar_DL.xi_r(m,k,nr)*d_dj_k_n';
            end
            term_2_temp = term_2_tempp + term_2_temp;
        else
            for g = 1:J
                term_3_temp = term_3_temp + P_dg_m*d_dj_k_n*radar_DL.xi_r(k,k,nr)*d_dj_k_n';
            end
        end
        term_1 = JB.'*Sigma_nr*term_1_temp+term_1;
        term_2 = eta_Bm_nr^2*term_2_temp + term_2;
        term_3 = eta_Bm_nr^2*term_3_temp + term_3;
    end
end
DLcomm.Cd_fixed = -term_1-term_2-term_3;
DLcomm.Bd_fixed = Bd;
end
