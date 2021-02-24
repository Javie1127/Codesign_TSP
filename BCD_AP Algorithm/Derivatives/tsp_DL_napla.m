function [fdcomm] = tsp_DL_napla(k, jj, fdcomm, radar, radar_comm)
% I = fdcomm.UL_num;
K = radar.codelength;
J = fdcomm.DL_num;
% napla_PiB_R_UL_k = cell(I,1);
% napla_PiB_R_DL_k = cell(J,1);
H_BB = fdcomm.BBchannel;
%% Ad
xi_UL_k = fdcomm.xi_UL(:,:,k);
xi_DL = fdcomm.xi_DL; 
term_1 = 0;
for g = 1:J
    H_Bg = fdcomm.DLchannels{g,1};
    term_1 = term_1 + H_Bg'*xi_DL{k,g}*H_Bg;
end
Ad = H_BB'*xi_UL_k*H_BB+term_1;
fdcomm.Ad_fixed = Ad;
%% Bd
d_DL_j = fdcomm.DLsymbols{jj,1};
n_Bm = radar_comm.n_Bm;
d_Bt_j_k = d_DL_j(:,1,k);
B_Bt_j_k = d_Bt_j_k*d_Bt_j_k';
d_Bm_j_k = d_DL_j(:,n_Bm,k);
B_Bm_j_k = d_Bm_j_k*d_Bm_j_k';
Nr = radar.Rx;

F_Bt_k = 0;
F_Bm_k = 0;
for nr = 1:Nr
    if isfield(radar_comm,'Sigma_Bt_Nr')
        Sigma_Bt_nr_k_k = radar_comm.Sigma_Bt_Nr{k,k,nr};
    else
        Sigma_Bt_nr_k_k = zeros(fdcomm.BSTx,fdcomm.BSTx);
    end
    F_Bt_k  = radar.xi_r(k,k,nr)*Sigma_Bt_nr_k_k+F_Bt_k;
    F_Bm_k  = radar.xi_r(k,k,nr)*radar_comm.Sigma_Bm_Nr{k,k,nr}+F_Bm_k;
end
fdcomm.B_Bm_j_k = B_Bm_j_k;
fdcomm.B_Bt_j_k = B_Bt_j_k;
fdcomm.F_Bt_k = F_Bt_k;
fdcomm.F_Bm_k = F_Bm_k;
%% Cd_fixed
JB = radar_comm.JB;
m_g_k_t = 0;
m_g_k_Bm = 0;
m_j_k_t = 0;
m_j_k_Bm = 0;
Delta_h_Bt = 0;
d_dj_k_t = d_DL_j(:,1,k);
d_dj_k_Bm = d_DL_j(:,n_Bm,k);
JH = radar_comm.JH;
for nr = 1:Nr
    Urnr = radar.WMMSE_RX{nr,1};
    urnr_k = Urnr(:,k);
    Wrnr = radar.WMMSE_weights{nr,1};
    m_g_k_t_nr = 0;
    m_g_k_Bm_nr = 0;
    m_j_k_t_nr = 0;
    m_j_k_Bm_nr = 0;
    Delta_h_Bt_nr = 0;
    for m = 1:K
        xi_r_nr_k = radar.xi_r(m,k,nr);
        if isfield(radar_comm,'Sigma_Bt_Nr')
            Sigma_Bt_nr_m_k = radar_comm.Sigma_Bt_Nr{m,k,nr};
        else
            Sigma_Bt_nr_m_k = zeros(fdcomm.BSTx,fdcomm.BSTx);
        end
        Sigma_Bm_nr_m_k = radar_comm.Sigma_Bm_Nr{m,k,nr};
        m_g_k_t_nr_m = 0;
        m_g_k_Bm_nr_m = 0;
        for g = 1:J
            d_dg = fdcomm.DLsymbols{g,1};
            if g == jj && m ~= k
                d_dj_m_Bm = d_dg(:,n_Bm,m);
                d_dj_m_Bt = d_dg(:,1,m);
                P_dj_m = fdcomm.DLprecoders{g,m};
                m_j_k_t_nr = m_j_k_t_nr+xi_r_nr_k*Sigma_Bt_nr_m_k*P_dj_m*d_dj_m_Bm*d_dj_k_Bm';
                m_j_k_Bm_nr = m_j_k_Bm_nr + xi_r_nr_k*Sigma_Bm_nr_m_k*P_dj_m*d_dj_m_Bt*d_dj_k_t';
            else
                d_dg_m_Bm = d_dg(:,n_Bm,m);
                d_dg_m_Bt = d_dg(:,1,m);
                P_dg_m = fdcomm.DLprecoders{g,m};
                m_g_k_t_nr_m = m_g_k_t_nr_m + xi_r_nr_k*Sigma_Bt_nr_m_k*P_dg_m*d_dg_m_Bt*d_dj_k_t';
                m_g_k_Bm_nr_m = m_g_k_Bm_nr_m + xi_r_nr_k*Sigma_Bm_nr_m_k*P_dg_m*d_dg_m_Bm*d_dj_k_Bm';
            end
        end
        m_g_k_t_nr = m_g_k_t_nr_m + m_g_k_t_nr;
        m_g_k_Bm_nr = m_g_k_Bm_nr + m_g_k_Bm_nr_m;
        if isfield(radar_comm,'Sigma_Bt_Nr')
            Sigma_Bt_nr_m_k = radar_comm.Sigma_Bt_Nr{m,k,nr};
        else
            Sigma_Bt_nr_m_k = zeros(fdcomm.BSTx,fdcomm.BSTx);
        end
        Delta_h_Bt_nr = Delta_h_Bt_nr+Sigma_Bt_nr_m_k*JB.'*JH{m}.'*Wrnr*urnr_k;
    end
    Delta_h_Bt = Delta_h_Bt + Delta_h_Bt_nr;
    m_g_k_t = m_g_k_t_nr + m_g_k_t;
    m_g_k_Bm = m_g_k_Bm + m_g_k_Bm_nr;
    m_j_k_Bm = m_j_k_Bm + m_j_k_Bm_nr;
    m_j_k_t = m_j_k_t + m_j_k_t_nr;
end
fdcomm.Cd_fixed = 2*(Delta_h_Bt-m_j_k_t-m_j_k_Bm-m_g_k_t-m_g_k_Bm);


end
