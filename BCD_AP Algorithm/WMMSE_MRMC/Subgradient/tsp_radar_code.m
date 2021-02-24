function [radar,cov] = tsp_radar_code(fdcomm, radar, radar_comm, cov,k)
% radar code matrix update
I = fdcomm.UL_num;
J = fdcomm.DL_num;
Mr = radar.Tx;
K = radar.codelength;
Nr = radar.Rx;
H_rB = radar_comm.radar2BSchannels;
%% Ar_k
xi_UL_k = fdcomm.xi_UL(:,:,k);
xi_DL = fdcomm.xi_DL; 
Ar_k = H_rB'*xi_UL_k*H_rB;
for g = 1:J
    Hrg = radar_comm.radar2DLchannnels{g};
    Ar_k = Ar_k + Hrg'*xi_DL{k,g}*Hrg;
end

%% Derivative of UL rate w.r.t a_k 
d_UL = fdcomm.ULstream_num; % number of data streams of the UL UEs
d_DL = fdcomm.DLstream_num; % number of data streams of the DL UEs
napla_ak_R_UL_sum = 0;
tilde_ar_k = radar.codematrix(k,:).';
for ii = 1 : I
    mu_iu_k = fdcomm.mu_UL(ii,k);
    R_in_iu_k = cov.in_UL{ii,k};
    H_iB = fdcomm.ULchannels{ii,1};
    P_iu_k = fdcomm.ULprecoders{ii,k};
    napla_ak_R_UL_sum = napla_ak_R_UL_sum - mu_iu_k*H_rB'/R_in_iu_k*H_iB*P_iu_k/...
        (eye(d_UL(ii))+(P_iu_k'*H_iB'/R_in_iu_k*H_iB*P_iu_k))*...
        P_iu_k'*H_iB'/R_in_iu_k*H_rB*tilde_ar_k;
end
napla_ak_R_DL_sum = 0;
for jj = 1:J
    mu_jd_k = fdcomm.mu_DL(jj,k);
    H_rj = radar_comm.radar2DLchannnels{jj};
    R_in_jd_k = cov.in_DL{jj,k};
    H_Bj = fdcomm.DLchannels{jj,1};
    P_jd_k = fdcomm.DLprecoders{jj,k};
    napla_ak_R_DL_sum = napla_ak_R_DL_sum-mu_jd_k*H_rj'/R_in_jd_k*H_Bj*P_jd_k/...
        (eye(d_DL(jj))+(P_jd_k'*H_Bj'/R_in_jd_k*H_Bj*P_jd_k))*...
        P_jd_k'*H_Bj'/R_in_jd_k*H_rj*tilde_ar_k;
end
%% Cr_k
Cr_k = napla_ak_R_DL_sum + napla_ak_R_UL_sum;
Jr = radar_comm.Jr; 
JH = radar_comm.JH;
Sigma_c = radar.Clutter_channel_cov;
Sigma_rt = radar.Sigma_h_tr;
Delta_h_rt = 0;
m_r = 0;
A = radar.codematrix;
F_r_k = 0;
for nr = 1:Nr
    Urnr = radar.WMMSE_RX{nr,1};
    urnr_k = Urnr(:,k);
    Wrnr = radar.WMMSE_weights{nr,1};
    Delta_h_rt_nr = 0;
    m_r_nr = 0;
    for m = 1:K
        Sigma_rt_nr_m_k = Sigma_rt(m,k,nr);
        xi_r_nr_k = radar.xi_r(m,k,nr);
        if m ~= k
            xi_r_nr_k = radar.xi_r(m,k,nr);
            m_r_nr = xi_r_nr_k*(Sigma_rt_nr_m_k+Sigma_c)*A(m,:).'+m_r_nr;
        end
        Delta_h_rt_nr = Delta_h_rt_nr+Sigma_rt_nr_m_k*Jr.'*JH{m}.'*Wrnr*urnr_k;
    end
    Delta_h_rt = Delta_h_rt + Delta_h_rt_nr;
    m_r = m_r + m_r_nr;
    Sigma_rt_nr_k_k = Sigma_rt(k,k,nr);
    F_r_k = xi_r_nr_k*(Sigma_rt_nr_k_k+Sigma_c)+F_r_k;
end
Cr_k = Cr_k-2*m_r+2*Delta_h_rt;
% a_k = sylvester(Ar_k,Br_k,Cr_k);
a_k = (kron(1,F_r_k\Ar_k)+kron(1,eye(Mr)))\Cr_k;
radar.codematrix(k,:) = a_k.'; 
end
