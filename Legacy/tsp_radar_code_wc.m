function [radar,cov] = tsp_radar_code_wc(fdcomm, radar, radar_comm, cov,k)
% radar code matrix update
I = fdcomm.UL_num;
J = fdcomm.DL_num;
Mr = radar.TX;
K = radar.codelength;
Nr = radar.RX;
%eta_B = radar_comm.Bmrchannelgains;
eta_rt = radar.channelgain;
H_rB = radar_comm.radar2BSchannels;
%% Ar_k
xi_UL_k = fdcomm.xi_UL(:,:,k);
xi_DL = fdcomm.xi_DL; 
Ar_k = H_rB'*xi_UL_k*H_rB;
for g = 1:J
    Hrg = radar_comm.radar2DLchannnels{g};
    Ar_k = Ar_k + Hrg'*xi_DL{k,g}*Hrg;
end
%% Br_k
Br_k = 0;
for nr = 1:Nr
    eta_rt_nr = eta_rt(nr);
    Br_k = Br_k + eta_rt_nr*radar.xi_r(k,k,nr);
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
Qr = radar.doppler;
Sigma = radar.Sigma;
for nr = 1:Nr
    Sigma_nr = Sigma(:,:,nr);
    Qr_nr_k = diag(Qr(:,k,nr));
    Cr_k_temp = 0;
    for m = 1:K
        if m ~= k
            st_nr_m = cov.S_tr(m,:,nr).';
            Cr_k_temp = Cr_k_temp + st_nr_m *radar.xi_r(k,k,nr);
        end
    end
    Cr_k = -Qr_nr_k*Sigma_nr*Cr_k_temp + Cr_k;
end
% a_k = sylvester(Ar_k,Br_k,Cr_k);
a_k = (kron(1,Ar_k)+kron(Br_k.',eye(Mr)))\Cr_k;
radar.codematrix(k,:) = a_k.'; 
end
