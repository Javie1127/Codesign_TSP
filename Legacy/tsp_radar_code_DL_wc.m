function [radar,cov] = tsp_radar_code_DL_wc(DLcomm, radar, radar_comm, cov,k)
% radar code matrix update
J = DLcomm.DL_num;
Mr = radar.TX;
K = radar.codelength;
Nr = radar.RX;
%eta_B = radar_comm.Bmrchannelgains;
eta_rt = radar.channelgain;
%% Ar_k
xi_DL = DLcomm.xi_DL; 
Ar_k = 0;
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
napla_ak_R_DL_sum = 0;
d_DL = DLcomm.DL_UE_Ant;
tilde_ar_k = radar.codematrix(k,:).';
for jj = 1:J
    mu_jd_k = DLcomm.mu_DL(jj,k);
    H_rj = radar_comm.radar2DLchannnels{jj};
    R_in_jd_k = cov.in_DL{jj,k};
    H_Bj = DLcomm.DLchannels{jj,1};
    P_jd_k = DLcomm.DLprecoders{jj,k};
    napla_ak_R_DL_sum = napla_ak_R_DL_sum-mu_jd_k*H_rj'/R_in_jd_k*H_Bj*P_jd_k/...
        (eye(d_DL(jj))+(P_jd_k'*H_Bj'/R_in_jd_k*H_Bj*P_jd_k))*...
        P_jd_k'*H_Bj'/R_in_jd_k*H_rj*tilde_ar_k;
end
%% Cr_k
Cr_k = napla_ak_R_DL_sum;
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
