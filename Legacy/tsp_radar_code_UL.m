function [radar_UL,cov_UL] = tsp_radar_code_UL(ULcomm, radar_UL, radar_comm_UL, cov_UL,k)
% radar code matrix update
I = ULcomm.UL_num;
Mr = radar_UL.TX;
K = radar_UL.codelength;
Nr = radar_UL.RX;
%eta_B = radar_comm.Bmrchannelgains;
eta_rt = radar_UL.channelgain;
H_rB = radar_comm_UL.radar2BSchannels;
%% Ar_k
xi_UL_k = ULcomm.xi_UL(:,:,k);
Ar_k = H_rB'*xi_UL_k*H_rB;
%% Br_k
Br_k = 0;
for nr = 1:Nr
    eta_rt_nr = eta_rt(nr);
    Br_k = Br_k + eta_rt_nr*radar_UL.xi_r(k,k,nr);
end
%% Derivative of UL rate w.r.t a_k 
d_UL = ULcomm.ULstream_num; % number of data streams of the UL UEs
napla_ak_R_UL_sum = 0;
tilde_ar_k = radar_UL.codematrix(k,:).';
for ii = 1 : I
    mu_iu_k = ULcomm.mu_UL(ii,k);
    R_in_iu_k = cov_UL.in_UL{ii,k};
    H_iB = ULcomm.ULchannels{ii,1};
    P_iu_k = ULcomm.ULprecoders{ii,k};
    napla_ak_R_UL_sum = napla_ak_R_UL_sum - mu_iu_k*H_rB'/R_in_iu_k*H_iB*P_iu_k/...
        (eye(d_UL(ii))+(P_iu_k'*H_iB'/R_in_iu_k*H_iB*P_iu_k))*...
        P_iu_k'*H_iB'/R_in_iu_k*H_rB*tilde_ar_k;
end
%% Cr_k
Cr_k = napla_ak_R_UL_sum;
Qr = radar_UL.doppler;
Jr = radar_comm_UL.Jr; 
Sigma = radar_UL.Sigma;
for nr = 1:Nr
    Sigma_nr = Sigma(:,:,nr);
    Qr_nr_k = diag(Qr(:,k,nr));
    Cr_k_temp = 0;
    for m = 1:K
        if m ~= k
            st_nr_m = cov_UL.S_tr(m,:,nr).';
            Cr_k_temp = Cr_k_temp + st_nr_m *radar_UL.xi_r(k,k,nr);
        end
    end
    Cr_k = -Qr_nr_k*Jr.'*Sigma_nr*Cr_k_temp + Cr_k;
end
% a_k = sylvester(Ar_k,Br_k,Cr_k);
a_k = (kron(1,Ar_k)+kron(Br_k.',eye(Mr)))\Cr_k;
radar_UL.codematrix(k,:) = a_k.'; 
end
