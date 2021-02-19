function [fdcomm, radar, radar_comm, cov] = tsp_ini_wc(radar,fdcomm,radar_comm, cov)
%%%%-----------------------------------
%% Precoding matrices initialization
I = fdcomm.UL_num;
J = fdcomm.DL_num;
K = radar.codelength;
Mr = radar.TX;
Nr = radar.RX;
Mc = fdcomm.BSTX;
%% Updating Sigma
M = Mc + Mr; 
Sigma = zeros(M,M,Nr);
eta_Btr = zeros(Nr,1);
radar_comm.Btrchannelgains = eta_Btr;
eta_rtr = radar.channelgain;
for nr = 1:Nr
    Sigma(:,:,nr) = blkdiag(eta_rtr(:,nr).*eye(Mr));
end
radar.Sigma = Sigma;
%% Initializing the covariance matrices
S_rtr = cov.S_rtr;
S_Btr = cov.S_Btr;
S_tr = cov.S_tr;
R_tr = zeros(K,K,Nr); %Covariance for the radar
for nr = 1:Nr
    eta_rtnr = eta_rtr(nr);
    eta_Btnr = eta_Btr(nr);
    Sigma_tnr = Sigma(:,:,nr);
    for m = 1:K 
        s_rt_nr_m = S_rtr(m,:,nr).';
        s_Bt_nr_m = S_Btr(m,:,nr).';
        s_tnr_m = S_tr(m,:,nr).';
        for l = 1:K
            s_rt_nr_l = S_rtr(l,:,nr).';
            s_Bt_nr_l = S_Btr(l,:,nr).';
            R_tr(m,l,nr) = eta_rtnr*s_rt_nr_l'*s_rt_nr_m+eta_Btnr*s_Bt_nr_l'*s_Bt_nr_m;
            s_tnr_l = S_tr(l,:,nr).';
            R_tr(m,l,nr) = trace(s_tnr_m*s_tnr_l'*Sigma_tnr);
        end
    end
    R_tr(:,:,nr) = R_tr(:,:,nr);
end
cov.target2radar = R_tr;
R_r = zeros(K,K,Nr);
R_in_r = cov.inr;
for nr = 1:Nr 
    R_r(:,:,nr) = R_in_r(:,:,nr)+R_tr(:,:,nr);
end
cov.total_r = R_r;
%% Initializing the MMSE matrices
fdcomm = tsp_Comm_MMSE(fdcomm,radar,cov);
radar = tsp_radar_MMSE(radar,cov);
%% Initializing the performance measures
radar = Xi_radar(radar);
for k = 1:K
    fdcomm = Xi_comm_k(fdcomm,k);
end
%% Alternating optimization
radar.iota_max = 8;% algorithm 3
fdcomm.tu_max = 7; % algorithm 1 max number of iterations to execute the UL subgradient method
fdcomm.td_max = 7; % algorithm 2 % max number of iterations to execute the DL suggradient method
%% Subgradient method Algorithm 1 & 2
fdcomm.lambda_UL = ones(I,K);
fdcomm.lambda_DL = ones(K,1);
fdcomm.mu_UL = ones(I,K);
fdcomm.mu_DL = ones(J,K);
end
