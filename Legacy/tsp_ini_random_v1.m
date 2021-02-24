function [fdcomm, radar, radar_comm, cov] = tsp_ini_random_v1(radar,fdcomm,radar_comm, snr_rtr)
%%%%-----------------------------------
%% Precoding matrices initialization
I = fdcomm.UL_num;
J = fdcomm.DL_num;
K = radar.codelength;
P_UL = fdcomm.ULpower;
P_UL_ini = cell(I,K); 
Mr = radar.TX;
Nr = radar.RX;
Nc = fdcomm.BSRX;
Mc = fdcomm.BSTX;
N_DL = fdcomm.DL_UE_Ant;
%% Updating radar channel matrices
M = Mc + Mr; 
eta_rtr = zeros(Mr,Nr);
Sigma = zeros(M,M,Nr);
eta_Btr = radar_comm.Btrchannelgains;
for nr = 1:Nr
    for mr = 1 : Mr
        eta_rtr(mr,nr) = db2pow(snr_rtr(mr,nr))*radar.noisepower;
        Sigma(:,:,nr) = blkdiag(eta_rtr(:,nr).*eye(Mr),eta_Btr(nr)*eye(Mc));
    end
end
radar.channelgain = eta_rtr;
radar.Sigma = Sigma;
%% radar_comm channels 
% radar TX - BS 
mu_r_BS = 0.5*ones(Mr,1);
kappa_BS = 0.5;
H_r_BS = zeros(Nc,Mr);
eta_rB = zeros(Mr,1);
for mr = 1 : Mr 
    eta_rB(mr) = sum(0.5*eta_rtr(mr,:));
    H_r_BS(:,mr) = sqrt(eta_rB(mr)/(kappa_BS+1)/2).*(randn(Nc,1)+1i*randn(Nc,1))+sqrt(kappa_BS/(kappa_BS+1))*mu_r_BS(mr)*ones(Nc,1); 
end
radar_comm.radar2BSchannels = H_r_BS;
radar_comm.r2Bchannelgains = eta_rB;
% radar to DL UEs
mu_DL = 0.2*ones(Mr,J);
kappa_DL = 0.3;
H_r_DL = cell(J,1);
eta_r_DL = zeros(Mr,J);
for jj = 1 : J
    H_r_j = zeros(N_DL(jj),Mr);
    for mr = 1:Mr
        eta_r_DL(mr,jj) = sum(0.7*eta_rtr(mr,:)); 
        H_r_j(:,mr) = sqrt(eta_r_DL(mr,jj)/(kappa_DL+1)/2)*(randn(N_DL(jj),1)+1i*randn(N_DL(jj),1))+sqrt(kappa_DL/(kappa_DL+1)*mu_DL(mr,jj)*ones(N_DL(jj),1));
    end
    H_r_DL{jj,1} = H_r_j;
end
radar_comm.radar2DLchannnels = H_r_DL;
radar_comm.radar2DLchannelgains = eta_r_DL;
%% Initialization approach 2 cholesky decompostion of the diagonal covariance matrix
d_UL = fdcomm.ULstream_num;
N_UL = fdcomm.UL_UE_Ant;
d_DL = fdcomm.DLstream_num;
for ii = 1:I
    P_UL_i  = P_UL(ii); %% Power level of the ith UL UE  
    d_ui    = d_UL(ii);
    N_ui    = N_UL(ii);
    Sigma_ui = [sqrt(P_UL_i/N_ui)*eye(d_ui);zeros(N_ui-d_ui,d_ui)];
    A1 = randn(N_ui) + 1i*randn(N_ui);
    [S_ui,~] = qr(A1);
    A2 = randn(d_ui) + 1i*randn(d_ui);
    [U_ui,~] = qr(A2);
    P_ui_ini = S_ui*Sigma_ui*U_ui';
    for k = 1 : K
        P_UL_ini{ii,k} = P_ui_ini;
    end
end
P_DL_ini = cell(J,K);
P_DL = fdcomm.DLpower;
for jj = 1:J
    d_dj = d_DL(jj);
    Sigma_dj = [sqrt(P_DL/J)*eye(d_dj);zeros(Mc-d_dj,d_dj)];
    A1 = randn(Mc) + 1i*randn(Mc);
    [S_dj,~] = qr(A1);
    A2 = randn(d_dj) + 1i*randn(d_dj);
    [U_dj,~] = qr(A2);
    P_dj_ini = S_dj*Sigma_dj*U_dj';
    for k = 1 : K
        P_DL_ini{jj,k} = P_dj_ini;
    end
end
Mr = radar.TX;
Pr = radar.Pr;
A_ini = zeros(K,Mr);
for mr = 1:Mr
    P_k = Pr(mr);
    A_ini(:,mr) = sqrt(P_k/K);
end
radar.codematrix = A_ini;
fdcomm.ULprecoders = P_UL_ini;
fdcomm.DLprecoders = P_DL_ini;
%% comm rate
fdcomm.R_DL = 0.3;
fdcomm.R_UL = 0.3;
%% Initializing the covariance matrices
cov = tsp_covmat(fdcomm,radar,radar_comm);
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
