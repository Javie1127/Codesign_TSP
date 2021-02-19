function [fdcomm, radar, radar_comm, cov] = tsp_ini_normal_rd_radar(radar,fdcomm,radar_comm, snr_rtr)
%%%%-----------------------------------
%% Precoding matrices initialization
I = fdcomm.UL_num;
J = fdcomm.DL_num;
K = radar.codelength;
d_UL = fdcomm.ULstream_num;
d_DL = fdcomm.DLstream_num;
P_UL = fdcomm.ULpower;
P_UL_ini = cell(I,K); 
H_UL = fdcomm.ULchannels;
H_DL = fdcomm.DLchannels;
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
%% Initialization approach 1 right sigular matrix 
for ii = 1:I
    P_UL_i  = P_UL(ii); %% Power level of the ith UL UE  
    HiB     = H_UL{ii}; %% UL channel matrix with dimension Mc*Ni
    d_ui    = d_UL(ii);
    [~,~,V] = svd(HiB);
    P_iB_ini = V(:,1:d_ui);
    [U_UL_i,s,V_UL_i] = svd(P_iB_ini);
    s(1:d_ui,:) = diag(sqrt(P_UL_i/d_ui)*ones(d_ui,1));
    %S_sort = diag(s_sort_nonzero);
    P_iB_ini = U_UL_i*s*V_UL_i';
    for k = 1 : K
        P_UL_ini{ii,k} = P_iB_ini;
    end
end
P_DL_ini = cell(J,K);
P_DL = fdcomm.DLpower;
for jj = 1:J
    HBj = H_DL{jj}; %% UL channel matrix with dimension Ni*Mc
    %Mc = size(HBj,2);
    d_dj = d_DL(jj);
    [~,~,V] = svd(HBj);
    P_dj_ini = V(:,1:d_dj);
    [U_DL_j,s,V_DL_j] = svd(P_dj_ini);
    s(1:d_dj,:) = diag(sqrt(P_DL/J/d_dj)*ones(d_dj,1));
    %HBj_H = HBj';
    %P_Bj_ini = HBj_H(:,1:d_Bj);
    %[U_DL_j,~,V_DL_j] = svd(P_Bj_ini);
    %s_DL = diag(S_DL);
    %[~,Idx] = sort(s_DL,'descend');
    %U_sort_j = U_DL(:,Idx);
    %V_sort_j = V_DL(:,Idx);
    %s_sort_nonzero_j = sqrt(P_DL/J/d_Bj)*ones(d_Bj,1);
    %S_sort_j = [diag(s_sort_nonzero_j);zeros(Mc-d_Bj,d_Bj)];
%     P_Bj_ini = U_DL_j*S_sort_j*V_DL_j';
    P_Bj_ini = U_DL_j*s*V_DL_j';
    for k = 1 : K
        P_DL_ini{jj,k} = P_Bj_ini;
    end
end
fdcomm.ULprecoders = P_UL_ini;
fdcomm.DLprecoders = P_DL_ini;
Mr = radar.TX;
Pr = radar.Pr;
%% random coded radar codes 
radar.codematrix = randomcode(Pr,Mr,K);
%% comm rate
fdcomm.R_DL = 0.5;
fdcomm.R_UL = 0.1;
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
