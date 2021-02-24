function [ULcomm, radar_UL, cov,radar_comm] = tsp_ini_normal_UL(radar_UL,ULcomm,radar_comm,radar_fd_ini)
%%%%-----------------------------------
%% Precoding matrices initialization
I = ULcomm.UL_num;
K = radar_UL.codelength;
d_UL = ULcomm.ULstream_num;
P_UL = ULcomm.ULpower;
P_UL_ini = cell(I,K); 
H_UL = ULcomm.ULchannels;
radar_UL = radar_fd_ini;
Nr = radar_UL.RX;
Nc = ULcomm.BSRX;
%% Update radar_channel
radar_UL.channelgain = radar_fd_ini.channelgain;
eta_rtr = radar_UL.channelgain;
Mr = radar_UL.TX;
Sigma = zeros(Mr,Mr,Nr);
for nr = 1:Nr
    Sigma(:,:,nr) = eta_rtr(:,nr).*eye(Mr);
end
radar_UL.Sigma = Sigma;
%% Update UL -> Radar channel info
% radar TX - BS 
mu_r_BS = 0.5*ones(Mr,1);
kappa_BS = 0.5;
H_r_BS = zeros(Nc,Mr);
eta_rB = zeros(Mr,1);
for mr = 1 : Mr 
    eta_rB(mr) = sum(0.7*eta_rtr(mr,:));
    H_r_BS(:,mr) = sqrt(eta_rB(mr)/(kappa_BS+1)/2).*(randn(Nc,1)+1i*randn(Nc,1))+sqrt(kappa_BS/(kappa_BS+1))*mu_r_BS(mr)*ones(Nc,1); 
end
radar_comm.radar2BSchannels = H_r_BS;
radar_comm.r2Bchannelgains = eta_rB;
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
ULcomm.ULprecoders = P_UL_ini;
%% radar waveform
radar_UL.codematrix = radar_fd_ini.codematrix;
%% comm rate
ULcomm.R_UL = 0.3;
%% Initializing the covariance matrices
cov = tsp_covmat_UL(ULcomm,radar_UL,radar_comm);
%% Initializing the MMSE matrices
ULcomm = tsp_Comm_MMSE_UL(ULcomm,radar_UL,cov);
radar_UL = tsp_radar_MMSE(radar_UL,cov);
%% Initializing the performance measures
radar_UL = Xi_radar(radar_UL);
for k = 1:K
    ULcomm = Xi_comm_k_UL(ULcomm,k);
end
%% Alternating optimization
radar_UL.iota_max = 8;% algorithm 3
ULcomm.tu_max = 7; % algorithm 1 max number of iterations to execute the UL subgradient method
%% Subgradient method Algorithm 1 & 2
ULcomm.lambda_UL = ones(I,K);
ULcomm.mu_UL = ones(I,K);
end
