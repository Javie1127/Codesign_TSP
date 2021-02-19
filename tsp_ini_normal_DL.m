function [DLcomm, radar_DL, cov_DL,radar_comm] = ...
    tsp_ini_normal_DL(radar_DL,DLcomm,radar_comm,radar_fd_ini)
%%%%-----------------------------------
%% Updating radar channel matrices
radar_DL.channelgain = radar_fd_ini.channelgain;
radar_DL.Sigma = radar_fd_ini.Sigma;
%% Precoding matrices initialization
J = DLcomm.DL_num;
K = radar_DL.codelength;
d_DL = DLcomm.DLstream_num;
H_DL = DLcomm.DLchannels;
Mr = radar_DL.TX;
%% Updating radar -> DL UEs
% radar to DL UEs
eta_rtr = radar_DL.channelgain;
mu_DL = 0.2*ones(Mr,J);
kappa_DL = 0.3;
H_r_DL = cell(J,1);
eta_r_DL = zeros(Mr,J);
N_DL = DLcomm.DL_UE_Ant;
for jj = 1 : J
    H_r_j = zeros(N_DL(jj),Mr);
    for mr = 1:Mr
        eta_r_DL(mr,jj) = sum(0.5*eta_rtr(mr,:)); 
        H_r_j(:,mr) = sqrt(eta_r_DL(mr,jj)/(kappa_DL+1)/2)*(randn(N_DL(jj),1)+1i*randn(N_DL(jj),1))+sqrt(kappa_DL/(kappa_DL+1)*mu_DL(mr,jj)*ones(N_DL(jj),1));
    end
    H_r_DL{jj,1} = H_r_j;
end
radar_comm.radar2DLchannnels = H_r_DL;
radar_comm.radar2DLchannelgains = eta_r_DL;

%% Initialization approach 1 right sigular matrix 
P_DL_ini = cell(J,K);
P_DL = DLcomm.DLpower;
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
DLcomm.DLprecoders = P_DL_ini;
%% Radar waveform
radar_DL.codematrix = radar_fd_ini.codematrix;
%% comm rate
DLcomm.R_DL = 0.3;
%% Initializing the covariance matrices
cov_DL = tsp_covmat_DL(DLcomm,radar_DL,radar_comm);
%% Initializing the MMSE matrices
DLcomm = tsp_Comm_MMSE_DL(DLcomm,radar_DL,cov_DL);
radar_DL = tsp_radar_MMSE(radar_DL,cov_DL);
%% Initializing the performance measures
radar_DL = Xi_radar(radar_DL);
for k = 1:K
    DLcomm = Xi_comm_k_DL(DLcomm,k);
end
%% Alternating optimization
radar_DL.iota_max = 8;% algorithm 3
DLcomm.td_max = 7; % algorithm 2 % max number of iterations to execute the DL suggradient method
%% Subgradient method Algorithm 1 & 2
DLcomm.lambda_DL = ones(K,1);
DLcomm.mu_DL = ones(J,K);

end
