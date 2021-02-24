function [DLcomm, radar_DL, radar_comm_DL] = tsp_parameters_DL(SNR_DL,radar_fd,DLcomm)
%%%%-----------------------------------
%------- SNR.SNR_rtr: Mr*Nr real matrix
%------- SNR.SNR_DL:
radar_DL = radar_fd;
Mr = radar_fd.TX;
Nr = radar_fd.RX;
sigma_0 = radar_fd.noisepower;
Mc = DLcomm.BSTX;
K = radar_fd.codelength;
N = radar_fd.PRI_num;
%% FD comm parameters
J= DLcomm.DL_num;
N_DL = 2*ones(J,1); % Number of the DL UE antennas
% number of stream for each DL UE d_DL(i) <= N_DL(i) 
d_DL = N_DL;
DLcomm.DL_UE_Ant = N_DL;
DLcomm.DLstream_num = d_DL;
%% FD Comm channels
% BS - DL_UEs
H_DL = cell(J,1); % J DL Channels
eta_DL = ones(J,1); %
for jj = 1:J
    eta_DL(jj) = db2pow(SNR_DL.BS_DL(jj))*sigma_0;
    H_DL{jj,1} = sqrt(eta_DL(jj)/2/N_DL(jj))*(randn(N_DL(jj),Mc)+1i*randn(N_DL(jj),Mc)); 
end
DLcomm.DLchannelgains = eta_DL;
DLcomm.DLchannels = H_DL;
%% FD comm symbols 
D_DL = cell(J,1);
for jj = 1:J
    d_j = d_DL(jj);
    frames = zeros(d_j, N, K);
    for kk = 1:K
        for nn = 1:N
            v = randn(d_j,2);
            v = bsxfun(@rdivide,v,sqrt(sum(v.^2,2)));
            frames(:,nn,kk) = complex(v(:,1),v(:,2));
        end
    end
    D_DL{jj,1} = frames;
end
DLcomm.DLsymbols = D_DL;
%% radar_comm channels 
% % BS to DL UEs
% mu_DL = 0.2*ones(Mr,J);
% kappa_DL = 0.3;
% H_r_DL = cell(J,1);
% eta_r_DL = zeros(Mr,J);
% for jj = 1 : J
%     H_r_j = zeros(N_DL(jj),Mr);
%     for mr = 1:Mr
%         eta_r_DL(mr,jj) = db2pow(SNR_DL.r_DL(mr,jj))*sigma_0; 
%         H_r_j(:,mr) = sqrt(eta_r_DL(mr,jj)/(kappa_DL+1)/2)*(randn(N_DL(jj),1)+1i*randn(N_DL(jj),1))+sqrt(kappa_DL/(kappa_DL+1)*mu_DL(mr,jj)*ones(N_DL(jj),1));
%     end
%     H_r_DL{jj,1} = H_r_j;
% end
% radar_comm_DL.radar2DLchannnels = H_r_DL;
% radar_comm_DL.radar2DLchannelgains = eta_r_DL;
% BS - target - radar RX and BS - multi-path - radar RX
eta_Btr = zeros(Nr,1);
H_Btr = zeros(Mc,Nr);
eta_Bmr = zeros(Nr,1);
H_Bmr = zeros(Mc,Nr);
% for nr = 1:Nr
%     % SNR.Btr : Nr * 1
%     eta_Btr(nr) = db2pow(SNR_DL.Btr(nr))*sigma_0;
%     eta_Bmr(nr) = db2pow(SNR_DL.Bmr(nr))*sigma_0;
%     H_Bmr(:,nr) = sqrt(eta_Bmr(nr)/2).*(randn(Mc,1)+1i*randn(Mc,1));
%     H_Btr(:,nr) = sqrt(eta_Btr(nr)/2).*(randn(Mc,1)+1i*randn(Mc,1));
%     Sigma(:,:,nr) = blkdiag(eta_rtr(:,nr).*eye(Mr),eta_Btr(nr)*eye(Mc));
% end
for nr = 1:Nr
    eta_Btr(nr) = 0.5*sum(eta_DL)/Nr;
    eta_Bmr(nr) = 0.7*sum(eta_DL)/Nr;
    H_Bmr(:,nr) = sqrt(eta_Bmr(nr)/2).*(randn(Mc,1)+1i*randn(Mc,1));
    H_Btr(:,nr) = sqrt(eta_Btr(nr)/2).*(randn(Mc,1)+1i*randn(Mc,1));
end
radar_comm_DL.Bmrchannels = H_Bmr;
radar_comm_DL.Bmrchannelgains = eta_Bmr;
radar_comm_DL.Btrchannels = H_Btr;
radar_comm_DL.Btrchannelgains = eta_Btr;
radar_comm_DL.Jr = [eye(Mr);zeros(Mc,Mr)];
radar_comm_DL.JB = [zeros(Mr,Mc);eye(Mc)];
%% radar_comm doppler 
kk = 0:K-1;
f_Bm_Nr = zeros(Nr,1);
Q_Btr = zeros(K,Nr);
Q_Bmr = zeros(K,Nr);
f_Bt_Nr = zeros(Nr,1);
for nr = 1:Nr
    f_Bt_nr = 0.15*(-1*rand(1)+0.5)+0.25; % Normalized Doppler frequency
    f_Bm_nr = 0.15*(-1*rand(1)+0.5)+0.25;
    Q_Btr(:,nr) = exp(1j*2*pi*f_Bt_nr.*kk'); % Doppler Domain Steering vector 
    Q_Bmr(:,nr) = exp(1j*2*pi*f_Bm_nr.*kk');
    f_Bt_Nr(nr) = f_Bt_nr;
    f_Bm_Nr(nr) = f_Bm_nr;
end
radar_comm_DL.f_Bt_Nr = f_Bt_Nr;
radar_comm_DL.f_Bm_Nr = f_Bm_Nr;
radar_comm_DL.Btrdoppler = Q_Btr;
radar_comm_DL.Bmrdoppler = Q_Bmr;
end
