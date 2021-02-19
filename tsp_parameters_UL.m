function [ULcomm, radar_UL, radar_comm_UL] = ...
    tsp_parameters_UL(SNR,radar_fd,ULcomm)
%%%%-----------------------------------
%------- SNR.SNR_rtr: Mr*Nr real matrix
%------- SNR.SNR_DL:
%% 
radar_UL = radar_fd;
Mr = radar_fd.TX;
Nr = radar_fd.RX;
%% FD comm parameters
I = ULcomm.UL_num;
N_UL = 2*ones(I,1); % Number of the UL UE antennas
% number of stream for each DL UE d_DL(i) <= N_DL(i) 
d_UL = N_UL;
ULcomm.UL_UE_Ant = N_UL;
ULcomm.ULstream_num = d_UL;
sigma_0 = radar_fd.noisepower;
Nc = ULcomm.BSRX;
%% FD Comm channels
% UL UEs - BS
eta_UL = zeros(I,1);
H_UL = cell(I,1);
for ii = 1:I
   eta_UL(ii) = db2pow(SNR.UL_BS(ii)) *sigma_0;
   H_UL{ii,1} = sqrt(eta_UL(ii)/2)*(randn(Nc,N_UL(ii))+1i*randn(Nc,N_UL(ii)));
end
ULcomm.ULchannelgains = eta_UL;
ULcomm.ULchannels = H_UL;

%% UL comm symbols 
D_UL = cell(I,1);
K = radar_fd.codelength;
N = radar_fd.PRI_num;
for ii = 1:I
    d_i = d_UL(ii);
    frames = zeros(d_i, N, K);
    for kk = 1:K
        for nn = 1:N
            v = randn(d_i, 2);
            v = bsxfun(@rdivide,v,sqrt(sum(v.^2,2)));
            frames(:,nn,kk) = complex(v(:,1),v(:,2));
        end
    end
    D_UL{ii,1} = frames;
end
ULcomm.ULsymbols = D_UL;
% %% radar_comm channels 
% % radar TX - BS 
% Mr = radar_co.TX;
% Nr = radar_co.RX;
% mu_r_BS = 0.5*ones(Mr,1);
% kappa_BS = 0.5;
% H_r_BS = zeros(Nc,Mr);
% eta_rB = zeros(Mr,1);
% for mr = 1 : Mr 
%     eta_rB(mr) = db2pow(SNR.r_B(mr))*sigma_0;
%     H_r_BS(:,mr) = sqrt(eta_rB(mr)/(kappa_BS+1)/2).*(randn(Nc,1)+1i*randn(Nc,1))+sqrt(kappa_BS/(kappa_BS+1))*mu_r_BS(mr)*ones(Nc,1); 
% end
% radar_comm_UL.radar2BSchannels = H_r_BS;
% radar_comm_UL.r2Bchannelgains = eta_rB;
radar_comm_UL.Jr = eye(Mr);
radar_comm_UL.JB = 0;
%% UL -> radar
eta_UL_r = zeros(I,Nr);
H_UL_r   = cell(I,Nr);
for nr = 1 : Nr 
    for ii = 1:I
        eta_UL_r(ii,nr) = 0.5*eta_UL(ii)/Nr;
        H_UL_r{ii,nr} = sqrt(eta_UL_r(ii,nr)/2)*(randn(N_UL(ii),1)+1i*randn(N_UL(ii),1));
    end
end
radar_comm_UL.UL2rchannelgains = eta_UL_r;
radar_comm_UL.UL2rchannles = H_UL_r;
%% radar_comm doppler 
kk = 0:K-1;
Q_Ir  = zeros(K,Nr,I);
for ii = 1:I
    for nr = 1:Nr
        f_i_nr = 0.15*(-1*rand(1)+0.5)+0.25;
        Q_Ir(:,nr,ii) = exp(1j*2*pi*f_i_nr.*kk');
    end
end
radar_comm_UL.ULrdoppler = Q_Ir;
end
