function [fdcomm, radar, radar_comm] = tsp_parameters_UL_varying(SNR,radar,fdcomm)
%%%%-----------------------------------
%------- SNR.SNR_rtr: Mr*Nr real matrix
%------- SNR.SNR_DL:
Mr = radar.TX;
Nr = radar.RX;
Mc = fdcomm.BSTX;
Nc = fdcomm.BSRX;
K = radar.codelength; % The length of the radar code; or the number of PRIs
N = 8;% Number of range cells
n = 4; % CUT index
radar.codelength = K;
radar.PRI_num = N;
radar.CUT_Idx = n;
sigma_0 = radar.noisepower;
%% radar channels
% % radar TX - target - radar RX
% eta_rtr = zeros(Mr,Nr);
% H_rtr = zeros(Mr,Nr);
% for nr = 1:Nr
%     for mr = 1 : Mr
%         eta_rtr(mr,nr) = db2pow(SNR.rtr(mr,nr))*sigma_0;
%         H_rtr(mr,nr) = sqrt(eta_rtr(mr,nr)/2)*(randn(1,1)+1i*randn(1,1));
%     end
% end
% radar.channel = H_rtr;
% radar.channelgain = eta_rtr;

%% radar Doppler               
Qr = zeros(Mr,K,Nr);
for nr = 1:Nr
    for k = 1 : K
        % the model of f_tnr is the same as ...
        % "An Information Theoretic Approach to Robust Constrained Code Design for MIMO Radars"
        f_mr_nr =  -1*rand(Mr,1)+0.5;  %[-0.5,0.5]
        f_mr_t_nr = 0.15*f_mr_nr+0.25; % Normalized Doppler frequency
        Qr(:,k,nr) = exp(1i*2*pi*(k-1).*f_mr_t_nr); % Doppler Domain Steering vector 
    end
end
radar.doppler = Qr; %radar temporal steering vector
%% clutter 
CNR = SNR.CNR;
R_C = zeros(K,K,Nr); 
rho = 0.5;
% Method 1 C is a gaussian rv (0,sigma_c);
for nr = 1:Nr
    sigma_c_nr = db2pow(CNR(nr))*sigma_0;
    R_C_nr = zeros(K,K);
    for ii = 1 : K 
        for jj = 1:K
            R_C_nr(ii,jj) = sigma_c_nr*rho^(abs(ii-jj));
        end
    end
    R_C(:,:,nr) = R_C_nr;
end
radar.cluttercov = R_C;
%% Radar PAR constraints
radar.gamma_r = 2*ones(Mr,1); % PAR level
%% FD comm parameters
I = fdcomm.UL_num;
J= fdcomm.DL_num;
N_UL = 2*ones(I,1); % Number of the UL UE antennas
N_DL = 2*ones(J,1); % Number of the DL UE antennas
% number of stream for each DL UE d_DL(i) <= N_DL(i) 
d_DL = N_DL;
d_UL = N_UL;
fdcomm.UL_UE_Ant = N_UL;
fdcomm.DL_UE_Ant = N_DL;
fdcomm.DLstream_num = d_DL;
fdcomm.ULstream_num = d_UL;
%% FD Comm channels

% BS - DL_UEs
H_DL = cell(J,1); % J DL Channels
eta_DL = ones(J,1); %
for jj = 1:J
    eta_DL(jj) = db2pow(SNR.BS_DL(jj))*sigma_0;
    H_DL{jj,1} = sqrt(eta_DL(jj)/2/N_DL(jj))*(randn(N_DL(jj),Mc)+1i*randn(N_DL(jj),Mc)); 
end
fdcomm.DLchannelgains = eta_DL;
fdcomm.DLchannels = H_DL;
% BS - BS
eta_BB = 0.8*sum(eta_DL);
H_BB = sqrt(eta_BB/2)*(randn(Nc,Mc)+1i*randn(Nc,Mc));
fdcomm.BBchannel = H_BB;
fdcomm.BBchannelgain = eta_BB;
%% DL -> radar
% BS - target - radar RX and BS - multi-path - radar RX
eta_Btr = zeros(Nr,1);
H_Btr = zeros(Mc,Nr);
eta_Bmr = zeros(Nr,1);
H_Bmr = zeros(Mc,Nr);
for nr = 1:Nr
    eta_Btr(nr) = 0.5*sum(eta_DL)/Nr;
    eta_Bmr(nr) = 0.7*sum(eta_DL)/Nr;
    H_Bmr(:,nr) = sqrt(eta_Bmr(nr)/2).*(randn(Mc,1)+1i*randn(Mc,1));
    H_Btr(:,nr) = sqrt(eta_Btr(nr)/2).*(randn(Mc,1)+1i*randn(Mc,1));
end
radar_comm.Bmrchannels = H_Bmr;
radar_comm.Bmrchannelgains = eta_Bmr;
radar_comm.Btrchannels = H_Btr;
radar_comm.Btrchannelgains = eta_Btr;
radar_comm.Jr = [eye(Mr);zeros(Mc,Mr)];
radar_comm.JB = [zeros(Mr,Mc);eye(Mc)];
%% Updating radar channel matrices
M = Mc + Mr; 
eta_rtr = zeros(Mr,Nr);
Sigma = zeros(M,M,Nr);
eta_Btr = radar_comm.Btrchannelgains;
for nr = 1:Nr
    for mr = 1 : Mr
        eta_rtr(mr,nr) = db2pow(SNR.rtr(mr,nr))*radar.noisepower;
        Sigma(:,:,nr) = blkdiag(eta_rtr(:,nr).*eye(Mr),eta_Btr(nr)*eye(Mc));
    end
end
radar.channelgain = eta_rtr;
radar.Sigma = Sigma;

%% radar ->comm channels 
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
%% FD comm symbols 
D_DL = cell(J,1);
D_UL = cell(I,1);
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
fdcomm.DLsymbols = D_DL;
fdcomm.ULsymbols = D_UL;
%% radar_comm doppler 
kk = 0:K-1;
f_Bm_Nr = zeros(Nr,1);
Q_Btr = zeros(K,Nr);
Q_Bmr = zeros(K,Nr);
Q_Ir  = zeros(K,Nr,I);
f_Bt_Nr = zeros(Nr,1);
for nr = 1:Nr
    f_Bt_nr = 0.15*(-1*rand(1)+0.5)+0.25; % Normalized Doppler frequency
    f_Bm_nr = 0.15*(-1*rand(1)+0.5)+0.25;
    Q_Btr(:,nr) = exp(1j*2*pi*f_Bt_nr.*kk'); % Doppler Domain Steering vector 
    Q_Bmr(:,nr) = exp(1j*2*pi*f_Bm_nr.*kk');
    f_Bt_Nr(nr) = f_Bt_nr;
    f_Bm_Nr(nr) = f_Bm_nr;
end
for ii = 1:I
    for nr = 1:Nr
        f_i_nr = 0.15*(-1*rand(1)+0.5)+0.25;
        Q_Ir(:,nr,ii) = exp(1j*2*pi*f_i_nr.*kk');
    end
end
radar_comm.f_Bt_Nr = f_Bt_Nr;
radar_comm.f_Bm_Nr = f_Bm_Nr;
radar_comm.Btrdoppler = Q_Btr;
radar_comm.Bmrdoppler = Q_Bmr;
radar_comm.ULrdoppler = Q_Ir;
end
