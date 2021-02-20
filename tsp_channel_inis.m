function [fdcomm, radar, radar_comm] = tsp_channel_inis(fdcomm,radar,radar_comm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Read parameters
Mr = radar.TX;
Nr = radar.RX;
Mc = fdcomm.BSTX;
Nc = fdcomm.BSRX;
K = radar.codelength; % Number of frames
N = radar.num_range_cell; % Number of symbols in a frame
%% radar Rx channels
H_rtr = zeros(Mr,K,Nr);
h_rtr = zeros(K*Mr,Nr);
% radar Tx-target-radar Rx
for nr = 1:Nr
    alpha_r_nr = sqrt(1/2)*(randn(Mr,1)+1i*randn(Mr,1));
    H_rt_nr = zeros(Mr,K);
    for k = 1 : K
        % the model of f_tnr is the same as ...
        % "An Information Theoretic Approach to Robust Constrained Code Design for MIMO Radars"
        f_mr_nr =  -1*rand(Mr,1)+0.5;  %[-0.5,0.5]
        f_mr_t_nr = 0.15*f_mr_nr+0.25; % Normalized Doppler frequency
        h_rt_nr_k = alpha_r_nr.*exp(1i*2*pi*(k-1).*f_mr_t_nr); % Doppler Domain Steering vector 
        H_rt_nr(:,k) = h_rt_nr_k;
    end
    H_rtr(:,:,nr) = H_rt_nr;
    h_rtr(:,nr) = reshape(H_rt_nr,[],1);
end
radar.channels.rtr = H_rtr;
radar.channels.tr = H_rtr;
radar.channels.Bmr = 0;
radar.channels.ur = 0;
if fdcomm.DL_num>0 % DL is enabled
    H_Bmr = zeros(Mc,K,Nr);
    for nr = 1:Nr
        alpha_Bm_nr = randn(Mc,1);
        H_Bm_nr = zeros(Mc,K);
        for k = 1 : K
            f_Bm_nr = 0.15*(-1*rand(1)+0.5)+0.25; % Normalized Doppler frequency
            h_Bm_nr_k = alpha_Bm_nr*exp(1i*2*pi*(k-1)*f_Bm_nr);
            H_Bm_nr(:,k) = h_Bm_nr_k;
        end
        H_Bmr(:,:,nr) = H_Bm_nr;
    end
    radar.channels.Bmr = H_Bmr;
    if radar_comm.radar_comm.isCollaborate
        H_Btr = zeros(Mc,K,Nr);
        H_tr = zeros(Mc+Mr,K,Nr);
        theta_Bt = fdcomm.theta_Bt;
        mc = 0:Mc-1;
        a_theta_Bt = exp(1i*pi*sin(theta_Bt).*mc); %steering vector
        for nr = 1:Nr
            alpha_Bt_nr = randn(1,1);
            H_Bt_nr = zeros(Mc,K);
            H_t_nr = zeros(Mc+Mr,K);
            for k = 1 : K
                f_Bt_nr = 0.15*(-1*rand(1)+0.5)+0.25; % Normalized Doppler frequency
                h_Bt_nr_k = alpha_Bt_nr*a_theta_Bt'*exp(1i*2*pi*(k-1)*f_Bt_nr);
                H_Bt_nr(:,k) = h_Bt_nr_k;
                h_rt_nr_k = H_rtr(:,k,nr);
                h_t_nr_k = [h_rt_nr_k;h_Bt_nr_k];
                H_t_nr(:,k) = h_t_nr_k;
            end
            H_Btr(:,:,nr) = H_Bt_nr;
            H_tr(:,:,nr) = h_t_nr_k;
        end
        radar.channels.tr = H_tr;
        radar.channels.Btr = H_Btr;
    else % no collaboration 
        radar.channels.tr = H_rtr;
        radar.channels.Btr = 0;
    end
    % Generating DL symbols
    J = fdcomm.DL_num;
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
    fdcomm.DLsymbols = D_DL;
    % Generating Channels BS-DL UEs
    % BS - DL_UEs
    H_DL = cell(J,1); % J DL Channels
    eta_DL = fdcomm.pathloss*ones(J,1); %
    for jj = 1:J
        %eta_DL(jj) = db2pow(SNR.BS_DL(jj))*sigma_0;
        H_DL{jj,1} = sqrt(eta_DL(jj)/2/N_DL(jj))*(randn(N_DL(jj),Mc)+1i*randn(N_DL(jj),Mc)); 
    end
    fdcomm.DLchannelgains = eta_DL;
    fdcomm.DLchannels = H_DL;
    % Generating radar-DL channels
    %mu_DL = 0.2*ones(Mr,J);
    %kappa_DL = 0.3;
    mu_DL = radar.Rician_direct;
    kappa_DL = radar.K_factor;
    H_r_DL = cell(J,1);
    eta_r_DL = zeros(Mr,J);
    for jj = 1 : J
        H_r_j = zeros(N_DL(jj),Mr);
        for mr = 1:Mr
            eta_r_DL(mr,jj) = db2pow(SNR.r_DL(mr,jj))*sigma_0; 
            H_r_j(:,mr) = sqrt(eta_r_DL(mr,jj)/(kappa_DL+1)/2)*(randn(N_DL(jj),1)+1i*randn(N_DL(jj),1))+sqrt(kappa_DL/(kappa_DL+1)*mu_DL(mr,jj)*ones(N_DL(jj),1));
        end
        H_r_DL{jj,1} = H_r_j;
    end
    radar_comm.radar2DLchannnels = H_r_DL;
    radar_comm.radar2DLchannelgains = eta_r_DL;
end
if fdcomm.UL_num>0 % UL Enabled
    I = fdcomm.UL_num;
    H_UL_r = cell(I,Nr);
    for nr = 1:Nr
        for ii = 1:I
            Nu = fdcomm.UL_UE_Ant;
            alpha_i_nr = sqrt(1/2)*(randn(Nu,1)+1i*randn(Nu,1));
            H_i_nr = zeros(Nu,K);
            for k = 1 : K
                f_i_nr = 0.15*(-1*rand(1)+0.5)+0.25; % Normalized Doppler frequency
                h_i_nr_k = alpha_i_nr*exp(1i*2*pi*(k-1)*f_i_nr);
                H_i_nr(:,k) = h_i_nr_k;
            end
            H_UL_r{ii,nr} = H_i_nr;
        end
    end
    radar.channels.ur = H_UL_r; %H_UL_r is a cell
    %Generating UL symbols
    D_UL = cell(I,1);
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
    fdcomm.ULsymbols = D_UL;
    % Generating UL channels: UL UEs - BS
    eta_UL = fdcomm.pathloss*ones(I,1);
    H_UL = cell(I,1);
    for ii = 1:I
        %eta_UL(ii) = db2pow(SNR.UL_BS(ii)) *sigma_0;
        H_UL{ii,1} = sqrt(eta_UL(ii)/2)*(randn(Nc,N_UL(ii))+1i*randn(Nc,N_UL(ii)));  
    end
    fdcomm.ULchannelgains = eta_UL;
    fdcomm.ULchannels = H_UL;
    % radar TX - BS 
    %mu_r_BS = 0.5*ones(Mr,1);
    %kappa_BS = 0.5;
    mu_r_BS = radar.Rician_direct;
    kappa_BS = radar.K_factor;
    H_r_BS = zeros(Nc,Mr);
    eta_rB = zeros(Mr,1);
    for mr = 1 : Mr 
        eta_rB(mr) = db2pow(SNR.r_B(mr))*sigma_0;
        H_r_BS(:,mr) = sqrt(eta_rB(mr)/(kappa_BS+1)/2).*(randn(Nc,1)+1i*randn(Nc,1))+sqrt(kappa_BS/(kappa_BS+1))*mu_r_BS(mr)*ones(Nc,1); 
    end
    radar_comm.radar2BSchannels = H_r_BS;
    radar_comm.r2Bchannelgains = eta_rB;
end
%% FD Comm
if fdcomm.UL_num>0 && fdcomm.UL_num>0
    % BS - BS
    K_B = fdcomm.K_factor ;
    eta_BB = 1/(1+K_B);
    H_BB = sqrt(eta_BB/2)*(randn(Nc,Mc)+1i*randn(Nc,Mc))+sqrt(k_B/(K_B+1))*ones(Nc,Mc);
    fdcomm.BBchannel = H_BB;
    fdcomm.BBchannelgain = eta_BB;
    % UL-DL
    eta_UL_DL = fdcomm.pathloss*ones(I,J);
    H_UL_DL = cell(I,J);
    for ii = 1:I
        for jj = 1:J
            %eta_UL_DL(ii,jj) = db2pow(SNR.UL_DL(ii,jj))*sigma_0;
            H_UL_DL{ii,jj} = sqrt(eta_UL_DL(ii,jj)/2)*(randn(N_DL(jj),N_UL(ii))+1i*randn(N_DL(jj),N_UL(ii)));
        end
    end
    fdcomm.ULDLchannels = H_UL_DL;
    fdcomm.ULDLchannelgains = eta_UL_DL;
end
end

