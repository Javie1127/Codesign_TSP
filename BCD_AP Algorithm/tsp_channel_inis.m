function [fdcomm, radar, radar_comm] = tsp_channel_inis(fdcomm,radar,radar_comm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Read parameters
Mr = radar.Tx;
Nr = radar.Rx;
Mc = fdcomm.BSTx;
Nc = fdcomm.BSRx;
K = radar.codelength; % Number of frames
N = radar.num_range_cell; % Number of symbols in a frame
M = radar.total_Tx;

%% Initiating channels
H_r_BS = zeros(Nc,Mr);
%% radar Rx channels
H_rtr = zeros(Mr,K,Nr);
h_rtr = zeros(K*Mr,Nr);
Qr = zeros(Mr,K,Nr); % store Doppler info
Alpha_r = zeros(Mr,Nr);
% radar Tx-target-radar Rx
%% Radar & Clutter
%sigma_c = radar.clutter_power;
rho = 0.5;
clutter_power = 10^(radar.CNR/10);
Sigma_c_Nr = cell(K,K,Nr);
% sigma_c = 10^((radar.CNR-radar.SNR)/10);
% Rho_c_r = zeros(Mr,Nr);
% for nr = 1:Nr
%    Rho_c_r(:,nr)=sqrt(sigma_c/2)*(randn(Mr,1)+1i*randn(Mr,1));
% end
% radar.Clutter_channel = Rho_c_r;
% radar.Clutter_channel_cov = eye(Mr)*sigma_c;
for nr = 1:Nr
    alpha_r_nr = sqrt(1/2)*(randn(Mr,1)+1i*randn(Mr,1));
    Alpha_r(:,nr) = alpha_r_nr;
    H_rt_nr = zeros(Mr,K);
    f_mr_nr =  -1*rand(Mr,1)+0.5;  %[-0.5,0.5]
    f_mr_t_nr = 0.15*f_mr_nr+0.25; % Normalized Doppler frequency
    for k = 1 : K
        Qr(:,k,nr) = exp(1i*2*pi*(k-1).*f_mr_t_nr); % Doppler Domain Steering vector 
        h_rt_nr_k = alpha_r_nr.*Qr(:,k,nr); % Doppler Domain Steering vector 
        H_rt_nr(:,k) = h_rt_nr_k;
    end
    H_rtr(:,:,nr) = H_rt_nr;
    h_rtr(:,nr) = reshape(H_rt_nr,[],1);
end
Sigma_rt_Nr = cell(K,K,Nr);
for nr = 1:Nr
    for m = 1:K
        qr_m_nr = Qr(:,m,nr);
        for l = 1:K
            qr_l_nr = Qr(:,l,nr);
            Sigma_rt_nr_m_l = diag(qr_m_nr.*conj(qr_l_nr));
            Sigma_rt_Nr{m,l,nr} = Sigma_rt_nr_m_l;
            Sigma_c_Nr{m,l,nr} = clutter_power*rho^(abs(m-l))*eye(Mr);
        end
    end
end
radar.Clutter_channel_cov = Sigma_c_Nr;
radar.Sigma_rt_Nr= Sigma_rt_Nr;
Sigma_t_Nr = Sigma_rt_Nr;
radar.channlegains = Alpha_r;
radar.doppler = Qr; 
radar.channel = H_rtr;
%radar.total_channel = H_rtr;
H_tr = H_rtr;
radar_comm.Bmr = 0;
radar_comm.ur = 0;
%----- DL-Radar----------------------------------
if fdcomm.DL_num>0 % DL is enabled
    %% Generating DL symbols
    J = fdcomm.DL_num;
    D_DL = cell(J,1);
    d_DL = fdcomm.DLstream_num;
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
    %% BS-radar Rx channels
    H_Bmr = zeros(Mc,K,Nr);
    Q_Bmr = zeros(K,Nr);
    Sigma_Bm_Nr = cell(K,K,Nr);
    alpha_Bm = zeros(Nr,Mc);
    for nr = 1:Nr
        alpha_Bm_nr = sqrt(1/2)*(randn(Mc,1)+1i*randn(Mc,1));
        alpha_Bm(nr,:) = alpha_Bm_nr.';
        H_Bm_nr = zeros(Mc,K);
        f_Bm_nr = 0.15*(-1*rand(1)+0.5)+0.25; % Normalized Doppler frequency
        for k = 1 : K
            Q_Bmr(k,nr) = exp(1i*2*pi*(k-1)*f_Bm_nr);
            h_Bm_nr_k = alpha_Bm_nr*Q_Bmr(k,nr);
            H_Bm_nr(:,k) = h_Bm_nr_k;
        end
        H_Bmr(:,:,nr) = H_Bm_nr;
    end
    for nr = 1:Nr
        for m = 1:K
            q_Bm_nr_m = Q_Bmr(m,nr);
            for l = 1:K
                q_Bm_nr_l = Q_Bmr(l,nr);
                Sigma_Bm_nr_m_l = q_Bm_nr_m.*conj(q_Bm_nr_l)*eye(Mc);
                Sigma_Bm_Nr{m,l,nr} = Sigma_Bm_nr_m_l;
            end
        end
    end
    radar_comm.Sigma_Bm_Nr = Sigma_Bm_Nr;
    radar_comm.Bmr = H_Bmr;
    radar_comm.doppler_Bmr = Q_Bmr;
    radar_comm.Bmr_channel_matrix = alpha_Bm;
    if radar_comm.isCollaborate
        H_Btr = zeros(Mc,K,Nr);
        H_tr = zeros(Mc+Mr,K,Nr);
        Q_Btr = zeros(K,Nr);
        theta_Bt = fdcomm.theta_BT;
        mc = 0:Mc-1;
        a_theta_Bt = exp(1i*2*pi*sin(theta_Bt).*mc); %steering vector
        for nr = 1:Nr
            alpha_Bt_nr = sqrt(1/2)*(randn(1,1)+1i*randn(1,1));
            H_Bt_nr = zeros(Mc,K);
            H_t_nr = zeros(Mc+Mr,K);
            f_Bt_nr = 0.15*(-1*rand(1)+0.5)+0.25; % Normalized Doppler frequency
            for k = 1 : K
                Q_Btr(k,nr) = exp(1i*2*pi*(k-1)*f_Bt_nr);
                h_Bt_nr_k = alpha_Bt_nr*a_theta_Bt'*Q_Btr(k,nr);
                H_Bt_nr(:,k) = h_Bt_nr_k;
                h_rt_nr_k = H_rtr(:,k,nr);
                h_t_nr_k = [h_rt_nr_k;h_Bt_nr_k];
                H_t_nr(:,k) = h_t_nr_k;
            end
            H_Btr(:,:,nr) = H_Bt_nr;
            H_tr(:,:,nr) = H_t_nr;
        end
        Sigma_Bt_Nr = cell(K,K,Nr);
        Sigma_t_Nr = cell(K,K,Nr);
        for nr = 1:Nr
            for m = 1:K
                q_Bt_nr_m = Q_Btr(m,nr);
                for l = 1:K
                    q_Bt_nr_l = Q_Btr(l,nr);
                    Sigma_Bt_nr_m_l = q_Bt_nr_m.*conj(q_Bt_nr_l)*eye(Mc);
                    Sigma_Bt_Nr{m,l,nr} = Sigma_Bt_nr_m_l;
                    Sigma_t_Nr{m,l,nr} = radar_comm.Jr *Sigma_rt_Nr{m,l,nr}*(radar_comm.Jr).'...
                        +radar_comm.JB*Sigma_Bt_Nr{m,l,nr}*radar_comm.JB.';
                end
            end
        end
        radar_comm.Sigma_Bt_Nr = Sigma_Bt_Nr;
        radar_comm.doppler_Btr = Q_Btr;
        radar_comm.Btr = H_Btr;
    end
    
    % Generating Channels BS-DL UEs
    % BS - DL_UEs
    H_DL = cell(J,1); % J DL Channels
    N_DL = fdcomm.DL_UE_Ant;
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
    mu_DL = radar.Rician_direct*ones(Mr,J);
    kappa_DL = radar.K_factor;
    H_r_DL = cell(J,1);
    eta_r_DL = 1/(1+kappa_DL)*ones(Mr,J);
    for jj = 1 : J
        H_r_j = zeros(N_DL(jj),Mr);
        for mr = 1:Mr
            %eta_r_DL(mr,jj) = db2pow(SNR.r_DL(mr,jj))*sigma_0; 
            H_r_j(:,mr) = sqrt(eta_r_DL(mr,jj)/(kappa_DL+1)/2)*(randn(N_DL(jj),1)+1i*randn(N_DL(jj),1))+sqrt(kappa_DL/(kappa_DL+1)*mu_DL(mr,jj)*ones(N_DL(jj),1));
        end
        H_r_DL{jj,1} = H_r_j;
    end
    radar_comm.radar2DLchannnels = H_r_DL;
    radar_comm.radar2DLchannelgains = eta_r_DL;
else
    
end
if fdcomm.UL_num>0 % UL Enabled
    I = fdcomm.UL_num;
    H_UL_r = cell(I,Nr);
    Q_UL_r = zeros(I,K,Nr);
    for nr = 1:Nr
        Nu = fdcomm.UL_UE_Ant;
        for ii = 1:I
            alpha_i_nr = sqrt(1/2)*(randn(Nu(ii),1)+1i*randn(Nu(ii),1));
            H_i_nr = zeros(Nu(ii),K);
            for k = 1 : K
                f_i_nr = 0.15*(-1*rand(1)+0.5)+0.25; % Normalized Doppler frequency
                Q_UL_r(ii,k,nr) = exp(1i*2*pi*(k-1)*f_i_nr);
                h_i_nr_k = alpha_i_nr*exp(1i*2*pi*(k-1)*f_i_nr);
                H_i_nr(:,k) = h_i_nr_k;
            end
            H_UL_r{ii,nr} = H_i_nr;
        end
    end
    Sigma_U_Nr = cell(K,K,I,Nr);
    for nr = 1:Nr
        for ii = 1:I
            for m = 1:K
                q_i_nr_m = Q_UL_r(ii,m,Nr);
                for l = 1:K
                    q_i_nr_l = Q_UL_r(ii,l,nr);
                    Sigma_i_mr_m_l = diag(q_i_nr_m.*conj(q_i_nr_l))*eye(Nu(ii));
                    Sigma_U_Nr{m,l,ii,nr} = Sigma_i_mr_m_l;
                end
            end
        end
    end
    radar_comm.Sigma_U_Nr = Sigma_U_Nr;
    radar_comm.doppler_ur = Q_UL_r;
    radar_comm.ur = H_UL_r; %H_UL_r is a cell
    %Generating UL symbols
    D_UL = cell(I,1);
    d_UL = fdcomm.ULstream_num;
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
    N_UL = fdcomm.UL_UE_Ant;
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
    mu_r_BS = radar.Rician_direct*ones(Mr,1);
    kappa_BS = radar.K_factor;
    %H_r_BS = zeros(Nc,Mr);
    eta_rB = 1/(1+kappa_BS)*ones(Mr,1);
    for mr = 1 : Mr 
        %eta_rB(mr) = db2pow(SNR.r_B(mr))*sigma_0;
        H_r_BS(:,mr) = sqrt(eta_rB(mr)/(kappa_BS+1)/2).*(randn(Nc,1)+1i*randn(Nc,1))+sqrt(kappa_BS/(kappa_BS+1))*mu_r_BS(mr)*ones(Nc,1); 
    end
    
    radar_comm.r2Bchannelgains = eta_rB;
end
%% FD Comm
if fdcomm.UL_num>0 && fdcomm.DL_num>0
    % BS - BS
    K_B = fdcomm.K_factor ;
    eta_BB = 1/(1+K_B);
    H_BB = sqrt(eta_BB/2)*(randn(Nc,Mc)+1i*randn(Nc,Mc))+sqrt(K_B/(K_B+1))*ones(Nc,Mc);
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

radar_comm.radar2BSchannels = H_r_BS;


radar.total_channel = H_tr;
radar.Sigma_t_Nr = Sigma_t_Nr;
Sigma_h_tr = zeros(K*M,K*M,Nr);
JH = radar_comm.JH;
for nr = 1:Nr
    for m = 1:K
        Jhm = JH{m,1};
        for l = 1:K
            Jhl = JH{l,1};
            Sigma_h_tr(:,:,nr) = Jhm*Sigma_t_Nr{m,l,nr}*Jhl.'+Sigma_h_tr(:,:,nr);
        end
    end
end
radar.Sigma_h_tr = Sigma_h_tr;
end

