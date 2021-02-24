% Pd vs \vu for different SNR_Btr values

%% Array Parameters
Mr = 6; % Number of radar TX antennas
Nr = 4; % Number of radar RX antennas

%% FD Comm parameters
Mc = 4; % Number of BS TX antennas
Nc = 4;% Number of BS RX antennas
I = 4; % Number of UL UEs
J = 2; % Number of DL UEs

%% Set the SNRs in dB
% SNR.rtr = randi([-5,5],Mr,Nr);
SNR_Btr = -20:5:10;
radar.TX = Mr;
radar.RX = Nr;
radar.noisepower = 0.01;
SNR.Btr = ones(Nr,1);
fdcomm.BSTX = Mc;
fdcomm.BSRX = Nc;
fdcomm.UL_num = I;
fdcomm.DL_num = J;
fdcomm.ULpower = ones(I,1);
fdcomm.DLpower = J;
SNR.rtr = 2*ones(Mr,Nr);
SNR.Bmr = 1*ones(Nr,1);
SNR.BB = 1;
SNR.BS_DL = 5*ones(J,1);
SNR.UL_BS = 2*ones(I,1);
SNR.r_B = 2*ones(Mr,1);
SNR.UL_r = 2*ones(I,Nr);
SNR.UL_DL = 1*ones(I,J);
SNR.r_DL = 1*ones(Mr,J);
radar.Pr = 1*ones(Mr,1);
% Clutter 
SNR.CNR = ones(Nr,1);
radar.ell_max = 7; % algorithm 5

% %% priority unequal weight
% fdcomm.alpha_UL = 0.1*ones(I,1);
% fdcomm.alpha_DL = 0.05*ones(J,1);
% radar.alpha_r = (1-0.1*I-0.05*J)/Nr*ones(Nr,1);
%% priority weight equal
fdcomm.alpha_UL = 1/(Nr+I+J)*ones(I,1);
fdcomm.alpha_DL = 1/(Nr+I+J)*ones(J,1);
radar.alpha_r = 1/(Nr+I+J)*ones(Nr,1);

%% simulation parameters
[fdcomm_para, radar_para, radar_comm_para] = tsp_parameters(SNR,radar,fdcomm);
eta_rtr = radar_para.channelgain;
M = Mc + Mr;
ell_max = radar.ell_max;
for snr_Btr = SNR_Btr
    %% Update the channel gain 
    eta_Btr = zeros(Nr,1);
    SNR.Btr = snr_Btr*ones(Nr,1); 
    Sigma = zeros(M,M,Nr);
    for nr = 1:Nr
        % SNR.Btr : Nr * 1
        eta_Btr(nr) = db2pow(SNR.Btr(nr))*radar.noisepower;
        Sigma(:,:,nr) = blkdiag(eta_rtr(nr)*eye(Mr),eta_Btr(nr)*eye(Mc));
    end
    radar_comm_para.Btrchannelgains = eta_Btr;
    radar_para.Sigma = Sigma;
    %% initialization 
    [fdcomm_ini, radar_ini, cov_temp] = tsp_ini_normal(radar_para,fdcomm_para,radar_comm_para);
    %% WMMSE alternating algorithm
    radar_op = radar_ini;
    fdcomm_op = fdcomm_ini;
    fdcomm = fdcomm_ini;
    radar = radar_ini; 
    I_total = zeros(ell_max,1);
    I_total_op = zeros(ell_max,1);
    I_max = sum(fdcomm_op.MI_UL(:))+sum(fdcomm_op.MI_DL(:))+ sum(radar_op.MI_radar);
    ell = 1;
    while ell <= ell_max
        [fdcomm, radar_ast] =...
            tsp_WMMSE_algorithm(fdcomm, radar, radar_comm_para, cov_temp);
        % Calculate A^star
        radar = tsp_Nearest_PAR(radar_ast);
        % Update covariance matrices
        cov_temp = tsp_covmat(fdcomm,radar,radar_comm_para);
        % Update linear receivers
        fdcomm = tsp_Comm_MMSE(fdcomm, radar, cov_temp);
        radar = tsp_radar_MMSE(radar, cov_temp);
        I_total(ell) = sum(fdcomm.MI_UL(:))+sum(fdcomm.MI_DL(:))+ sum(radar.MI_radar);
        if I_total(ell) > I_max
            I_total_op(ell) = I_total(ell);
            I_max = I_total(ell);
            radar_op = radar;
            fdcomm_op = fdcomm;
        else
            I_total_op(ell) = I_max;
        end
        ell = ell+1;
    end
    fdcomm_op.I_total = I_total;
    fdcomm_op.I_total_op = I_total_op;
    cov_op = tsp_covmat(fdcomm_op,radar_op,radar_comm_para);
    file_name = ['uniformweight_BTR_', num2str(snr_Btr),'dB_co.mat'];
    save(file_name);
end

 
