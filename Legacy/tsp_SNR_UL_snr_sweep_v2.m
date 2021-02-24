% Pd vs \vu for different SNR_Btr values

%% Array Parameters
Mr = 4; % Number of radar TX antennas
Nr = 4; % Number of radar RX antennas

%% FD Comm parameters
Mc = 4; % Number of BS TX antennas
Nc = 4;% Number of BS RX antennas
I = 2; % Number of UL UEs
J = 2; % Number of DL UEs
K = 8;
radar.codelength = K;
%% Set the SNRs in dB
% SNR.rtr = randi([-5,5],Mr,Nr);
SNR.Btr = 2*ones(Nr,1);
radar.TX = Mr;
radar.RX = Nr;
sigma_0 = 0.01;
radar.noisepower = sigma_0;
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
% SNR.UL_r = 2*ones(I,Nr);
SNR.UL_DL = 1*ones(I,J);
SNR.r_DL = 1*ones(Mr,J);
radar.Pr = 1*ones(Mr,1);
SNR.CNR = ones(Nr,1);
radar.ell_max = 20; % algorithm 5
% %% priority unequal weight
% fdcomm.alpha_UL = 0.1*ones(I,1);
% fdcomm.alpha_DL = 0.05*ones(J,1);
% radar.alpha_r = (1-0.1*I-0.05*J)/Nr*ones(Nr,1);
%% priority weight equal
fdcomm.alpha_UL = 1/(Nr+I+J)*ones(I,1);
fdcomm.alpha_DL = 1/(Nr+I+J)*ones(J,1);
radar.alpha_r = 1/(Nr+I+J)*ones(Nr,1);
%% simulation parameters
%eta_rtr = radar_para.channelgain;
ell_max = radar.ell_max;
%% Parameters to save
UL_SNR_sweep = -20:5:20;
% radar_save = cell(length(UL_radar_SNR_sweep),1);
% fdcomm_save = cell(length(UL_radar_SNR_sweep),1);
% cov_save = cell(length(UL_radar_SNR_sweep),1);
%% Loop through CNRs 
M = 15;
[fdcomm_para, radar_para, radar_comm_para] = tsp_parameters_UL_varying(SNR,radar,fdcomm);
I_total = zeros(length(UL_SNR_sweep),1);
I_total_UL = zeros(length(UL_SNR_sweep),1);
I_total_DL = zeros(length(UL_SNR_sweep),1);
I_total_FD = zeros(length(UL_SNR_sweep),1);
I_total_radar = zeros(length(UL_SNR_sweep),1);
I_total_UL_avg = zeros(length(UL_SNR_sweep),1);
I_total_DL_avg = zeros(length(UL_SNR_sweep),1);
I_total_FD_avg = zeros(length(UL_SNR_sweep),1);
I_total_radar_avg = zeros(length(UL_SNR_sweep),1);
for ii = 1:length(UL_SNR_sweep)
    %% random initialization 
    snr_UL = UL_SNR_sweep(ii)*ones(I,1);
    radar_para_ii = radar_para;
    fdcomm_para_ii = fdcomm_para;
    radar_comm_para_ii = radar_comm_para;
    [fdcomm_ini_ii,radar_ini_ii,radar_comm_ini_ii] =...
        tsp_ini_random_UL_sweep(radar_para_ii,fdcomm_para_ii,radar_comm_para_ii, snr_UL);
    I_total_ii = zeros(M,1);
    I_total_radar_ii = zeros(M,1);
    I_total_FD_ii = zeros(M,1);
    I_total_DL_ii = zeros(M,1);
    I_total_UL_ii = zeros(M,1);
    parfor m = 1:M
        [fdcomm_ii_m,radar_ii_m,cov_ii_m] = tsp_ini_random_precoders(radar_ini_ii,fdcomm_ini_ii,radar_comm_ini_ii);
        [fdcomm_op_ii_m, radar_op_ii_m, cov_op_ii_m,~] =...
        tsp_altermating_projection_v1(fdcomm_ii_m, radar_ii_m,radar_comm_ini_ii,cov_ii_m);
        I_total_ii(m) = fdcomm_op_ii_m.I_total_op(end);
        I_total_DL_ii(m) = fdcomm_op_ii_m.I_DL_op(end);
        I_total_UL_ii(m) = fdcomm_op_ii_m.I_UL_op(end);
        I_total_FD_ii(m) =I_total_DL_ii(m) + I_total_UL_ii(m);
        I_total_radar_ii(m) = fdcomm_op_ii_m.I_radar_op(end);
    end
    
    I_total(ii) = mean(I_total_ii);
    I_total_radar(ii) = mean(I_total_radar_ii);
    I_total_UL(ii) = mean(I_total_UL_ii);
    I_total_DL(ii) = mean(I_total_DL_ii);
    I_total_FD(ii) = mean(I_total_FD_ii);
end
save('tsp_UL_SNR_new.mat')
 
