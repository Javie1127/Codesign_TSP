
%% Array Parameters
Mr = 4; % Number of radar TX antennas
Nr = 4; % Number of radar RX antennas

%% FD Comm parameters
Mc = 4; % Number of BS TX antennas
Nc = 4;% Number of BS RX antennas
I = 2; % Number of UL UEs
J = 2; % Number of DL UEs
radar.codelength = 8;
%% Set the SNRs in dB
% SNR.rtr = randi([-5,5],Mr,Nr);
radar.TX = Mr;
radar.RX = Nr;
radar.noisepower = 0.01;
SNR.Btr = 2*ones(Nr,1);
fdcomm.BSTX = Mc;
fdcomm.BSRX = Nc;
fdcomm.UL_num = I;
fdcomm.DL_num = J;
fdcomm.ULpower = ones(I,1);
fdcomm.DLpower = J;
SNR.rtr = 4*ones(Mr,Nr);
SNR.Bmr = 1*ones(Nr,1);
SNR.BB = 1;
SNR.BS_DL = 2*ones(J,1);
SNR.UL_BS = 2*ones(I,1);
SNR.r_B = 2*ones(Mr,1);
SNR.UL_r = 2*ones(I,Nr);
SNR.UL_DL = 1*ones(I,J);
SNR.r_DL = 1*ones(Mr,J);
radar.Pr = 1*ones(Mr,1);
% Clutter 
SNR.CNR = ones(Nr,1);
radar.ell_max = 5; % algorithm 5

% %% priority unequal weight
% fdcomm.alpha_UL = 0.1*ones(I,1);
% fdcomm.alpha_DL = 0.05*ones(J,1);
% radar.alpha_r = (1-0.1*I-0.05*J)/Nr*ones(Nr,1);
%% priority weight equal
fdcomm.alpha_UL = 1/(Nr+I+J)*ones(I,1);
fdcomm.alpha_DL = 1/(Nr+I+J)*ones(J,1);
radar.alpha_r = 1/(Nr+I+J)*ones(Nr,1);

dmode = 'radar';
SNR_rtr = -5:5:5;
radar_co_save = cell(length(SNR_rtr),1);
fdcomm_co_save = cell(length(SNR_rtr),1);
cov_co_save = cell(length(SNR_rtr),1);
radar_co_wc_save = cell(length(SNR_rtr),1);
fdcomm_co_wc_save = cell(length(SNR_rtr),1);
cov_co_wc_save = cell(length(SNR_rtr),1);
radar_uc_save = cell(length(SNR_rtr),1);
fdcomm_uc_save = cell(length(SNR_rtr),1);
cov_uc_save = cell(length(SNR_rtr),1);
radar_uc_wc_save = cell(length(SNR_rtr),1);
fdcomm_uc_wc_save = cell(length(SNR_rtr),1);
cov_uc_wc_save = cell(length(SNR_rtr),1);
radar_rc_save = cell(length(SNR_rtr),1);
fdcomm_rc_save = cell(length(SNR_rtr),1);
cov_rc_save =  cell(length(SNR_rtr),1);
radar_rc_wc_save = cell(length(SNR_rtr),1);
fdcomm_rc_wc_save = cell(length(SNR_rtr),1);
cov_rc_wc_save =  cell(length(SNR_rtr),1);
%% simulation parameters
[fdcomm_para_co, radar_para_co, radar_comm_para_co] = ...
    tsp_parameters_v1(SNR,radar,fdcomm);
parfor ii = 1:length(SNR_rtr)
    snr_rtr = SNR_rtr(ii)*ones(Mr,Nr);
    %% Update the channel gain
    radar_comm_para_co_ii = radar_comm_para_co;
    radar_para_co_ii = radar_para_co;
    fdcomm_para_co_ii = fdcomm_para_co;
    %% initialization 
    [fdcomm_co_ini_ii, radar_co_ini_ii, radar_comm_co_ini_ii,cov_co_ini_ii] =...
        tsp_ini_normal_v1(radar_para_co_ii,fdcomm_para_co_ii,radar_comm_para_co_ii,snr_rtr);
    %% co-design with co-op
    [fdcomm_co_op_ii, radar_co_op_ii, cov_co_op_ii,radar_co_ini_ii] =...
        tsp_altermating_projection_v1(fdcomm_co_ini_ii, radar_co_ini_ii, radar_comm_co_ini_ii,cov_co_ini_ii);
    fdcomm_co_save{ii} = fdcomm_co_op_ii;
    radar_co_save{ii} = radar_co_op_ii;
    cov_co_save{ii} = cov_co_op_ii;
    %% co-design without co-op
    [fdcomm_co_wc_ini_ii, radar_co_wc_ini_ii, radar_comm_co_wc_ini_ii,cov_co_wc_ini_ii] =...
        tsp_ini_normal_radar_snr_wc(radar_co_ini_ii,fdcomm_co_ini_ii,radar_comm_co_ini_ii,snr_rtr);
    [fdcomm_co_wc_op_ii, radar_co_wc_op_ii, cov_co_wc_op_ii,~] =...
        tsp_altermating_projection_wc(fdcomm_co_wc_ini_ii, radar_co_wc_ini_ii, radar_comm_co_wc_ini_ii,cov_co_wc_ini_ii);
    fdcomm_co_wc_save{ii} = fdcomm_co_wc_op_ii;
    radar_co_wc_save{ii} = radar_co_wc_op_ii;
    cov_co_wc_save{ii} = cov_co_wc_op_ii;
    %% uncoded with co-op
    [fdcomm_uc_op_ii, radar_uc_op_ii, cov_uc_op_ii,~] =...
        tsp_altermating_projection_wo_radar(fdcomm_co_ini_ii, radar_co_ini_ii, radar_comm_co_ini_ii,cov_co_ini_ii);
    fdcomm_uc_save{ii} = fdcomm_uc_op_ii;
    radar_uc_save{ii} = radar_uc_op_ii;
    cov_uc_save{ii} = cov_uc_op_ii;
    %% uncoded without co-op
%     [fdcomm_uc_wc_ini_ii, radar_uc_wc_ini_ii, radar_comm_uc_wc_ini_ii,cov_uc_wc_ini_ii] =...
%        tsp_ini_wc(radar_co_ini_ii,fdcomm_co_ini_ii,radar_comm_co_ini_ii,cov_co_ini_ii);
     [fdcomm_uc_wc_op_ii, radar_uc_wc_op_ii, cov_uc_wc_op_ii,~] =...
        tsp_altermating_projection_wo_radar_wc(fdcomm_co_wc_ini_ii, radar_co_wc_ini_ii, radar_comm_co_wc_ini_ii,cov_co_wc_ini_ii);
    fdcomm_uc_wc_save{ii} = fdcomm_uc_wc_op_ii;
    radar_uc_wc_save{ii} = radar_uc_wc_op_ii;
    cov_uc_wc_save{ii} = cov_uc_wc_op_ii;
    %% randomly coded with co-op
    [fdcomm_rc_ini_ii, radar_rc_ini_ii, radar_comm_rc_ini_ii,cov_rc_ini_ii] =...
        tsp_ini_normal_rd_radar(radar_para_co_ii,fdcomm_para_co_ii,radar_comm_para_co_ii,snr_rtr);
    [fdcomm_rc_op_ii, radar_rc_op_ii, cov_rc_op_ii,radar_rc_ini_ii] =...
        tsp_altermating_projection_wo_radar(fdcomm_rc_ini_ii, radar_rc_ini_ii, radar_comm_rc_ini_ii,cov_rc_ini_ii);
    fdcomm_rc_save{ii} = fdcomm_rc_op_ii;
    radar_rc_save{ii} = radar_rc_op_ii;
    cov_rc_save{ii} = cov_rc_op_ii;
    %% randomly coded without co-op
    [fdcomm_rc_wc_ini_ii, radar_rc_wc_ini_ii, radar_comm_rc_wc_ini_ii,cov_rc_wc_ini_ii]...
        = tsp_ini_rd_radar_snr_wc(radar_rc_ini_ii,fdcomm_rc_ini_ii,radar_comm_rc_ini_ii,snr_rtr);
    [fdcomm_rc_wc_op_ii, radar_rc_wc_op_ii, cov_rc_wc_op_ii,~] =...
        tsp_altermating_projection_wo_radar_wc(fdcomm_rc_wc_ini_ii, radar_rc_wc_ini_ii, radar_comm_rc_wc_ini_ii,cov_rc_wc_ini_ii);
    fdcomm_rc_wc_save{ii} = fdcomm_rc_wc_op_ii;
    radar_rc_wc_save{ii} = radar_rc_wc_op_ii;
    cov_rc_wc_save{ii} = cov_rc_wc_op_ii;
end
save('tsp_radar_snr_comparison_v5.mat');
% function savetofile(data,fullfilename)
%     save(fullfilename,'data');
% end
 
