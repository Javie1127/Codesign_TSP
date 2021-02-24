%%% Convergence analysis for equal weights

%% Array Parameters
Mr = 4; % Number of radar TX antennas
Nr = 4; % Number of radar RX antennas
radar.TX = Mr;
radar.RX = Nr;      
radar.noisepower = 0.01;
K = 8;
radar.codelength = K;
%% FD Comm parameters
Mc = 4; % Number of BS TX antennas
Nc = 4;% Number of BS RX antennas
I = 2; % Number of UL UEs
J = 2; % Number of DL UEs
fdcomm.BSTX = Mc;
fdcomm.BSRX = Nc;
fdcomm.UL_num = I;
fdcomm.DL_num = J;
fdcomm.ULpower = ones(I,1);
fdcomm.DLpower = J;
%% Set the SNRs in dB
%SNR.rtr = randi([-5,5],Mr,Nr);
SNR.rtr = 5*ones(Mr,Nr);
SNR.Bmr = 1*ones(Nr,1);
SNR.Btr = 1*ones(Nr,1); 
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
radar.ell_max = 30; % algorithm 5
M = 1;
%% simulation parameters
[fdcomm_para, radar_para, radar_comm_para] = tsp_parameters_v1(SNR,radar,fdcomm);
SNR_rtr = -5:5:10;
I_total_ew_ni = zeros(radar.ell_max,length(SNR_rtr));
I_UL_ew_ni = zeros(radar.ell_max,length(SNR_rtr));
I_DL_ew_ni = zeros(radar.ell_max,length(SNR_rtr));
I_radar_ew_ni = zeros(radar.ell_max,length(SNR_rtr));
I_total_ew_ri = zeros(radar.ell_max,length(SNR_rtr));
I_UL_ew_ri = zeros(radar.ell_max,length(SNR_rtr));
I_DL_ew_ri = zeros(radar.ell_max,length(SNR_rtr));
I_radar_ew_ri = zeros(radar.ell_max,length(SNR_rtr));
% I_total_uew_ni = zeros(radar.ell_max,length(SNR_rtr));
% I_UL_uew_ni = zeros(radar.ell_max,length(SNR_rtr));
% I_DL_uew_ni = zeros(radar.ell_max,length(SNR_rtr));
% I_radar_uew_ni = zeros(radar.ell_max,length(SNR_rtr));
% I_total_uew_ri = zeros(radar.ell_max,length(SNR_rtr));
% I_UL_uew_ri = zeros(radar.ell_max,length(SNR_rtr));
% I_DL_uew_ri = zeros(radar.ell_max,length(SNR_rtr));
% I_radar_uew_ri = zeros(radar.ell_max,length(SNR_rtr));
% equal weights
fdcomm_ew_alpha_UL = 1/(Nr+I+J)*ones(I,1);
fdcomm_ew_alpha_DL = 1/(Nr+I+J)*ones(J,1);
radar_ew_alpha_r = 1/(Nr+I+J)*ones(Nr,1);
% % unequal weights 
% fdcomm_uew_alpha_UL = 0.1*ones(I,1);
% fdcomm_uew_alpha_DL = 0.05*ones(J,1);
% radar_uew_alpha_r = (1-0.1*I-0.05*J)/Nr*ones(Nr,1);
parfor ii = 1:length(SNR_rtr)
    snr_rtr = SNR_rtr(ii)*ones(Mr,Nr);
    [fdcomm_ii, radar_ii, radar_comm_ii] =...
        tsp_parameter_radarsnr(radar_para,fdcomm_para,radar_comm_para,snr_rtr);
    % equal weights normal initialization
    [fdcomm_ni_ii, radar_ni_ii,cov_ni_ii] =...
        tsp_ini_normal_v2(radar_ii,fdcomm_ii,radar_comm_ii);
    fdcomm_ew_ni_ii = fdcomm_ni_ii;
    radar_ew_ni_ii = radar_ni_ii;
    radar_comm_ew_ni_ii = radar_comm_ii;
    fdcomm_ew_ni_ii.alpha_UL = fdcomm_ew_alpha_UL;
    fdcomm_ew_ni_ii.alpha_DL = fdcomm_ew_alpha_DL;
    radar_ew_ni_ii.alpha_r = radar_ew_alpha_r;
    cov_ew_ni_ii = cov_ni_ii;
    % update mmse
    [fdcomm_ew_ni_ii,radar_ew_ni_ii] = tsp_ini_mmse(radar_ew_ni_ii,fdcomm_ew_ni_ii,cov_ni_ii);
    [fdcomm_ew_ni_op_ii, radar_ew_ni_op_ii, cov_ew_ni_op_ii,~] =...
        tsp_altermating_projection_v1(fdcomm_ew_ni_ii, radar_ew_ni_ii, radar_comm_ew_ni_ii,cov_ew_ni_ii);
    I_total_ew_ni(:,ii) = fdcomm_ew_ni_op_ii.I_total_op;
    I_UL_ew_ni(:,ii) = fdcomm_ew_ni_op_ii.I_UL_op;
    I_DL_ew_ni(:,ii) = fdcomm_ew_ni_op_ii.I_DL_op;
    I_radar_ew_ni(:,ii) = fdcomm_ew_ni_op_ii.I_radar_op;
%     % unequall weights normal initialization
%     fdcomm_uew_ni_ii = fdcomm_ni_ii;
%     radar_uew_ni_ii = radar_ni_ii;
%     radar_comm_uew_ni_ii = radar_comm_ii;
%     cov_uew_ni_ii = cov_ni_ii;
%     % update mmse with new weights
%     [fdcomm_uew_ni_ii,radar_uew_ni_ii] = tsp_ini_mmse(radar_uew_ni_ii,fdcomm_uew_ni_ii,cov_uew_ni_ii);
%     fdcomm_uew_ni_ii.alpha_UL = fdcomm_uew_alpha_UL;
%     fdcomm_uew_ni_ii.alpha_DL = fdcomm_uew_alpha_DL;
%     radar_uew_ni_ii.alpha_r = radar_uew_alpha_r;
%     [fdcomm_uew_ni_op_ii, radar_uew_ni_op_ii, cov_uew_ni_op_ii,~] =...
%         tsp_altermating_projection_v1(fdcomm_uew_ni_ii, radar_uew_ni_ii, radar_comm_uew_ni_ii,cov_uew_ni_ii);
%     I_total_uew_ni(:,ii) = fdcomm_uew_ni_op_ii.I_total_op;
%     I_UL_uew_ni(:,ii) = fdcomm_uew_ni_op_ii.I_UL_op;
%     I_DL_uew_ni(:,ii) = fdcomm_uew_ni_op_ii.I_DL_op;
%     I_radar_uew_ni(:,ii) = fdcomm_uew_ni_op_ii.I_radar_op;
%     
%     equal weights random initializations
    I_total_ew_ri_ii = zeros(radar_para.ell_max,M);
    I_UL_ew_ri_ii = zeros(radar_para.ell_max,M);
    I_DL_ew_ri_ii = zeros(radar_para.ell_max,M);
    I_radar_ew_ri_ii = zeros(radar_para.ell_max,M);
    % unequal weights
    I_total_uew_ri_ii = zeros(radar_para.ell_max,M);
    I_UL_uew_ri_ii = zeros(radar_para.ell_max,M);
    I_DL_uew_ri_ii = zeros(radar_para.ell_max,M);
    I_radar_uew_ri_ii = zeros(radar_para.ell_max,M);
    for m = 1:M
        % Random Initialization
        [fdcomm_ri_ii_m, radar_ri_ii_m,cov_ri_ii_m] =...
        tsp_ini_random_v2(radar_ii,fdcomm_ii,radar_comm_ii);
        % equal weights
        fdcomm_ew_ri_ii_m = fdcomm_ri_ii_m;
        radar_ew_ri_ii_m = radar_ri_ii_m;
        radar_comm_ew_ri_ii_m = radar_comm_ii;
        fdcomm_ew_ri_ii_m.alpha_DL = fdcomm_ew_alpha_DL;
        fdcomm_ew_ri_ii_m.alpha_UL = fdcomm_ew_alpha_UL;
        radar_ew_ri_ii_m.alpha_r = radar_ew_alpha_r;
        cov_ew_ri_ii_m = cov_ri_ii_m;
        [fdcomm_ew_ri_ii_m,radar_ew_ri_ii_m] = tsp_ini_mmse(radar_ew_ri_ii_m,fdcomm_ew_ri_ii_m,cov_ew_ri_ii_m);
        [fdcomm_ew_ri_op_ii_m, radar_ew_ri_op_m, cov_ew_ri_op_m,~] =...
            tsp_altermating_projection_v1(fdcomm_ew_ri_ii_m, radar_ew_ri_ii_m,radar_comm_ew_ri_ii_m,cov_ew_ri_ii_m);
%         fdcomm_m_op_ew_ri.I_total = I_total_m_ew_ri;
%         fdcomm_m_op_ew_ri.I_total_op = I_total_m_op_ew_ri;
%         cov_op = tsp_covmat(fdcomm_op,radar_op,radar_comm_para);
        I_total_ew_ri_ii(:,m) = fdcomm_ew_ri_op_ii_m.I_total_op;
        I_UL_ew_ri_ii(:,m) = fdcomm_ew_ri_op_ii_m.I_UL_op;
        I_DL_ew_ri_ii(:,m) = fdcomm_ew_ri_op_ii_m.I_DL_op;
        I_radar_ew_ri_ii(:,m) = fdcomm_ew_ri_op_ii_m.I_radar_op;
%         % unequal weights
%         fdcomm_uew_ri_ii_m = fdcomm_ri_ii_m;
%         radar_uew_ri_ii_m = radar_ri_ii_m;
%         radar_comm_uew_ri_ii_m = radar_comm_ii;
%         fdcomm_uew_ri_ii_m.alpha_DL = fdcomm_uew_alpha_DL;
%         fdcomm_uew_ri_ii_m.alpha_UL = fdcomm_uew_alpha_UL;
%         radar_uew_ri_ii_m.alpha_r = radar_uew_alpha_r;
%         cov_uew_ri_ii_m = cov_ri_ii_m;
%         [fdcomm_uew_ri_ii_m,radar_uew_ri_ii_m] = tsp_ini_mmse(radar_uew_ri_ii_m,fdcomm_uew_ri_ii_m,cov_uew_ri_ii_m);
%         [fdcomm_uew_ri_op_ii_m, radar_uew_ri_op_m, cov_uew_ri_op_m,~] =...
%             tsp_altermating_projection_v1(fdcomm_uew_ri_ii_m, radar_uew_ri_ii_m,radar_comm_uew_ri_ii_m,cov_uew_ri_ii_m);
% %         fdcomm_m_op_ew_ri.I_total = I_total_m_ew_ri;
% %         fdcomm_m_op_ew_ri.I_total_op = I_total_m_op_ew_ri;
% %         cov_op = tsp_covmat(fdcomm_op,radar_op,radar_comm_para);
%         I_total_uew_ri_ii(:,m) = fdcomm_uew_ri_op_ii_m.I_total_op;
%         I_UL_uew_ri_ii(:,m) = fdcomm_uew_ri_op_ii_m.I_UL_op;
%         I_DL_uew_ri_ii(:,m) = fdcomm_uew_ri_op_ii_m.I_DL_op;
%         I_radar_uew_ri_ii(:,m) = fdcomm_uew_ri_op_ii_m.I_radar_op;
    end
    I_total_ew_ri(:,ii) = mean(I_total_ew_ri_ii,2);
    I_UL_ew_ri(:,ii) = mean(I_UL_ew_ri_ii,2);
    I_DL_ew_ri(:,ii) = mean(I_DL_ew_ri_ii,2);
    I_radar_ew_ri(:,ii) = mean(I_radar_ew_ri_ii,2);
    
%     I_total_uew_ri(:,ii) = mean(I_total_uew_ri_ii,2);
%     I_UL_uew_ri(:,ii) = mean(I_UL_uew_ri_ii,2);
%     I_DL_uew_ri(:,ii) = mean(I_DL_uew_ri_ii,2);
%     I_radar_uew_ri(:,ii) = mean(I_radar_uew_ri_ii,2);
    

    
end
save('tsp_convergence_equalweights.mat')


