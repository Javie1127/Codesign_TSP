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
SNR.rtr = 0*ones(Mr,Nr);
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
%% priority weight equal
fdcomm.alpha_UL = 1/(Nr+I+J)*ones(I,1);
fdcomm.alpha_DL = 1/(Nr+I+J)*ones(J,1);
radar.alpha_r = 1/(Nr+I+J)*ones(Nr,1);

%% simulation parameters
[fdcomm_para_ew, radar_para_ew, radar_comm_para] = tsp_parameters(SNR,radar,fdcomm);
%% random initialization
M = 20;
I_total_ew_ri = zeros(radar.ell_max,M);
I_UL_ew_ri = zeros(radar.ell_max,M);
I_DL_ew_ri = zeros(radar.ell_max,M);
I_radar_ew_ri = zeros(radar.ell_max,M);
%snr_rtr = 0*ones(Mr,Nr);
parfor m = 1:M
    fdcomm_para_ew_m = fdcomm_para_ew;
    radar_para_ew_m = radar_para_ew;
    radar_comm_para_m = radar_comm_para;
    %% Random Initialization
    %[fdcomm_ew_ri_ini_m, radar_ew_ri_ini_m, cov_ew_ri_ini_m] = tsp_ini_random_v1(radar_para_ew_m,fdcomm_para_ew_m,radar_comm_para_m,snr_rtr);
    [fdcomm_ew_ri_ini_m, radar_ew_ri_ini_m, cov_ew_ri_ini_m] = tsp_ini_random(radar_para_ew_m,fdcomm_para_ew_m,radar_comm_para_m);
    [fdcomm_ew_ri_op_ii_m, radar_ew_ri_op_m, cov_ew_ri_op_m,~] =...
        tsp_altermating_projection_v1(fdcomm_ew_ri_ini_m, radar_ew_ri_ini_m,radar_comm_para_m,cov_ew_ri_ini_m);
    %fdcomm_m_op_ew_ri.I_total = I_total_m_ew_ri;
    %fdcomm_m_op_ew_ri.I_total_op = I_total_m_op_ew_ri;
    %cov_op = tsp_covmat(fdcomm_op,radar_op,radar_comm_para);
    I_total_ew_ri(:,m) = fdcomm_ew_ri_op_ii_m.I_total_op;
    I_UL_ew_ri(:,m) = fdcomm_ew_ri_op_ii_m.I_UL_op;
    I_DL_ew_ri(:,m) = fdcomm_ew_ri_op_ii_m.I_DL_op;
    I_radar_ew_ri(:,m) = fdcomm_ew_ri_op_ii_m.I_radar_op;
end
I_total_avg_ew_ri = mean(I_total_ew_ri,2);
I_UL_avg_ew_ri = mean(I_UL_ew_ri,2);
I_DL_avg_ew_ri = mean(I_DL_ew_ri,2);
I_radar_avg_ew_ri = mean(I_radar_ew_ri,2);
%% normal initialization equal weights
[fdcomm_ini_ew_ni, radar_ini_ew_ni, cov_ew_ni] = tsp_ini_normal(radar_para_ew,fdcomm_para_ew,radar_comm_para);
[fdcomm_ew_ni_op, radar_ew_ni_op, cov_ew_ni_op,~] =...
        tsp_altermating_projection_v1(fdcomm_ini_ew_ni, radar_ini_ew_ni,radar_comm_para,cov_ew_ni);
I_total_op_ew_ni = fdcomm_ew_ni_op.I_total_op;
I_UL_op_ew_ni = fdcomm_ew_ni_op.I_UL_op;
I_DL_op_ew_ni = fdcomm_ew_ni_op.I_DL_op;
I_radar_op_ew_ni = fdcomm_ew_ni_op.I_radar_op;
%%------------------------- unequal weights--------------------------------
fdcomm_para_uew = fdcomm_para_ew;
radar_para_uew = radar_para_ew;
%% priority unequal weight
fdcomm_para_uew.alpha_UL = 0.1*ones(I,1);
fdcomm_para_uew.alpha_DL = 0.05*ones(J,1);
radar_para_uew.alpha_r = (1-0.1*I-0.05*J)/Nr*ones(Nr,1);
%% random initialization
I_total_uew_ri = zeros(radar.ell_max,M);
I_UL_uew_ri = zeros(radar.ell_max,M);
I_DL_uew_ri = zeros(radar.ell_max,M);
I_radar_uew_ri = zeros(radar.ell_max,M);
parfor m = 1:M
    %snr_rtr = 5*ones(Mr,Nr);
    fdcomm_para_uew_m = fdcomm_para_uew;
    radar_para_uew_m = radar_para_uew;
    radar_comm_para_m = radar_comm_para;
    %% Random Initialization
    [fdcomm_uew_ri_ini_m, radar_uew_ri_ini_m, cov_uew_ri_ini_m] = tsp_ini_random(radar_para_uew_m,fdcomm_para_uew_m,radar_comm_para_m);
    [fdcomm_uew_ri_op_ii_m, radar_uew_ri_op_m, cov_uew_ri_op_m,~] =...
        tsp_altermating_projection_v1(fdcomm_uew_ri_ini_m, radar_uew_ri_ini_m,radar_comm_para_m,cov_uew_ri_ini_m);
    %fdcomm_m_op_ew_ri.I_total = I_total_m_ew_ri;
    %fdcomm_m_op_ew_ri.I_total_op = I_total_m_op_ew_ri;
    %cov_op = tsp_covmat(fdcomm_op,radar_op,radar_comm_para);
    I_total_uew_ri(:,m) = fdcomm_uew_ri_op_ii_m.I_total_op;
    I_UL_uew_ri(:,m) = fdcomm_uew_ri_op_ii_m.I_UL_op;
    I_DL_uew_ri(:,m) = fdcomm_uew_ri_op_ii_m.I_DL_op;
    I_radar_uew_ri(:,m) = fdcomm_uew_ri_op_ii_m.I_radar_op;
end
I_total_avg_uew_ri = mean(I_total_uew_ri,2);
I_UL_avg_uew_ri = mean(I_UL_uew_ri,2);
I_DL_avg_uew_ri = mean(I_DL_uew_ri,2);
I_radar_avg_uew_ri = mean(I_radar_uew_ri,2);
%% normal initialization
[fdcomm_uew_ni_ini, radar_uew_ni_ini, cov_uew_ni_ini] = tsp_ini_normal(radar_para_uew,fdcomm_para_uew,radar_comm_para);
[fdcomm_uew_ni, radar_uew_ni, cov_uew_ni,~] =...
        tsp_altermating_projection_v1(fdcomm_uew_ni_ini, radar_uew_ni_ini, radar_comm_para,cov_uew_ni_ini);
I_total_op_uew_ni = fdcomm_uew_ni.I_total_op;
I_UL_op_uew_ni = fdcomm_uew_ni.I_UL_op;
I_DL_op_uew_ni = fdcomm_uew_ni.I_DL_op;
I_radar_op_uew_ni = fdcomm_uew_ni.I_radar_op;
save('tsp_convergence_0dB.mat')


